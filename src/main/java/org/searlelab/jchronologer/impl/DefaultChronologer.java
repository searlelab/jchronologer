package org.searlelab.jchronologer.impl;

import java.util.ArrayList;
import java.util.List;
import org.searlelab.jchronologer.api.AcceptedPrediction;
import org.searlelab.jchronologer.api.Chronologer;
import org.searlelab.jchronologer.api.ChronologerOptions;
import org.searlelab.jchronologer.api.PredictionResult;
import org.searlelab.jchronologer.api.RejectedPrediction;
import org.searlelab.jchronologer.inference.BatchPredictor;
import org.searlelab.jchronologer.preprocessing.ChronologerPreprocessor;
import org.searlelab.jchronologer.preprocessing.PreprocessingMetadataLoader;
import org.searlelab.jchronologer.preprocessing.PreprocessingOutcome;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Default {@link Chronologer} implementation that orchestrates preprocessing and batched
 * inference.
 *
 * <p>This class is intentionally stateful: model and metadata resources are loaded once at
 * construction time, then reused across repeated {@link #predict(List)} calls until
 * {@link #close()}.
 */
public final class DefaultChronologer implements Chronologer {

    private static final Logger LOGGER = LoggerFactory.getLogger(DefaultChronologer.class);

    private final ChronologerOptions options;
    private final ChronologerPreprocessor preprocessor;
    private final BatchPredictor batchPredictor;

    public DefaultChronologer(ChronologerOptions options) {
        this.options = options;
        PreprocessingMetadataLoader.CompiledPreprocessingMetadata metadata =
                PreprocessingMetadataLoader.loadFromClasspath(options.getPreprocessingResource());
        this.preprocessor = new ChronologerPreprocessor(metadata);
        this.batchPredictor = new BatchPredictor(options.getModelResource());
        LOGGER.info("Loaded Chronologer model from classpath resource {}", options.getModelResource());
    }

    /**
     * Runs end-to-end Chronologer inference for a batch of peptide sequences.
     *
     * <p>Each row is first preprocessed independently; accepted rows are scored in configured
     * batch chunks and rejected rows are returned with diagnostics.
     */
    @Override
    public PredictionResult predict(List<String> peptideModSeqs) {
        List<AcceptedDraft> acceptedDrafts = new ArrayList<>();
        List<AcceptedPrediction> accepted = new ArrayList<>();
        List<RejectedPrediction> rejected = new ArrayList<>();

        for (int rowIndex = 0; rowIndex < peptideModSeqs.size(); rowIndex++) {
            String peptide = peptideModSeqs.get(rowIndex);
            PreprocessingOutcome outcome = preprocessor.preprocess(peptide);
            if (outcome.isAccepted()) {
                acceptedDrafts.add(new AcceptedDraft(
                        rowIndex,
                        peptide,
                        outcome.getPatchedPeptideModSeq(),
                        outcome.getCodedPeptideSeq(),
                        outcome.getTokenArray()));
            } else {
                rejected.add(new RejectedPrediction(
                        rowIndex,
                        peptide,
                        outcome.getPatchedPeptideModSeq(),
                        outcome.getRejectionReason(),
                        outcome.getErrorDetail()));
            }
        }

        int batchSize = options.getBatchSize();
        for (int i = 0; i < acceptedDrafts.size(); i += batchSize) {
            int end = Math.min(i + batchSize, acceptedDrafts.size());
            long[][] batchTokens = new long[end - i][];
            for (int j = i; j < end; j++) {
                batchTokens[j - i] = acceptedDrafts.get(j).tokenArray;
            }

            float[] batchPredictions = batchPredictor.predict(batchTokens);
            for (int j = i; j < end; j++) {
                AcceptedDraft draft = acceptedDrafts.get(j);
                accepted.add(new AcceptedPrediction(
                        draft.rowIndex,
                        draft.peptideModSeq,
                        draft.patchedPeptideModSeq,
                        draft.codedPeptideSeq,
                        draft.tokenArray,
                        batchPredictions[j - i]));
            }
        }

        return new PredictionResult(accepted, rejected);
    }

    @Override
    public void close() {
        batchPredictor.close();
    }

    /**
     * Internal accepted-row payload retained between preprocessing and model scoring.
     */
    private static final class AcceptedDraft {
        private final int rowIndex;
        private final String peptideModSeq;
        private final String patchedPeptideModSeq;
        private final String codedPeptideSeq;
        private final long[] tokenArray;

        private AcceptedDraft(
                int rowIndex,
                String peptideModSeq,
                String patchedPeptideModSeq,
                String codedPeptideSeq,
                long[] tokenArray) {
            this.rowIndex = rowIndex;
            this.peptideModSeq = peptideModSeq;
            this.patchedPeptideModSeq = patchedPeptideModSeq;
            this.codedPeptideSeq = codedPeptideSeq;
            this.tokenArray = tokenArray;
        }
    }
}

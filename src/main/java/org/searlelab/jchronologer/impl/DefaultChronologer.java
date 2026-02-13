package org.searlelab.jchronologer.impl;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;
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
    private final ExecutorService inferenceExecutor;

    public DefaultChronologer(ChronologerOptions options) {
        this.options = options;
        PreprocessingMetadataLoader.CompiledPreprocessingMetadata metadata =
                PreprocessingMetadataLoader.loadFromClasspath(options.getPreprocessingResource());
        this.preprocessor = new ChronologerPreprocessor(metadata);
        this.batchPredictor = new BatchPredictor(options.getModelResource());
        this.inferenceExecutor = options.getInferenceThreads() > 1
                ? newInferenceExecutor(options.getInferenceThreads())
                : null;
        LOGGER.info(
                "Loaded Chronologer model from classpath resource {} (batchSize={}, inferenceThreads={})",
                options.getModelResource(),
                options.getBatchSize(),
                options.getInferenceThreads());
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

        float[] predictionsByAcceptedIndex = scoreAcceptedDrafts(acceptedDrafts);
        for (int i = 0; i < acceptedDrafts.size(); i++) {
            AcceptedDraft draft = acceptedDrafts.get(i);
            accepted.add(new AcceptedPrediction(
                    draft.rowIndex,
                    draft.peptideModSeq,
                    draft.patchedPeptideModSeq,
                    draft.codedPeptideSeq,
                    draft.tokenArray,
                    predictionsByAcceptedIndex[i]));
        }

        return new PredictionResult(accepted, rejected);
    }

    @Override
    public void close() {
        shutdownInferenceExecutor();
        batchPredictor.close();
    }

    private float[] scoreAcceptedDrafts(List<AcceptedDraft> acceptedDrafts) {
        float[] predictionsByAcceptedIndex = new float[acceptedDrafts.size()];
        if (acceptedDrafts.isEmpty()) {
            return predictionsByAcceptedIndex;
        }

        List<BatchRequest> batchRequests = partitionIntoBatches(acceptedDrafts, options.getBatchSize());
        if (inferenceExecutor == null || batchRequests.size() == 1) {
            for (BatchRequest request : batchRequests) {
                float[] batchPredictions = batchPredictor.predict(request.tokenBatch);
                System.arraycopy(
                        batchPredictions,
                        0,
                        predictionsByAcceptedIndex,
                        request.startIndex,
                        batchPredictions.length);
            }
            return predictionsByAcceptedIndex;
        }

        List<Future<BatchResult>> futures = new ArrayList<>(batchRequests.size());
        for (BatchRequest request : batchRequests) {
            futures.add(inferenceExecutor.submit(
                    () -> new BatchResult(request.startIndex, batchPredictor.predict(request.tokenBatch))));
        }

        try {
            for (Future<BatchResult> future : futures) {
                BatchResult batchResult = future.get();
                System.arraycopy(
                        batchResult.predictions,
                        0,
                        predictionsByAcceptedIndex,
                        batchResult.startIndex,
                        batchResult.predictions.length);
            }
            return predictionsByAcceptedIndex;
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            futures.forEach(future -> future.cancel(true));
            throw new IllegalStateException("Chronologer inference was interrupted.", e);
        } catch (ExecutionException e) {
            futures.forEach(future -> future.cancel(true));
            Throwable cause = e.getCause() == null ? e : e.getCause();
            throw new IllegalStateException("Failed to run Chronologer inference.", cause);
        }
    }

    private static List<BatchRequest> partitionIntoBatches(List<AcceptedDraft> acceptedDrafts, int batchSize) {
        List<BatchRequest> batchRequests = new ArrayList<>();
        for (int i = 0; i < acceptedDrafts.size(); i += batchSize) {
            int end = Math.min(i + batchSize, acceptedDrafts.size());
            long[][] batchTokens = new long[end - i][];
            for (int j = i; j < end; j++) {
                batchTokens[j - i] = acceptedDrafts.get(j).tokenArray;
            }
            batchRequests.add(new BatchRequest(i, batchTokens));
        }
        return batchRequests;
    }

    private static ExecutorService newInferenceExecutor(int threads) {
        AtomicInteger threadCounter = new AtomicInteger(1);
        return Executors.newFixedThreadPool(threads, runnable -> {
            Thread thread = new Thread(runnable);
            thread.setName("chronologer-inference-" + threadCounter.getAndIncrement());
            thread.setDaemon(true);
            return thread;
        });
    }

    private void shutdownInferenceExecutor() {
        if (inferenceExecutor == null) {
            return;
        }
        inferenceExecutor.shutdown();
        try {
            if (!inferenceExecutor.awaitTermination(30, TimeUnit.SECONDS)) {
                inferenceExecutor.shutdownNow();
            }
        } catch (InterruptedException e) {
            inferenceExecutor.shutdownNow();
            Thread.currentThread().interrupt();
        }
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

    private static final class BatchRequest {
        private final int startIndex;
        private final long[][] tokenBatch;

        private BatchRequest(int startIndex, long[][] tokenBatch) {
            this.startIndex = startIndex;
            this.tokenBatch = tokenBatch;
        }
    }

    private static final class BatchResult {
        private final int startIndex;
        private final float[] predictions;

        private BatchResult(int startIndex, float[] predictions) {
            this.startIndex = startIndex;
            this.predictions = predictions;
        }
    }
}

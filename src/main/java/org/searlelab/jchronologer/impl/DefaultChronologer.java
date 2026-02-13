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
 * construction time, then reused across repeated inference calls until {@link #close()}.
 */
public final class DefaultChronologer implements Chronologer {

    private static final Logger LOGGER = LoggerFactory.getLogger(DefaultChronologer.class);
    private static final String DJL_PLATFORM_LOG_LEVEL_PROPERTY = "org.slf4j.simpleLogger.log.ai.djl.util.Platform";
    private static final String DJL_PT_ENGINE_LOG_LEVEL_PROPERTY =
            "org.slf4j.simpleLogger.log.ai.djl.pytorch.engine.PtEngine";

    private final ChronologerOptions options;
    private final ChronologerPreprocessor preprocessor;
    private final BatchPredictor batchPredictor;
    private final ExecutorService inferenceExecutor;
    private volatile boolean closed;

    public DefaultChronologer(ChronologerOptions options) {
        long initStart = System.nanoTime();
        this.options = options;

        configureDjlLogging(options.isVerboseLogging());

        LOGGER.info(
                "Initializing Chronologer (modelResource={}, preprocessingResource={}, batchSize={}, inferenceThreads={})",
                options.getModelResource(),
                options.getPreprocessingResource(),
                options.getBatchSize(),
                options.getInferenceThreads());

        long metadataLoadStart = System.nanoTime();
        PreprocessingMetadataLoader.CompiledPreprocessingMetadata metadata =
                PreprocessingMetadataLoader.loadFromClasspath(options.getPreprocessingResource());
        long metadataLoadMillis = elapsedMillis(metadataLoadStart);
        logVerbose(
                "Loaded preprocessing metadata resource {} in {} ms",
                options.getPreprocessingResource(),
                metadataLoadMillis);

        long preprocessorCreateStart = System.nanoTime();
        this.preprocessor = new ChronologerPreprocessor(metadata);
        long preprocessorCreateMillis = elapsedMillis(preprocessorCreateStart);
        logVerbose("Constructed preprocessor in {} ms", preprocessorCreateMillis);

        long predictorInitStart = System.nanoTime();
        this.batchPredictor = new BatchPredictor(options.getModelResource(), options.isVerboseLogging());
        long predictorInitMillis = elapsedMillis(predictorInitStart);

        long executorInitStart = System.nanoTime();
        this.inferenceExecutor = options.getInferenceThreads() > 1
                ? newInferenceExecutor(options.getInferenceThreads())
                : null;
        long executorInitMillis = elapsedMillis(executorInitStart);

        logVerbose(
                "Loaded Chronologer model from classpath resource {} (batchSize={}, inferenceThreads={})",
                options.getModelResource(),
                options.getBatchSize(),
                options.getInferenceThreads());
        LOGGER.info(
                "Chronologer initialization complete in {} ms (metadataLoadMs={}, preprocessorInitMs={}, modelInitMs={}, executorInitMs={})",
                elapsedMillis(initStart),
                metadataLoadMillis,
                preprocessorCreateMillis,
                predictorInitMillis,
                executorInitMillis);
    }

    @Override
    public void init() {
        if (closed) {
            throw new IllegalStateException("Chronologer has been closed.");
        }
    }

    /**
     * Runs end-to-end Chronologer inference for a batch of peptide sequences.
     *
     * <p>Each row is first preprocessed independently; accepted rows are scored in configured
     * batch chunks and rejected rows are returned with diagnostics.
     */
    @Override
    public PredictionResult predict(List<String> peptideModSeqs) {
        init();

        long predictStart = System.nanoTime();
        List<AcceptedDraft> acceptedDrafts = new ArrayList<>();
        List<AcceptedPrediction> accepted = new ArrayList<>();
        List<RejectedPrediction> rejected = new ArrayList<>();

        long preprocessingStart = System.nanoTime();
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
        long preprocessingMillis = elapsedMillis(preprocessingStart);

        long inferenceStart = System.nanoTime();
        float[] predictionsByAcceptedIndex =
                scoreAcceptedDrafts(acceptedDrafts, batchPredictor, inferenceExecutor);
        long inferenceMillis = elapsedMillis(inferenceStart);

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

        long totalMillis = elapsedMillis(predictStart);
        LOGGER.info(
                "Predict completed: inputRows={}, accepted={}, rejected={}, preprocessMs={}, inferenceMs={}, totalMs={}",
                peptideModSeqs.size(),
                accepted.size(),
                rejected.size(),
                preprocessingMillis,
                inferenceMillis,
                totalMillis);

        return new PredictionResult(accepted, rejected);
    }

    @Override
    public void close() {
        if (closed) {
            return;
        }
        closed = true;
        shutdownInferenceExecutor(inferenceExecutor);
        batchPredictor.close();
    }

    private float[] scoreAcceptedDrafts(
            List<AcceptedDraft> acceptedDrafts,
            BatchPredictor predictor,
            ExecutorService executor) {
        float[] predictionsByAcceptedIndex = new float[acceptedDrafts.size()];
        if (acceptedDrafts.isEmpty()) {
            return predictionsByAcceptedIndex;
        }

        List<BatchRequest> batchRequests = partitionIntoBatches(acceptedDrafts, options.getBatchSize());
        if (executor == null || batchRequests.size() == 1) {
            for (BatchRequest request : batchRequests) {
                float[] batchPredictions = predictor.predict(request.tokenBatch);
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
            futures.add(executor.submit(
                    () -> new BatchResult(request.startIndex, predictor.predict(request.tokenBatch))));
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

    private static void shutdownInferenceExecutor(ExecutorService executor) {
        if (executor == null) {
            return;
        }
        executor.shutdown();
        try {
            if (!executor.awaitTermination(30, TimeUnit.SECONDS)) {
                executor.shutdownNow();
            }
        } catch (InterruptedException e) {
            executor.shutdownNow();
            Thread.currentThread().interrupt();
        }
    }

    private static long elapsedMillis(long startNanos) {
        return (System.nanoTime() - startNanos) / 1_000_000L;
    }

    private static void configureDjlLogging(boolean verbose) {
        if (verbose) {
            System.setProperty(DJL_PLATFORM_LOG_LEVEL_PROPERTY, "info");
            System.setProperty(DJL_PT_ENGINE_LOG_LEVEL_PROPERTY, "info");
        } else {
            System.setProperty(DJL_PLATFORM_LOG_LEVEL_PROPERTY, "warn");
            System.setProperty(DJL_PT_ENGINE_LOG_LEVEL_PROPERTY, "warn");
        }
    }

    private void logVerbose(String message, Object... args) {
        if (options.isVerboseLogging()) {
            LOGGER.info(message, args);
        } else {
            LOGGER.debug(message, args);
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

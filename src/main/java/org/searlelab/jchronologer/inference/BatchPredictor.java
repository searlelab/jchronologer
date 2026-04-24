package org.searlelab.jchronologer.inference;

import ai.djl.Model;
import ai.djl.MalformedModelException;
import ai.djl.engine.EngineException;
import ai.djl.inference.Predictor;
import ai.djl.ndarray.NDArray;
import ai.djl.ndarray.NDList;
import ai.djl.ndarray.NDManager;
import ai.djl.translate.NoopTranslator;
import ai.djl.translate.TranslateException;
import java.io.IOException;
import java.nio.file.Path;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import org.searlelab.jchronologer.util.ResourceUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Thin DJL wrapper for running Chronologer TorchScript inference on token batches.
 *
 * <p>The model is loaded once from a classpath resource copied to a temporary file, then reused
 * for subsequent prediction calls.
 */
public final class BatchPredictor implements AutoCloseable {

    private static final Logger LOGGER = LoggerFactory.getLogger(BatchPredictor.class);

    private final Model model;
    private final Set<Predictor<NDList, NDList>> predictors;
    private final ThreadLocal<Predictor<NDList, NDList>> threadLocalPredictor;
    private final boolean verboseLogging;

    public BatchPredictor(String modelResource) {
        this(modelResource, false);
    }

    public BatchPredictor(String modelResource, boolean verboseLogging) {
        this.verboseLogging = verboseLogging;
        long initStart = System.nanoTime();
        logVerbose("Starting BatchPredictor initialization for model resource {}", modelResource);

        long copyStart = System.nanoTime();
        Path modelPath = ResourceUtils.copyClasspathResourceToTempFile(modelResource, ".torchscript.pt");
        logVerbose(
                "Copied model resource {} to temporary file {} in {} ms",
                modelResource,
                modelPath,
                elapsedMillis(copyStart));

        LOGGER.info(
                "Found matching platform from current JVM: {}-{}",
                System.getProperty("os.name"),
                System.getProperty("os.arch"));

        try {
            TorchModelLoader modelLoader = new TorchModelLoader();
            long modelInitStart = System.nanoTime();
            TorchModelLoader.LoadedModel loadedModel = modelLoader.load("Chronologer", modelPath);
            this.model = loadedModel.model();
            logVerbose(
                    "Loaded TorchScript model {} on {} in {} ms",
                    modelPath.getFileName(),
                    loadedModel.device(),
                    elapsedMillis(modelInitStart));

            this.predictors = ConcurrentHashMap.newKeySet();
            this.threadLocalPredictor = ThreadLocal.withInitial(() -> {
                Predictor<NDList, NDList> predictor = model.newPredictor(new NoopTranslator());
                predictors.add(predictor);
                return predictor;
            });

            long predictorInitStart = System.nanoTime();
            Predictor<NDList, NDList> initialPredictor = model.newPredictor(new NoopTranslator());
            predictors.add(initialPredictor);
            threadLocalPredictor.set(initialPredictor);
            logVerbose("Created initial predictor in {} ms", elapsedMillis(predictorInitStart));

            logVerbose("BatchPredictor initialization complete in {} ms", elapsedMillis(initStart));
        } catch (IOException | MalformedModelException | EngineException e) {
            throw new IllegalStateException("Failed to load Chronologer model from resource: " + modelResource, e);
        }
    }

    /**
     * Executes model inference for a batch of tokenized peptides.
     *
     * @param tokenBatch token matrix with shape {@code [batch, maxLen + 2]}
     * @return one {@code Pred_HI} value per input row
     */
    public float[] predict(long[][] tokenBatch) {
        try (NDManager manager = model.getNDManager().newSubManager()) {
            Predictor<NDList, NDList> predictor = threadLocalPredictor.get();
            NDArray input = manager.create(tokenBatch);
            NDList output = predictor.predict(new NDList(input));
            return output.singletonOrThrow().toFloatArray();
        } catch (TranslateException e) {
            throw new IllegalStateException("Failed to run Chronologer inference.", e);
        }
    }

    @Override
    public void close() {
        for (Predictor<NDList, NDList> predictor : predictors) {
            predictor.close();
        }
        model.close();
    }

    private static long elapsedMillis(long startNanos) {
        return (System.nanoTime() - startNanos) / 1_000_000L;
    }

    private void logVerbose(String message, Object... args) {
        if (verboseLogging) {
            LOGGER.info(message, args);
        } else {
            LOGGER.debug(message, args);
        }
    }
}

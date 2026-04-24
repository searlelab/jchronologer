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

/**
 * DJL wrapper for Cartographer TorchScript inference with 3 inputs:
 * tokens, precursor charge one-hot, and NCE.
 */
public final class CartographerBatchPredictor implements AutoCloseable {

    private final Model model;
    private final Set<Predictor<NDList, NDList>> predictors;
    private final ThreadLocal<Predictor<NDList, NDList>> threadLocalPredictor;

    public CartographerBatchPredictor(String modelResource) {
        Path modelPath = ResourceUtils.copyClasspathResourceToTempFile(modelResource, ".torchscript.pt");
        try {
            TorchModelLoader.LoadedModel loadedModel = new TorchModelLoader().load("Cartographer", modelPath);
            this.model = loadedModel.model();

            this.predictors = ConcurrentHashMap.newKeySet();
            this.threadLocalPredictor = ThreadLocal.withInitial(() -> {
                Predictor<NDList, NDList> predictor = model.newPredictor(new NoopTranslator());
                predictors.add(predictor);
                return predictor;
            });

            Predictor<NDList, NDList> initialPredictor = model.newPredictor(new NoopTranslator());
            predictors.add(initialPredictor);
            threadLocalPredictor.set(initialPredictor);
        } catch (IOException | MalformedModelException | EngineException e) {
            throw new IllegalStateException("Failed to load Cartographer model from resource: " + modelResource, e);
        }
    }

    /**
     * Executes one Cartographer batch.
     */
    public float[][] predict(long[][] tokenBatch, float[][] chargeBatch, float[][] nceBatch) {
        if (tokenBatch.length != chargeBatch.length || tokenBatch.length != nceBatch.length) {
            throw new IllegalArgumentException("Cartographer batch inputs must have matching batch size.");
        }
        if (tokenBatch.length == 0) {
            return new float[0][CartographerSpectrumDecoder.VECTOR_LENGTH];
        }

        try (NDManager manager = model.getNDManager().newSubManager()) {
            Predictor<NDList, NDList> predictor = threadLocalPredictor.get();
            NDArray tokens = manager.create(tokenBatch);
            NDArray charges = manager.create(chargeBatch);
            NDArray nce = manager.create(nceBatch);

            NDList output = predictor.predict(new NDList(tokens, charges, nce));
            return reshapeOutput(output.singletonOrThrow().toFloatArray(), tokenBatch.length, "Cartographer");
        } catch (TranslateException e) {
            throw new IllegalStateException("Failed to run Cartographer inference.", e);
        }
    }

    /**
     * Executes one two-input batch (tokens + charge one-hot), used by Sculptor.
     */
    public float[][] predict(long[][] tokenBatch, float[][] chargeBatch) {
        if (tokenBatch.length != chargeBatch.length) {
            throw new IllegalArgumentException("Cartographer batch inputs must have matching batch size.");
        }
        if (tokenBatch.length == 0) {
            return new float[0][1];
        }

        try (NDManager manager = model.getNDManager().newSubManager()) {
            Predictor<NDList, NDList> predictor = threadLocalPredictor.get();
            NDArray tokens = manager.create(tokenBatch);
            NDArray charges = manager.create(chargeBatch);

            NDList output = predictor.predict(new NDList(tokens, charges));
            return reshapeOutput(output.singletonOrThrow().toFloatArray(), tokenBatch.length, "Sculptor");
        } catch (TranslateException e) {
            throw new IllegalStateException("Failed to run Sculptor inference.", e);
        }
    }

    private static float[][] reshapeOutput(float[] flattened, int batchSize, String modelName) {
        int width = flattened.length / batchSize;
        if (width * batchSize != flattened.length) {
            throw new IllegalStateException("Unexpected " + modelName + " output shape.");
        }

        float[][] reshaped = new float[batchSize][width];
        for (int i = 0; i < batchSize; i++) {
            System.arraycopy(flattened, i * width, reshaped[i], 0, width);
        }
        return reshaped;
    }

    @Override
    public void close() {
        for (Predictor<NDList, NDList> predictor : predictors) {
            predictor.close();
        }
        model.close();
    }
}

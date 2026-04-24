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
 * DJL wrapper for Electrician TorchScript inference with 1 input: tokens.
 */
public final class ElectricianBatchPredictor implements AutoCloseable {

    private static final int CHARGE_STATE_COUNT = 6;

    private final Model model;
    private final Set<Predictor<NDList, NDList>> predictors;
    private final ThreadLocal<Predictor<NDList, NDList>> threadLocalPredictor;

    public ElectricianBatchPredictor(String modelResource) {
        Path modelPath = ResourceUtils.copyClasspathResourceToTempFile(modelResource, ".torchscript.pt");
        try {
            TorchModelLoader.LoadedModel loadedModel = new TorchModelLoader().load("Electrician", modelPath);
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
            throw new IllegalStateException("Failed to load Electrician model from resource: " + modelResource, e);
        }
    }

    /**
     * Executes one Electrician batch.
     */
    public float[][] predict(long[][] tokenBatch) {
        if (tokenBatch.length == 0) {
            return new float[0][CHARGE_STATE_COUNT];
        }

        try (NDManager manager = model.getNDManager().newSubManager()) {
            Predictor<NDList, NDList> predictor = threadLocalPredictor.get();
            NDArray tokens = manager.create(tokenBatch);

            NDList output = predictor.predict(new NDList(tokens));
            NDArray out = output.singletonOrThrow();
            float[] flattened = out.toFloatArray();

            int batchSize = tokenBatch.length;
            int width = flattened.length / batchSize;
            if (width * batchSize != flattened.length) {
                throw new IllegalStateException("Unexpected Electrician output shape.");
            }

            float[][] reshaped = new float[batchSize][width];
            for (int i = 0; i < batchSize; i++) {
                System.arraycopy(flattened, i * width, reshaped[i], 0, width);
            }
            return reshaped;
        } catch (TranslateException e) {
            throw new IllegalStateException("Failed to run Electrician inference.", e);
        }
    }

    @Override
    public void close() {
        for (Predictor<NDList, NDList> predictor : predictors) {
            predictor.close();
        }
        model.close();
    }
}

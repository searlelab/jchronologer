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
 * DJL wrapper for Scout TorchScript inference with 3 inputs and 4 outputs.
 */
public final class ScoutBatchPredictor implements AutoCloseable {

    private final Model model;
    private final Set<Predictor<NDList, NDList>> predictors;
    private final ThreadLocal<Predictor<NDList, NDList>> threadLocalPredictor;

    public ScoutBatchPredictor(String modelResource) {
        Path modelPath = ResourceUtils.copyClasspathResourceToTempFile(modelResource, ".torchscript.pt");
        try {
        	this.model = Model.newInstance("scout", "PyTorch");
            this.model.load(modelPath.getParent(), modelPath.getFileName().toString());
            
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
            throw new IllegalStateException("Failed to load Scout model from resource: " + modelResource, e);
        }
    }

    public ScoutBatchPredictionResult predict(long[][] tokenBatch, float[][] chargeBatch, float[][] nceBatch) {
        if (tokenBatch.length != chargeBatch.length || tokenBatch.length != nceBatch.length) {
            throw new IllegalArgumentException("Scout batch inputs must have matching batch size.");
        }
        if (tokenBatch.length == 0) {
            return new ScoutBatchPredictionResult(new float[0][0], new float[0][0], new float[0][0], new float[0][0]);
        }

        try (NDManager manager = model.getNDManager().newSubManager()) {
            Predictor<NDList, NDList> predictor = threadLocalPredictor.get();
            NDArray tokens = manager.create(tokenBatch);
            NDArray charges = manager.create(chargeBatch);
            NDArray nce = manager.create(nceBatch);
            NDList output = predictor.predict(new NDList(tokens, charges, nce));
            if (output.size() != 4) {
                throw new IllegalStateException("Unexpected Scout output count: " + output.size() + " expected=4");
            }
            return new ScoutBatchPredictionResult(
                    reshapeOutput(output.get(0).toFloatArray(), tokenBatch.length, "Scout ms2"),
                    reshapeOutput(output.get(1).toFloatArray(), tokenBatch.length, "Scout irt"),
                    reshapeOutput(output.get(2).toFloatArray(), tokenBatch.length, "Scout ccs"),
                    reshapeOutput(output.get(3).toFloatArray(), tokenBatch.length, "Scout charge_dist"));
        } catch (TranslateException e) {
            throw new IllegalStateException("Failed to run Scout inference.", e);
        }
    }

    private static float[][] reshapeOutput(float[] flattened, int batchSize, String outputName) {
        int width = flattened.length / batchSize;
        if (width * batchSize != flattened.length) {
            throw new IllegalStateException("Unexpected " + outputName + " shape.");
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

    public static final class ScoutBatchPredictionResult {
        private final float[][] ms2;
        private final float[][] irt;
        private final float[][] ccs;
        private final float[][] chargeDist;

        public ScoutBatchPredictionResult(float[][] ms2, float[][] irt, float[][] ccs, float[][] chargeDist) {
            this.ms2 = ms2;
            this.irt = irt;
            this.ccs = ccs;
            this.chargeDist = chargeDist;
        }

        public float[][] getMs2() {
            return ms2;
        }

        public float[][] getIrt() {
            return irt;
        }

        public float[][] getCcs() {
            return ccs;
        }

        public float[][] getChargeDist() {
            return chargeDist;
        }
    }
}

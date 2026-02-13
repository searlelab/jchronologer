package org.searlelab.jchronologer.inference;

import ai.djl.Model;
import ai.djl.MalformedModelException;
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
 * Thin DJL wrapper for running Chronologer TorchScript inference on token batches.
 *
 * <p>The model is loaded once from a classpath resource copied to a temporary file, then reused
 * for subsequent prediction calls.
 */
public final class BatchPredictor implements AutoCloseable {

    private final Model model;
    private final Set<Predictor<NDList, NDList>> predictors;
    private final ThreadLocal<Predictor<NDList, NDList>> threadLocalPredictor;

    public BatchPredictor(String modelResource) {
        Path modelPath = ResourceUtils.copyClasspathResourceToTempFile(modelResource, ".torchscript.pt");
        try {
            this.model = Model.newInstance("chronologer", "PyTorch");
            this.model.load(modelPath.getParent(), modelPath.getFileName().toString());
            this.predictors = ConcurrentHashMap.newKeySet();
            this.threadLocalPredictor = ThreadLocal.withInitial(() -> {
                Predictor<NDList, NDList> predictor = model.newPredictor(new NoopTranslator());
                predictors.add(predictor);
                return predictor;
            });
        } catch (IOException | MalformedModelException e) {
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
}

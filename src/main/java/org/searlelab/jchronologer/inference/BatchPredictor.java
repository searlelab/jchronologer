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
import org.searlelab.jchronologer.util.ResourceUtils;

public final class BatchPredictor implements AutoCloseable {

    private final Model model;
    private final Predictor<NDList, NDList> predictor;

    public BatchPredictor(String modelResource) {
        Path modelPath = ResourceUtils.copyClasspathResourceToTempFile(modelResource, ".torchscript.pt");
        try {
            this.model = Model.newInstance("chronologer", "PyTorch");
            this.model.load(modelPath.getParent(), modelPath.getFileName().toString());
            this.predictor = model.newPredictor(new NoopTranslator());
        } catch (IOException | MalformedModelException e) {
            throw new IllegalStateException("Failed to load Chronologer model from resource: " + modelResource, e);
        }
    }

    public float[] predict(long[][] tokenBatch) {
        try (NDManager manager = model.getNDManager().newSubManager()) {
            NDArray input = manager.create(tokenBatch);
            NDList output = predictor.predict(new NDList(input));
            return output.singletonOrThrow().toFloatArray();
        } catch (TranslateException e) {
            throw new IllegalStateException("Failed to run Chronologer inference.", e);
        }
    }

    @Override
    public void close() {
        predictor.close();
        model.close();
    }
}

package org.searlelab.jchronologer.inference;

import ai.djl.Device;
import ai.djl.MalformedModelException;
import ai.djl.Model;
import ai.djl.engine.Engine;
import java.io.IOException;
import java.nio.file.Path;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Loads TorchScript models with automatic CUDA detection and CPU fallback.
 */
final class TorchModelLoader {

    private static final Logger LOGGER = LoggerFactory.getLogger(TorchModelLoader.class);
    private static final String PYTORCH_ENGINE = "PyTorch";

    private final EngineSupplier engineSupplier;
    private final ModelFactory modelFactory;

    TorchModelLoader() {
        this(
                () -> Engine.getEngine(PYTORCH_ENGINE),
                (modelName, device) -> Model.newInstance(modelName, device, PYTORCH_ENGINE));
    }

    TorchModelLoader(EngineSupplier engineSupplier, ModelFactory modelFactory) {
        this.engineSupplier = engineSupplier;
        this.modelFactory = modelFactory;
    }

    LoadedModel load(String logicalModelName, Path modelPath) throws IOException, MalformedModelException {
        Engine engine = engineSupplier.get();
        int gpuCount = resolveGpuCount(logicalModelName, engine);
        if (gpuCount > 0) {
            try {
                LoadedModel loadedModel = loadOnDevice(logicalModelName, modelPath, Device.gpu(), true);
                LOGGER.info("Loaded {} model on device {}", logicalModelName, loadedModel.device());
                return loadedModel;
            } catch (RuntimeException | UnsatisfiedLinkError | IOException | MalformedModelException e) {
                LOGGER.warn(
                        "CUDA initialization failed for {}; falling back to CPU: {}: {}",
                        logicalModelName,
                        e.getClass().getSimpleName(),
                        e.getMessage());
            }
        }

        LoadedModel loadedModel = loadOnDevice(logicalModelName, modelPath, Device.cpu(), false);
        LOGGER.info("Loaded {} model on device {}", logicalModelName, loadedModel.device());
        return loadedModel;
    }

    private int resolveGpuCount(String logicalModelName, Engine engine) {
        try {
            return engine.getGpuCount();
        } catch (RuntimeException | UnsatisfiedLinkError e) {
            LOGGER.warn(
                    "CUDA probe failed for {}; falling back to CPU: {}: {}",
                    logicalModelName,
                    e.getClass().getSimpleName(),
                    e.getMessage());
            return 0;
        }
    }

    private LoadedModel loadOnDevice(String logicalModelName, Path modelPath, Device device, boolean cudaEnabled)
            throws IOException, MalformedModelException {
        Model model = modelFactory.newModel(logicalModelName, device);
        try {
            model.load(modelPath.getParent(), modelPath.getFileName().toString());
            return new LoadedModel(model, device, cudaEnabled);
        } catch (IOException | MalformedModelException | RuntimeException e) {
            model.close();
            throw e;
        }
    }

    @FunctionalInterface
    interface EngineSupplier {
        Engine get();
    }

    @FunctionalInterface
    interface ModelFactory {
        Model newModel(String logicalModelName, Device device);
    }

    record LoadedModel(Model model, Device device, boolean cudaEnabled) {
    }
}

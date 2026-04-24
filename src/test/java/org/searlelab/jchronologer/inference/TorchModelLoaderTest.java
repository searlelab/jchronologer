package org.searlelab.jchronologer.inference;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import ai.djl.Device;
import ai.djl.Model;
import ai.djl.engine.Engine;
import ai.djl.ndarray.NDManager;
import java.io.IOException;
import java.lang.reflect.Proxy;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import org.junit.jupiter.api.Test;

class TorchModelLoaderTest {

    @Test
    void loadSelectsCpuWhenNoGpuIsAvailable() throws Exception {
        List<Device> devices = new ArrayList<>();
        TorchModelLoader loader = new TorchModelLoader(
                () -> new TestEngine(0),
                (modelName, device) -> {
                    devices.add(device);
                    return newModelProxy(false, false);
                });

        TorchModelLoader.LoadedModel loadedModel = loader.load("Chronologer", tempModelPath());

        assertEquals(List.of(Device.cpu()), devices);
        assertEquals(Device.cpu(), loadedModel.device());
        assertFalse(loadedModel.cudaEnabled());
        loadedModel.model().close();
    }

    @Test
    void loadSelectsGpuWhenAvailable() throws Exception {
        List<Device> devices = new ArrayList<>();
        TorchModelLoader loader = new TorchModelLoader(
                () -> new TestEngine(1),
                (modelName, device) -> {
                    devices.add(device);
                    return newModelProxy(false, false);
                });

        TorchModelLoader.LoadedModel loadedModel = loader.load("Chronologer", tempModelPath());

        assertEquals(1, devices.size());
        assertTrue(devices.get(0).isGpu());
        assertTrue(loadedModel.device().isGpu());
        assertTrue(loadedModel.cudaEnabled());
        loadedModel.model().close();
    }

    @Test
    void loadFallsBackToCpuWhenGpuModelLoadFails() throws Exception {
        List<Device> devices = new ArrayList<>();
        TorchModelLoader loader = new TorchModelLoader(
                () -> new TestEngine(1),
                (modelName, device) -> {
                    devices.add(device);
                    return newModelProxy(device.isGpu(), false);
                });

        TorchModelLoader.LoadedModel loadedModel = loader.load("Chronologer", tempModelPath());

        assertEquals(2, devices.size());
        assertTrue(devices.get(0).isGpu());
        assertEquals(Device.cpu(), devices.get(1));
        assertEquals(Device.cpu(), loadedModel.device());
        assertFalse(loadedModel.cudaEnabled());
        loadedModel.model().close();
    }

    @Test
    void loadFallsBackToCpuWhenGpuProbeFails() throws Exception {
        List<Device> devices = new ArrayList<>();
        TorchModelLoader loader = new TorchModelLoader(
                () -> new TestEngine(new IllegalStateException("gpu probe failed")),
                (modelName, device) -> {
                    devices.add(device);
                    return newModelProxy(false, false);
                });

        TorchModelLoader.LoadedModel loadedModel = loader.load("Chronologer", tempModelPath());

        assertEquals(List.of(Device.cpu()), devices);
        assertEquals(Device.cpu(), loadedModel.device());
        assertFalse(loadedModel.cudaEnabled());
        loadedModel.model().close();
    }

    @Test
    void loadRethrowsCpuFailures() {
        TorchModelLoader loader = new TorchModelLoader(
                () -> new TestEngine(0),
                (modelName, device) -> newModelProxy(true, false));

        IllegalStateException error = assertThrows(
                IllegalStateException.class,
                () -> loader.load("Chronologer", tempModelPath()));

        assertEquals("simulated load failure", error.getMessage());
    }

    private static Path tempModelPath() throws IOException {
        return Files.createTempFile("torch-model-loader", ".pt");
    }

    private static Model newModelProxy(boolean failLoad, boolean failClose) {
        return (Model) Proxy.newProxyInstance(
                Model.class.getClassLoader(),
                new Class<?>[] {Model.class},
                (proxy, method, args) -> {
                    switch (method.getName()) {
                        case "load":
                            if (failLoad) {
                                throw new IllegalStateException("simulated load failure");
                            }
                            return null;
                        case "close":
                            if (failClose) {
                                throw new IllegalStateException("simulated close failure");
                            }
                            return null;
                        case "toString":
                            return "test-model";
                        case "hashCode":
                            return System.identityHashCode(proxy);
                        case "equals":
                            return proxy == args[0];
                        default:
                            throw new UnsupportedOperationException("Unexpected method: " + method.getName());
                    }
                });
    }

    private static final class TestEngine extends Engine {
        private final Integer gpuCount;
        private final RuntimeException gpuCountFailure;

        private TestEngine(int gpuCount) {
            this.gpuCount = gpuCount;
            this.gpuCountFailure = null;
        }

        private TestEngine(RuntimeException gpuCountFailure) {
            this.gpuCount = null;
            this.gpuCountFailure = gpuCountFailure;
        }

        @Override
        public int getGpuCount() {
            if (gpuCountFailure != null) {
                throw gpuCountFailure;
            }
            return gpuCount;
        }

        @Override
        public Engine getAlternativeEngine() {
            return null;
        }

        @Override
        public String getEngineName() {
            return "PyTorch";
        }

        @Override
        public int getRank() {
            return 0;
        }

        @Override
        public String getVersion() {
            return "test";
        }

        @Override
        public boolean hasCapability(String capability) {
            return false;
        }

        @Override
        public Model newModel(String name, Device device) {
            throw new UnsupportedOperationException();
        }

        @Override
        public NDManager newBaseManager() {
            throw new UnsupportedOperationException();
        }

        @Override
        public NDManager newBaseManager(Device device) {
            throw new UnsupportedOperationException();
        }
    }
}

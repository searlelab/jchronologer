package org.searlelab.jchronologer.inference;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import org.junit.jupiter.api.Test;
import org.searlelab.jchronologer.api.ChronologerLibraryOptions;

class ElectricianBatchPredictorTest {

    @Test
    void predictReturnsExpectedOutputShape() {
        try (ElectricianBatchPredictor predictor = new ElectricianBatchPredictor(
                ChronologerLibraryOptions.DEFAULT_ELECTRICIAN_MODEL_RESOURCE)) {
            long[][] tokens = new long[][] {new long[33], new long[33], new long[33]};

            float[][] output = predictor.predict(tokens);

            assertEquals(3, output.length);
            assertEquals(6, output[0].length);
            assertEquals(6, output[1].length);
            assertEquals(6, output[2].length);
        }
    }

    @Test
    void predictWithEmptyBatchReturnsEmptyOutput() {
        try (ElectricianBatchPredictor predictor = new ElectricianBatchPredictor(
                ChronologerLibraryOptions.DEFAULT_ELECTRICIAN_MODEL_RESOURCE)) {
            float[][] output = predictor.predict(new long[][] {});
            assertEquals(0, output.length);
        }
    }

    @Test
    void constructorRejectsMissingModelResource() {
        assertThrows(
                IllegalArgumentException.class,
                () -> new ElectricianBatchPredictor("models/not_a_real_electrician_model.pt"));
    }

    @Test
    void constructorWrapsMalformedModelLoadFailures() {
        IllegalStateException error = assertThrows(
                IllegalStateException.class,
                () -> new ElectricianBatchPredictor("models/invalid_cartographer_model.torchscript.pt"));
        assertTrue(error.getMessage().contains("Failed to load Electrician model"));
    }

    @Test
    void predictInitializesThreadLocalPredictorInWorkerThread() throws Exception {
        try (ElectricianBatchPredictor predictor = new ElectricianBatchPredictor(
                ChronologerLibraryOptions.DEFAULT_ELECTRICIAN_MODEL_RESOURCE)) {
            ExecutorService executor = Executors.newSingleThreadExecutor();
            try {
                Future<float[][]> outputFuture = executor.submit(() -> predictor.predict(new long[][] {new long[33]}));
                float[][] output = outputFuture.get();
                assertEquals(1, output.length);
                assertEquals(6, output[0].length);
            } catch (InterruptedException e) {
                Thread.currentThread().interrupt();
                throw e;
            } catch (ExecutionException e) {
                throw new AssertionError("Worker thread prediction failed", e.getCause());
            } finally {
                executor.shutdownNow();
            }
        }
    }
}

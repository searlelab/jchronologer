package org.searlelab.jchronologer.inference;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import org.junit.jupiter.api.Test;
import org.searlelab.jchronologer.api.ChronologerOptions;

class ScoutBatchPredictorTest {

    @Test
    void predictReturnsExpectedOutputShapes() {
        try (ScoutBatchPredictor predictor = new ScoutBatchPredictor(ChronologerOptions.DEFAULT_SCOUT_MODEL_RESOURCE)) {
            long[][] tokens = new long[][] {new long[33], new long[33]};
            float[][] charge = new float[][] {
                new float[] {0f, 1f, 0f, 0f, 0f, 0f},
                new float[] {0f, 0f, 1f, 0f, 0f, 0f}
            };
            float[][] nce = new float[][] {
                new float[] {0.3f},
                new float[] {0.35f}
            };

            ScoutBatchPredictor.ScoutBatchPredictionResult output = predictor.predict(tokens, charge, nce);

            assertEquals(2, output.getMs2().length);
            assertEquals(116, output.getMs2()[0].length);
            assertEquals(1, output.getIrt()[0].length);
            assertEquals(1, output.getCcs()[0].length);
            assertEquals(6, output.getChargeDist()[0].length);
        }
    }

    @Test
    void predictRejectsMismatchedBatchSizes() {
        try (ScoutBatchPredictor predictor = new ScoutBatchPredictor(ChronologerOptions.DEFAULT_SCOUT_MODEL_RESOURCE)) {
            IllegalArgumentException error = assertThrows(
                    IllegalArgumentException.class,
                    () -> predictor.predict(new long[][] {new long[33]}, new float[][] {}, new float[][] {new float[] {0.3f}}));
            assertTrue(error.getMessage().contains("matching batch size"));
        }
    }

    @Test
    void predictWithEmptyBatchReturnsEmptyOutput() {
        try (ScoutBatchPredictor predictor = new ScoutBatchPredictor(ChronologerOptions.DEFAULT_SCOUT_MODEL_RESOURCE)) {
            ScoutBatchPredictor.ScoutBatchPredictionResult output =
                    predictor.predict(new long[][] {}, new float[][] {}, new float[][] {});
            assertEquals(0, output.getMs2().length);
            assertEquals(0, output.getChargeDist().length);
        }
    }

    @Test
    void constructorRejectsMissingModelResource() {
        assertThrows(IllegalArgumentException.class, () -> new ScoutBatchPredictor("models/not_a_real_scout_model.pt"));
    }

    @Test
    void predictInitializesThreadLocalPredictorInWorkerThread() throws Exception {
        try (ScoutBatchPredictor predictor = new ScoutBatchPredictor(ChronologerOptions.DEFAULT_SCOUT_MODEL_RESOURCE)) {
            ExecutorService executor = Executors.newSingleThreadExecutor();
            try {
                Future<ScoutBatchPredictor.ScoutBatchPredictionResult> outputFuture = executor.submit(() -> predictor.predict(
                        new long[][] {new long[33]},
                        new float[][] {new float[] {0f, 1f, 0f, 0f, 0f, 0f}},
                        new float[][] {new float[] {0.3f}}));
                ScoutBatchPredictor.ScoutBatchPredictionResult output = outputFuture.get();
                assertEquals(1, output.getMs2().length);
                assertEquals(6, output.getChargeDist()[0].length);
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

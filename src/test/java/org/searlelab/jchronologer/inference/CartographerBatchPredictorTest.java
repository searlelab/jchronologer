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

class CartographerBatchPredictorTest {

    @Test
    void predictReturnsExpectedOutputShape() {
        try (CartographerBatchPredictor predictor = new CartographerBatchPredictor(
        		ChronologerLibraryOptions.DEFAULT_CARTOGRAPHER_MODEL_RESOURCE)) {
            long[][] tokens = new long[][] {new long[33], new long[33]};
            float[][] charge = new float[][] {
                new float[] {0f, 1f, 0f, 0f, 0f, 0f},
                new float[] {0f, 0f, 1f, 0f, 0f, 0f}
            };
            float[][] nce = new float[][] {
                new float[] {0.3f},
                new float[] {0.35f}
            };

            float[][] output = predictor.predict(tokens, charge, nce);

            assertEquals(2, output.length);
            assertEquals(CartographerSpectrumDecoder.VECTOR_LENGTH, output[0].length);
            assertEquals(CartographerSpectrumDecoder.VECTOR_LENGTH, output[1].length);
        }
    }

    @Test
    void predictTwoInputReturnsExpectedOutputShapeForSculptor() {
        try (CartographerBatchPredictor predictor = new CartographerBatchPredictor(
                ChronologerLibraryOptions.DEFAULT_SCULPTOR_MODEL_RESOURCE)) {
            long[][] tokens = new long[][] {new long[52], new long[52]};
            float[][] charge = new float[][] {
                new float[] {0f, 1f, 0f, 0f, 0f, 0f},
                new float[] {0f, 0f, 1f, 0f, 0f, 0f}
            };

            float[][] output = predictor.predict(tokens, charge);

            assertEquals(2, output.length);
            assertEquals(1, output[0].length);
            assertEquals(1, output[1].length);
        }
    }

    @Test
    void predictRejectsMismatchedBatchSizes() {
        try (CartographerBatchPredictor predictor = new CartographerBatchPredictor(
                ChronologerLibraryOptions.DEFAULT_CARTOGRAPHER_MODEL_RESOURCE)) {
            IllegalArgumentException error = assertThrows(
                    IllegalArgumentException.class,
                    () -> predictor.predict(
                            new long[][] {new long[33]},
                            new float[][] {},
                            new float[][] {new float[] {0.3f}}));
            assertTrue(error.getMessage().contains("matching batch size"));
        }
    }

    @Test
    void predictWithEmptyBatchReturnsEmptyOutput() {
        try (CartographerBatchPredictor predictor = new CartographerBatchPredictor(
                ChronologerLibraryOptions.DEFAULT_CARTOGRAPHER_MODEL_RESOURCE)) {
            float[][] output = predictor.predict(new long[][] {}, new float[][] {}, new float[][] {});
            assertEquals(0, output.length);
        }
    }

    @Test
    void constructorRejectsMissingModelResource() {
        assertThrows(
                IllegalArgumentException.class,
                () -> new CartographerBatchPredictor("models/not_a_real_cartographer_model.pt"));
    }

    @Test
    void constructorWrapsMalformedModelLoadFailures() {
        IllegalStateException error = assertThrows(
                IllegalStateException.class,
                () -> new CartographerBatchPredictor("models/invalid_cartographer_model.torchscript.pt"));
        assertTrue(error.getMessage().contains("Failed to load Cartographer model"));
    }

    @Test
    void predictWrapsTranslateFailures() {
        try (CartographerBatchPredictor predictor = new CartographerBatchPredictor(
                ChronologerLibraryOptions.DEFAULT_CARTOGRAPHER_MODEL_RESOURCE)) {
            IllegalStateException error = assertThrows(
                    IllegalStateException.class,
                    () -> predictor.predict(
                            new long[][] {new long[33]},
                            new float[][] {new float[] {1f, 0f, 0f, 0f, 0f}},
                            new float[][] {new float[] {0.3f}}));
            assertTrue(error.getMessage().contains("Failed to run Cartographer inference"));
        }
    }

    @Test
    void predictInitializesThreadLocalPredictorInWorkerThread() throws Exception {
        try (CartographerBatchPredictor predictor = new CartographerBatchPredictor(
                ChronologerLibraryOptions.DEFAULT_CARTOGRAPHER_MODEL_RESOURCE)) {
            ExecutorService executor = Executors.newSingleThreadExecutor();
            try {
                Future<float[][]> outputFuture = executor.submit(() -> predictor.predict(
                        new long[][] {new long[33]},
                        new float[][] {new float[] {0f, 1f, 0f, 0f, 0f, 0f}},
                        new float[][] {new float[] {0.3f}}));
                float[][] output = outputFuture.get();
                assertEquals(1, output.length);
                assertEquals(CartographerSpectrumDecoder.VECTOR_LENGTH, output[0].length);
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

package org.searlelab.jchronologer.inference;

import static org.junit.jupiter.api.Assertions.assertEquals;

import org.junit.jupiter.api.Test;

class CartographerBatchPredictorTest {

    @Test
    void predictReturnsExpectedOutputShape() {
        try (CartographerBatchPredictor predictor = new CartographerBatchPredictor(
                "models/Cartographer_20260214222413.torchscript.pt")) {
            long[][] tokens = new long[][] {new long[52], new long[52]};
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
}

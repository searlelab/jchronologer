package org.searlelab.jchronologer.inference;

import static org.junit.jupiter.api.Assertions.assertEquals;

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
}

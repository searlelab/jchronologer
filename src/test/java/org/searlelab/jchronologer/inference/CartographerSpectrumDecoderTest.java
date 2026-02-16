package org.searlelab.jchronologer.inference;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;
import org.searlelab.jchronologer.preprocessing.PeptideSequenceConverter;
import org.searlelab.jchronologer.preprocessing.PeptideSequenceConverter.ParsedUnimodSequence;

class CartographerSpectrumDecoderTest {

    @Test
    void decodeFiltersNegativeLowAndImpossibleIons() {
        ParsedUnimodSequence peptide = PeptideSequenceConverter.parseNormalizedUnimod("[]-PEPTIDE-[]");

        float[] vector = new float[CartographerSpectrumDecoder.VECTOR_LENGTH];
        for (int i = 0; i < vector.length; i++) {
            vector[i] = -1.0f;
        }

        vector[0] = 0.5f;    // 1+y1 keep
        vector[1] = 0.009f;  // 2+y1 drop (<0.01)
        vector[3] = 0.3f;    // 1+b1 keep
        vector[5] = 0.4f;    // 3+b1 drop (precursor charge 2)
        vector[(20 - 1) * 6] = 0.8f; // 1+y20 drop (impossible for 7-mer)

        CartographerSpectrumDecoder.DecodedSpectrum decoded = CartographerSpectrumDecoder.decode(
                peptide,
                (byte) 2,
                vector,
                0.01f);

        assertEquals(2, decoded.getIonTypeArray().length);
        assertEquals("1+y1", decoded.getIonTypeArray()[0]);
        assertEquals("1+b1", decoded.getIonTypeArray()[1]);

        assertEquals(2, decoded.getMassArray().length);
        assertEquals(2, decoded.getIntensityArray().length);
        assertTrue(decoded.getMassArray()[0] > 0.0);
        assertTrue(decoded.getMassArray()[1] > 0.0);
        assertTrue(decoded.getIntensityArray()[0] >= 0.01f);
        assertTrue(decoded.getIntensityArray()[1] >= 0.01f);
    }
}

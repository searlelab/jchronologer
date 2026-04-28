package org.searlelab.jchronologer.inference;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;
import org.searlelab.jchronologer.preprocessing.PeptideSequenceConverter;

class ScoutSpectrumDecoderTest {

    @Test
    void rejectsUnexpectedVectorLength() {
        assertThrows(IllegalArgumentException.class, () -> ScoutSpectrumDecoder.decode(
                PeptideSequenceConverter.parseNormalizedUnimod("[]-PEPTIDEK-[]"),
                (byte) 2,
                new float[10],
                0.01f));
    }

    @Test
    void decodesOnlyOneAndTwoChargeFragments() {
        float[] vector = new float[ScoutSpectrumDecoder.VECTOR_LENGTH];
        vector[0] = 0.5f;
        vector[1] = 0.4f;
        vector[2] = 0.3f;
        vector[3] = 0.2f;

        var decoded = ScoutSpectrumDecoder.decode(
                PeptideSequenceConverter.parseNormalizedUnimod("[]-PEPTIDEK-[]"),
                (byte) 2,
                vector,
                0.01f);

        assertEquals(4, decoded.getIonTypeArray().length);
        assertArrayEquals(new String[] {"1+y1", "2+y1", "1+b1", "2+b1"}, decoded.getIonTypeArray());
        for (String ionType : decoded.getIonTypeArray()) {
            assertTrue(ionType.startsWith("1+") || ionType.startsWith("2+"));
        }
    }
}

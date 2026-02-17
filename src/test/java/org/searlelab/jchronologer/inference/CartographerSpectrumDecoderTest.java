package org.searlelab.jchronologer.inference;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.Arrays;
import java.util.LinkedHashSet;
import java.util.Set;
import org.junit.jupiter.api.Test;
import org.searlelab.jchronologer.preprocessing.PeptideSequenceConverter;
import org.searlelab.jchronologer.preprocessing.PeptideSequenceConverter.ParsedUnimodSequence;

class CartographerSpectrumDecoderTest {

    /**
     * Prosit ion-interleaved index: (ionNumber - 1) * 6 + channel.
     *
     * <pre>
     *   channel 0 = y+1, channel 1 = y+2, channel 2 = y+3,
     *   channel 3 = b+1, channel 4 = b+2, channel 5 = b+3
     * </pre>
     */
    private static int idx(int channel, int ionNumber) {
        return (ionNumber - 1) * 6 + channel;
    }

    private static float[] emptyVector() {
        float[] v = new float[CartographerSpectrumDecoder.VECTOR_LENGTH];
        Arrays.fill(v, -1.0f);
        return v;
    }

    // ── Original smoke test ──────────────────────────────────────────────

    @Test
    void decodeFiltersNegativeLowAndImpossibleIons() {
        ParsedUnimodSequence peptide = PeptideSequenceConverter.parseNormalizedUnimod("[]-PEPTIDE-[]");

        float[] vector = emptyVector();

        // Prosit ion-interleaved layout (stride 6):
        //   index = (ionNumber - 1) * 6 + channelOffset
        vector[idx(0, 1)] = 0.5f;    // y1+1  keep
        vector[idx(1, 1)] = 0.009f;  // y1+2  drop (<0.01)
        vector[idx(3, 1)] = 0.3f;    // b1+1  keep
        vector[idx(5, 1)] = 0.4f;    // b1+3  drop (precursor charge 2)
        vector[idx(0, 20)] = 0.8f;   // y20+1 drop (impossible for 7-mer)

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

    // ── Prosit index-to-ion-type mapping ────────────────────────────────

    @Test
    void allSixChannelsMapToCorrectIonTypes() {
        // 20-mer, charge 3 → all 6 channels can produce ions
        ParsedUnimodSequence peptide = PeptideSequenceConverter.parseNormalizedUnimod(
                "[]-ACDEFGHIKLMNPQRSTVWY-[]");
        float[] vector = emptyVector();

        // Place one ion in each channel at ion number 5
        vector[idx(0, 5)] = 1.0f;  // y5+1
        vector[idx(1, 5)] = 1.0f;  // y5+2
        vector[idx(2, 5)] = 1.0f;  // y5+3
        vector[idx(3, 5)] = 1.0f;  // b5+1
        vector[idx(4, 5)] = 1.0f;  // b5+2
        vector[idx(5, 5)] = 1.0f;  // b5+3

        CartographerSpectrumDecoder.DecodedSpectrum decoded = CartographerSpectrumDecoder.decode(
                peptide, (byte) 3, vector, 0.0f);

        Set<String> ions = new LinkedHashSet<>(Arrays.asList(decoded.getIonTypeArray()));
        assertEquals(6, ions.size());
        assertTrue(ions.contains("1+y5"));
        assertTrue(ions.contains("2+y5"));
        assertTrue(ions.contains("3+y5"));
        assertTrue(ions.contains("1+b5"));
        assertTrue(ions.contains("2+b5"));
        assertTrue(ions.contains("3+b5"));
    }

    // ── Ion-number boundary masking ─────────────────────────────────────

    @Test
    void ionNumberBoundaryLastValidIonIncluded() {
        // 10-mer → max valid ion number is 9
        ParsedUnimodSequence peptide = PeptideSequenceConverter.parseNormalizedUnimod(
                "[]-ACDEFGHIKL-[]");

        float[] vector = emptyVector();
        vector[idx(0, 9)] = 1.0f;  // y9+1 → last valid

        CartographerSpectrumDecoder.DecodedSpectrum decoded = CartographerSpectrumDecoder.decode(
                peptide, (byte) 2, vector, 0.0f);

        assertEquals(1, decoded.getIonTypeArray().length);
        assertEquals("1+y9", decoded.getIonTypeArray()[0]);
    }

    @Test
    void ionNumberBoundaryAtPeptideLengthMasked() {
        // 10-mer → ion number 10 is impossible (needs peptideLength > ionNumber)
        ParsedUnimodSequence peptide = PeptideSequenceConverter.parseNormalizedUnimod(
                "[]-ACDEFGHIKL-[]");

        float[] vector = emptyVector();
        vector[idx(0, 10)] = 1.0f;  // y10+1 → impossible
        vector[idx(3, 10)] = 1.0f;  // b10+1 → impossible

        CartographerSpectrumDecoder.DecodedSpectrum decoded = CartographerSpectrumDecoder.decode(
                peptide, (byte) 2, vector, 0.0f);

        assertEquals(0, decoded.getIonTypeArray().length);
    }

    @Test
    void ionNumberBoundaryBeyondPeptideLengthMasked() {
        // 10-mer → ion numbers 15, 20, 29 all impossible
        ParsedUnimodSequence peptide = PeptideSequenceConverter.parseNormalizedUnimod(
                "[]-ACDEFGHIKL-[]");

        float[] vector = emptyVector();
        vector[idx(0, 15)] = 1.0f;  // y15+1
        vector[idx(0, 20)] = 1.0f;  // y20+1
        vector[idx(0, 29)] = 1.0f;  // y29+1
        vector[idx(3, 15)] = 1.0f;  // b15+1
        vector[idx(3, 20)] = 1.0f;  // b20+1
        vector[idx(3, 29)] = 1.0f;  // b29+1

        CartographerSpectrumDecoder.DecodedSpectrum decoded = CartographerSpectrumDecoder.decode(
                peptide, (byte) 3, vector, 0.0f);

        assertEquals(0, decoded.getIonTypeArray().length);
    }

    @Test
    void ionNumberBoundaryBothSidesOfCutoff() {
        // 15-mer → ion 14 valid, ion 15 masked
        ParsedUnimodSequence peptide = PeptideSequenceConverter.parseNormalizedUnimod(
                "[]-ACDEFGHIKLMNPQR-[]");

        float[] vector = emptyVector();
        // y ions
        vector[idx(0, 14)] = 0.9f;  // y14+1 → valid (last)
        vector[idx(0, 15)] = 0.9f;  // y15+1 → masked
        // b ions
        vector[idx(3, 14)] = 0.8f;  // b14+1 → valid (last)
        vector[idx(3, 15)] = 0.8f;  // b15+1 → masked

        CartographerSpectrumDecoder.DecodedSpectrum decoded = CartographerSpectrumDecoder.decode(
                peptide, (byte) 2, vector, 0.0f);

        assertEquals(2, decoded.getIonTypeArray().length);
        Set<String> ions = new LinkedHashSet<>(Arrays.asList(decoded.getIonTypeArray()));
        assertTrue(ions.contains("1+y14"));
        assertTrue(ions.contains("1+b14"));
    }

    // ── Charge masking ──────────────────────────────────────────────────

    @Test
    void chargeOnePrecursorOnlyCharge1FragmentsSurvive() {
        ParsedUnimodSequence peptide = PeptideSequenceConverter.parseNormalizedUnimod(
                "[]-ACDEFGHIKL-[]");

        float[] vector = emptyVector();
        // Set ion 3 in all 6 channels
        for (int ch = 0; ch < 6; ch++) {
            vector[idx(ch, 3)] = 1.0f;
        }

        CartographerSpectrumDecoder.DecodedSpectrum decoded = CartographerSpectrumDecoder.decode(
                peptide, (byte) 1, vector, 0.0f);

        // Only charge 1 fragments pass: y3+1 and b3+1
        assertEquals(2, decoded.getIonTypeArray().length);
        Set<String> ions = new LinkedHashSet<>(Arrays.asList(decoded.getIonTypeArray()));
        assertTrue(ions.contains("1+y3"));
        assertTrue(ions.contains("1+b3"));
    }

    @Test
    void chargeTwoPrecursorCharge1And2FragmentsSurvive() {
        ParsedUnimodSequence peptide = PeptideSequenceConverter.parseNormalizedUnimod(
                "[]-ACDEFGHIKL-[]");

        float[] vector = emptyVector();
        for (int ch = 0; ch < 6; ch++) {
            vector[idx(ch, 3)] = 1.0f;
        }

        CartographerSpectrumDecoder.DecodedSpectrum decoded = CartographerSpectrumDecoder.decode(
                peptide, (byte) 2, vector, 0.0f);

        // Charge 1 and 2 pass: y3+1, b3+1, y3+2, b3+2
        assertEquals(4, decoded.getIonTypeArray().length);
        Set<String> ions = new LinkedHashSet<>(Arrays.asList(decoded.getIonTypeArray()));
        assertTrue(ions.contains("1+y3"));
        assertTrue(ions.contains("1+b3"));
        assertTrue(ions.contains("2+y3"));
        assertTrue(ions.contains("2+b3"));
    }

    @Test
    void chargeThreePrecursorAllFragmentChargesSurvive() {
        ParsedUnimodSequence peptide = PeptideSequenceConverter.parseNormalizedUnimod(
                "[]-ACDEFGHIKL-[]");

        float[] vector = emptyVector();
        for (int ch = 0; ch < 6; ch++) {
            vector[idx(ch, 3)] = 1.0f;
        }

        CartographerSpectrumDecoder.DecodedSpectrum decoded = CartographerSpectrumDecoder.decode(
                peptide, (byte) 3, vector, 0.0f);

        // All 6 pass
        assertEquals(6, decoded.getIonTypeArray().length);
        Set<String> ions = new LinkedHashSet<>(Arrays.asList(decoded.getIonTypeArray()));
        assertTrue(ions.contains("1+y3"));
        assertTrue(ions.contains("2+y3"));
        assertTrue(ions.contains("3+y3"));
        assertTrue(ions.contains("1+b3"));
        assertTrue(ions.contains("2+b3"));
        assertTrue(ions.contains("3+b3"));
    }

    // ── Full vector: total ion count ────────────────────────────────────

    @Test
    void fullVectorTotalIonCountMatchesExpected() {
        // 10-mer, charge 2 → valid ions 1..9, charges 1..2, y+b
        // Expected: 9 * 2 * 2 = 36
        ParsedUnimodSequence peptide = PeptideSequenceConverter.parseNormalizedUnimod(
                "[]-ACDEFGHIKL-[]");

        float[] vector = new float[CartographerSpectrumDecoder.VECTOR_LENGTH];
        Arrays.fill(vector, 1.0f);

        CartographerSpectrumDecoder.DecodedSpectrum decoded = CartographerSpectrumDecoder.decode(
                peptide, (byte) 2, vector, 0.0f);

        assertEquals(36, decoded.getIonTypeArray().length);
        assertEquals(36, decoded.getMassArray().length);
        assertEquals(36, decoded.getIntensityArray().length);
    }

    @Test
    void fullVectorCharge3TotalIonCount() {
        // 10-mer, charge 3 → 9 * 3 * 2 = 54
        ParsedUnimodSequence peptide = PeptideSequenceConverter.parseNormalizedUnimod(
                "[]-ACDEFGHIKL-[]");

        float[] vector = new float[CartographerSpectrumDecoder.VECTOR_LENGTH];
        Arrays.fill(vector, 1.0f);

        CartographerSpectrumDecoder.DecodedSpectrum decoded = CartographerSpectrumDecoder.decode(
                peptide, (byte) 3, vector, 0.0f);

        assertEquals(54, decoded.getIonTypeArray().length);
    }

    @Test
    void fullVectorCharge1TotalIonCount() {
        // 10-mer, charge 1 → 9 * 1 * 2 = 18
        ParsedUnimodSequence peptide = PeptideSequenceConverter.parseNormalizedUnimod(
                "[]-ACDEFGHIKL-[]");

        float[] vector = new float[CartographerSpectrumDecoder.VECTOR_LENGTH];
        Arrays.fill(vector, 1.0f);

        CartographerSpectrumDecoder.DecodedSpectrum decoded = CartographerSpectrumDecoder.decode(
                peptide, (byte) 1, vector, 0.0f);

        assertEquals(18, decoded.getIonTypeArray().length);
    }

    @Test
    void fullVectorShortPeptide() {
        // 6-mer, charge 2 → 5 * 2 * 2 = 20
        ParsedUnimodSequence peptide = PeptideSequenceConverter.parseNormalizedUnimod(
                "[]-PEPTIK-[]");

        float[] vector = new float[CartographerSpectrumDecoder.VECTOR_LENGTH];
        Arrays.fill(vector, 1.0f);

        CartographerSpectrumDecoder.DecodedSpectrum decoded = CartographerSpectrumDecoder.decode(
                peptide, (byte) 2, vector, 0.0f);

        assertEquals(20, decoded.getIonTypeArray().length);
    }

    @Test
    void fullVectorLongPeptide() {
        // 25-mer, charge 3 → 24 * 3 * 2 = 144
        ParsedUnimodSequence peptide = PeptideSequenceConverter.parseNormalizedUnimod(
                "[]-ACDEFGHIKLMNPQRSTVWYACDEF-[]");

        float[] vector = new float[CartographerSpectrumDecoder.VECTOR_LENGTH];
        Arrays.fill(vector, 1.0f);

        CartographerSpectrumDecoder.DecodedSpectrum decoded = CartographerSpectrumDecoder.decode(
                peptide, (byte) 3, vector, 0.0f);

        assertEquals(144, decoded.getIonTypeArray().length);
    }

    // ── Intensity threshold edge cases ──────────────────────────────────

    @Test
    void intensityExactlyAtThresholdIncluded() {
        ParsedUnimodSequence peptide = PeptideSequenceConverter.parseNormalizedUnimod(
                "[]-ACDEFGHIKL-[]");

        float[] vector = emptyVector();
        vector[idx(0, 1)] = 0.01f;  // exactly at threshold

        CartographerSpectrumDecoder.DecodedSpectrum decoded = CartographerSpectrumDecoder.decode(
                peptide, (byte) 2, vector, 0.01f);

        assertEquals(1, decoded.getIonTypeArray().length);
        assertEquals("1+y1", decoded.getIonTypeArray()[0]);
    }

    @Test
    void intensityJustBelowThresholdExcluded() {
        ParsedUnimodSequence peptide = PeptideSequenceConverter.parseNormalizedUnimod(
                "[]-ACDEFGHIKL-[]");

        float[] vector = emptyVector();
        vector[idx(0, 1)] = 0.0099f;  // just below 0.01

        CartographerSpectrumDecoder.DecodedSpectrum decoded = CartographerSpectrumDecoder.decode(
                peptide, (byte) 2, vector, 0.01f);

        assertEquals(0, decoded.getIonTypeArray().length);
    }

    @Test
    void zeroIntensityWithZeroThresholdIncluded() {
        ParsedUnimodSequence peptide = PeptideSequenceConverter.parseNormalizedUnimod(
                "[]-ACDEFGHIKL-[]");

        float[] vector = emptyVector();
        vector[idx(0, 1)] = 0.0f;

        CartographerSpectrumDecoder.DecodedSpectrum decoded = CartographerSpectrumDecoder.decode(
                peptide, (byte) 2, vector, 0.0f);

        assertEquals(1, decoded.getIonTypeArray().length);
        assertEquals("1+y1", decoded.getIonTypeArray()[0]);
    }

    @Test
    void nanAndInfinityExcluded() {
        ParsedUnimodSequence peptide = PeptideSequenceConverter.parseNormalizedUnimod(
                "[]-ACDEFGHIKL-[]");

        float[] vector = emptyVector();
        vector[idx(0, 1)] = Float.NaN;
        vector[idx(0, 2)] = Float.POSITIVE_INFINITY;
        vector[idx(0, 3)] = Float.NEGATIVE_INFINITY;
        vector[idx(3, 1)] = 0.5f;  // only valid ion

        CartographerSpectrumDecoder.DecodedSpectrum decoded = CartographerSpectrumDecoder.decode(
                peptide, (byte) 2, vector, 0.0f);

        assertEquals(1, decoded.getIonTypeArray().length);
        assertEquals("1+b1", decoded.getIonTypeArray()[0]);
    }

    // ── Combined masking: charge + ion number + intensity ───────────────

    @Test
    void combinedMaskingAllConditionsExercised() {
        // 8-mer PEPTIDER, charge 2 → valid ions 1..7, charges 1..2
        ParsedUnimodSequence peptide = PeptideSequenceConverter.parseNormalizedUnimod(
                "[]-PEPTIDER-[]");

        float[] vector = emptyVector();

        // Valid: y3+1 (ion 3, charge 1, 8-mer, good intensity)
        vector[idx(0, 3)] = 0.7f;

        // Valid: b5+2 (ion 5, charge 2, 8-mer, good intensity)
        vector[idx(4, 5)] = 0.6f;

        // Masked by charge: y3+3 (charge 3 > precursor charge 2)
        vector[idx(2, 3)] = 0.9f;

        // Masked by ion number: b8+1 (ion 8 >= peptide length 8)
        vector[idx(3, 8)] = 0.9f;

        // Masked by ion number: y10+2 (ion 10 >= peptide length 8)
        vector[idx(1, 10)] = 0.9f;

        // Masked by intensity threshold
        vector[idx(0, 1)] = 0.005f;

        // Masked by negative intensity
        vector[idx(3, 1)] = -0.1f;

        CartographerSpectrumDecoder.DecodedSpectrum decoded = CartographerSpectrumDecoder.decode(
                peptide, (byte) 2, vector, 0.01f);

        assertEquals(2, decoded.getIonTypeArray().length);
        Set<String> ions = new LinkedHashSet<>(Arrays.asList(decoded.getIonTypeArray()));
        assertTrue(ions.contains("1+y3"));
        assertTrue(ions.contains("2+b5"));
    }

    // ── Verify every valid ion appears exactly once ─────────────────────

    @Test
    void noIonTypesDuplicated() {
        // 10-mer, charge 3 → 54 ions; fill all with 1.0
        ParsedUnimodSequence peptide = PeptideSequenceConverter.parseNormalizedUnimod(
                "[]-ACDEFGHIKL-[]");

        float[] vector = new float[CartographerSpectrumDecoder.VECTOR_LENGTH];
        Arrays.fill(vector, 1.0f);

        CartographerSpectrumDecoder.DecodedSpectrum decoded = CartographerSpectrumDecoder.decode(
                peptide, (byte) 3, vector, 0.0f);

        String[] ionTypes = decoded.getIonTypeArray();
        Set<String> unique = new LinkedHashSet<>(Arrays.asList(ionTypes));
        assertEquals(ionTypes.length, unique.size(), "Duplicate ion types in decoded spectrum");
    }

    @Test
    void everyExpectedIonPresentInFullVector() {
        // 10-mer, charge 2 → verify every expected ion string is present
        ParsedUnimodSequence peptide = PeptideSequenceConverter.parseNormalizedUnimod(
                "[]-ACDEFGHIKL-[]");

        float[] vector = new float[CartographerSpectrumDecoder.VECTOR_LENGTH];
        Arrays.fill(vector, 1.0f);

        CartographerSpectrumDecoder.DecodedSpectrum decoded = CartographerSpectrumDecoder.decode(
                peptide, (byte) 2, vector, 0.0f);

        Set<String> ions = new LinkedHashSet<>(Arrays.asList(decoded.getIonTypeArray()));
        for (int ionNum = 1; ionNum <= 9; ionNum++) {
            for (int charge = 1; charge <= 2; charge++) {
                String yLabel = charge + "+y" + ionNum;
                String bLabel = charge + "+b" + ionNum;
                assertTrue(ions.contains(yLabel), "Missing expected ion: " + yLabel);
                assertTrue(ions.contains(bLabel), "Missing expected ion: " + bLabel);
            }
        }
        // Verify no ion number >= peptideLength leaked through
        for (int ionNum = 10; ionNum <= 29; ionNum++) {
            for (int charge = 1; charge <= 3; charge++) {
                String yLabel = charge + "+y" + ionNum;
                String bLabel = charge + "+b" + ionNum;
                assertTrue(!ions.contains(yLabel), "Impossible ion leaked: " + yLabel);
                assertTrue(!ions.contains(bLabel), "Impossible ion leaked: " + bLabel);
            }
        }
    }
}

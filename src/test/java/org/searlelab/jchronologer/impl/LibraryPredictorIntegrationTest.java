package org.searlelab.jchronologer.impl;

import static org.junit.jupiter.api.Assertions.assertDoesNotThrow;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.junit.jupiter.api.Test;
import org.searlelab.jchronologer.api.ChronologerLibraryEntry;
import org.searlelab.jchronologer.api.ChronologerLibraryPredictor;
import org.searlelab.jchronologer.api.LibraryPredictionRequest;
import org.searlelab.jchronologer.api.PrecursorCondition;
import org.searlelab.jchronologer.preprocessing.PeptideSequenceConverter;

class LibraryPredictorIntegrationTest {
    private static final Pattern ION_TYPE_PATTERN = Pattern.compile("([1-3])\\+([yb])(\\d+)");

    @Test
    void predictTASEFDSAIAQDK() {
        LibraryPredictionRequest request = new LibraryPredictionRequest(
                "TASEFDSAIAQDK",
                List.of(new PrecursorCondition((byte) 2, 27.0)));

        try (ChronologerLibraryPredictor predictor = ChronologerFactory.createLibraryPredictorDefault()) {
            List<ChronologerLibraryEntry> entries = predictor.predict(List.of(request));
            assertEquals(1, entries.size());

            ChronologerLibraryEntry entry = entries.get(0);
            assertEquals("[]-TASEFDSAIAQDK-[]", entry.getUnimodPeptideSequence());
            assertEquals("TASEFDSAIAQDK", entry.getPeptideModSeq());
            assertEquals((byte) 2, entry.getPrecursorCharge());
            assertTrue(entry.getPrecursorMz() > 0.0);
            assertTrue(entry.getRetentionTimeInSeconds() > 0.0f);
            assertEntryHasValidSpectrum(entry, 5);
            assertTrue(countAbove(entry.getIntensityArray(), 0.10f) >= 2);
        }
    }

    @Test
    void massAndUnimodModifiedPeptidesProduceEquivalentEntries() {
        LibraryPredictionRequest massEncoded = new LibraryPredictionRequest(
                "K[458.325864]PGLAITFAK[229.162932]",
                List.of(new PrecursorCondition((byte) 2, 30.0)));
        LibraryPredictionRequest unimodEncoded = new LibraryPredictionRequest(
                "[UNIMOD:737]-K[UNIMOD:737]PGLAITFAK[UNIMOD:737]-[]",
                List.of(new PrecursorCondition((byte) 2, 30.0)));

        try (ChronologerLibraryPredictor predictor = ChronologerFactory.createLibraryPredictorDefault()) {
            List<ChronologerLibraryEntry> entries = predictor.predict(List.of(massEncoded, unimodEncoded));
            assertEquals(2, entries.size());

            ChronologerLibraryEntry first = entries.get(0);
            ChronologerLibraryEntry second = entries.get(1);
            String canonicalUnimod = "[UNIMOD:737]-K[UNIMOD:737]PGLAITFAK[UNIMOD:737]-[]";
            String canonicalMass = "[+229.162932]K[+229.162932]PGLAITFAK[+229.162932]";

            assertEquals(canonicalUnimod, first.getUnimodPeptideSequence());
            assertEquals(canonicalUnimod, second.getUnimodPeptideSequence());
            assertEquals(canonicalMass, first.getPeptideModSeq());
            assertEquals(canonicalMass, second.getPeptideModSeq());
            assertEquals(first.getPrecursorCharge(), second.getPrecursorCharge());
            assertEquals(first.getPrecursorNce(), second.getPrecursorNce());
            assertEquals(first.getPrecursorMz(), second.getPrecursorMz(), 1e-9);
            assertEquals(first.getRetentionTimeInSeconds(), second.getRetentionTimeInSeconds(), 1e-5f);

            assertEntryHasValidSpectrum(first, 5);
            assertEntryHasValidSpectrum(second, 5);

            List<String> firstTopIons = topIonTypes(first, 8);
            List<String> secondTopIons = topIonTypes(second, 8);
            assertTrue(overlapCount(firstTopIons, secondTopIons) >= 6,
                    "Expected strong top-ion overlap between equivalent modified-peptide encodings.");
        }
    }

    @Test
    void predictsMultipleChargeNcePairsAndFiltersMissingIons() {
        LibraryPredictionRequest tmtRequest = new LibraryPredictionRequest(
                "K[458.325864]PGLAITFAK[229.162932]",
                List.of(
                        new PrecursorCondition((byte) 2, 30.0),
                        new PrecursorCondition((byte) 3, 35.0)));
        LibraryPredictionRequest phosphoOxRequest = new LibraryPredictionRequest(
                "ACDM[+15.994915]STY[+79.966331]K[+42.010565]",
                List.of(new PrecursorCondition((byte) 2, 27.0)));

        try (ChronologerLibraryPredictor predictor = ChronologerFactory.createLibraryPredictorDefault()) {
            List<ChronologerLibraryEntry> entries = predictor.predict(List.of(tmtRequest, phosphoOxRequest));
            assertEquals(3, entries.size());

            ChronologerLibraryEntry tmtCharge2 = entries.get(0);
            ChronologerLibraryEntry tmtCharge3 = entries.get(1);
            ChronologerLibraryEntry phosphoOx = entries.get(2);

            assertEquals("[UNIMOD:737]-K[UNIMOD:737]PGLAITFAK[UNIMOD:737]-[]", tmtCharge2.getUnimodPeptideSequence());
            assertEquals("[UNIMOD:737]-K[UNIMOD:737]PGLAITFAK[UNIMOD:737]-[]", tmtCharge3.getUnimodPeptideSequence());
            assertEquals((byte) 2, tmtCharge2.getPrecursorCharge());
            assertEquals((byte) 3, tmtCharge3.getPrecursorCharge());
            assertEquals(tmtCharge2.getRetentionTimeInSeconds(), tmtCharge3.getRetentionTimeInSeconds(), 1e-5f);

            assertTrue(phosphoOx.getUnimodPeptideSequence().contains("UNIMOD:35"));
            assertTrue(phosphoOx.getUnimodPeptideSequence().contains("UNIMOD:21"));
            assertTrue(phosphoOx.getUnimodPeptideSequence().contains("UNIMOD:1"));

            for (ChronologerLibraryEntry entry : entries) {
                assertEntryHasValidSpectrum(entry, 5);
                assertTrue(countAbove(entry.getIntensityArray(), 0.10f) >= 2,
                        "Expected >=2 high-intensity ions for " + entry.getUnimodPeptideSequence());
            }
        }
    }

    @Test
    void predictorValidatesRequestsAndLifecycle() {
        try (ChronologerLibraryPredictor predictor = ChronologerFactory.createLibraryPredictorDefault()) {
            assertEquals(0, predictor.predict(null).size());
            assertEquals(0, predictor.predict(List.of()).size());

            List<LibraryPredictionRequest> nullRequestList = new ArrayList<>();
            nullRequestList.add(null);
            assertThrows(IllegalArgumentException.class, () -> predictor.predict(nullRequestList));
            assertThrows(IllegalArgumentException.class, () -> predictor.predict(List.of(
                    new LibraryPredictionRequest("PEPTIDEK", List.of()))));
            assertThrows(IllegalArgumentException.class, () -> predictor.predict(List.of(
                    new LibraryPredictionRequest("PEPTIDEK", Arrays.asList((PrecursorCondition) null)))));
            assertThrows(IllegalArgumentException.class, () -> predictor.predict(List.of(
                    new LibraryPredictionRequest("PEPTIDEK", List.of(new PrecursorCondition((byte) 0, 30.0))))));
            assertThrows(IllegalArgumentException.class, () -> predictor.predict(List.of(
                    new LibraryPredictionRequest("PEPTIDEK", List.of(new PrecursorCondition((byte) 7, 30.0))))));
            assertThrows(IllegalArgumentException.class, () -> predictor.predict(List.of(
                    new LibraryPredictionRequest("PEPTIDEK", List.of(new PrecursorCondition((byte) 2, Double.NaN))))));
            assertThrows(IllegalArgumentException.class, () -> predictor.predict(List.of(
                    new LibraryPredictionRequest("PEPTIDEK", List.of(new PrecursorCondition((byte) 2, 9.99))))));
            assertThrows(IllegalArgumentException.class, () -> predictor.predict(List.of(
                    new LibraryPredictionRequest("PEPTIDEK", List.of(new PrecursorCondition((byte) 2, 61.0))))));
        }

        ChronologerLibraryPredictor closedPredictor = ChronologerFactory.createLibraryPredictorDefault();
        closedPredictor.close();
        assertDoesNotThrow(closedPredictor::close);
        assertThrows(IllegalStateException.class, closedPredictor::init);
        assertThrows(IllegalStateException.class, () -> closedPredictor.predict(List.of(
                new LibraryPredictionRequest("PEPTIDEK", List.of(new PrecursorCondition((byte) 2, 30.0))))));
    }

    private static void assertEntryHasValidSpectrum(ChronologerLibraryEntry entry, int minimumIonCount) {
        double[] masses = entry.getMassArray();
        float[] intensities = entry.getIntensityArray();
        String[] ionTypes = entry.getIonTypeArray();
        int peptideLength = PeptideSequenceConverter.parseNormalizedUnimod(entry.getUnimodPeptideSequence())
                .getResidues()
                .length();

        assertEquals(masses.length, intensities.length);
        assertEquals(masses.length, ionTypes.length);
        assertTrue(masses.length >= minimumIonCount);

        for (int i = 0; i < ionTypes.length; i++) {
            assertTrue(masses[i] > 0.0, "Expected positive fragment m/z for " + ionTypes[i]);
            assertTrue(intensities[i] >= 0.01f, "Expected filtered intensity for " + ionTypes[i]);

            Matcher matcher = ION_TYPE_PATTERN.matcher(ionTypes[i]);
            assertTrue(matcher.matches(), "Unexpected ion type label: " + ionTypes[i]);
            int fragmentCharge = Integer.parseInt(matcher.group(1));
            int ionNumber = Integer.parseInt(matcher.group(3));
            assertTrue(fragmentCharge <= entry.getPrecursorCharge(), "Fragment charge exceeds precursor charge.");
            assertTrue(ionNumber > 0 && ionNumber < peptideLength, "Ion number out of peptide bounds.");
        }
    }

    private static List<String> topIonTypes(ChronologerLibraryEntry entry, int count) {
        float[] intensities = entry.getIntensityArray();
        String[] ionTypes = entry.getIonTypeArray();
        Integer[] indices = new Integer[intensities.length];
        for (int i = 0; i < intensities.length; i++) {
            indices[i] = i;
        }
        Arrays.sort(indices, (left, right) -> Float.compare(intensities[right], intensities[left]));

        List<String> top = new ArrayList<>();
        int limit = Math.min(count, indices.length);
        for (int i = 0; i < limit; i++) {
            top.add(ionTypes[indices[i]]);
        }
        return top;
    }

    private static int overlapCount(List<String> left, List<String> right) {
        int overlap = 0;
        for (String value : left) {
            if (right.contains(value)) {
                overlap++;
            }
        }
        return overlap;
    }

    private static int countAbove(float[] values, float threshold) {
        int count = 0;
        for (float value : values) {
            if (value > threshold) {
                count++;
            }
        }
        return count;
    }
}

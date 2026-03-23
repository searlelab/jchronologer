package org.searlelab.jchronologer.impl;

import static org.junit.jupiter.api.Assertions.assertDoesNotThrow;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.junit.jupiter.api.Test;
import org.searlelab.jchronologer.api.ChronologerLibraryEntry;
import org.searlelab.jchronologer.api.ChronologerLibraryOptions;
import org.searlelab.jchronologer.api.ChronologerLibraryPredictor;
import org.searlelab.jchronologer.api.LibraryPredictionRequest;
import org.searlelab.jchronologer.api.PrecursorCondition;
import org.searlelab.jchronologer.inference.ElectricianBatchPredictor;
import org.searlelab.jchronologer.preprocessing.ChronologerPreprocessor;
import org.searlelab.jchronologer.preprocessing.PeptideSequenceConverter;
import org.searlelab.jchronologer.preprocessing.PreprocessingMetadataLoader;
import org.searlelab.jchronologer.preprocessing.PreprocessingOutcome;

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
            assertTrue(entry.getCCS().isPresent());
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
            assertTrue(first.getCCS().isEmpty());
            assertTrue(second.getCCS().isEmpty());

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
            assertTrue(tmtCharge2.getCCS().isEmpty());
            assertTrue(tmtCharge3.getCCS().isEmpty());

            assertTrue(phosphoOx.getUnimodPeptideSequence().contains("UNIMOD:35"));
            assertTrue(phosphoOx.getUnimodPeptideSequence().contains("UNIMOD:21"));
            assertTrue(phosphoOx.getUnimodPeptideSequence().contains("UNIMOD:1"));
            assertTrue(phosphoOx.getCCS().isPresent());

            for (ChronologerLibraryEntry entry : entries) {
                assertEntryHasValidSpectrum(entry, 5);
                assertTrue(countAbove(entry.getIntensityArray(), 0.10f) >= 2,
                        "Expected >=2 high-intensity ions for " + entry.getUnimodPeptideSequence());
            }
        }
    }

    @Test
    void automaticChargeSelectionExpandsToAllChargesAtZeroThreshold() {
        LibraryPredictionRequest automaticRequest = new LibraryPredictionRequest("TASEFDSAIAQDK", 30.0, 0.0);

        try (ChronologerLibraryPredictor predictor = ChronologerFactory.createLibraryPredictorDefault()) {
            List<ChronologerLibraryEntry> entries = predictor.predict(List.of(automaticRequest));
            assertEquals(6, entries.size());

            for (int i = 0; i < entries.size(); i++) {
                ChronologerLibraryEntry entry = entries.get(i);
                assertEquals((byte) (i + 1), entry.getPrecursorCharge());
                assertEquals(30.0, entry.getPrecursorNce());
                assertEquals("[]-TASEFDSAIAQDK-[]", entry.getUnimodPeptideSequence());
                assertTrue(entry.getCCS().isPresent());
                assertEntryHasValidSpectrum(entry, 5);
            }
        }
    }

    @Test
    void automaticChargeSelectionSkipsWhenNoChargeMeetsThreshold() {
        double strictThreshold = probabilityJustAboveMaximumChargeProbability("TASEFDSAIAQDK");
        LibraryPredictionRequest automaticRequest = new LibraryPredictionRequest("TASEFDSAIAQDK", 30.0, strictThreshold);

        try (ChronologerLibraryPredictor predictor = ChronologerFactory.createLibraryPredictorDefault()) {
            List<ChronologerLibraryEntry> entries = predictor.predict(List.of(automaticRequest));
            assertEquals(0, entries.size());
        }
    }

    @Test
    void ccsPredictionsAreDenormalizedToPhysicalScale() {
        List<String> peptides = List.of(
                "[]-HC[UNIMOD:4]VDPAVIAAIISR-[]",
                "[]-TLLISSLSPALPAEHLEDR-[]",
                "[]-TPIGSFLGSLS-[]",
                "[]-ANAEKTSGSNVKIVKVKKE-[]",
                "[UNIMOD:737]-K[UNIMOD:737]PGLAITFAK[UNIMOD:737]-[]");
        List<LibraryPredictionRequest> requests = new ArrayList<>();
        for (String peptide : peptides) {
            requests.add(new LibraryPredictionRequest(peptide, 33.0, 0.0));
        }

        Map<String, float[]> expectedByPeptide = new LinkedHashMap<>();
        expectedByPeptide.put(peptides.get(0), new float[] {389.715f, 435.296f, 503.681f, 620.384f, 738.567f, 890.815f});
        expectedByPeptide.put(peptides.get(1), new float[] {456.165f, 501.746f, 570.131f, 686.834f, 805.017f, 957.265f});
        expectedByPeptide.put(peptides.get(2), new float[] {320.342f, 365.924f, 434.309f, 551.011f, 669.194f, 821.443f});
        expectedByPeptide.put(peptides.get(3), new float[] {459.073f, 504.655f, 573.040f, 689.743f, 807.925f, 960.174f});

        try (ChronologerLibraryPredictor predictor = ChronologerFactory.createLibraryPredictorDefault()) {
            List<ChronologerLibraryEntry> entries = predictor.predict(requests);
            assertEquals(30, entries.size());

            Map<String, List<ChronologerLibraryEntry>> entriesByPeptide = new LinkedHashMap<>();
            for (ChronologerLibraryEntry entry : entries) {
                entriesByPeptide.computeIfAbsent(entry.getUnimodPeptideSequence(), key -> new ArrayList<>()).add(entry);
            }

            for (Map.Entry<String, float[]> expected : expectedByPeptide.entrySet()) {
                List<ChronologerLibraryEntry> peptideEntries = entriesByPeptide.get(expected.getKey());
                assertEquals(6, peptideEntries.size(), "Expected six charge states for " + expected.getKey());
                for (int i = 0; i < peptideEntries.size(); i++) {
                    ChronologerLibraryEntry entry = peptideEntries.get(i);
                    assertEquals((byte) (i + 1), entry.getPrecursorCharge());
                    assertTrue(entry.getCCS().isPresent(), "Expected CCS for " + expected.getKey());
                    assertEquals(expected.getValue()[i], entry.getCCS().get(), 0.02f);
                }
            }

            List<ChronologerLibraryEntry> unsupported = entriesByPeptide.get(peptides.get(4));
            assertEquals(6, unsupported.size(), "Expected six charge states for " + peptides.get(4));
            for (ChronologerLibraryEntry entry : unsupported) {
                assertTrue(entry.getCCS().isEmpty(), "Expected unsupported CCS for " + peptides.get(4));
            }
        }
    }

    @Test
    void automaticChargeSelectionMatchesExpectedChargesAndDistributions() {
        List<String> peptides = List.of(
                "[]-HC[UNIMOD:4]VDPAVIAAIISR-[]",
                "[]-TLLISSLSPALPAEHLEDR-[]",
                "[]-TPIGSFLGSLS-[]",
                "[]-ANAEKTSGSNVKIVKVKKE-[]");
        List<LibraryPredictionRequest> requests = List.of(
                new LibraryPredictionRequest(peptides.get(0), 33.0, 0.05),
                new LibraryPredictionRequest(peptides.get(1), 33.0, 0.05),
                new LibraryPredictionRequest(peptides.get(2), 33.0, 0.05),
                new LibraryPredictionRequest(peptides.get(3), 33.0, 0.05));

        try (ChronologerLibraryPredictor predictor = ChronologerFactory.createLibraryPredictorDefault()) {
            List<ChronologerLibraryEntry> entries = predictor.predict(requests);
            assertEquals(6, entries.size());
            for (int i=0; i<entries.size(); i++) {
				System.out.println(i+") "+entries.get(i).getPrecursorCharge());
			}

            assertEquals(peptides.get(0), entries.get(0).getUnimodPeptideSequence());
            assertEquals((byte) 2, entries.get(0).getPrecursorCharge());
            assertEquals(peptides.get(0), entries.get(1).getUnimodPeptideSequence());
            assertEquals((byte) 3, entries.get(1).getPrecursorCharge());

            assertEquals(peptides.get(1), entries.get(2).getUnimodPeptideSequence());
            assertEquals((byte) 2, entries.get(2).getPrecursorCharge());
            assertEquals(peptides.get(1), entries.get(3).getUnimodPeptideSequence());
            assertEquals((byte) 3, entries.get(3).getPrecursorCharge());

            assertEquals(peptides.get(2), entries.get(4).getUnimodPeptideSequence());
            assertEquals((byte) 1, entries.get(4).getPrecursorCharge());

            assertEquals(peptides.get(3), entries.get(5).getUnimodPeptideSequence());
            assertEquals((byte) 4, entries.get(5).getPrecursorCharge());

            for (ChronologerLibraryEntry entry : entries) {
                assertEquals(33.0, entry.getPrecursorNce());
                assertTrue(entry.getCCS().isPresent());
                assertEntryHasValidSpectrum(entry, 5);
            }
        }

        Map<String, float[]> actual = electricianChargeDistributions(peptides);
        Map<String, byte[]> expected = new LinkedHashMap<>();
        // 0=0, 1=below 1%, 2=1-10%, 3=10-33%, 4=>33%
        expected.put(peptides.get(0), new byte[] {0, 4, 4, 0, 0, 0});
        expected.put(peptides.get(1), new byte[] {0, 2, 4, 0, 0, 0});
        expected.put(peptides.get(2), new byte[] {4, 1, 0, 0, 0, 0});
        expected.put(peptides.get(3), new byte[] {0, 0, 2, 4, 1, 0});

        for (Map.Entry<String, byte[]> entry : expected.entrySet()) {
            float[] observed = actual.get(entry.getKey());
            byte[] expectedDistribution = entry.getValue();
            assertEquals(6, observed.length);
            for (int i = 0; i < expectedDistribution.length; i++) {
            	boolean accept = switch (expectedDistribution[i]) {
					case 0 -> observed[i]<=0.0001f;
					case 1 -> observed[i]>0.0001f&&observed[i]<=0.01f;
					case 2 -> observed[i]>0.01f&&observed[i]<=0.1f;
					case 3 -> observed[i]>0.1f&&observed[i]<=0.33f;
					case 4 -> observed[i]>0.33f;
					default -> false;
				};
                assertTrue(accept,
                        "Unexpected charge probability for " + entry.getKey() + " at z=" + (i + 1)+" is "+observed[i]);
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
            assertThrows(IllegalArgumentException.class, () -> predictor.predict(List.of(
                    new LibraryPredictionRequest("PEPTIDEK", Double.NaN, 0.1))));
            assertThrows(IllegalArgumentException.class, () -> predictor.predict(List.of(
                    new LibraryPredictionRequest("PEPTIDEK", 9.99, 0.1))));
            assertThrows(IllegalArgumentException.class, () -> predictor.predict(List.of(
                    new LibraryPredictionRequest("PEPTIDEK", 61.0, 0.1))));
            assertThrows(IllegalArgumentException.class, () -> predictor.predict(List.of(
                    new LibraryPredictionRequest("PEPTIDEK", 30.0, Double.NaN))));
            assertThrows(IllegalArgumentException.class, () -> predictor.predict(List.of(
                    new LibraryPredictionRequest("PEPTIDEK", 30.0, -0.0001))));
            assertThrows(IllegalArgumentException.class, () -> predictor.predict(List.of(
                    new LibraryPredictionRequest("PEPTIDEK", 30.0, 1.0001))));
        }

        ChronologerLibraryPredictor closedPredictor = ChronologerFactory.createLibraryPredictorDefault();
        closedPredictor.close();
        assertDoesNotThrow(closedPredictor::close);
        assertThrows(IllegalStateException.class, closedPredictor::init);
        assertThrows(IllegalStateException.class, () -> closedPredictor.predict(List.of(
                new LibraryPredictionRequest("PEPTIDEK", List.of(new PrecursorCondition((byte) 2, 30.0))))));
    }

    @Test
    void ccsPredictionCanBeDisabledWithoutLoadingSculptorResources() {
        ChronologerLibraryOptions options = ChronologerLibraryOptions.builder()
                .ccsPredictionEnabled(false)
                .sculptorModelResource("models/missing_sculptor_model.torchscript.pt")
                .sculptorPreprocessingResource("models/missing_sculptor_preprocessing.json")
                .build();

        LibraryPredictionRequest request = new LibraryPredictionRequest(
                "TASEFDSAIAQDK",
                List.of(new PrecursorCondition((byte) 2, 27.0)));

        try (ChronologerLibraryPredictor predictor = ChronologerFactory.createLibraryPredictor(options)) {
            List<ChronologerLibraryEntry> entries = predictor.predict(List.of(request));
            assertEquals(1, entries.size());
            assertTrue(entries.get(0).getCCS().isEmpty());
            assertEntryHasValidSpectrum(entries.get(0), 5);
        }
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

    private static double probabilityJustAboveMaximumChargeProbability(String peptide) {
        ChronologerPreprocessor preprocessor = new ChronologerPreprocessor(
                PreprocessingMetadataLoader.loadFromClasspath(
                        ChronologerLibraryOptions.DEFAULT_ELECTRICIAN_PREPROCESSING_RESOURCE));
        PreprocessingOutcome outcome = preprocessor.preprocess(peptide);
        assertTrue(outcome.isAccepted(), "Failed to preprocess test peptide: " + peptide);

        try (ElectricianBatchPredictor predictor = new ElectricianBatchPredictor(
                ChronologerLibraryOptions.DEFAULT_ELECTRICIAN_MODEL_RESOURCE)) {
            float[][] distributions = predictor.predict(new long[][] {outcome.getTokenArray()});
            assertEquals(1, distributions.length);
            assertEquals(6, distributions[0].length);

            float[] normalized = normalizeDistribution(distributions[0]);
            float max = 0.0f;
            for (float value : normalized) {
                if (value > max) {
                    max = value;
                }
            }
            double strictThreshold = Math.nextUp((double) max);
            assertTrue(
                    strictThreshold <= 1.0 && strictThreshold > max,
                    "Could not construct strict threshold above max charge probability.");
            return strictThreshold;
        }
    }

    private static float[] normalizeDistribution(float[] rawDistribution) {
        float[] normalized = new float[rawDistribution.length];
        double sum = 0.0;
        for (int i = 0; i < rawDistribution.length; i++) {
            float value = rawDistribution[i];
            if (Float.isFinite(value) && value > 0.0f) {
                normalized[i] = value;
                sum += value;
            }
        }
        if (sum <= 0.0) {
            return normalized;
        }
        for (int i = 0; i < normalized.length; i++) {
            normalized[i] = (float) (normalized[i] / sum);
        }
        return normalized;
    }

    private static Map<String, float[]> electricianChargeDistributions(List<String> unimodPeptides) {
        ChronologerPreprocessor preprocessor = new ChronologerPreprocessor(
                PreprocessingMetadataLoader.loadFromClasspath(
                        ChronologerLibraryOptions.DEFAULT_ELECTRICIAN_PREPROCESSING_RESOURCE));
        long[][] tokenBatch = new long[unimodPeptides.size()][];
        for (int i = 0; i < unimodPeptides.size(); i++) {
            String massEncoded = PeptideSequenceConverter.unimodToMassEncoded(unimodPeptides.get(i));
            PreprocessingOutcome outcome = preprocessor.preprocess(massEncoded);
            assertTrue(outcome.isAccepted(), "Failed to preprocess test peptide: " + unimodPeptides.get(i));
            tokenBatch[i] = outcome.getTokenArray();
        }

        try (ElectricianBatchPredictor predictor = new ElectricianBatchPredictor(
                ChronologerLibraryOptions.DEFAULT_ELECTRICIAN_MODEL_RESOURCE)) {
            float[][] distributions = predictor.predict(tokenBatch);
            assertEquals(unimodPeptides.size(), distributions.length);

            LinkedHashMap<String, float[]> byPeptide = new LinkedHashMap<>();
            for (int i = 0; i < unimodPeptides.size(); i++) {
                byPeptide.put(unimodPeptides.get(i), normalizeDistribution(distributions[i]));
            }
            return byPeptide;
        }
    }
}

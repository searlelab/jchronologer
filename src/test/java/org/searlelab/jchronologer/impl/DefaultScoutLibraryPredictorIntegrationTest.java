package org.searlelab.jchronologer.impl;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertDoesNotThrow;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.lang.reflect.Field;
import java.util.List;
import java.util.concurrent.ExecutorService;
import org.junit.jupiter.api.Test;
import org.searlelab.jchronologer.api.ChronologerLibraryEntry;
import org.searlelab.jchronologer.api.ChronologerLibraryOptions;
import org.searlelab.jchronologer.api.ChronologerLibraryPredictor;
import org.searlelab.jchronologer.api.LibraryPredictionRequest;
import org.searlelab.jchronologer.api.PrecursorCondition;

class DefaultScoutLibraryPredictorIntegrationTest {

    @Test
    void explicitPredictionReturnsCcsAndNoThreePlusIons() {
        LibraryPredictionRequest request = new LibraryPredictionRequest(
                "TASEFDSAIAQDK",
                List.of(new PrecursorCondition((byte) 2, 27.0)));

        try (ChronologerLibraryPredictor predictor = ChronologerFactory.createFastLibraryPredictorDefault()) {
            List<ChronologerLibraryEntry> entries = predictor.predict(List.of(request));
            assertEquals(1, entries.size());
            ChronologerLibraryEntry entry = entries.get(0);
            assertTrue(entry.getPrecursorMz() > 0.0);
            assertTrue(entry.getRetentionTimeInSeconds() > 0.0f);
            assertTrue(entry.getCCS().isPresent());
            assertTrue(entry.getCCS().get() > 100.0f);
            for (String ionType : entry.getIonTypeArray()) {
                assertTrue(!ionType.startsWith("3+"));
            }
        }
    }

    @Test
    void automaticChargeSelectionUsesScoutOwnedExpansion() {
        LibraryPredictionRequest automaticRequest = new LibraryPredictionRequest("TASEFDSAIAQDK", 30.0, 0.0);

        try (ChronologerLibraryPredictor predictor = ChronologerFactory.createFastLibraryPredictorDefault()) {
            List<ChronologerLibraryEntry> entries = predictor.predict(List.of(automaticRequest));
            assertTrue(!entries.isEmpty());
            byte previousCharge = 0;
            for (ChronologerLibraryEntry entry : entries) {
                assertTrue(entry.getCCS().isPresent());
                assertTrue(entry.getChargeProbability().isPresent());
                assertTrue(entry.getChargeProbability().orElseThrow() >= 0.0f);
                assertTrue(entry.getPrecursorCharge() >= previousCharge);
                previousCharge = entry.getPrecursorCharge();
            }
        }
    }

    @Test
    void automaticChargeSelectionRespectsMinimumProbabilityThreshold() {
        double threshold = 0.20;
        LibraryPredictionRequest automaticRequest = new LibraryPredictionRequest("TASEFDSAIAQDK", 30.0, threshold);

        try (ChronologerLibraryPredictor predictor = ChronologerFactory.createFastLibraryPredictorDefault()) {
            List<ChronologerLibraryEntry> entries = predictor.predict(List.of(automaticRequest));
            assertTrue(!entries.isEmpty());
            for (ChronologerLibraryEntry entry : entries) {
                assertTrue(entry.getChargeProbability().orElseThrow() >= threshold);
            }
        }
    }

    @Test
    void ccsPredictionFlagIsIgnoredForFastMode() {
        ChronologerLibraryOptions options = ChronologerLibraryOptions.builder()
                .ccsPredictionEnabled(false)
                .build();
        LibraryPredictionRequest request = new LibraryPredictionRequest(
                "TASEFDSAIAQDK",
                List.of(new PrecursorCondition((byte) 2, 27.0)));

        try (ChronologerLibraryPredictor predictor = ChronologerFactory.createFastLibraryPredictor(options)) {
            List<ChronologerLibraryEntry> entries = predictor.predict(List.of(request));
            assertEquals(1, entries.size());
            assertTrue(entries.get(0).getCCS().isPresent());
        }
    }

    @Test
    void parallelScoutPredictionMatchesSingleThreadOutput() {
        List<LibraryPredictionRequest> requests = List.of(
                new LibraryPredictionRequest(
                        "TASEFDSAIAQDK",
                        List.of(new PrecursorCondition((byte) 2, 27.0), new PrecursorCondition((byte) 3, 31.0))),
                new LibraryPredictionRequest(
                        "[]-HC[UNIMOD:4]VDPAVIAAIISR-[]",
                        List.of(new PrecursorCondition((byte) 2, 33.0), new PrecursorCondition((byte) 3, 33.0))));
        ChronologerLibraryOptions serialOptions = ChronologerLibraryOptions.builder()
                .inferenceThreads(1)
                .batchSize(2)
                .build();
        ChronologerLibraryOptions parallelOptions = ChronologerLibraryOptions.builder()
                .inferenceThreads(4)
                .batchSize(2)
                .build();

        List<ChronologerLibraryEntry> serialEntries;
        List<ChronologerLibraryEntry> parallelEntries;
        try (ChronologerLibraryPredictor serial = ChronologerFactory.createFastLibraryPredictor(serialOptions);
                ChronologerLibraryPredictor parallel = ChronologerFactory.createFastLibraryPredictor(parallelOptions)) {
            serialEntries = serial.predict(requests);
            parallelEntries = parallel.predict(requests);
        }

        assertEquals(serialEntries.size(), parallelEntries.size());
        for (int i = 0; i < serialEntries.size(); i++) {
            ChronologerLibraryEntry expected = serialEntries.get(i);
            ChronologerLibraryEntry observed = parallelEntries.get(i);
            assertEquals(expected.getUnimodPeptideSequence(), observed.getUnimodPeptideSequence());
            assertEquals(expected.getPrecursorCharge(), observed.getPrecursorCharge());
            assertEquals(expected.getPrecursorNce(), observed.getPrecursorNce(), 1e-9);
            assertEquals(expected.getPrecursorMz(), observed.getPrecursorMz(), 1e-9);
            assertEquals(expected.getRetentionTimeInSeconds(), observed.getRetentionTimeInSeconds(), 1e-4f);
            assertEquals(expected.getCCS().orElseThrow(), observed.getCCS().orElseThrow(), 1e-3f);
            assertEquals(expected.getChargeProbability(), observed.getChargeProbability());
            assertArrayEquals(expected.getIonTypeArray(), observed.getIonTypeArray());
            assertArrayEquals(expected.getIntensityArray(), observed.getIntensityArray(), 1e-5f);
        }
    }

    @Test
    void predictorValidatesLifecycleAndCanCloseTwice() throws Exception {
        DefaultScoutLibraryPredictor predictor = new DefaultScoutLibraryPredictor(ChronologerLibraryOptions.builder()
                .inferenceThreads(3)
                .batchSize(2)
                .build());
        ExecutorService executor = extractScoutExecutor(predictor);
        predictor.close();
        assertTrue(executor.isShutdown());
        assertDoesNotThrow(predictor::close);
        assertThrows(IllegalStateException.class, predictor::init);
    }

    private static ExecutorService extractScoutExecutor(DefaultScoutLibraryPredictor predictor) throws Exception {
        Field field = DefaultScoutLibraryPredictor.class.getDeclaredField("scoutExecutor");
        field.setAccessible(true);
        return (ExecutorService) field.get(predictor);
    }
}

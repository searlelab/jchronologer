package org.searlelab.jchronologer.impl;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertDoesNotThrow;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutorService;
import org.junit.jupiter.api.Test;
import org.searlelab.jchronologer.api.ChronologerLibraryEntry;
import org.searlelab.jchronologer.api.ChronologerLibraryOptions;
import org.searlelab.jchronologer.api.ChronologerLibraryPredictor;
import org.searlelab.jchronologer.api.LibraryPredictionRequest;
import org.searlelab.jchronologer.api.PrecursorCondition;

class DefaultChronologerLibraryPredictorParallelTest {

    @Test
    void parallelCartographerPredictionMatchesSingleThreadOutput() {
        List<LibraryPredictionRequest> requests = buildExplicitRequests();

        ChronologerLibraryOptions serialOptions = ChronologerLibraryOptions.builder()
                .inferenceThreads(1)
                .batchSize(2048)
                .cartographerBatchSize(3)
                .build();
        ChronologerLibraryOptions parallelOptions = ChronologerLibraryOptions.builder()
                .inferenceThreads(4)
                .batchSize(2048)
                .cartographerBatchSize(3)
                .build();

        List<ChronologerLibraryEntry> serialEntries;
        List<ChronologerLibraryEntry> parallelEntries;
        try (ChronologerLibraryPredictor serial = ChronologerFactory.createLibraryPredictor(serialOptions);
                ChronologerLibraryPredictor parallel = ChronologerFactory.createLibraryPredictor(parallelOptions)) {
            serialEntries = serial.predict(requests);
            parallelEntries = parallel.predict(requests);
        }

        assertEquals(serialEntries.size(), parallelEntries.size());
        for (int i = 0; i < serialEntries.size(); i++) {
            assertEquivalentEntries(serialEntries.get(i), parallelEntries.get(i));
        }
    }

    @Test
    void preservesEntryOrderingAcrossParallelBatches() {
        List<LibraryPredictionRequest> requests = List.of(
                new LibraryPredictionRequest(
                        "[]-TASEFDSAIAQDK-[]",
                        List.of(
                                new PrecursorCondition((byte) 2, 27.0),
                                new PrecursorCondition((byte) 3, 31.0))),
                new LibraryPredictionRequest(
                        "[]-ACDEFGHIK-[]",
                        List.of(
                                new PrecursorCondition((byte) 2, 30.0),
                                new PrecursorCondition((byte) 3, 34.0))),
                new LibraryPredictionRequest(
                        "[]-MPEPTIDER-[]",
                        List.of(
                                new PrecursorCondition((byte) 2, 29.0),
                                new PrecursorCondition((byte) 3, 33.0))));

        List<String> expected = List.of(
                "[]-TASEFDSAIAQDK-[]|2",
                "[]-TASEFDSAIAQDK-[]|3",
                "[]-ACDEFGHIK-[]|2",
                "[]-ACDEFGHIK-[]|3",
                "[]-MPEPTIDER-[]|2",
                "[]-MPEPTIDER-[]|3");

        ChronologerLibraryOptions parallelOptions = ChronologerLibraryOptions.builder()
                .inferenceThreads(4)
                .batchSize(2048)
                .cartographerBatchSize(2)
                .build();

        List<ChronologerLibraryEntry> entries;
        try (ChronologerLibraryPredictor predictor = ChronologerFactory.createLibraryPredictor(parallelOptions)) {
            entries = predictor.predict(requests);
        }

        assertEquals(expected.size(), entries.size());
        for (int i = 0; i < entries.size(); i++) {
            ChronologerLibraryEntry entry = entries.get(i);
            assertEquals(expected.get(i), entry.getUnimodPeptideSequence() + "|" + entry.getPrecursorCharge());
        }
    }

    @Test
    void closeShutsDownCartographerExecutorAndRemainsIdempotent() throws Exception {
        ChronologerLibraryOptions options = ChronologerLibraryOptions.builder()
                .inferenceThreads(3)
                .batchSize(2048)
                .cartographerBatchSize(3)
                .build();

        DefaultChronologerLibraryPredictor predictor = new DefaultChronologerLibraryPredictor(options);
        ExecutorService executor = extractCartographerExecutor(predictor);
        assertNotNull(executor);
        assertFalse(executor.isShutdown());

        predictor.close();
        assertTrue(executor.isShutdown());
        assertDoesNotThrow(predictor::close);
    }

    private static ExecutorService extractCartographerExecutor(DefaultChronologerLibraryPredictor predictor) throws Exception {
        Field field = DefaultChronologerLibraryPredictor.class.getDeclaredField("cartographerExecutor");
        field.setAccessible(true);
        return (ExecutorService) field.get(predictor);
    }

    private static List<LibraryPredictionRequest> buildExplicitRequests() {
        List<LibraryPredictionRequest> requests = new ArrayList<>();
        requests.add(new LibraryPredictionRequest(
                "TASEFDSAIAQDK",
                List.of(new PrecursorCondition((byte) 2, 27.0), new PrecursorCondition((byte) 3, 31.0))));
        requests.add(new LibraryPredictionRequest(
                "ACDM[+15.994915]STY[+79.966331]K[+42.010565]",
                List.of(new PrecursorCondition((byte) 2, 30.0), new PrecursorCondition((byte) 3, 34.0))));
        requests.add(new LibraryPredictionRequest(
                "K[458.325864]PGLAITFAK[229.162932]",
                List.of(new PrecursorCondition((byte) 2, 30.0), new PrecursorCondition((byte) 3, 35.0))));
        requests.add(new LibraryPredictionRequest(
                "[]-HC[UNIMOD:4]VDPAVIAAIISR-[]",
                List.of(new PrecursorCondition((byte) 2, 33.0), new PrecursorCondition((byte) 3, 33.0))));
        requests.add(new LibraryPredictionRequest(
                "[]-TPIGSFLGSLS-[]",
                List.of(new PrecursorCondition((byte) 1, 33.0), new PrecursorCondition((byte) 2, 33.0))));
        return requests;
    }

    private static void assertEquivalentEntries(ChronologerLibraryEntry expected, ChronologerLibraryEntry observed) {
        assertEquals(expected.getUnimodPeptideSequence(), observed.getUnimodPeptideSequence());
        assertEquals(expected.getPrecursorCharge(), observed.getPrecursorCharge());
        assertEquals(expected.getPrecursorNce(), observed.getPrecursorNce(), 1e-9);
        assertEquals(expected.getPrecursorMz(), observed.getPrecursorMz(), 1e-9);
        assertEquals(expected.getRetentionTimeInSeconds(), observed.getRetentionTimeInSeconds(), 1e-5f);
        assertArrayEquals(expected.getMassArray(), observed.getMassArray(), 1e-6);
        assertArrayEquals(expected.getIntensityArray(), observed.getIntensityArray(), 1e-6f);
        assertArrayEquals(expected.getIonTypeArray(), observed.getIonTypeArray());
    }
}

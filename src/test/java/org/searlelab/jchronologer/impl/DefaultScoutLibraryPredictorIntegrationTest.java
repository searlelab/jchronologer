package org.searlelab.jchronologer.impl;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertDoesNotThrow;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.List;
import java.util.Optional;
import java.util.concurrent.ExecutorService;
import org.junit.jupiter.api.Test;
import org.searlelab.jchronologer.api.ChronologerLibraryEntry;
import org.searlelab.jchronologer.api.ChronologerLibraryOptions;
import org.searlelab.jchronologer.api.ChronologerLibraryPredictor;
import org.searlelab.jchronologer.api.LibraryPredictionRequest;
import org.searlelab.jchronologer.api.PrecursorCondition;
import org.searlelab.jchronologer.inference.ScoutSpectrumDecoder;
import org.searlelab.jchronologer.preprocessing.PeptideSequenceConverter;

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
    void electricianAutomaticChargeSelectionProducesProbabilities() {
        ChronologerLibraryOptions options = ChronologerLibraryOptions.builder()
                .batchSize(2)
                .build();
        LibraryPredictionRequest automaticRequest = new LibraryPredictionRequest("TASEFDSAIAQDK", 30.0, 0.0);

        try (DefaultScoutLibraryPredictor predictor = new DefaultScoutLibraryPredictor(
                options,
                DefaultScoutLibraryPredictor.AutomaticChargeSelectionMode.ELECTRICIAN_SELECTION)) {
            List<ChronologerLibraryEntry> entries = predictor.predict(List.of(automaticRequest));
            assertFalse(entries.isEmpty());
            for (ChronologerLibraryEntry entry : entries) {
                assertTrue(entry.getChargeProbability().isPresent());
                assertTrue(entry.getChargeProbability().orElseThrow() >= 0.0f);
            }
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

    @Test
    void predictRejectsNullEmptyAndInvalidRequests() {
        try (DefaultScoutLibraryPredictor predictor = new DefaultScoutLibraryPredictor(
                ChronologerLibraryOptions.builder().build())) {
            assertEquals(List.of(), predictor.predict(null));
            assertEquals(List.of(), predictor.predict(List.of()));

            List<LibraryPredictionRequest> requestsWithNull = new ArrayList<>();
            requestsWithNull.add(null);
            IllegalArgumentException nullRequest = assertThrows(
                    IllegalArgumentException.class,
                    () -> predictor.predict(requestsWithNull));
            assertEquals("Request at index 0 is null.", nullRequest.getMessage());

            IllegalArgumentException missingConditions = assertThrows(
                    IllegalArgumentException.class,
                    () -> predictor.predict(List.of(new LibraryPredictionRequest("TASEFDSAIAQDK", List.of()))));
            assertEquals("Request at index 0 has no precursor conditions.", missingConditions.getMessage());
        }
    }

    @Test
    void predictRejectsUnsupportedChargesAndAutomaticProbability() {
        try (DefaultScoutLibraryPredictor predictor = new DefaultScoutLibraryPredictor(
                ChronologerLibraryOptions.builder().build())) {
            IllegalArgumentException badCharge = assertThrows(
                    IllegalArgumentException.class,
                    () -> predictor.predict(List.of(
                            new LibraryPredictionRequest(
                                    "TASEFDSAIAQDK",
                                    List.of(new PrecursorCondition((byte) 0, 27.0))))));
            assertTrue(badCharge.getMessage().contains("outside supported range"));

            IllegalArgumentException badProbability = assertThrows(
                    IllegalArgumentException.class,
                    () -> predictor.predict(List.of(new LibraryPredictionRequest("TASEFDSAIAQDK", 30.0, 1.5))));
            assertTrue(badProbability.getMessage().contains("outside the supported range 0.0-1.0"));

            IllegalArgumentException badNce = assertThrows(
                    IllegalArgumentException.class,
                    () -> predictor.predict(List.of(new LibraryPredictionRequest("TASEFDSAIAQDK", 9.0, 0.1))));
            assertTrue(badNce.getMessage().contains("outside the supported range"));
        }
    }

    @Test
    void helperMethodsNormalizeDistributionsAndFormatDiagnostics() throws Exception {
        try (DefaultScoutLibraryPredictor predictor = new DefaultScoutLibraryPredictor(
                ChronologerLibraryOptions.builder().build())) {
            float[] normalized = (float[]) invoke(
                    predictor,
                    "normalizeChargeDistribution",
                    new Class<?>[] {float[].class},
                    new float[] {Float.NaN, -1.0f, 2.0f, 6.0f});
            assertArrayEquals(new float[] {0.0f, 0.0f, 0.25f, 0.75f}, normalized, 1e-6f);

            float[] allInvalid = (float[]) invoke(
                    predictor,
                    "normalizeChargeDistribution",
                    new Class<?>[] {float[].class},
                    new float[] {Float.NaN, -1.0f});
            assertArrayEquals(new float[] {0.0f, 0.0f}, allInvalid, 1e-6f);

            String diagnostic = (String) invoke(
                    predictor,
                    "formatMissedSeedDiagnostic",
                    new Class<?>[] {String.class, byte.class, double.class, float[].class},
                    "PEPTIDE",
                    (byte) 3,
                    0.20,
                    new float[] {0.1f, 0.2f, 0.7f});
            assertTrue(diagnostic.contains("peptide=PEPTIDE"));
            assertTrue(diagnostic.contains("seedCharge=3"));
            assertTrue(diagnostic.contains("1=0.1000"));
            assertTrue(diagnostic.contains("3=0.7000"));
        }
    }

    @Test
    void helperMethodsEstimateSeedChargeAndElapsedMillis() throws Exception {
        try (DefaultScoutLibraryPredictor predictor = new DefaultScoutLibraryPredictor(
                ChronologerLibraryOptions.builder().build())) {
            Object parsedShort = PeptideSequenceConverter.parseNormalizedUnimod("[]-AK-[]");
            Object parsedLong = PeptideSequenceConverter.parseNormalizedUnimod("[]-AKRKRKAAAKR-[]");

            int shortCharge = (int) invoke(
                    predictor,
                    "selectSeedCharge",
                    new Class<?>[] {parsedShort.getClass()},
                    parsedShort);
            int longCharge = (int) invoke(
                    predictor,
                    "selectSeedCharge",
                    new Class<?>[] {parsedLong.getClass()},
                    parsedLong);

            assertTrue(shortCharge >= 2);
            assertTrue(longCharge >= shortCharge);

            long elapsed = (long) invokeStatic(
                    DefaultScoutLibraryPredictor.class,
                    "elapsedMillis",
                    new Class<?>[] {long.class},
                    2_500_000L);
            assertEquals(2L, elapsed);
        }
    }

    @Test
    void libraryEntryOmitsNonFiniteCcs() throws Exception {
        try (DefaultScoutLibraryPredictor predictor = new DefaultScoutLibraryPredictor(
                ChronologerLibraryOptions.builder().build())) {
            Object parsed = PeptideSequenceConverter.parseNormalizedUnimod("[]-TASEFDSAIAQDK-[]");
            Object job = buildExplicitPredictionJob(parsed, (byte) 2, 27.0);
            Object record = buildScoutPredictionRecord(
                    new float[ScoutSpectrumDecoder.VECTOR_LENGTH],
                    new float[] {0.0f},
                    new float[] {Float.NaN},
                    new float[] {1.0f, 0.0f, 0.0f, 0.0f, 0.0f});

            ChronologerLibraryEntry entry = (ChronologerLibraryEntry) invoke(
                    predictor,
                    "toLibraryEntry",
                    new Class<?>[] {job.getClass(), record.getClass()},
                    job,
                    record);

            assertFalse(entry.getCCS().isPresent());
            assertEquals(Optional.empty(), entry.getChargeProbability());
            assertEquals("[]-TASEFDSAIAQDK-[]", entry.getUnimodPeptideSequence());
            assertEquals(2, entry.getPrecursorCharge());
            assertTrue(entry.getRetentionTimeInSeconds() > 0.0f);
            assertEquals(0, entry.getMassArray().length);
            assertEquals(0, entry.getIntensityArray().length);
            assertEquals(0, entry.getIonTypeArray().length);
        }
    }

    private static ExecutorService extractScoutExecutor(DefaultScoutLibraryPredictor predictor) throws Exception {
        Field field = DefaultScoutLibraryPredictor.class.getDeclaredField("scoutExecutor");
        field.setAccessible(true);
        return (ExecutorService) field.get(predictor);
    }

    private static Object buildExplicitPredictionJob(Object parsed, byte charge, double nce) throws Exception {
        Class<?> jobClass = Class.forName("org.searlelab.jchronologer.impl.DefaultScoutLibraryPredictor$PredictionJob");
        Method explicit = jobClass.getDeclaredMethod(
                "explicit",
                String.class,
                parsed.getClass(),
                long[].class,
                byte.class,
                double.class);
        explicit.setAccessible(true);
        return explicit.invoke(null, "[]-TASEFDSAIAQDK-[]", parsed, new long[] {1L, 2L, 3L}, charge, nce);
    }

    private static Object buildScoutPredictionRecord(float[] ms2, float[] irt, float[] ccs, float[] chargeDist) throws Exception {
        Class<?> recordClass = Class.forName("org.searlelab.jchronologer.impl.DefaultScoutLibraryPredictor$ScoutPredictionRecord");
        var constructor = recordClass.getDeclaredConstructor(float[].class, float[].class, float[].class, float[].class);
        constructor.setAccessible(true);
        return constructor.newInstance(ms2, irt, ccs, chargeDist);
    }

    private static Object invoke(Object target, String methodName, Class<?>[] parameterTypes, Object... args) throws Exception {
        Method method = target.getClass().getDeclaredMethod(methodName, parameterTypes);
        method.setAccessible(true);
        try {
            return method.invoke(target, args);
        } catch (InvocationTargetException e) {
            Throwable cause = e.getCause();
            if (cause instanceof Exception exception) {
                throw exception;
            }
            throw e;
        }
    }

    private static Object invokeStatic(Class<?> type, String methodName, Class<?>[] parameterTypes, Object... args) throws Exception {
        Method method = type.getDeclaredMethod(methodName, parameterTypes);
        method.setAccessible(true);
        return method.invoke(null, args);
    }
}

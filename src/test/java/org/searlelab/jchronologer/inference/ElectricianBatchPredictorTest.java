package org.searlelab.jchronologer.inference;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;
import java.util.Random;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import org.junit.jupiter.api.Test;
import org.searlelab.jchronologer.api.ChronologerOptions;
import org.searlelab.jchronologer.preprocessing.ChronologerPreprocessor;
import org.searlelab.jchronologer.preprocessing.PreprocessingMetadataLoader;
import org.searlelab.jchronologer.preprocessing.PreprocessingOutcome;

class ElectricianBatchPredictorTest {
    private static final int RANDOM_PEPTIDE_COUNT = 500;
    private static final int RANDOM_MIN_LENGTH = 7;
    private static final int RANDOM_MAX_LENGTH = 31;
    private static final int WARMUP_ROUNDS = 1;
    private static final int MEASURE_ATTEMPTS = 10;
    private static final int INFERENCE_REPEATS_PER_ATTEMPT = 1;
    private static final double EXPECTED_MIN_SPEEDUP = 1.5d;
    private static final char[] CANONICAL_AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY".toCharArray();
    private static final long RANDOM_SEED = 23_837_401L;


    @Test
    void predictReturnsExpectedOutputShape() {
        try (ElectricianBatchPredictor predictor = new ElectricianBatchPredictor(
        		ChronologerOptions.DEFAULT_ELECTRICIAN_MODEL_RESOURCE)) {
            long[][] tokens = new long[][] {new long[33], new long[33], new long[33]};

            float[][] output = predictor.predict(tokens);

            assertEquals(3, output.length);
            assertEquals(6, output[0].length);
            assertEquals(6, output[1].length);
            assertEquals(6, output[2].length);
        }
    }

    @Test
    void predictWithEmptyBatchReturnsEmptyOutput() {
        try (ElectricianBatchPredictor predictor = new ElectricianBatchPredictor(
        		ChronologerOptions.DEFAULT_ELECTRICIAN_MODEL_RESOURCE)) {
            float[][] output = predictor.predict(new long[][] {});
            assertEquals(0, output.length);
        }
    }

    @Test
    void constructorRejectsMissingModelResource() {
        assertThrows(
                IllegalArgumentException.class,
                () -> new ElectricianBatchPredictor("models/not_a_real_electrician_model.pt"));
    }

    @Test
    void constructorWrapsMalformedModelLoadFailures() {
        IllegalStateException error = assertThrows(
                IllegalStateException.class,
                () -> new ElectricianBatchPredictor("models/invalid_cartographer_model.torchscript.pt"));
        assertTrue(error.getMessage().contains("Failed to load Electrician model"));
    }

    @Test
    void predictInitializesThreadLocalPredictorInWorkerThread() throws Exception {
        try (ElectricianBatchPredictor predictor = new ElectricianBatchPredictor(
        		ChronologerOptions.DEFAULT_ELECTRICIAN_MODEL_RESOURCE)) {
            ExecutorService executor = Executors.newSingleThreadExecutor();
            try {
                Future<float[][]> outputFuture = executor.submit(() -> predictor.predict(new long[][] {new long[33]}));
                float[][] output = outputFuture.get();
                assertEquals(1, output.length);
                assertEquals(6, output[0].length);
            } catch (InterruptedException e) {
                Thread.currentThread().interrupt();
                throw e;
            } catch (ExecutionException e) {
                throw new AssertionError("Worker thread prediction failed", e.getCause());
            } finally {
                executor.shutdownNow();
            }
        }
    }

    @Test
    void electricianIsInTwoXBallparkComparedToChronologerOnRandomPeptides() {
        long testStartNanos = System.nanoTime();
        List<String> randomPeptides = generateRandomPeptides();
        ChronologerPreprocessor chronologerPreprocessor = new ChronologerPreprocessor(
                PreprocessingMetadataLoader.loadFromClasspath(ChronologerOptions.DEFAULT_PREPROCESSING_RESOURCE));
        ChronologerPreprocessor electricianPreprocessor = new ChronologerPreprocessor(
                PreprocessingMetadataLoader.loadFromClasspath(
                		ChronologerOptions.DEFAULT_ELECTRICIAN_PREPROCESSING_RESOURCE));

        long[][] chronologerTokens = tokenize(randomPeptides, chronologerPreprocessor);
        long[][] electricianTokens = tokenize(randomPeptides, electricianPreprocessor);

        try (BatchPredictor chronologerPredictor = new BatchPredictor(ChronologerOptions.DEFAULT_MODEL_RESOURCE);
                ElectricianBatchPredictor electricianPredictor = new ElectricianBatchPredictor(
                		ChronologerOptions.DEFAULT_ELECTRICIAN_MODEL_RESOURCE)) {
            warmUp(chronologerPredictor, chronologerTokens, electricianPredictor, electricianTokens);

            double[] chronologerPerBatchNanos = new double[MEASURE_ATTEMPTS];
            double[] electricianPerBatchNanos = new double[MEASURE_ATTEMPTS];
            for (int i = 0; i < MEASURE_ATTEMPTS; i++) {
                if ((i & 1) == 0) {
                    long chronologerNanos = runChronologerInferenceRepeats(
                            chronologerPredictor,
                            chronologerTokens,
                            INFERENCE_REPEATS_PER_ATTEMPT);
                    long electricianNanos = runElectricianInferenceRepeats(
                            electricianPredictor,
                            electricianTokens,
                            INFERENCE_REPEATS_PER_ATTEMPT);
                    chronologerPerBatchNanos[i] = chronologerNanos / (double) INFERENCE_REPEATS_PER_ATTEMPT;
                    electricianPerBatchNanos[i] = electricianNanos / (double) INFERENCE_REPEATS_PER_ATTEMPT;
                } else {
                    long electricianNanos = runElectricianInferenceRepeats(
                            electricianPredictor,
                            electricianTokens,
                            INFERENCE_REPEATS_PER_ATTEMPT);
                    long chronologerNanos = runChronologerInferenceRepeats(
                            chronologerPredictor,
                            chronologerTokens,
                            INFERENCE_REPEATS_PER_ATTEMPT);
                    chronologerPerBatchNanos[i] = chronologerNanos / (double) INFERENCE_REPEATS_PER_ATTEMPT;
                    electricianPerBatchNanos[i] = electricianNanos / (double) INFERENCE_REPEATS_PER_ATTEMPT;
                }
            }

            double chronologerMedianNanos = median(chronologerPerBatchNanos);
            double electricianMedianNanos = median(electricianPerBatchNanos);
            double speedup = chronologerMedianNanos / electricianMedianNanos;
            double testElapsedMillis = (System.nanoTime() - testStartNanos) / 1_000_000.0d;

            System.out.printf(
                    Locale.US,
                    "Electrician benchmark: speedup=%.2fx, chronologerMedian=%.3f ms, electricianMedian=%.3f ms, testElapsed=%.3f ms%n",
                    speedup,
                    chronologerMedianNanos / 1_000_000.0d,
                    electricianMedianNanos / 1_000_000.0d,
                    testElapsedMillis);

            assertTrue(
                    speedup >= EXPECTED_MIN_SPEEDUP,
                    String.format(
                            Locale.US,
                            "Expected Electrician to be in the ~2x speedup ballpark (>= %.2fx), but measured %.2fx (Chronologer median=%.3f ms, Electrician median=%.3f ms).",
                            EXPECTED_MIN_SPEEDUP,
                            speedup,
                            chronologerMedianNanos / 1_000_000.0d,
                            electricianMedianNanos / 1_000_000.0d));
        }
    }

    private static void warmUp(
            BatchPredictor chronologerPredictor,
            long[][] chronologerTokens,
            ElectricianBatchPredictor electricianPredictor,
            long[][] electricianTokens) {
        for (int i = 0; i < WARMUP_ROUNDS; i++) {
            assertEquals(RANDOM_PEPTIDE_COUNT, chronologerPredictor.predict(chronologerTokens).length);
            assertEquals(RANDOM_PEPTIDE_COUNT, electricianPredictor.predict(electricianTokens).length);
        }
    }

    private static long runChronologerInferenceRepeats(
            BatchPredictor chronologerPredictor,
            long[][] chronologerTokens,
            int repeats) {
        long start = System.nanoTime();
        for (int i = 0; i < repeats; i++) {
            float[] output = chronologerPredictor.predict(chronologerTokens);
            assertEquals(RANDOM_PEPTIDE_COUNT, output.length);
        }
        return System.nanoTime() - start;
    }

    private static long runElectricianInferenceRepeats(
            ElectricianBatchPredictor electricianPredictor,
            long[][] electricianTokens,
            int repeats) {
        long start = System.nanoTime();
        for (int i = 0; i < repeats; i++) {
            float[][] output = electricianPredictor.predict(electricianTokens);
            assertEquals(RANDOM_PEPTIDE_COUNT, output.length);
        }
        return System.nanoTime() - start;
    }

    private static List<String> generateRandomPeptides() {
        Random random = new Random(RANDOM_SEED);
        List<String> peptides = new ArrayList<>(RANDOM_PEPTIDE_COUNT);
        for (int i = 0; i < RANDOM_PEPTIDE_COUNT; i++) {
            int length = RANDOM_MIN_LENGTH + random.nextInt(RANDOM_MAX_LENGTH - RANDOM_MIN_LENGTH + 1);
            StringBuilder peptide = new StringBuilder(length);
            for (int position = 0; position < length; position++) {
                peptide.append(CANONICAL_AMINO_ACIDS[random.nextInt(CANONICAL_AMINO_ACIDS.length)]);
            }
            peptides.add(peptide.toString());
        }
        return peptides;
    }

    private static long[][] tokenize(List<String> peptides, ChronologerPreprocessor preprocessor) {
        long[][] tokens = new long[peptides.size()][];
        for (int i = 0; i < peptides.size(); i++) {
            String peptide = peptides.get(i);
            PreprocessingOutcome outcome = preprocessor.preprocess(peptide);
            assertTrue(outcome.isAccepted(), "Failed to preprocess random peptide: " + peptide);
            tokens[i] = outcome.getTokenArray();
        }
        return tokens;
    }

    private static double median(double[] values) {
        double[] copy = Arrays.copyOf(values, values.length);
        Arrays.sort(copy);
        int mid = copy.length / 2;
        if ((copy.length & 1) == 0) {
            return (copy[mid - 1] + copy[mid]) / 2.0d;
        }
        return copy[mid];
    }
}

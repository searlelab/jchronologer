package org.searlelab.jchronologer.smoke;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.junit.jupiter.api.Test;
import org.searlelab.jchronologer.api.AcceptedPrediction;
import org.searlelab.jchronologer.api.Chronologer;
import org.searlelab.jchronologer.api.PredictionResult;
import org.searlelab.jchronologer.impl.ChronologerFactory;
import org.searlelab.jchronologer.util.TsvTable;

class PeptidesOutSmokeTest {

    private static final String PEPTIDES_OUT_RESOURCE = "data/demo/peptides_out.txt";
    private static final double MAX_ABS_ERROR = 0.001d;

    @Test
    void predictionsMatchPythonChronologerWithinTolerance() throws Exception {
        long testStartNanos = System.nanoTime();
        TsvTable table = TsvTable.read(copyResourceToTemp(PEPTIDES_OUT_RESOURCE));
        long loadEndNanos = System.nanoTime();
        int peptideColumn = table.columnIndex("PeptideModSeq");
        int expectedColumn = table.columnIndex("Pred_HI");
        assertTrue(peptideColumn >= 0, "Missing PeptideModSeq column in peptides_out resource.");
        assertTrue(expectedColumn >= 0, "Missing Pred_HI column in peptides_out resource.");

        List<String> peptides = new ArrayList<>(table.getRows().size());
        List<Double> expectedPredictions = new ArrayList<>(table.getRows().size());
        for (String[] row : table.getRows()) {
            peptides.add(row[peptideColumn]);
            expectedPredictions.add(Double.parseDouble(row[expectedColumn]));
        }

        PredictionResult result;
        long predictStartNanos = System.nanoTime();
        try (Chronologer chronologer = ChronologerFactory.createDefault()) {
            result = chronologer.predict(peptides);
        }
        long predictEndNanos = System.nanoTime();

        assertEquals(
                peptides.size(),
                result.getAcceptedCount() + result.getRejectedCount(),
                "Accepted and rejected counts should cover every input row.");
        assertEquals(peptides.size(), result.getAcceptedCount(), "Expected all peptides_out rows to be accepted.");
        assertEquals(0, result.getRejectedCount(), "Expected no rejected rows in peptides_out.");

        Map<Integer, AcceptedPrediction> acceptedByRow = new HashMap<>();
        for (AcceptedPrediction accepted : result.getAccepted()) {
            acceptedByRow.put(accepted.getRowIndex(), accepted);
        }

        double maxAbsError = 0.0d;
        double sumAbsError = 0.0d;
        System.out.println("RowIndex\tPeptideModSeq\tExpectedPred_HI\tActualPred_HI\tAbsError");
        for (int rowIndex = 0; rowIndex < peptides.size(); rowIndex++) {
            AcceptedPrediction accepted = acceptedByRow.get(rowIndex);
            assertNotNull(accepted, "Missing accepted prediction for row " + rowIndex);

            double expected = expectedPredictions.get(rowIndex);
            double actual = accepted.getPredHi();
            double absError = Math.abs(actual - expected);
            sumAbsError += absError;
            if (absError > maxAbsError) {
                maxAbsError = absError;
            }
            System.out.printf(
                    "%d\t%s\t%.8f\t%.8f\t%.8f%n",
                    rowIndex,
                    peptides.get(rowIndex),
                    expected,
                    actual,
                    absError);
            assertTrue(
                    absError <= MAX_ABS_ERROR,
                    "Prediction mismatch at row " + rowIndex + " for peptide " + peptides.get(rowIndex)
                            + ": expected " + expected + ", actual " + actual + ", absError=" + absError);
        }

        long compareEndNanos = System.nanoTime();
        double loadMillis = (loadEndNanos - testStartNanos) / 1_000_000.0d;
        double predictMillis = (predictEndNanos - predictStartNanos) / 1_000_000.0d;
        double compareMillis = (compareEndNanos - predictEndNanos) / 1_000_000.0d;
        double totalMillis = (compareEndNanos - testStartNanos) / 1_000_000.0d;
        double meanAbsError = peptides.isEmpty() ? 0.0d : sumAbsError / peptides.size();
        System.out.printf(
                "TimingSummary\tload_ms=%.3f\tpredict_ms=%.3f\tcompare_ms=%.3f\ttotal_ms=%.3f%n",
                loadMillis,
                predictMillis,
                compareMillis,
                totalMillis);
        System.out.printf(
                "ErrorSummary\trows=%d\tmax_abs_error=%.8f\tmean_abs_error=%.8f\ttolerance=%.8f%n",
                peptides.size(),
                maxAbsError,
                meanAbsError,
                MAX_ABS_ERROR);
    }

    private static Path copyResourceToTemp(String resource) throws Exception {
        try (InputStream stream = Thread.currentThread().getContextClassLoader().getResourceAsStream(resource)) {
            assertTrue(stream != null, "Missing test resource: " + resource);
            Path tempFile = Files.createTempFile("jchronologer-peptides-out", ".tsv");
            Files.copy(stream, tempFile, StandardCopyOption.REPLACE_EXISTING);
            return tempFile;
        }
    }
}

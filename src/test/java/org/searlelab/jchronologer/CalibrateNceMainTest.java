package org.searlelab.jchronologer;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertIterableEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
import org.junit.jupiter.api.Test;

class CalibrateNceMainTest {

    @Test
    void cosineScoreIsOneForExactMatchedSpectra() {
        double score = CalibrateNceMain.scorePredictedIonCosine(
                new double[] {100.0, 200.0},
                new float[] {1.0f, 2.0f},
                new double[] {100.0, 200.0},
                new float[] {1.0f, 2.0f},
                10.0);

        assertEquals(1.0, score, 1e-9);
    }

    @Test
    void cosineScoreReturnsZeroWhenNoPredictedPeaksMatch() {
        double score = CalibrateNceMain.scorePredictedIonCosine(
                new double[] {100.0, 200.0},
                new float[] {1.0f, 2.0f},
                new double[] {300.0, 400.0},
                new float[] {5.0f, 6.0f},
                10.0);

        assertEquals(0.0, score, 1e-9);
    }

    @Test
    void cosineScoreUsesNearestObservedPeakWithinTolerance() {
        double score = CalibrateNceMain.scorePredictedIonCosine(
                new double[] {100.0, 200.0},
                new float[] {1.0f, 1.0f},
                new double[] {100.0008, 100.0002, 200.0001},
                new float[] {20.0f, 5.0f, 1.0f},
                10.0);

        double expected = 6.0 / (Math.sqrt(2.0) * Math.sqrt(26.0));
        assertEquals(expected, score, 1e-9);
    }

    @Test
    void selectionHonorsChargeThenMzCoverageThenIntensityFill() {
        List<CalibrateNceMain.CalibrationRow> rows = List.of(
                row("pepA", 2, 450.0, 50.0),
                row("pepB", 2, 520.0, 40.0),
                row("pepC", 3, 470.0, 30.0),
                row("pepD", 3, 620.0, 20.0),
                row("pepE", 4, 720.0, 10.0));

        List<CalibrateNceMain.CalibrationRow> selected = CalibrateNceMain.selectRows(rows, 4, 1);

        assertEquals(5, selected.size());
        assertTrue(selected.stream().anyMatch(row -> row.unimodPeptideSequence().equals("pepA")));
        assertTrue(selected.stream().anyMatch(row -> row.unimodPeptideSequence().equals("pepB")));
        assertTrue(selected.stream().anyMatch(row -> row.unimodPeptideSequence().equals("pepC")));
        assertTrue(selected.stream().anyMatch(row -> row.unimodPeptideSequence().equals("pepE")));
        assertTrue(selected.stream().anyMatch(row -> row.unimodPeptideSequence().equals("pepD")));
    }

    @Test
    void selectionAllowsSoftCapToExpandForCoverage() {
        List<CalibrateNceMain.CalibrationRow> rows = List.of(
                row("pepA", 2, 450.0, 50.0),
                row("pepB", 3, 470.0, 40.0),
                row("pepC", 4, 490.0, 30.0));

        List<CalibrateNceMain.CalibrationRow> selected = CalibrateNceMain.selectRows(rows, 2, 2);

        assertEquals(3, selected.size());
        assertIterableEquals(List.of("pepA", "pepB", "pepC"), selected.stream().map(
                CalibrateNceMain.CalibrationRow::unimodPeptideSequence).toList());
    }

    @Test
    void percentileSummaryReturnsExpectedBoxplotStatistics() {
        CalibrateNceMain.SummaryStats summary = CalibrateNceMain.summarize(List.of(10.0, 20.0, 30.0, 40.0, 50.0));

        assertEquals(5, summary.count());
        assertEquals(12.0, summary.p05(), 1e-9);
        assertEquals(20.0, summary.p25(), 1e-9);
        assertEquals(30.0, summary.p50(), 1e-9);
        assertEquals(40.0, summary.p75(), 1e-9);
        assertEquals(48.0, summary.p95(), 1e-9);
    }

    @Test
    void histogramBinsAreDeterministic() {
        List<CalibrateNceMain.HistogramBin> histogram = CalibrateNceMain.buildHistogram(List.of(32.1, 32.9, 33.0, 34.2));

        assertEquals(3, histogram.size());
        assertEquals(2, histogram.get(0).count());
        assertEquals(1, histogram.get(1).count());
        assertEquals(1, histogram.get(2).count());
    }

    @Test
    void searchUsesCachedEvaluationsAndPrefersLowerNceOnTies() {
        AtomicInteger evaluations = new AtomicInteger();
        List<String> phases = new ArrayList<>();

        CalibrateNceMain.SearchResult result = CalibrateNceMain.searchBestNce(
                33.0,
                nce -> new CalibrateNceMain.Evaluation(
                        nce,
                        nce == 33.0 || nce == 33.75 ? 1.0 : 0.5,
                        List.of(1.0),
                        List.of(row("pepA", 2, 450.0, 10.0)),
                        evaluations.incrementAndGet(),
                        nce,
                        java.util.Optional.empty(),
                        java.util.Optional.empty()),
                (phase, evaluation) -> phases.add(phase + ":" + evaluation.nce()));

        assertEquals(33.0, result.bestEvaluation().nce(), 1e-9);
        assertEquals(result.evaluations().size(), evaluations.get());
        assertTrue(phases.stream().anyMatch(phase -> phase.startsWith("seed:")));
        assertTrue(phases.stream().anyMatch(phase -> phase.startsWith("refine:")));
    }

    @Test
    void endianTolerantDecodingHandlesBigEndianExternalDlibs() {
        double[] decodedMasses = CalibrateNceMain.decodeMassArray(encodeDoubles(new double[] {213.087, 221.548, 294.181}, ByteOrder.BIG_ENDIAN));
        float[] decodedIntensities = CalibrateNceMain.decodeIntensityArray(
                encodeFloats(new float[] {360102.97f, 370290.2f, 378825.3f}, ByteOrder.BIG_ENDIAN));

        assertEquals(213.087, decodedMasses[0], 1e-6);
        assertEquals(221.548, decodedMasses[1], 1e-6);
        assertEquals(294.181, decodedMasses[2], 1e-6);
        assertEquals(360102.97f, decodedIntensities[0], 1e-2f);
        assertEquals(370290.2f, decodedIntensities[1], 1e-2f);
        assertEquals(378825.3f, decodedIntensities[2], 1e-2f);
    }

    private static CalibrateNceMain.CalibrationRow row(String peptide, int charge, double precursorMz, double intensity) {
        double binStart = CalibrateNceMain.mzBinStart(precursorMz, 100.0);
        return new CalibrateNceMain.CalibrationRow(
                peptide,
                (byte) charge,
                precursorMz,
                new double[] {100.0, 200.0},
                new float[] {(float) intensity, 1.0f},
                intensity,
                binStart,
                CalibrateNceMain.formatMzBinLabel(binStart, 100.0));
    }

    private static byte[] encodeDoubles(double[] values, ByteOrder byteOrder) {
        ByteBuffer buffer = ByteBuffer.allocate(values.length * Double.BYTES).order(byteOrder);
        for (double value : values) {
            buffer.putDouble(value);
        }
        return buffer.array();
    }

    private static byte[] encodeFloats(float[] values, ByteOrder byteOrder) {
        ByteBuffer buffer = ByteBuffer.allocate(values.length * Float.BYTES).order(byteOrder);
        for (float value : values) {
            buffer.putFloat(value);
        }
        return buffer.array();
    }
}

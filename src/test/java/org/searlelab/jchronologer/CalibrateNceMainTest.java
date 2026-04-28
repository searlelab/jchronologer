package org.searlelab.jchronologer;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertIterableEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
import org.junit.jupiter.api.Test;
import org.searlelab.jchronologer.dlib.DlibDatabase;
import org.searlelab.jchronologer.dlib.DlibEntryRecord;
import org.searlelab.jchronologer.dlib.DlibMetadata;

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

    @Test
    void parseArgsAcceptsAliasesAndValidatesValues() {
        CalibrateNceMain.CliArgs cliArgs = CalibrateNceMain.parseArgs(new String[] {
                "--start-nce", "31.5",
                "--ppm_tolerance", "15.0",
                "--max-entries", "50",
                "--group_target", "4",
                "--mz-bin-width", "75.0",
                "--batch_size", "8",
                "--inference-threads", "2",
                "--verbose",
                "input.dlib"
        });

        assertEquals(Path.of("input.dlib"), cliArgs.input());
        assertEquals(31.5, cliArgs.startNce(), 1e-9);
        assertEquals(15.0, cliArgs.ppmTolerance(), 1e-9);
        assertEquals(50, cliArgs.maxEntries());
        assertEquals(4, cliArgs.groupTarget());
        assertEquals(75.0, cliArgs.mzBinWidth(), 1e-9);
        assertEquals(8, cliArgs.batchSize());
        assertEquals(2, cliArgs.inferenceThreads());
        assertTrue(cliArgs.verbose());
        assertFalse(cliArgs.help());
    }

    @Test
    void parseArgsRejectsMissingInputAndOutOfRangeOptions() {
        IllegalArgumentException missingInput = assertThrows(
                IllegalArgumentException.class,
                () -> CalibrateNceMain.parseArgs(new String[] {"--verbose"}));
        assertEquals("Expected input DLIB path.", missingInput.getMessage());

        IllegalArgumentException badStartNce = assertThrows(
                IllegalArgumentException.class,
                () -> CalibrateNceMain.parseArgs(new String[] {"--start-nce", "61.0", "input.dlib"}));
        assertTrue(badStartNce.getMessage().contains("Start NCE must be a finite number in the range 10.0-60.0."));

        IllegalArgumentException badThreads = assertThrows(
                IllegalArgumentException.class,
                () -> CalibrateNceMain.parseArgs(new String[] {"--inference-threads", "0", "input.dlib"}));
        assertEquals("Inference threads must be positive.", badThreads.getMessage());
    }

    @Test
    void runPrintsHelpAndReturnsZero() {
        ByteArrayOutputStream stdout = new ByteArrayOutputStream();
        ByteArrayOutputStream stderr = new ByteArrayOutputStream();

        int exitCode = CalibrateNceMain.run(
                new String[] {"--help"},
                new PrintStream(stdout),
                new PrintStream(stderr));

        assertEquals(0, exitCode);
        assertTrue(stdout.toString().contains("Usage: calibrate-nce"));
        assertEquals("", stderr.toString());
    }

    @Test
    void runReportsArgumentErrors() {
        ByteArrayOutputStream stdout = new ByteArrayOutputStream();
        ByteArrayOutputStream stderr = new ByteArrayOutputStream();

        int exitCode = CalibrateNceMain.run(
                new String[] {"--unknown"},
                new PrintStream(stdout),
                new PrintStream(stderr));

        assertEquals(2, exitCode);
        assertTrue(stderr.toString().contains("Unknown option: --unknown"));
        assertTrue(stderr.toString().contains("Usage: calibrate-nce"));
    }

    @Test
    void readCalibrationRowsKeepsHighestIntensityDuplicateAndSkipsInvalidRows() throws Exception {
        Path dlib = Files.createTempFile("calibrate-rows", ".dlib");
        try (DlibDatabase database = new DlibDatabase(dlib, DlibMetadata.defaults())) {
            database.writeBatch(List.of(
                    dlibEntry("[]-VATVSLPR-[]", 2, 450.0, new double[] {100.0, 200.0}, new float[] {10.0f, 5.0f}),
                    dlibEntry("[]-VATVSLPR-[]", 2, 450.0, new double[] {100.0, 200.0}, new float[] {50.0f, 25.0f}),
                    dlibEntry("[]-TASEFDSAIAQDK-[]", 3, 520.0, new double[] {150.0}, new float[] {3.0f}),
                    dlibEntry("[]-BAD-[]", 2, 600.0, new double[] {}, new float[] {})),
                    List.of());
        }

        List<CalibrateNceMain.CalibrationRow> rows = CalibrateNceMain.readCalibrationRows(dlib, 100.0);

        assertEquals(2, rows.size());
        assertEquals("[]-VATVSLPR-[]", rows.get(0).unimodPeptideSequence());
        assertTrue(rows.get(0).summedObservedIntensity() > 0.0);
        assertEquals(2, rows.get(0).observedMasses().length);
        assertEquals(2, rows.get(0).observedIntensities().length);
        assertEquals("500.0-599.9999", rows.get(1).mzBinLabel());
    }

    @Test
    void distributionReportPrefersLowerNceOnPerRowTies() {
        List<CalibrateNceMain.CalibrationRow> rows = List.of(
                row("pepA", 2, 450.0, 10.0),
                row("pepB", 3, 550.0, 10.0));
        CalibrateNceMain.Evaluation lower = new CalibrateNceMain.Evaluation(
                32.0,
                0.9,
                List.of(0.7, 0.8),
                rows,
                1,
                32.0,
                java.util.Optional.empty(),
                java.util.Optional.of(34.0));
        CalibrateNceMain.Evaluation higher = new CalibrateNceMain.Evaluation(
                34.0,
                0.9,
                List.of(0.7, 0.8),
                rows,
                2,
                32.0,
                java.util.Optional.of(32.0),
                java.util.Optional.empty());

        CalibrateNceMain.DistributionReport report = CalibrateNceMain.buildDistributionReport(List.of(higher, lower), lower);

        assertEquals(32.0, report.overall().p50(), 1e-9);
        assertEquals(32.0, report.byCharge().get((byte) 2).p50(), 1e-9);
        assertEquals(32.0, report.byMzBin().get("400.0-499.9999").p50(), 1e-9);
    }

    @Test
    void percentileHandlesEdgePercentiles() {
        List<Double> values = List.of(10.0, 20.0, 30.0);

        assertEquals(10.0, CalibrateNceMain.percentile(values, 0.0), 1e-9);
        assertEquals(30.0, CalibrateNceMain.percentile(values, 1.0), 1e-9);
        assertEquals(20.0, CalibrateNceMain.percentile(values, 0.5), 1e-9);
    }

    @Test
    void integrationRunPrintsExpectedDistributionSections() throws Exception {
        Path dlib = Files.createTempFile("calibrate-main-int", ".dlib");
        try (DlibDatabase database = new DlibDatabase(dlib, DlibMetadata.defaults())) {
            database.writeBatch(List.of(
                    dlibEntry("[]-TASEFDSAIAQDK-[]", 2, 450.0, new double[] {100.0, 200.0}, new float[] {10.0f, 5.0f}),
                    dlibEntry("[]-VATVSLPR-[]", 2, 520.0, new double[] {120.0, 240.0}, new float[] {9.0f, 4.0f})),
                    List.of());
        }

        ByteArrayOutputStream stdout = new ByteArrayOutputStream();
        ByteArrayOutputStream stderr = new ByteArrayOutputStream();
        int exitCode = CalibrateNceMain.run(
                new String[] {dlib.toString(), "--start_nce", "33.0", "--max_entries", "2"},
                new PrintStream(stdout, true, StandardCharsets.UTF_8),
                new PrintStream(stderr, true, StandardCharsets.UTF_8));

        String output = stdout.toString(StandardCharsets.UTF_8);
        assertEquals(0, exitCode);
        assertEquals("", stderr.toString(StandardCharsets.UTF_8));
        assertTrue(output.contains("phase=seed"));
        assertTrue(output.contains("Global calibration summary"));
        assertTrue(output.contains("Histogram of best-sampled NCE values"));
    }

    @Test
    void validationHelpersRejectInvalidInputs() {
        assertThrows(IllegalArgumentException.class, () -> CalibrateNceMain.selectRows(List.of(row("pep", 2, 400.0, 1.0)), 0, 1));
        assertThrows(IllegalArgumentException.class, () -> CalibrateNceMain.selectRows(List.of(row("pep", 2, 400.0, 1.0)), 1, 0));
        assertThrows(IllegalArgumentException.class, () -> CalibrateNceMain.scorePredictedIonCosine(
                new double[] {100.0},
                new float[] {},
                new double[] {100.0},
                new float[] {1.0f},
                10.0));
        assertThrows(IllegalArgumentException.class, () -> CalibrateNceMain.summarize(List.of()));
        assertEquals(List.of(), CalibrateNceMain.buildHistogram(List.of()));
    }

    @Test
    void searchWrapperAndComparisonHelpersBehaveAsExpected() {
        CalibrateNceMain.SearchResult result = CalibrateNceMain.searchBestNce(
                10.0,
                nce -> new CalibrateNceMain.Evaluation(
                        nce,
                        -Math.abs(nce - 10.0),
                        List.of(1.0),
                        List.of(row("pepA", 2, 450.0, 10.0)),
                        1,
                        nce,
                        java.util.Optional.empty(),
                        java.util.Optional.empty()));

        assertEquals(10.0, result.bestEvaluation().nce(), 1e-9);

        CalibrateNceMain.Evaluation better = new CalibrateNceMain.Evaluation(
                30.0,
                0.9,
                List.of(1.0),
                List.of(row("pepA", 2, 450.0, 10.0)),
                1,
                30.0,
                java.util.Optional.empty(),
                java.util.Optional.empty());
        CalibrateNceMain.Evaluation worse = new CalibrateNceMain.Evaluation(
                31.0,
                0.8,
                List.of(1.0),
                List.of(row("pepA", 2, 450.0, 10.0)),
                2,
                30.0,
                java.util.Optional.empty(),
                java.util.Optional.empty());
        assertTrue(CalibrateNceMain.compareEvaluationsAscending(better, worse) > 0);
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

    private static DlibEntryRecord dlibEntry(
            String peptideModSeq,
            int charge,
            double precursorMz,
            double[] massArray,
            float[] intensityArray) {
        return new DlibEntryRecord(
                precursorMz,
                charge,
                peptideModSeq,
                peptideModSeq.replace("[]-", "").replace("-[]", "").replace("[UNIMOD:4]", ""),
                1,
                100.0f,
                0.0f,
                massArray,
                intensityArray,
                java.util.Optional.empty(),
                "calibrate-test");
    }
}

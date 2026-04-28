package org.searlelab.jchronologer;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import org.searlelab.jchronologer.api.ChronologerLibraryEntry;
import org.searlelab.jchronologer.api.ChronologerLibraryPredictor;
import org.searlelab.jchronologer.api.LibraryPredictionRequest;
import org.searlelab.jchronologer.api.PrecursorCondition;
import org.searlelab.jchronologer.dlib.DlibDatabase;
import org.searlelab.jchronologer.dlib.DlibEntryRecord;
import org.searlelab.jchronologer.dlib.DlibMetadata;
import org.searlelab.jchronologer.impl.ChronologerFactory;
import org.searlelab.jchronologer.preprocessing.PeptideSequenceConverter;
import org.junit.jupiter.api.Test;

class BatchCalibrateNceMainTest {

    @Test
    void parsesAcquiredNceFromFileName() {
        assertEquals(
                33.0,
                BatchCalibrateNceMain.parseAcquiredNce("20260403_no10_tne_hela_NCE33_OrbiDDA_01.dia.dlib"),
                1e-9);
        assertEquals(
                27.5,
                BatchCalibrateNceMain.parseAcquiredNce("sample_NCE27.5_AstDDA_01.dlib"),
                1e-9);
        assertNull(BatchCalibrateNceMain.parseAcquiredNce("sample_without_nce.dlib"));
    }

    @Test
    void findsOnlyMatchingDlibs() throws Exception {
        Path dir = Files.createTempDirectory("batch-calibrate-nce");
        Files.writeString(dir.resolve("a.dlib"), "");
        Files.writeString(dir.resolve("b.txt"), "");
        Files.writeString(dir.resolve("c.dia.dlib"), "");

        List<Path> matches = BatchCalibrateNceMain.findDlibs(dir, "*.dlib");

        assertEquals(2, matches.size());
        assertTrue(matches.get(0).getFileName().toString().endsWith(".dlib"));
        assertTrue(matches.get(1).getFileName().toString().endsWith(".dlib"));
    }

    @Test
    void parseArgsAcceptsLongAndUnderscoreFlags() throws Exception {
        Path dir = Files.createTempDirectory("batch-calibrate-args");

        BatchCalibrateNceMain.CliArgs cliArgs = BatchCalibrateNceMain.parseArgs(new String[] {
                "--glob", "*_NCE*.dlib",
                "--start-offset", "4.5",
                "--fallback_start_nce", "31.0",
                "--ppm-tolerance", "12.0",
                "--max_entries", "25",
                "--group-target", "3",
                "--mz_bin_width", "50.0",
                "--batch-size", "9",
                "--inference_threads", "2",
                "--verbose",
                dir.toString()
        });

        assertEquals(dir, cliArgs.inputDir());
        assertEquals("*_NCE*.dlib", cliArgs.glob());
        assertEquals(4.5, cliArgs.startOffset(), 1e-9);
        assertEquals(31.0, cliArgs.fallbackStartNce(), 1e-9);
        assertEquals(12.0, cliArgs.ppmTolerance(), 1e-9);
        assertEquals(25, cliArgs.maxEntries());
        assertEquals(3, cliArgs.groupTarget());
        assertEquals(50.0, cliArgs.mzBinWidth(), 1e-9);
        assertEquals(9, cliArgs.batchSize());
        assertEquals(2, cliArgs.inferenceThreads());
        assertTrue(cliArgs.verbose());
        assertFalse(cliArgs.help());
    }

    @Test
    void parseArgsRejectsInvalidInputs() {
        Path missingDir = Path.of("missing-batch-calibrate-dir");

        IllegalArgumentException missingDirError = org.junit.jupiter.api.Assertions.assertThrows(
                IllegalArgumentException.class,
                () -> BatchCalibrateNceMain.parseArgs(new String[] {missingDir.toString()}));
        assertTrue(missingDirError.getMessage().contains("Input directory must exist"));

        IllegalArgumentException badFallbackError = org.junit.jupiter.api.Assertions.assertThrows(
                IllegalArgumentException.class,
                () -> BatchCalibrateNceMain.parseArgs(new String[] {"--fallback-start-nce", "9.0"}));
        assertTrue(badFallbackError.getMessage().contains("Fallback start NCE must be in the range 10.0-60.0."));

        IllegalArgumentException badNumericError = org.junit.jupiter.api.Assertions.assertThrows(
                IllegalArgumentException.class,
                () -> BatchCalibrateNceMain.parseArgs(new String[] {"--batch-size", "0"}));
        assertTrue(badNumericError.getMessage().contains("Batch calibration numeric options must be positive."));
    }

    @Test
    void runPrintsHelpAndReturnsZero() {
        ByteArrayOutputStream stdout = new ByteArrayOutputStream();
        ByteArrayOutputStream stderr = new ByteArrayOutputStream();

        int exitCode = BatchCalibrateNceMain.run(
                new String[] {"--help"},
                new PrintStream(stdout),
                new PrintStream(stderr));

        assertEquals(0, exitCode);
        assertTrue(stdout.toString().contains("Usage: batch-calibrate-nce"));
        assertEquals("", stderr.toString());
    }

    @Test
    void runReportsMissingDlibsAsFailure() throws Exception {
        Path dir = Files.createTempDirectory("batch-calibrate-empty");
        ByteArrayOutputStream stdout = new ByteArrayOutputStream();
        ByteArrayOutputStream stderr = new ByteArrayOutputStream();

        int exitCode = BatchCalibrateNceMain.run(
                new String[] {dir.toString()},
                new PrintStream(stdout),
                new PrintStream(stderr));

        assertEquals(1, exitCode);
        assertTrue(stderr.toString().contains("No matching DLIB files found"));
    }

    @Test
    void runProcessesMultipleDlibsAndPrintsBatchSummary() throws Exception {
        Path dir = Files.createTempDirectory("batch-calibrate-run");
        writePredictedLibrary(dir.resolve("sample_NCE27_A.dlib"), List.of(
                new LibraryPredictionRequest("[]-TASEFDSAIAQDK-[]", List.of(new PrecursorCondition((byte) 2, 33.0))),
                new LibraryPredictionRequest("[]-VATVSLPR-[]", List.of(new PrecursorCondition((byte) 2, 33.0)))));
        writePredictedLibrary(dir.resolve("sample_without_nce.dlib"), List.of(
                new LibraryPredictionRequest("[]-HC[UNIMOD:4]VDPAVIAAIISR-[]", List.of(new PrecursorCondition((byte) 2, 33.0))),
                new LibraryPredictionRequest("[]-TASEFDSAIAQDK-[]", List.of(new PrecursorCondition((byte) 3, 33.0)))));

        ByteArrayOutputStream stdout = new ByteArrayOutputStream();
        ByteArrayOutputStream stderr = new ByteArrayOutputStream();

        int exitCode = BatchCalibrateNceMain.run(
                new String[] {"--max_entries", "2", "--group_target", "1", dir.toString()},
                new PrintStream(stdout, true, StandardCharsets.UTF_8),
                new PrintStream(stderr, true, StandardCharsets.UTF_8));

        String output = stdout.toString(StandardCharsets.UTF_8);
        assertEquals(0, exitCode);
        assertEquals("", stderr.toString(StandardCharsets.UTF_8));
        assertTrue(output.contains("== sample_NCE27_A.dlib acquired_nce=27.0000 start_nce=33.0000 =="));
        assertTrue(output.contains("== sample_without_nce.dlib acquired_nce=NA start_nce=33.0000 =="));
        assertTrue(output.contains("completed sample_NCE27_A.dlib"));
        assertTrue(output.contains("charge_median_best_nce"));
        assertTrue(output.contains("mz_bin_median_best_nce"));
        assertTrue(output.contains("Batch summary"));
        assertTrue(output.contains("EstimatedNCE"));
    }

    @Test
    void pathMatcherSupportsGlobWildcards() {
        BatchCalibrateNceMain.PathMatcher matcher = new BatchCalibrateNceMain.PathMatcher("*_NCE??.*.dlib");

        assertTrue(matcher.matches("sample_NCE33.a.dlib"));
        assertFalse(matcher.matches("sample_NCE3.a.dlib"));
        assertFalse(matcher.matches("sample_NCE333.a.dlib"));
    }

    @Test
    void parseArgsHelpersRejectMissingAndInvalidValues() {
        IllegalArgumentException missingValue = org.junit.jupiter.api.Assertions.assertThrows(
                IllegalArgumentException.class,
                () -> BatchCalibrateNceMain.parseArgs(new String[] {"--glob"}));
        assertTrue(missingValue.getMessage().contains("Missing value for option: --glob"));

        IllegalArgumentException badDouble = org.junit.jupiter.api.Assertions.assertThrows(
                IllegalArgumentException.class,
                () -> BatchCalibrateNceMain.parseArgs(new String[] {"--start-offset", "abc"}));
        assertTrue(badDouble.getMessage().contains("Invalid numeric value for option: --start-offset"));

        IllegalArgumentException badInt = org.junit.jupiter.api.Assertions.assertThrows(
                IllegalArgumentException.class,
                () -> BatchCalibrateNceMain.parseArgs(new String[] {"--max-entries", "abc"}));
        assertTrue(badInt.getMessage().contains("Invalid integer value for option: --max-entries"));
    }

    private static void writePredictedLibrary(Path output, List<LibraryPredictionRequest> requests) throws Exception {
        try (ChronologerLibraryPredictor predictor = ChronologerFactory.createLibraryPredictorDefault();
                DlibDatabase database = new DlibDatabase(output, DlibMetadata.defaults())) {
            List<ChronologerLibraryEntry> entries = predictor.predict(requests);
            List<DlibEntryRecord> dlibEntries = new ArrayList<>();
            for (ChronologerLibraryEntry entry : entries) {
                String peptideSeq = PeptideSequenceConverter.parseNormalizedUnimod(entry.getUnimodPeptideSequence()).getResidues();
                dlibEntries.add(new DlibEntryRecord(
                        entry.getPrecursorMz(),
                        entry.getPrecursorCharge(),
                        entry.getPeptideModSeq(),
                        peptideSeq,
                        1,
                        entry.getRetentionTimeInSeconds(),
                        0.0f,
                        entry.getMassArray(),
                        entry.getIntensityArray(),
                        entry.getCCS(),
                        "batch-test"));
            }
            database.writeBatch(dlibEntries, List.of());
        }
    }
}

package org.searlelab.jchronologer;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import org.junit.jupiter.api.Test;
import org.searlelab.jchronologer.api.ChronologerLibraryEntry;
import org.searlelab.jchronologer.api.ChronologerLibraryPredictor;
import org.searlelab.jchronologer.api.LibraryPredictionRequest;
import org.searlelab.jchronologer.api.PrecursorCondition;
import org.searlelab.jchronologer.dlib.DlibDatabase;
import org.searlelab.jchronologer.dlib.DlibEntryRecord;
import org.searlelab.jchronologer.dlib.DlibMetadata;
import org.searlelab.jchronologer.impl.ChronologerFactory;
import org.searlelab.jchronologer.preprocessing.PeptideSequenceConverter;

class CalibrateNceMainIntegrationTest {

    @Test
    void calibrateMainPrintsProgressAndDistributionSections() throws Exception {
        Path dlib = Files.createTempFile("calibrate-nce-main", ".dlib");
        writePredictedLibrary(dlib);

        RunResult result = runMain(dlib.toString(), "--start_nce", "33.0", "--max_entries", "2");

        assertEquals(0, result.code);
        assertTrue(result.stdout.contains("phase=seed"));
        assertTrue(result.stdout.contains("Global calibration summary"));
        assertTrue(result.stdout.contains("best_nce=33.0000"));
        assertTrue(result.stdout.contains("Overall best-sampled NCE distribution"));
        assertTrue(result.stdout.contains("Best-sampled NCE distribution by charge"));
        assertTrue(result.stdout.contains("Best-sampled NCE distribution by precursor m/z bin"));
        assertTrue(result.stdout.contains("Histogram of best-sampled NCE values"));
    }

    @Test
    void invalidDlibReturnsFailure() throws Exception {
        Path dlib = Files.createTempFile("calibrate-nce-empty", ".dlib");
        try (DlibDatabase ignored = new DlibDatabase(dlib, DlibMetadata.defaults())) {
            // empty DLIB
        }

        RunResult result = runMain(dlib.toString());

        assertEquals(1, result.code);
        assertTrue(result.stderr.contains("No valid calibration rows were found"));
    }

    private static void writePredictedLibrary(Path output) throws Exception {
        List<LibraryPredictionRequest> requests = List.of(
                new LibraryPredictionRequest(
                        "[]-TASEFDSAIAQDK-[]",
                        List.of(new PrecursorCondition((byte) 2, 33.0))),
                new LibraryPredictionRequest(
                        "[]-VATVSLPR-[]",
                        List.of(new PrecursorCondition((byte) 2, 33.0))));

        try (ChronologerLibraryPredictor predictor = ChronologerFactory.createLibraryPredictorDefault();
                DlibDatabase database = new DlibDatabase(output, DlibMetadata.defaults())) {
            List<ChronologerLibraryEntry> entries = predictor.predict(requests);
            database.writeBatch(entries.stream().map(CalibrateNceMainIntegrationTest::toDlibEntry).toList(), List.of());
        }
    }

    private static DlibEntryRecord toDlibEntry(ChronologerLibraryEntry entry) {
        String peptideSeq = PeptideSequenceConverter.parseNormalizedUnimod(entry.getUnimodPeptideSequence()).getResidues();
        return new DlibEntryRecord(
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
                "integration-test");
    }

    private static RunResult runMain(String... args) {
        ByteArrayOutputStream stdoutBytes = new ByteArrayOutputStream();
        ByteArrayOutputStream stderrBytes = new ByteArrayOutputStream();
        PrintStream stdout = new PrintStream(stdoutBytes, true, StandardCharsets.UTF_8);
        PrintStream stderr = new PrintStream(stderrBytes, true, StandardCharsets.UTF_8);
        int code = CalibrateNceMain.run(args, stdout, stderr);
        return new RunResult(
                code,
                stdoutBytes.toString(StandardCharsets.UTF_8),
                stderrBytes.toString(StandardCharsets.UTF_8));
    }

    private record RunResult(int code, String stdout, String stderr) {
    }
}

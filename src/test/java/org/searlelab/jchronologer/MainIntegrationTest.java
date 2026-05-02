package org.searlelab.jchronologer;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.List;
import org.junit.jupiter.api.Test;
import org.searlelab.jchronologer.api.ChronologerLibraryEntry;
import org.searlelab.jchronologer.api.ChronologerLibraryPredictor;
import org.searlelab.jchronologer.api.LibraryPredictionRequest;
import org.searlelab.jchronologer.dlib.DlibCodec;
import org.searlelab.jchronologer.impl.ChronologerFactory;

class MainIntegrationTest {

    @Test
    void plainTextInputWritesDlibWithEntriesMetadataAndProteinMappings() throws Exception {
        Path input = Files.createTempFile("jchronologer-main-plain", ".txt");
        Path fasta = Files.createTempFile("jchronologer-main-proteins", ".fasta");
        Path output = Files.createTempFile("jchronologer-main", ".dlib");
        Files.writeString(input, "VATVSLPR\nTASEFDSAIAQDK\nNOPEPTIDE\n", StandardCharsets.UTF_8);
        Files.writeString(
                fasta,
                ">sp|P1|protein1\nMAVATVSLPRGKTASEFDSAIAQDKL\n>sp|P2|protein2\nPEPTIDER\n",
                StandardCharsets.UTF_8);

        RunResult result = runMain(input.toString(), fasta.toString(), output.toString(), "--batch_size", "2");

        assertEquals(0, result.code);
        assertTrue(result.stderr.contains("Dropped 1 peptides for invalid input."));
        try (Connection connection = DriverManager.getConnection("jdbc:sqlite:" + output.toAbsolutePath());
                Statement statement = connection.createStatement()) {
            ResultSet entryCount = statement.executeQuery("select count(*) from entries");
            assertTrue(entryCount.next());
            assertTrue(entryCount.getInt(1) > 0);

            ResultSet metadataCount = statement.executeQuery("select count(*) from metadata");
            assertTrue(metadataCount.next());
            assertEquals(4, metadataCount.getInt(1));

            ResultSet proteinCount = statement.executeQuery("select count(*) from peptidetoprotein");
            assertTrue(proteinCount.next());
            assertTrue(proteinCount.getInt(1) >= 2);
        }
    }

    @Test
    void generatedLibraryRoundTripsThroughCalibrationReader() throws Exception {
        Path input = Files.createTempFile("jchronologer-main-roundtrip", ".txt");
        Path fasta = Files.createTempFile("jchronologer-main-roundtrip", ".fasta");
        Path output = Files.createTempFile("jchronologer-main-roundtrip", ".dlib");
        Files.writeString(input, "VATVSLPR\n", StandardCharsets.UTF_8);
        Files.writeString(fasta, ">P1\nMAVATVSLPRGK\n", StandardCharsets.UTF_8);

        RunResult result = runMain(input.toString(), fasta.toString(), output.toString());

        assertEquals(0, result.code);
        List<CalibrateNceMain.CalibrationRow> rows = CalibrateNceMain.readCalibrationRows(output, 100.0);
        assertFalse(rows.isEmpty());

        CalibrateNceMain.CalibrationRow vatvslpr = rows.stream()
                .filter(row -> "[]-VATVSLPR-[]".equals(row.unimodPeptideSequence()))
                .findFirst()
                .orElseThrow();
        assertEquals(2, vatvslpr.precursorCharge());
        assertTrue(vatvslpr.observedMasses().length > 0);
        assertEquals(vatvslpr.observedMasses().length, vatvslpr.observedIntensities().length);
        for (double mass : vatvslpr.observedMasses()) {
            assertTrue(Double.isFinite(mass));
            assertTrue(mass > 0.0);
            assertTrue(mass < 4000.0);
        }
        boolean hasPositiveIntensity = false;
        for (float intensity : vatvslpr.observedIntensities()) {
            assertTrue(Float.isFinite(intensity));
            assertTrue(intensity >= 0.0f);
            if (intensity > 0.0f) {
                hasPositiveIntensity = true;
            }
        }
        assertTrue(hasPositiveIntensity);
    }

    @Test
    void tsvInputSupportsCustomPeptideColumnAndPersistsScoreFromChargeProbability() throws Exception {
        Path input = Files.createTempFile("jchronologer-main-tsv", ".tsv");
        Path fasta = Files.createTempFile("jchronologer-main-fasta", ".fasta");
        Path output = Files.createTempFile("jchronologer-main-out", ".dlib");
        Files.writeString(input, "Seq\tId\nTASEFDSAIAQDK\t1\n", StandardCharsets.UTF_8);
        Files.writeString(fasta, ">P1\nMTASEFDSAIAQDKAA\n", StandardCharsets.UTF_8);

        RunResult result = runMain(
                input.toString(),
                fasta.toString(),
                output.toString(),
                "--peptide_column",
                "Seq");
        assertEquals(0, result.code);

        try (ChronologerLibraryPredictor predictor = ChronologerFactory.createLibraryPredictorDefault();
                Connection connection = DriverManager.getConnection("jdbc:sqlite:" + output.toAbsolutePath());
                Statement statement = connection.createStatement()) {
            List<ChronologerLibraryEntry> expected = predictor.predict(List.of(
                    new LibraryPredictionRequest("[]-TASEFDSAIAQDK-[]", 33.0, 0.01)));
            ChronologerLibraryEntry firstExpected = expected.get(0);

            ResultSet rs = statement.executeQuery(
                    "select Score, MassEncodedLength, MassArray, IntensityEncodedLength, IntensityArray "
                            + "from entries order by PrecursorCharge limit 1");
            assertTrue(rs.next());
            float observedScore = rs.getFloat(1);
            assertEquals(1.0f - firstExpected.getChargeProbability().orElseThrow(), observedScore, 1e-5f);

            byte[] massBlob = DlibCodec.decompress(rs.getBytes(3));
            byte[] intensityBlob = DlibCodec.decompress(rs.getBytes(5));
            assertEquals(rs.getInt(2), massBlob.length);
            assertEquals(rs.getInt(4), intensityBlob.length);
            assertArrayEquals(firstExpected.getMassArray(), DlibCodec.decodeDoubles(massBlob), 1e-9);
            assertArrayEquals(firstExpected.getIntensityArray(), DlibCodec.decodeFloats(intensityBlob), 1e-6f);
        }
    }

    @Test
    void entriesWithoutProteinMatchesAreDropped() throws Exception {
        Path input = Files.createTempFile("jchronologer-main-drop", ".txt");
        Path fasta = Files.createTempFile("jchronologer-main-drop", ".fasta");
        Path output = Files.createTempFile("jchronologer-main-drop", ".dlib");
        Files.writeString(input, "VATVSLPR\nTASEFDSAIAQDK\n", StandardCharsets.UTF_8);
        Files.writeString(fasta, ">P1\nVATVSLPR\n", StandardCharsets.UTF_8);

        RunResult result = runMain(input.toString(), fasta.toString(), output.toString());
        assertEquals(0, result.code);
        assertTrue(result.stderr.contains("Dropped 1 peptides and"));

        try (Connection connection = DriverManager.getConnection("jdbc:sqlite:" + output.toAbsolutePath());
                Statement statement = connection.createStatement()) {
            ResultSet rs = statement.executeQuery("select distinct PeptideSeq from entries order by PeptideSeq");
            List<String> peptideSeqs = new ArrayList<>();
            while (rs.next()) {
                peptideSeqs.add(rs.getString(1));
            }
            assertEquals(List.of("VATVSLPR"), peptideSeqs);
        }
    }

    @Test
    void noEntriesWrittenDeletesOutputFileAndReturnsFailure() throws Exception {
        Path input = Files.createTempFile("jchronologer-main-no-entries", ".txt");
        Path fasta = Files.createTempFile("jchronologer-main-no-entries", ".fasta");
        Path output = Files.createTempFile("jchronologer-main-no-entries", ".dlib");
        Files.writeString(input, "SHORT\n", StandardCharsets.UTF_8);
        Files.writeString(fasta, ">P1\nSHORT\n", StandardCharsets.UTF_8);

        RunResult result = runMain(input.toString(), fasta.toString(), output.toString());
        assertEquals(1, result.code);
        assertTrue(result.stderr.contains("No library entries were written."));
        assertEquals(false, Files.exists(output));
    }

    private static RunResult runMain(String... args) {
        ByteArrayOutputStream stdoutBytes = new ByteArrayOutputStream();
        ByteArrayOutputStream stderrBytes = new ByteArrayOutputStream();
        PrintStream stdout = new PrintStream(stdoutBytes, true, StandardCharsets.UTF_8);
        PrintStream stderr = new PrintStream(stderrBytes, true, StandardCharsets.UTF_8);
        int code = Main.run(args, stdout, stderr);
        return new RunResult(
                code,
                stdoutBytes.toString(StandardCharsets.UTF_8),
                stderrBytes.toString(StandardCharsets.UTF_8));
    }

    private static final class RunResult {
        private final int code;
        private final String stdout;
        private final String stderr;

        private RunResult(int code, String stdout, String stderr) {
            this.code = code;
            this.stdout = stdout;
            this.stderr = stderr;
        }
    }
}

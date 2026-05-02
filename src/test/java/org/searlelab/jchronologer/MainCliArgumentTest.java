package org.searlelab.jchronologer;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import org.junit.jupiter.api.Test;

class MainCliArgumentTest {

    @Test
    void helpLongOptionReturnsZeroAndPrintsUsage() {
        RunResult result = runMain("--help");
        assertEquals(0, result.code);
        assertTrue(result.stdout.contains("Usage: jchronologer"));
        assertTrue(result.stdout.contains("<proteins.fasta> <output.dlib>"));
        assertTrue(result.stdout.contains("--fast_mode"));
    }

    @Test
    void noArgsReturnsErrorCodeAndUsageOnStderr() {
        RunResult result = runMain();
        assertEquals(2, result.code);
        assertTrue(result.stderr.contains("Missing required input peptide file."));
    }

    @Test
    void missingFastaOrOutputReturnsError() throws Exception {
        Path input = Files.createTempFile("main-cli-input", ".txt");
        Files.writeString(input, "VATVSLPR\n", StandardCharsets.UTF_8);

        RunResult result = runMain(input.toString());
        assertEquals(2, result.code);
        assertTrue(result.stderr.contains("Expected input peptides, FASTA, and output DLIB."));
    }

    @Test
    void tooManyPositionalArgumentsReturnsError() throws Exception {
        Path input = Files.createTempFile("main-cli-input", ".txt");
        Path fasta = Files.createTempFile("main-cli-fasta", ".fasta");
        Files.writeString(input, "VATVSLPR\n", StandardCharsets.UTF_8);
        Files.writeString(fasta, ">P1\nVATVSLPR\n", StandardCharsets.UTF_8);

        RunResult result = runMain(input.toString(), fasta.toString(), "a.dlib", "b.dlib");
        assertEquals(2, result.code);
        assertTrue(result.stderr.contains("Too many positional arguments."));
    }

    @Test
    void invalidMinimumChargeProbabilityReturnsError() throws Exception {
        Path input = Files.createTempFile("main-cli-input", ".txt");
        Path fasta = Files.createTempFile("main-cli-fasta", ".fasta");
        Path output = Files.createTempFile("main-cli-output", ".dlib");
        Files.writeString(input, "VATVSLPR\n", StandardCharsets.UTF_8);
        Files.writeString(fasta, ">P1\nVATVSLPR\n", StandardCharsets.UTF_8);

        RunResult result = runMain(
                input.toString(),
                fasta.toString(),
                output.toString(),
                "--min_charge_probability",
                "1.1");
        assertEquals(2, result.code);
        assertTrue(result.stderr.contains("Minimum charge probability must be a finite number in the range 0.0-1.0."));
    }

    @Test
    void invalidNceReturnsError() throws Exception {
        Path input = Files.createTempFile("main-cli-input", ".txt");
        Path fasta = Files.createTempFile("main-cli-fasta", ".fasta");
        Path output = Files.createTempFile("main-cli-output", ".dlib");
        Files.writeString(input, "VATVSLPR\n", StandardCharsets.UTF_8);
        Files.writeString(fasta, ">P1\nVATVSLPR\n", StandardCharsets.UTF_8);

        RunResult result = runMain(
                input.toString(),
                fasta.toString(),
                output.toString(),
                "--nce",
                "9.0");
        assertEquals(2, result.code);
        assertTrue(result.stderr.contains("NCE must be a finite number in the range 10-60."));
    }

    @Test
    void invalidBatchSizeReturnsError() throws Exception {
        Path input = Files.createTempFile("main-cli-input", ".txt");
        Path fasta = Files.createTempFile("main-cli-fasta", ".fasta");
        Path output = Files.createTempFile("main-cli-output", ".dlib");
        Files.writeString(input, "VATVSLPR\n", StandardCharsets.UTF_8);
        Files.writeString(fasta, ">P1\nVATVSLPR\n", StandardCharsets.UTF_8);

        RunResult result = runMain(
                input.toString(),
                fasta.toString(),
                output.toString(),
                "--batch_size",
                "0");
        assertEquals(2, result.code);
        assertTrue(result.stderr.contains("Batch size must be positive."));
    }

    @Test
    void fastModeFlagIsAccepted() throws Exception {
        Path input = Files.createTempFile("main-cli-input", ".txt");
        Path fasta = Files.createTempFile("main-cli-fasta", ".fasta");
        Path output = Files.createTempFile("main-cli-output", ".dlib");
        Files.writeString(input, "VATVSLPR\n", StandardCharsets.UTF_8);
        Files.writeString(fasta, ">P1\nVATVSLPR\n", StandardCharsets.UTF_8);

        RunResult result = runMain(
                input.toString(),
                fasta.toString(),
                output.toString(),
                "--fast_mode");
        assertTrue(result.code == 0 || result.code == 1);
        assertTrue(!result.stderr.contains("Unknown option: --fast_mode"));
    }

    @Test
    void emptyInputFileReturnsPredictionFailure() throws Exception {
        Path input = Files.createTempFile("main-cli-empty", ".txt");
        Path fasta = Files.createTempFile("main-cli-fasta", ".fasta");
        Path output = Files.createTempFile("main-cli-output", ".dlib");
        Files.writeString(fasta, ">P1\nVATVSLPR\n", StandardCharsets.UTF_8);

        RunResult result = runMain(input.toString(), fasta.toString(), output.toString());
        assertEquals(1, result.code);
        assertTrue(result.stderr.contains("DLIB generation failed: Input file is empty"));
    }

    @Test
    void missingPeptideColumnInTsvReturnsPredictionFailure() throws Exception {
        Path input = Files.createTempFile("main-cli-input", ".tsv");
        Path fasta = Files.createTempFile("main-cli-fasta", ".fasta");
        Path output = Files.createTempFile("main-cli-output", ".dlib");
        Files.writeString(input, "Other\tId\nVATVSLPR\t1\n", StandardCharsets.UTF_8);
        Files.writeString(fasta, ">P1\nVATVSLPR\n", StandardCharsets.UTF_8);

        RunResult result = runMain(input.toString(), fasta.toString(), output.toString());
        assertEquals(1, result.code);
        assertTrue(result.stderr.contains("Input TSV is missing peptide column: PeptideModSeq"));
    }

    @Test
    void malformedFastaReturnsPredictionFailure() throws Exception {
        Path input = Files.createTempFile("main-cli-input", ".txt");
        Path fasta = Files.createTempFile("main-cli-fasta", ".fasta");
        Path output = Files.createTempFile("main-cli-output", ".dlib");
        Files.writeString(input, "VATVSLPR\n", StandardCharsets.UTF_8);
        Files.writeString(fasta, "VATVSLPR\n", StandardCharsets.UTF_8);

        RunResult result = runMain(input.toString(), fasta.toString(), output.toString());
        assertEquals(1, result.code);
        assertTrue(result.stderr.contains("DLIB generation failed: Malformed FASTA"));
    }

    @Test
    void missingInputFileReturnsError() {
        RunResult result = runMain("missing.txt", "missing.fasta", tempOutputPath().toString());
        assertEquals(2, result.code);
        assertTrue(result.stderr.contains("Input peptide file does not exist"));
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

    private static Path tempOutputPath() {
        try {
            return Files.createTempFile("main-cli-output", ".dlib");
        } catch (Exception e) {
            throw new IllegalStateException("Failed to create temp output path for test.", e);
        }
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

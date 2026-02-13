package org.searlelab.jchronologer.cli;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import org.junit.jupiter.api.Test;

class CliMainArgumentTest {

    @Test
    void noArgsPrintsUsageAndReturnsOne() {
        RunResult result = runCli();
        assertEquals(1, result.code);
        assertTrue(result.stdout.contains("Usage:"));
    }

    @Test
    void helpPrintsUsageAndReturnsZero() {
        RunResult result = runCli("--help");
        assertEquals(0, result.code);
        assertTrue(result.stdout.contains("jchronologer predict"));
    }

    @Test
    void unknownCommandReturnsTwo() {
        RunResult result = runCli("unknown");
        assertEquals(2, result.code);
        assertTrue(result.stderr.contains("Unknown command: unknown"));
    }

    @Test
    void missingInputOptionReturnsTwo() {
        RunResult result = runCli("predict", "--output", "out.tsv");
        assertEquals(2, result.code);
        assertTrue(result.stderr.contains("Missing required option: --input"));
    }

    @Test
    void missingOutputOptionReturnsTwo() {
        RunResult result = runCli("predict", "--input", "in.tsv");
        assertEquals(2, result.code);
        assertTrue(result.stderr.contains("Missing required option: --output"));
    }

    @Test
    void missingValueForInputOptionReturnsTwo() {
        RunResult result = runCli("predict", "--input");
        assertEquals(2, result.code);
        assertTrue(result.stderr.contains("Missing value for option: --input"));
    }

    @Test
    void unknownOptionReturnsTwo() {
        RunResult result = runCli("predict", "--input", "in.tsv", "--output", "out.tsv", "--bad");
        assertEquals(2, result.code);
        assertTrue(result.stderr.contains("Unknown option: --bad"));
    }

    @Test
    void invalidBatchSizeReturnsTwo() {
        RunResult result = runCli("predict", "--input", "in.tsv", "--output", "out.tsv", "--batch-size", "0");
        assertEquals(2, result.code);
        assertTrue(result.stderr.contains("Batch size must be positive."));
    }

    @Test
    void predictHelpOptionReturnsTwoWithGuidance() {
        RunResult result = runCli("predict", "--help");
        assertEquals(2, result.code);
        assertTrue(result.stderr.contains("Use `jchronologer --help` to print usage."));
    }

    @Test
    void predictionFailurePathReturnsOne() throws IOException {
        Path input = Files.createTempFile("jchronologer-cli-no-col", ".tsv");
        Path output = Files.createTempFile("jchronologer-cli-out", ".tsv");
        Files.writeString(input, "Other\nVATVSLPR\n", StandardCharsets.UTF_8);

        RunResult result = runCli("predict", "--input", input.toString(), "--output", output.toString());
        assertEquals(1, result.code);
        assertTrue(result.stderr.contains("Prediction failed: Input TSV is missing peptide column: PeptideModSeq"));
    }

    private static RunResult runCli(String... args) {
        PrintStream originalOut = System.out;
        PrintStream originalErr = System.err;
        ByteArrayOutputStream stdoutBytes = new ByteArrayOutputStream();
        ByteArrayOutputStream stderrBytes = new ByteArrayOutputStream();
        try {
            System.setOut(new PrintStream(stdoutBytes, true, StandardCharsets.UTF_8));
            System.setErr(new PrintStream(stderrBytes, true, StandardCharsets.UTF_8));
            int code = Main.run(args);
            return new RunResult(
                    code,
                    stdoutBytes.toString(StandardCharsets.UTF_8),
                    stderrBytes.toString(StandardCharsets.UTF_8));
        } finally {
            System.setOut(originalOut);
            System.setErr(originalErr);
        }
    }

    private record RunResult(int code, String stdout, String stderr) {
    }
}

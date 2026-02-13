package org.searlelab.jchronologer;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import org.junit.jupiter.api.Test;
import org.searlelab.jchronologer.util.TsvTable;

class MainCliArgumentTest {

    @Test
    void helpLongOptionReturnsZeroAndPrintsUsage() {
        RunResult result = runMain("--help");
        assertEquals(0, result.code);
        assertTrue(result.stdout.contains("Usage: jchronologer"));
    }

    @Test
    void helpShortOptionReturnsZeroAndPrintsUsage() {
        RunResult result = runMain("-h");
        assertEquals(0, result.code);
        assertTrue(result.stdout.contains("Predict RTs using Chronologer."));
    }

    @Test
    void noArgsReturnsErrorCodeAndUsageOnStderr() {
        RunResult result = runMain();
        assertEquals(2, result.code);
        assertTrue(result.stderr.contains("Missing required input file."));
        assertTrue(result.stderr.contains("Usage: jchronologer"));
    }

    @Test
    void unknownOptionReturnsErrorCodeAndUsageOnStderr() {
        RunResult result = runMain("--unknown");
        assertEquals(2, result.code);
        assertTrue(result.stderr.contains("Unknown option: --unknown"));
    }

    @Test
    void tooManyPositionalArgumentsReturnsError() {
        RunResult result = runMain("a.tsv", "b.tsv", "c.tsv");
        assertEquals(2, result.code);
        assertTrue(result.stderr.contains("Too many positional arguments."));
    }

    @Test
    void missingBatchSizeValueReturnsError() {
        RunResult result = runMain("input.tsv", "--batch_size");
        assertEquals(2, result.code);
        assertTrue(result.stderr.contains("Missing value for option: --batch_size"));
    }

    @Test
    void missingPeptideColumnValueReturnsError() {
        RunResult result = runMain("input.tsv", "--peptide_column");
        assertEquals(2, result.code);
        assertTrue(result.stderr.contains("Missing value for option: --peptide_column"));
    }

    @Test
    void invalidBatchSizeReturnsError() {
        RunResult result = runMain("input.tsv", "--batch_size", "0");
        assertEquals(2, result.code);
        assertTrue(result.stderr.contains("Batch size must be positive."));
    }

    @Test
    void emptyInputFileReturnsPredictionFailure() throws IOException {
        Path input = Files.createTempFile("jchronologer-main-empty", ".txt");
        RunResult result = runMain(input.toString());
        assertEquals(1, result.code);
        assertTrue(result.stderr.contains("Prediction failed: Input file is empty"));
    }

    @Test
    void missingPeptideColumnInTsvReturnsPredictionFailure() throws IOException {
        Path input = Files.createTempFile("jchronologer-main-no-col", ".tsv");
        Files.writeString(input, "Other\tId\nVATVSLPR\t1\n", StandardCharsets.UTF_8);
        Path output = Files.createTempFile("jchronologer-main-out", ".tsv");
        RunResult result = runMain(input.toString(), output.toString());
        assertEquals(1, result.code);
        assertTrue(result.stderr.contains("Input TSV is missing peptide column: PeptideModSeq"));
    }

    @Test
    void customPeptideColumnInTsvIsSupported() throws IOException {
        Path input = Files.createTempFile("jchronologer-main-custom-col", ".tsv");
        Files.writeString(
                input,
                "Seq\tId\nVATVSLPR\t1\n[42.010565]ACDEFGHIK\t2\n",
                StandardCharsets.UTF_8);
        Path output = Files.createTempFile("jchronologer-main-custom-col-out", ".tsv");

        RunResult result = runMain(input.toString(), output.toString(), "--peptide_column", "Seq");
        assertEquals(0, result.code);

        TsvTable table = TsvTable.read(output);
        assertEquals(3, table.getHeaders().size());
        assertEquals("Pred_HI", table.getHeaders().get(2));
        assertEquals(1, table.getRows().size());
        assertEquals("VATVSLPR", table.getRows().get(0)[0]);
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

    private record RunResult(int code, String stdout, String stderr) {
    }
}

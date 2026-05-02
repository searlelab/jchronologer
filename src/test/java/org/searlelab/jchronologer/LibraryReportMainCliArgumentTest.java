package org.searlelab.jchronologer;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import org.junit.jupiter.api.Test;

class LibraryReportMainCliArgumentTest {

    @Test
    void helpLongOptionReturnsZeroAndPrintsUsage() {
        RunResult result = runMain("--help");
        assertEquals(0, result.code);
        assertTrue(result.stdout.contains("Usage: library-report"));
    }

    @Test
    void helpShortOptionReturnsZeroAndPrintsUsage() {
        RunResult result = runMain("-h");
        assertEquals(0, result.code);
        assertTrue(result.stdout.contains("Generate a quick console report"));
    }

    @Test
    void unknownOptionReturnsErrorCodeAndUsageOnStderr() {
        RunResult result = runMain("--unknown");
        assertEquals(2, result.code);
        assertTrue(result.stderr.contains("Unknown option: --unknown"));
        assertTrue(result.stderr.contains("Usage: library-report"));
    }

    @Test
    void missingNceValueReturnsErrorCode() {
        RunResult result = runMain("--nce");
        assertEquals(2, result.code);
        assertTrue(result.stderr.contains("Missing value for option: --nce"));
    }

    @Test
    void outOfRangeNceReturnsErrorCode() {
        RunResult result = runMain("--nce", "9.99");
        assertEquals(2, result.code);
        assertTrue(result.stderr.contains("NCE must be a finite number in the range"));
    }

    @Test
    void outOfRangeMinimumChargeProbabilityReturnsErrorCode() {
        RunResult result = runMain("--minimum_charge_probability", "1.1");
        assertEquals(2, result.code);
        assertTrue(result.stderr.contains("Minimum charge probability must be a finite number in the range 0.0-1.0."));
    }

    @Test
    void invalidMinimumChargeProbabilityValueReturnsErrorCode() {
        RunResult result = runMain("--minimum_charge_probability", "abc");
        assertEquals(2, result.code);
        assertTrue(result.stderr.contains("Invalid numeric value for option: --minimum_charge_probability"));
    }

    @Test
    void strictUnimodValidationProducesErrorRow() {
        RunResult result = runMain("TASEFDSAIAQDK");
        assertEquals(0, result.code);
        assertTrue(result.stdout.contains("Peptide sequence (formatted with mods)"));
        assertTrue(result.stdout.contains("TASEFDSAIAQDK"));
        assertTrue(result.stdout.contains("ERROR"));
    }

    @Test
    void ionFormattingConvertsRawCartographerLabels() {
        assertEquals("y4", LibraryReportMain.formatIonTypeForDisplay("1+y4"));
        assertEquals("y4++", LibraryReportMain.formatIonTypeForDisplay("2+y4"));
        assertEquals("b6+++", LibraryReportMain.formatIonTypeForDisplay("3+b6"));
        assertEquals("precursor", LibraryReportMain.formatIonTypeForDisplay("precursor"));
    }

    private static RunResult runMain(String... args) {
        ByteArrayOutputStream stdoutBytes = new ByteArrayOutputStream();
        ByteArrayOutputStream stderrBytes = new ByteArrayOutputStream();
        PrintStream stdout = new PrintStream(stdoutBytes, true, StandardCharsets.UTF_8);
        PrintStream stderr = new PrintStream(stderrBytes, true, StandardCharsets.UTF_8);
        int code = LibraryReportMain.run(args, stdout, stderr);
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

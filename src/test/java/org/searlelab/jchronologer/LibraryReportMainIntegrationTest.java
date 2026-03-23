package org.searlelab.jchronologer;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.regex.Pattern;
import org.junit.jupiter.api.Test;
import org.searlelab.jchronologer.api.ChronologerLibraryOptions;
import org.searlelab.jchronologer.inference.ElectricianBatchPredictor;
import org.searlelab.jchronologer.preprocessing.ChronologerPreprocessor;
import org.searlelab.jchronologer.preprocessing.PeptideSequenceConverter;
import org.searlelab.jchronologer.preprocessing.PreprocessingMetadataLoader;
import org.searlelab.jchronologer.preprocessing.PreprocessingOutcome;

class LibraryReportMainIntegrationTest {

    private static final List<String> DEFAULT_PEPTIDES = List.of(
            "[]-HC[UNIMOD:4]VDPAVIAAIISR-[]",
            "[]-TLLISSLSPALPAEHLEDR-[]",
            "[]-TPIGSFLGSLS-[]",
            "[]-ANAEKTSGSNVKIVKVKKE-[]",
            "[UNIMOD:737]-K[UNIMOD:737]PGLAITFAK[UNIMOD:737]-[]");

    private static final Pattern FORMATTED_ION_PATTERN = Pattern.compile("[yb]\\d+(?:\\+{2,3})?");

    @Test
    void noArgsUsesDefaultPeptides() {
        RunResult result = runMain();
        assertEquals(0, result.code);

        List<String[]> rows = parseDataRows(result.stdout);
        assertFalse(rows.isEmpty());

        Set<String> peptidesInOutput = new HashSet<>();
        for (String[] row : rows) {
            peptidesInOutput.add(row[0]);
        }

        for (String peptide : DEFAULT_PEPTIDES) {
            assertTrue(peptidesInOutput.contains(peptide), "Missing peptide in output: " + peptide);
        }
    }

    @Test
    void emitsOneRowPerChargeStateWhenThresholdIsZero() {
        String peptide = "[]-HC[UNIMOD:4]VDPAVIAAIISR-[]";
        RunResult result = runMain("--minimum_charge_probability", "0.0", peptide);
        assertEquals(0, result.code);

        List<String[]> rows = parseDataRows(result.stdout);
        assertEquals(6, rows.size());

        Set<Integer> charges = new HashSet<>();
        for (String[] row : rows) {
            assertEquals(peptide, row[0]);
            charges.add(Integer.parseInt(row[2]));
        }
        assertEquals(Set.of(1, 2, 3, 4, 5, 6), charges);
    }

    @Test
    void topIonFormattingUsesSuffixStyle() {
        String peptide = "[]-HC[UNIMOD:4]VDPAVIAAIISR-[]";
        RunResult result = runMain("--minimum_charge_probability", "0.0", peptide);
        assertEquals(0, result.code);

        List<String[]> rows = parseDataRows(result.stdout);
        assertFalse(rows.isEmpty());

        for (String[] row : rows) {
            String topIons = row[4];
            assertFalse(topIons.contains("1+"));
            assertFalse(topIons.contains("2+"));
            assertFalse(topIons.contains("3+"));

            String[] tokens = topIons.split(", ");
            for (String token : tokens) {
                assertTrue(FORMATTED_ION_PATTERN.matcher(token).matches(), "Unexpected ion token: " + token);
            }
        }
    }

    @Test
    void emitsNoChargeRowWhenNoChargePassesThreshold() {
        String peptide = "[]-TASEFDSAIAQDK-[]";
        double strictThreshold = probabilityJustAboveMaximumChargeProbability(peptide);

        RunResult result = runMain(
                "--minimum_charge_probability",
                Double.toString(strictThreshold),
                peptide);
        assertEquals(0, result.code);

        List<String[]> rows = parseDataRows(result.stdout);
        assertEquals(1, rows.size());
        assertEquals(peptide, rows.get(0)[0]);
        assertEquals("NO_CHARGE", rows.get(0)[1]);
        assertEquals("NO_CHARGE", rows.get(0)[2]);
        assertEquals("NO_CHARGE", rows.get(0)[3]);
        assertEquals("NO_CHARGE", rows.get(0)[4]);
    }

    @Test
    void continuesAfterErrorRow() {
        String invalid = "TASEFDSAIAQDK";
        String valid = "[]-TPIGSFLGSLS-[]";

        RunResult result = runMain(invalid, valid);
        assertEquals(0, result.code);

        List<String[]> rows = parseDataRows(result.stdout);
        assertFalse(rows.isEmpty());

        boolean foundError = false;
        boolean foundValid = false;
        for (String[] row : rows) {
            if (row[0].equals(invalid) && row[1].equals("ERROR")) {
                foundError = true;
            }
            if (row[0].equals(valid) && !row[2].equals("ERROR") && !row[2].equals("NO_CHARGE")) {
                foundValid = true;
            }
        }

        assertTrue(foundError, "Expected an ERROR row for invalid peptide.");
        assertTrue(foundValid, "Expected at least one valid prediction row.");
    }

    private static double probabilityJustAboveMaximumChargeProbability(String unimodPeptide) {
        ChronologerPreprocessor preprocessor = new ChronologerPreprocessor(
                PreprocessingMetadataLoader.loadFromClasspath(
                        ChronologerLibraryOptions.DEFAULT_ELECTRICIAN_PREPROCESSING_RESOURCE));

        String massEncoded = PeptideSequenceConverter.unimodToMassEncoded(unimodPeptide);
        PreprocessingOutcome outcome = preprocessor.preprocess(massEncoded);
        assertTrue(outcome.isAccepted(), "Failed to preprocess test peptide: " + unimodPeptide);

        try (ElectricianBatchPredictor predictor = new ElectricianBatchPredictor(
                ChronologerLibraryOptions.DEFAULT_ELECTRICIAN_MODEL_RESOURCE)) {
            float[][] distributions = predictor.predict(new long[][] {outcome.getTokenArray()});
            assertEquals(1, distributions.length);
            assertEquals(6, distributions[0].length);

            float[] normalized = normalizeDistribution(distributions[0]);
            float max = 0.0f;
            for (float value : normalized) {
                if (value > max) {
                    max = value;
                }
            }

            double strictThreshold = Math.nextUp((double) max);
            assertTrue(
                    strictThreshold <= 1.0 && strictThreshold > max,
                    "Could not construct strict threshold above max charge probability.");
            return strictThreshold;
        }
    }

    private static float[] normalizeDistribution(float[] rawDistribution) {
        float[] normalized = new float[rawDistribution.length];
        double sum = 0.0;
        for (int i = 0; i < rawDistribution.length; i++) {
            float value = rawDistribution[i];
            if (Float.isFinite(value) && value > 0.0f) {
                normalized[i] = value;
                sum += value;
            }
        }
        if (sum <= 0.0) {
            return normalized;
        }
        for (int i = 0; i < normalized.length; i++) {
            normalized[i] = (float) (normalized[i] / sum);
        }
        return normalized;
    }

    private static List<String[]> parseDataRows(String output) {
        String[] lines = output.split("\\R");
        assertTrue(lines.length >= 2, "Expected at least header and separator lines.");

        List<String[]> rows = new ArrayList<>();
        for (int i = 2; i < lines.length; i++) {
            String line = lines[i];
            if (line.isBlank()) {
                continue;
            }
            String[] parts = line.split("\\|", -1);
            assertEquals(5, parts.length, "Unexpected number of table columns in line: " + line);
            for (int j = 0; j < parts.length; j++) {
                parts[j] = parts[j].trim();
            }
            rows.add(parts);
        }
        return rows;
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

    private record RunResult(int code, String stdout, String stderr) {
    }
}

package org.searlelab.jchronologer;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.ByteArrayOutputStream;
import java.io.InputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import org.junit.jupiter.api.Test;
import org.searlelab.jchronologer.api.AcceptedPrediction;
import org.searlelab.jchronologer.api.Chronologer;
import org.searlelab.jchronologer.api.PredictionResult;
import org.searlelab.jchronologer.impl.ChronologerFactory;
import org.searlelab.jchronologer.preprocessing.PeptideSequenceConverter;
import org.searlelab.jchronologer.util.TsvTable;

class UnimodPeptidesOutSmokeTest {

    private static final String UNIMOD_PEPTIDES_OUT_RESOURCE = "data/demo/unimod_peptides_out.txt";
    private static final double UNIMOD_MASS_MATCH_EPSILON = 1e-5;

    @Test
    void mainHandlesUnimodInputConsistentlyWithMassCodedChronologerPath() throws Exception {
        Path input = copyResourceToTemp(UNIMOD_PEPTIDES_OUT_RESOURCE);
        TsvTable inputTable = TsvTable.read(input);
        int peptideColumn = inputTable.columnIndex("PeptideModSeq");
        assertTrue(peptideColumn >= 0, "Missing PeptideModSeq column in unimod_peptides_out resource.");

        List<String> converted = new ArrayList<>(inputTable.getRows().size());
        for (String[] row : inputTable.getRows()) {
            converted.add(toChronologerInputSequence(row[peptideColumn]));
        }

        PredictionResult baseline;
        try (Chronologer chronologer = ChronologerFactory.createDefault()) {
            baseline = chronologer.predict(converted);
        }
        assertEquals(inputTable.getRows().size(), baseline.getAcceptedCount(), "Expected all unimod rows to be accepted.");
        assertEquals(0, baseline.getRejectedCount(), "Expected no rejected unimod rows.");

        Path output = Files.createTempFile("jchronologer-unimod-main-out", ".tsv");
        ByteArrayOutputStream stdoutBytes = new ByteArrayOutputStream();
        ByteArrayOutputStream stderrBytes = new ByteArrayOutputStream();
        PrintStream stdout = new PrintStream(stdoutBytes, true, StandardCharsets.UTF_8);
        PrintStream stderr = new PrintStream(stderrBytes, true, StandardCharsets.UTF_8);
        int code = Main.run(new String[] {input.toString(), output.toString()}, stdout, stderr);
        assertEquals(0, code, "CLI failed: " + stderrBytes.toString(StandardCharsets.UTF_8));

        TsvTable outputTable = TsvTable.read(output);
        int outputPeptideColumn = outputTable.columnIndex("PeptideModSeq");
        int outputPredictedColumn = outputTable.getHeaders().size() - 1;
        assertTrue(outputPeptideColumn >= 0, "Missing PeptideModSeq column in output.");
        assertEquals("Pred_HI", outputTable.getHeaders().get(outputPredictedColumn));
        assertEquals(
                baseline.getAcceptedCount(),
                outputTable.getRows().size(),
                "CLI output should contain exactly all accepted rows.");

        Map<Integer, AcceptedPrediction> acceptedByRow = baseline.getAcceptedByRowIndex();
        int outputRowIndex = 0;
        for (int rowIndex = 0; rowIndex < inputTable.getRows().size(); rowIndex++) {
            AcceptedPrediction accepted = acceptedByRow.get(rowIndex);
            assertNotNull(accepted, "Missing accepted prediction for row " + rowIndex);
            String[] sourceRow = inputTable.getRows().get(rowIndex);
            String[] outputRow = outputTable.getRows().get(outputRowIndex++);

            assertEquals(sourceRow[peptideColumn], outputRow[outputPeptideColumn]);
            assertEquals(accepted.getPredHi(), Double.parseDouble(outputRow[outputPredictedColumn]), 1e-6);
        }
        assertEquals(outputTable.getRows().size(), outputRowIndex);
    }

    private static String toChronologerInputSequence(String peptide) {
        try {
            String unimod = PeptideSequenceConverter.normalizeToUnimod(peptide, UNIMOD_MASS_MATCH_EPSILON);
            return PeptideSequenceConverter.unimodToMassEncoded(unimod);
        } catch (RuntimeException e) {
            return peptide;
        }
    }

    private static Path copyResourceToTemp(String resource) throws Exception {
        try (InputStream stream = Thread.currentThread().getContextClassLoader().getResourceAsStream(resource)) {
            assertNotNull(stream, "Missing test resource: " + resource);
            Path tempFile = Files.createTempFile("jchronologer-unimod-peptides-out", ".tsv");
            Files.copy(stream, tempFile, StandardCopyOption.REPLACE_EXISTING);
            return tempFile;
        }
    }
}

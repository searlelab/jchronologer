package org.searlelab.jchronologer.smoke;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.junit.jupiter.api.Test;
import org.searlelab.jchronologer.api.AcceptedPrediction;
import org.searlelab.jchronologer.api.Chronologer;
import org.searlelab.jchronologer.api.PredictionResult;
import org.searlelab.jchronologer.api.RejectedPrediction;
import org.searlelab.jchronologer.impl.ChronologerFactory;
import org.searlelab.jchronologer.util.TsvTable;

class DemoPeptidesSmokeTest {

    private static final String DEMO_PEPTIDES_RESOURCE = "data/demo/demo_peptides.txt";
    private static final String DEMO_RT_RESOURCE = "data/demo/demo_rt.txt";
    private static final double MAX_ABS_ERROR = 0.001d;

    @Test
    void predictsDemoPeptidesAndPrintsSequenceAndPrediction() throws Exception {
        TsvTable table = TsvTable.read(copyResourceToTemp(DEMO_PEPTIDES_RESOURCE));
        int peptideColumn = table.columnIndex("PeptideModSeq");
        assertTrue(peptideColumn >= 0, "Missing PeptideModSeq column in demo file.");

        List<String> peptides = new ArrayList<>();
        for (String[] row : table.getRows()) {
            peptides.add(row[peptideColumn]);
        }

        Map<String, Double> expectedByPeptide = loadExpectedPredictions();
        assertEquals(peptides.size(), expectedByPeptide.size(), "Input peptide count should match expected table.");

        PredictionResult result;
        try (Chronologer chronologer = ChronologerFactory.createDefault()) {
            result = chronologer.predict(peptides);
        }

        Map<Integer, AcceptedPrediction> acceptedByRow = result.getAcceptedByRowIndex();
        Map<Integer, RejectedPrediction> rejectedByRow = result.getRejectedByRowIndex();

        System.out.println("PeptideModSeq\tPred_HI");
        int compared = 0;
        for (int i = 0; i < peptides.size(); i++) {
            String peptide = peptides.get(i);
            AcceptedPrediction accepted = acceptedByRow.get(i);
            if (accepted != null) {
                Double expected = expectedByPeptide.get(peptide);
                float calculated=accepted.getPredHi();
				double delta=expected-calculated;
                System.out.printf("%s\t%.8f\t%.8f%n", peptide, calculated, delta);
                assertNotNull(expected, "Missing expected value for peptide: " + peptide);
                double absError = Math.abs(calculated - expected);
                assertTrue(
                        absError <= MAX_ABS_ERROR,
                        "Prediction mismatch for peptide " + peptide + ": expected " + expected
                                + ", actual " + calculated + ", absError=" + absError);
                compared++;
            } else {
                RejectedPrediction rejected = rejectedByRow.get(i);
                String reason = rejected == null ? "UNKNOWN_REJECTION" : rejected.getRejectionReason().name();
                System.out.printf("%s\tREJECTED(%s)%n", peptide, reason);
            }
        }

        assertEquals(peptides.size(), result.getAcceptedCount() + result.getRejectedCount());
        assertEquals(peptides.size(), result.getAcceptedCount(), "Expected all demo peptides to be accepted.");
        assertEquals(peptides.size(), compared, "Expected values should be compared for every peptide.");
    }

    private static Path copyResourceToTemp(String resource) throws Exception {
        try (InputStream stream = Thread.currentThread().getContextClassLoader().getResourceAsStream(resource)) {
            assertTrue(stream != null, "Missing test resource: " + resource);
            Path tempFile = Files.createTempFile("jchronologer-demo-peptides", ".tsv");
            Files.copy(stream, tempFile, StandardCopyOption.REPLACE_EXISTING);
            return tempFile;
        }
    }

    private static Map<String, Double> loadExpectedPredictions() throws Exception {
        TsvTable rtTable = TsvTable.read(copyResourceToTemp(DEMO_RT_RESOURCE));
        int peptideColumn = rtTable.columnIndex("PeptideModSeq");
        int predColumn = rtTable.columnIndex("Pred_HI");
        assertTrue(peptideColumn >= 0, "Missing PeptideModSeq column in demo_rt resource.");
        assertTrue(predColumn >= 0, "Missing Pred_HI column in demo_rt resource.");

        Map<String, Double> expectedByPeptide = new HashMap<>();
        for (String[] row : rtTable.getRows()) {
            String peptide = row[peptideColumn];
            double pred = Double.parseDouble(row[predColumn]);
            expectedByPeptide.put(peptide, pred);
        }
        return expectedByPeptide;
    }
}

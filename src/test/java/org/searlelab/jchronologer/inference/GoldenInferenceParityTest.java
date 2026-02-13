package org.searlelab.jchronologer.inference;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import org.junit.jupiter.api.Test;
import org.searlelab.jchronologer.api.AcceptedPrediction;
import org.searlelab.jchronologer.api.Chronologer;
import org.searlelab.jchronologer.api.PredictionResult;
import org.searlelab.jchronologer.api.RejectedPrediction;
import org.searlelab.jchronologer.impl.ChronologerFactory;
import org.searlelab.jchronologer.util.TsvTable;

class GoldenInferenceParityTest {

    private static final String GOLDEN_RESOURCE = "data/golden/chronologer_parity_cases.golden.json";
    private static final String PARITY_CASES_RESOURCE = "data/golden/parity_cases.tsv";

    @Test
    void predictionsMatchGoldenFixtureWithinTolerance() throws IOException {
        List<String> peptides = loadInputPeptides();
        JsonNode payload = readGoldenJson();

        PredictionResult result;
        try (Chronologer chronologer = ChronologerFactory.createDefault()) {
            result = chronologer.predict(peptides);
        }

        assertEquals(payload.get("accepted").size(), result.getAcceptedCount());
        assertEquals(payload.get("rejected").size(), result.getRejectedCount());

        Map<Integer, AcceptedPrediction> acceptedByRow = result.getAcceptedByRowIndex();

        for (JsonNode acceptedNode : payload.get("accepted")) {
            int rowIndex = acceptedNode.get("RowIndex").asInt();
            float expected = (float) acceptedNode.get("Pred_HI").asDouble();
            AcceptedPrediction actual = acceptedByRow.get(rowIndex);
            assertNotNull(actual, "Missing accepted prediction for row " + rowIndex);

            float diff = Math.abs(actual.getPredHi() - expected);
            float rel = expected == 0f ? diff : diff / Math.abs(expected);
            assertTrue(diff <= 1e-5f, "Absolute diff too large at row " + rowIndex + ": " + diff);
            assertTrue(rel <= 1e-5f, "Relative diff too large at row " + rowIndex + ": " + rel);
        }

        Map<Integer, RejectedPrediction> rejectedByRow = result.getRejectedByRowIndex();

        for (JsonNode rejectedNode : payload.get("rejected")) {
            int rowIndex = rejectedNode.get("RowIndex").asInt();
            RejectedPrediction actual = rejectedByRow.get(rowIndex);
            assertNotNull(actual, "Missing rejected prediction for row " + rowIndex);
            assertEquals(rejectedNode.get("RejectionReason").asText(), actual.getRejectionReason().name());
        }
    }

    private static List<String> loadInputPeptides() throws IOException {
        try (InputStream stream = Thread.currentThread()
                .getContextClassLoader()
                .getResourceAsStream(PARITY_CASES_RESOURCE)) {
            if (stream == null) {
                throw new IllegalStateException("Missing resource: " + PARITY_CASES_RESOURCE);
            }
            java.nio.file.Path tempFile = java.nio.file.Files.createTempFile("jchronologer-parity", ".tsv");
            java.nio.file.Files.copy(stream, tempFile, java.nio.file.StandardCopyOption.REPLACE_EXISTING);
            TsvTable table = TsvTable.read(tempFile);
            int idx = table.columnIndex("PeptideModSeq");
            if (idx < 0) {
                throw new IllegalStateException("Missing PeptideModSeq column in parity cases.");
            }
            List<String> peptides = new ArrayList<>();
            for (String[] row : table.getRows()) {
                peptides.add(row[idx]);
            }
            return peptides;
        }
    }

    private static JsonNode readGoldenJson() throws IOException {
        ObjectMapper mapper = new ObjectMapper();
        try (InputStream stream = Thread.currentThread().getContextClassLoader().getResourceAsStream(GOLDEN_RESOURCE)) {
            if (stream == null) {
                throw new IllegalStateException("Missing resource: " + GOLDEN_RESOURCE);
            }
            return mapper.readTree(stream);
        }
    }
}

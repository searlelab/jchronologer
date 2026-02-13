package org.searlelab.jchronologer.preprocessing;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;
import org.junit.jupiter.api.Test;
import org.searlelab.jchronologer.api.RejectionReason;
import org.searlelab.jchronologer.preprocessing.PreprocessingMetadataLoader.CompiledPreprocessingMetadata;

class PreprocessingParityTest {

    private static final String GOLDEN_RESOURCE = "data/golden/chronologer_parity_cases.golden.json";

    @Test
    void quantifierNormalizationMatchesPythonStylePatterns() {
        String normalized = PreprocessingMetadataLoader.normalizePythonQuantifiers("M\\[\\+15\\.99.{,6}\\]");
        assertEquals("M\\[\\+15\\.99.{0,6}\\]", normalized);
    }

    @Test
    void preprocessingMatchesGoldenAcceptedAndRejectedRecords() throws IOException {
        CompiledPreprocessingMetadata metadata =
                PreprocessingMetadataLoader.loadFromClasspath("models/Chronologer_20220601193755.preprocessing.json");
        ChronologerPreprocessor preprocessor = new ChronologerPreprocessor(metadata);
        JsonNode payload = readGoldenJson();

        for (JsonNode acceptedNode : payload.get("accepted")) {
            String peptide = acceptedNode.get("PeptideModSeq").asText();
            PreprocessingOutcome outcome = preprocessor.preprocess(peptide);
            assertTrue(outcome.isAccepted(), "Expected accepted peptide: " + peptide);
            assertEquals(acceptedNode.get("PatchedPeptideModSeq").asText(), outcome.getPatchedPeptideModSeq());
            assertEquals(acceptedNode.get("CodedPeptideSeq").asText(), outcome.getCodedPeptideSeq());
            assertArrayEquals(toLongArray(acceptedNode.get("TokenArray")), outcome.getTokenArray());
        }

        for (JsonNode rejectedNode : payload.get("rejected")) {
            String peptide = rejectedNode.get("PeptideModSeq").asText();
            PreprocessingOutcome outcome = preprocessor.preprocess(peptide);
            assertFalse(outcome.isAccepted(), "Expected rejected peptide: " + peptide);
            assertEquals(
                    RejectionReason.valueOf(rejectedNode.get("RejectionReason").asText()),
                    outcome.getRejectionReason());
            if (outcome.getRejectionReason() == RejectionReason.TOKENIZATION_ERROR) {
                assertNotNull(outcome.getErrorDetail());
            }
        }
    }

    @Test
    void shortNtermAnnotationIsRejectedWithTokenizationError() {
        CompiledPreprocessingMetadata metadata =
                PreprocessingMetadataLoader.loadFromClasspath("models/Chronologer_20220601193755.preprocessing.json");
        ChronologerPreprocessor preprocessor = new ChronologerPreprocessor(metadata);

        PreprocessingOutcome outcome = preprocessor.preprocess("[+1]ACDEFGHIK");

        assertFalse(outcome.isAccepted());
        assertEquals(RejectionReason.TOKENIZATION_ERROR, outcome.getRejectionReason());
        assertNotNull(outcome.getErrorDetail());
        assertTrue(outcome.getErrorDetail().contains("N-term annotation is too short"));
    }

    @Test
    void unknownNtermAnnotationIsRejectedWithTokenizationError() {
        CompiledPreprocessingMetadata metadata =
                PreprocessingMetadataLoader.loadFromClasspath("models/Chronologer_20220601193755.preprocessing.json");
        ChronologerPreprocessor preprocessor = new ChronologerPreprocessor(metadata);

        PreprocessingOutcome outcome = preprocessor.preprocess("[+99.999999]ACDEFGHIK");

        assertFalse(outcome.isAccepted());
        assertEquals(RejectionReason.TOKENIZATION_ERROR, outcome.getRejectionReason());
        assertNotNull(outcome.getErrorDetail());
        assertTrue(outcome.getErrorDetail().contains("Unknown N-term key"));
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

    private static long[] toLongArray(JsonNode node) {
        List<Long> values = new ArrayList<>();
        for (JsonNode item : node) {
            values.add(item.asLong());
        }
        long[] result = new long[values.size()];
        for (int i = 0; i < values.size(); i++) {
            result[i] = values.get(i);
        }
        return result;
    }
}

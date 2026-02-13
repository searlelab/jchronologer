package org.searlelab.jchronologer.api;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotSame;
import static org.junit.jupiter.api.Assertions.assertThrows;

import java.util.List;
import java.util.Map;
import org.junit.jupiter.api.Test;

class ApiModelTest {

    @Test
    void acceptedPredictionDefensivelyCopiesTokenArray() {
        long[] input = new long[] {1L, 2L, 3L};
        AcceptedPrediction prediction = new AcceptedPrediction(
                7,
                "PEPTIDE",
                "PEPTIDE",
                "_PEPTIDE_",
                input,
                12.34f);

        input[0] = 99L;
        long[] fromGetter = prediction.getTokenArray();
        fromGetter[1] = 88L;

        assertEquals(7, prediction.getRowIndex());
        assertEquals("PEPTIDE", prediction.getPeptideModSeq());
        assertEquals("PEPTIDE", prediction.getPatchedPeptideModSeq());
        assertEquals("_PEPTIDE_", prediction.getCodedPeptideSeq());
        assertEquals(12.34f, prediction.getPredHi());
        assertNotSame(input, fromGetter);
        assertArrayEquals(new long[] {1L, 2L, 3L}, prediction.getTokenArray());
    }

    @Test
    void rejectedPredictionGettersReturnConstructorValues() {
        RejectedPrediction rejected = new RejectedPrediction(
                5,
                "BADSEQ",
                "BADSEQ",
                RejectionReason.TOKENIZATION_ERROR,
                "Invalid token");

        assertEquals(5, rejected.getRowIndex());
        assertEquals("BADSEQ", rejected.getPeptideModSeq());
        assertEquals("BADSEQ", rejected.getPatchedPeptideModSeq());
        assertEquals(RejectionReason.TOKENIZATION_ERROR, rejected.getRejectionReason());
        assertEquals("Invalid token", rejected.getErrorDetail());
    }

    @Test
    void chronologerOptionsBuilderSupportsDefaultsAndCustomValues() {
        ChronologerOptions defaults = ChronologerOptions.builder().build();
        assertEquals(ChronologerOptions.DEFAULT_MODEL_RESOURCE, defaults.getModelResource());
        assertEquals(ChronologerOptions.DEFAULT_PREPROCESSING_RESOURCE, defaults.getPreprocessingResource());
        assertEquals(ChronologerOptions.DEFAULT_BATCH_SIZE, defaults.getBatchSize());
        assertEquals(ChronologerOptions.DEFAULT_INFERENCE_THREADS, defaults.getInferenceThreads());

        ChronologerOptions custom = ChronologerOptions.builder()
                .modelResource("models/custom.pt")
                .preprocessingResource("models/custom.json")
                .batchSize(128)
                .inferenceThreads(2)
                .build();
        assertEquals("models/custom.pt", custom.getModelResource());
        assertEquals("models/custom.json", custom.getPreprocessingResource());
        assertEquals(128, custom.getBatchSize());
        assertEquals(2, custom.getInferenceThreads());
    }

    @Test
    void chronologerOptionsBuilderRejectsInvalidValues() {
        assertThrows(IllegalArgumentException.class, () -> ChronologerOptions.builder()
                .modelResource(null)
                .build());
        assertThrows(IllegalArgumentException.class, () -> ChronologerOptions.builder()
                .modelResource(" ")
                .build());
        assertThrows(IllegalArgumentException.class, () -> ChronologerOptions.builder()
                .preprocessingResource(null)
                .build());
        assertThrows(IllegalArgumentException.class, () -> ChronologerOptions.builder()
                .preprocessingResource("")
                .build());
        assertThrows(IllegalArgumentException.class, () -> ChronologerOptions.builder()
                .batchSize(0)
                .build());
        assertThrows(IllegalArgumentException.class, () -> ChronologerOptions.builder()
                .inferenceThreads(0)
                .build());
    }

    @Test
    void predictionResultRowIndexViewsExposeAcceptedAndRejectedRows() {
        AcceptedPrediction acceptedFirst = new AcceptedPrediction(
                3,
                "AAA",
                "AAA",
                "_AAA_",
                new long[] {1L, 2L},
                1.0f);
        AcceptedPrediction acceptedSecond = new AcceptedPrediction(
                3,
                "BBB",
                "BBB",
                "_BBB_",
                new long[] {3L, 4L},
                2.0f);
        RejectedPrediction rejected = new RejectedPrediction(
                5,
                "BAD",
                "BAD",
                RejectionReason.TOKENIZATION_ERROR,
                "bad token");

        PredictionResult result = new PredictionResult(
                List.of(acceptedFirst, acceptedSecond),
                List.of(rejected));

        Map<Integer, AcceptedPrediction> acceptedByRow = result.getAcceptedByRowIndex();
        Map<Integer, RejectedPrediction> rejectedByRow = result.getRejectedByRowIndex();

        assertEquals(1, acceptedByRow.size());
        assertEquals("BBB", acceptedByRow.get(3).getPeptideModSeq());
        assertEquals(2.0f, acceptedByRow.get(3).getPredHi());
        assertEquals(1, rejectedByRow.size());
        assertEquals(RejectionReason.TOKENIZATION_ERROR, rejectedByRow.get(5).getRejectionReason());
        assertThrows(UnsupportedOperationException.class, () -> acceptedByRow.put(7, acceptedFirst));
        assertThrows(UnsupportedOperationException.class, () -> rejectedByRow.clear());
    }
}

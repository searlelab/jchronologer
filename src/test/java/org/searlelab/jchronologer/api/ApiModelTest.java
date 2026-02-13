package org.searlelab.jchronologer.api;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotSame;
import static org.junit.jupiter.api.Assertions.assertThrows;

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

        ChronologerOptions custom = ChronologerOptions.builder()
                .modelResource("models/custom.pt")
                .preprocessingResource("models/custom.json")
                .batchSize(128)
                .build();
        assertEquals("models/custom.pt", custom.getModelResource());
        assertEquals("models/custom.json", custom.getPreprocessingResource());
        assertEquals(128, custom.getBatchSize());
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
    }
}

package org.searlelab.jchronologer.inference;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import org.junit.jupiter.api.Test;
import org.searlelab.jchronologer.api.AcceptedPrediction;
import org.searlelab.jchronologer.api.Chronologer;
import org.searlelab.jchronologer.api.ChronologerOptions;
import org.searlelab.jchronologer.api.PredictionResult;
import org.searlelab.jchronologer.api.RejectedPrediction;
import org.searlelab.jchronologer.impl.ChronologerFactory;

class ParallelInferenceParityTest {

    @Test
    void parallelBatchInferenceMatchesSingleThreadOutput() {
        List<String> peptides = buildInputPeptides();

        ChronologerOptions serialOptions = ChronologerOptions.builder()
                .batchSize(4)
                .inferenceThreads(1)
                .build();

        ChronologerOptions parallelOptions = ChronologerOptions.builder()
                .batchSize(4)
                .inferenceThreads(2)
                .build();

        PredictionResult serialResult;
        PredictionResult parallelResult;
        try (Chronologer serial = ChronologerFactory.create(serialOptions);
                Chronologer parallel = ChronologerFactory.create(parallelOptions)) {
            serialResult = serial.predict(peptides);
            parallelResult = parallel.predict(peptides);
        }

        assertEquals(serialResult.getAcceptedCount(), parallelResult.getAcceptedCount());
        assertEquals(serialResult.getRejectedCount(), parallelResult.getRejectedCount());

        Map<Integer, AcceptedPrediction> serialAccepted = serialResult.getAcceptedByRowIndex();
        Map<Integer, AcceptedPrediction> parallelAccepted = parallelResult.getAcceptedByRowIndex();
        assertEquals(serialAccepted.size(), parallelAccepted.size());

        for (Map.Entry<Integer, AcceptedPrediction> entry : serialAccepted.entrySet()) {
            int rowIndex = entry.getKey();
            AcceptedPrediction serialPrediction = entry.getValue();
            AcceptedPrediction parallelPrediction = parallelAccepted.get(rowIndex);
            assertNotNull(parallelPrediction, "Missing accepted prediction for row " + rowIndex);
            assertEquals(serialPrediction.getPeptideModSeq(), parallelPrediction.getPeptideModSeq());
            assertEquals(serialPrediction.getPatchedPeptideModSeq(), parallelPrediction.getPatchedPeptideModSeq());
            assertEquals(serialPrediction.getCodedPeptideSeq(), parallelPrediction.getCodedPeptideSeq());
            assertArrayEquals(serialPrediction.getTokenArray(), parallelPrediction.getTokenArray());
            assertEquals(serialPrediction.getPredHi(), parallelPrediction.getPredHi(), 1e-5f);
        }

        Map<Integer, RejectedPrediction> serialRejected = serialResult.getRejectedByRowIndex();
        Map<Integer, RejectedPrediction> parallelRejected = parallelResult.getRejectedByRowIndex();
        assertEquals(serialRejected.size(), parallelRejected.size());
        for (Map.Entry<Integer, RejectedPrediction> entry : serialRejected.entrySet()) {
            int rowIndex = entry.getKey();
            RejectedPrediction serialRejection = entry.getValue();
            RejectedPrediction parallelRejection = parallelRejected.get(rowIndex);
            assertNotNull(parallelRejection, "Missing rejected prediction for row " + rowIndex);
            assertEquals(serialRejection.getRejectionReason(), parallelRejection.getRejectionReason());
            assertEquals(serialRejection.getErrorDetail(), parallelRejection.getErrorDetail());
        }
    }

    private static List<String> buildInputPeptides() {
        List<String> peptides = new ArrayList<>();
        for (int i = 0; i < 24; i++) {
            peptides.add("VATVSLPR");
        }
        for (int i = 0; i < 8; i++) {
            peptides.add("[42.010565]ACDEFGHIK");
        }
        assertTrue(peptides.size() > 8);
        return peptides;
    }
}

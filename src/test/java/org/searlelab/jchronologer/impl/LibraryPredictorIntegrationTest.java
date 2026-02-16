package org.searlelab.jchronologer.impl;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;
import org.junit.jupiter.api.Test;
import org.searlelab.jchronologer.api.ChronologerLibraryEntry;
import org.searlelab.jchronologer.api.ChronologerLibraryPredictor;
import org.searlelab.jchronologer.api.LibraryPredictionRequest;
import org.searlelab.jchronologer.api.PrecursorCondition;

class LibraryPredictorIntegrationTest {

    @Test
    void predictsMultipleChargeNcePairsAndFiltersMissingIons() {
        LibraryPredictionRequest request = new LibraryPredictionRequest(
                "K[458.325864]PGLAITFAK[229.162932]",
                List.of(
                        new PrecursorCondition((byte) 2, 0.30),
                        new PrecursorCondition((byte) 3, 0.35)));

        try (ChronologerLibraryPredictor predictor = ChronologerFactory.createLibraryPredictorDefault()) {
            List<ChronologerLibraryEntry> entries = predictor.predict(List.of(request));

            assertEquals(2, entries.size());
            assertEquals("[UNIMOD:737]-K[UNIMOD:737]PGLAITFAK[UNIMOD:737]-[]", entries.get(0).getUnimodPeptideSequence());
            assertEquals("[UNIMOD:737]-K[UNIMOD:737]PGLAITFAK[UNIMOD:737]-[]", entries.get(1).getUnimodPeptideSequence());
            assertEquals((byte) 2, entries.get(0).getPrecursorCharge());
            assertEquals((byte) 3, entries.get(1).getPrecursorCharge());

            assertEquals(entries.get(0).getRetentionTimeInSeconds(), entries.get(1).getRetentionTimeInSeconds(), 1e-5f);
            assertTrue(entries.get(0).getPrecursorMz() > 0.0);
            assertTrue(entries.get(1).getPrecursorMz() > 0.0);

            for (ChronologerLibraryEntry entry : entries) {
                assertEquals(entry.getMassArray().length, entry.getIntensityArray().length);
                assertEquals(entry.getMassArray().length, entry.getIonTypeArray().length);
                for (float intensity : entry.getIntensityArray()) {
                    assertTrue(intensity >= 0.01f);
                }
                for (double mz : entry.getMassArray()) {
                    assertTrue(mz > 0.0);
                }
            }
        }
    }

}

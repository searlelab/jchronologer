package org.searlelab.jchronologer.impl;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.List;
import org.junit.jupiter.api.Test;
import org.searlelab.jchronologer.api.ChronologerLibraryEntry;
import org.searlelab.jchronologer.api.ChronologerLibraryPredictor;
import org.searlelab.jchronologer.api.LibraryPredictionRequest;
import org.searlelab.jchronologer.api.PrecursorCondition;

class LibraryPredictorIntegrationTest {
	@Test
	void predictTASEFDSAIAQDK() {
		ChronologerLibraryPredictor chronologer = ChronologerFactory.createLibraryPredictorDefault();

        List<LibraryPredictionRequest> peptides = new ArrayList<>();
        List<PrecursorCondition> conditions=new ArrayList<>();
        conditions.add(new PrecursorCondition((byte)2, 27.0));
        peptides.add(new LibraryPredictionRequest("TASEFDSAIAQDK", conditions));

        List<ChronologerLibraryEntry> result = chronologer.predict(peptides);

        for (ChronologerLibraryEntry entry : result) {
        	assertEquals("TASEFDSAIAQDK", entry.getPeptideModSeq());
        	assertEquals((byte)2, entry.getPrecursorCharge());
        	assertEquals(691.8253461193799, entry.getPrecursorMz(), 0.0000001);
        	assertEquals(8.757546, entry.getRetentionTimeInSeconds()/60f, 0.0001);
        	
			assertTrue(entry.getMassArray().length>1);
			
			for (int i=0; i<entry.getMassArray().length; i++) {
				System.out.println("\t"+entry.getMassArray()[i]+"\t"+entry.getIntensityArray()[i]);
			}
		}
	}

    @Test
    void predictsMultipleChargeNcePairsAndFiltersMissingIons() {
        LibraryPredictionRequest request = new LibraryPredictionRequest(
                "K[458.325864]PGLAITFAK[229.162932]",
                List.of(
                        new PrecursorCondition((byte) 2, 30.0),
                        new PrecursorCondition((byte) 3, 35.0)));

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

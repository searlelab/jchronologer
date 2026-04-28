package org.searlelab.jchronologer.impl;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import org.junit.jupiter.api.Test;
import org.searlelab.jchronologer.api.ChronologerLibraryEntry;
import org.searlelab.jchronologer.api.ChronologerLibraryOptions;
import org.searlelab.jchronologer.api.ChronologerLibraryPredictor;
import org.searlelab.jchronologer.api.LibraryPredictionRequest;
import org.searlelab.jchronologer.util.TsvTable;

class ScoutPredictionBenchmarkTest {

    private static final Path DEMO_PEPTIDES = Path.of("src/test/resources/data/demo/crap_peptides.txt");
    private static final double DEFAULT_NCE = 33.0;
    private static final double DEFAULT_MIN_CHARGE_PROBABILITY = 0.01;

    @Test
    void scoutPredictionTimeBeatsLegacyPredictorOnDemoPeptides() throws Exception {
        TsvTable table = TsvTable.read(DEMO_PEPTIDES);
        int peptideColumn = table.columnIndex("PeptideModSeq");
        assertTrue(peptideColumn >= 0, "Missing PeptideModSeq column in demo peptide list.");

        List<LibraryPredictionRequest> requests = new ArrayList<>(table.getRows().size());
        for (String[] row : table.getRows()) {
            requests.add(new LibraryPredictionRequest(
                    row[peptideColumn],
                    DEFAULT_NCE,
                    DEFAULT_MIN_CHARGE_PROBABILITY));
        }

        ChronologerLibraryOptions options = ChronologerLibraryOptions.builder()
                .batchSize(2048)
                .cartographerBatchSize(2048)
                .inferenceThreads(1)
                .build();

        try (ChronologerLibraryPredictor legacy = ChronologerFactory.createLibraryPredictor(options);
                ChronologerLibraryPredictor scout = ChronologerFactory.createFastLibraryPredictor(options);
                ChronologerLibraryPredictor electricianSelectedScout = new DefaultScoutLibraryPredictor(
                        options,
                        DefaultScoutLibraryPredictor.AutomaticChargeSelectionMode.ELECTRICIAN_SELECTION)) {
            long legacyStartNanos = System.nanoTime();
            List<ChronologerLibraryEntry> legacyEntries = legacy.predict(requests);
            long legacyElapsedNanos = System.nanoTime() - legacyStartNanos;

            long scoutStartNanos = System.nanoTime();
            List<ChronologerLibraryEntry> scoutEntries = scout.predict(requests);
            long scoutElapsedNanos = System.nanoTime() - scoutStartNanos;

            long electricianSelectedScoutStartNanos = System.nanoTime();
            List<ChronologerLibraryEntry> electricianSelectedScoutEntries = electricianSelectedScout.predict(requests);
            long electricianSelectedScoutElapsedNanos = System.nanoTime() - electricianSelectedScoutStartNanos;

            assertTrue(!legacyEntries.isEmpty(), "Legacy predictor returned no entries.");
            assertTrue(!scoutEntries.isEmpty(), "Scout predictor returned no entries.");
            assertTrue(!electricianSelectedScoutEntries.isEmpty(),
                    "Electrician-selected Scout predictor returned no entries.");

            double legacyMillis = legacyElapsedNanos / 1_000_000.0;
            double scoutMillis = scoutElapsedNanos / 1_000_000.0;
            double electricianSelectedScoutMillis = electricianSelectedScoutElapsedNanos / 1_000_000.0;
            double speedupRatio = legacyMillis / scoutMillis;
            double electricianSelectedScoutSpeedupRatio = legacyMillis / electricianSelectedScoutMillis;

            System.out.printf(
                    "Scout benchmark on %d requests: slow_mode=%.3f fast_mode=%.3f ratio=%.3f electrician_selected_fast_mode=%.3f electrician_selected_ratio=%.3f%n",
                    requests.size(),
                    legacyMillis,
                    scoutMillis,
                    speedupRatio,
                    electricianSelectedScoutMillis,
                    electricianSelectedScoutSpeedupRatio);

            assertTrue(speedupRatio > 1.0,
                    "Expected Scout prediction time to beat legacy predictor, but ratio was " + speedupRatio);
        }
    }
}

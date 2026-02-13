package org.searlelab.jchronologer.api;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

import java.util.List;
import org.junit.jupiter.api.Test;
import org.searlelab.jchronologer.impl.ChronologerFactory;

class ChronologerLifecycleTest {

    @Test
    void initIsIdempotentAndSupportsMultiplePredictRounds() {
        ChronologerOptions options = ChronologerOptions.builder()
                .batchSize(4)
                .inferenceThreads(1)
                .build();

        try (Chronologer chronologer = ChronologerFactory.create(options)) {
            chronologer.init();
            chronologer.init();

            PredictionResult first = chronologer.predict(List.of(
                    "VATVSLPR",
                    "[42.010565]ACDEFGHIK"));
            PredictionResult second = chronologer.predict(List.of("LGEHNIDVLEGNEQFINAAK"));

            assertEquals(1, first.getAcceptedCount());
            assertEquals(1, first.getRejectedCount());
            assertEquals(1, second.getAcceptedCount());
            assertEquals(0, second.getRejectedCount());
        }
    }

    @Test
    void initAfterCloseThrows() {
        Chronologer chronologer = ChronologerFactory.createDefault();
        chronologer.close();
        assertThrows(IllegalStateException.class, chronologer::init);
    }
}

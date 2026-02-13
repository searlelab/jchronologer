package org.searlelab.jchronologer.inference;

import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;

class BatchPredictorTest {

    @Test
    void constructorRejectsNonModelResource() {
        RuntimeException error = assertThrows(
                RuntimeException.class,
                () -> new BatchPredictor("data/demo/demo_peptides.txt"));
        assertTrue(!error.getMessage().isBlank());
    }
}

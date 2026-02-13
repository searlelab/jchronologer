package org.searlelab.jchronologer.preprocessing;

import static org.junit.jupiter.api.Assertions.assertDoesNotThrow;
import static org.junit.jupiter.api.Assertions.assertNotNull;

import java.util.Random;
import org.junit.jupiter.api.Test;
import org.searlelab.jchronologer.preprocessing.PreprocessingMetadataLoader.CompiledPreprocessingMetadata;

class PreprocessingFuzzTest {

    private static final String ALPHABET = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789[]()+-._ \t";

    @Test
    void preprocessHandlesMalformedInputsWithoutThrowing() {
        CompiledPreprocessingMetadata metadata =
                PreprocessingMetadataLoader.loadFromClasspath("models/Chronologer_20220601193755.preprocessing.json");
        ChronologerPreprocessor preprocessor = new ChronologerPreprocessor(metadata);
        Random random = new Random(20260213L);

        for (int i = 0; i < 4000; i++) {
            final int caseIndex = i;
            final String candidate = (i % 41 == 0) ? null : randomString(random, random.nextInt(80));
            PreprocessingOutcome outcome = assertDoesNotThrow(
                    () -> preprocessor.preprocess(candidate),
                    () -> "Unexpected preprocess exception for input index " + caseIndex + ": " + candidate);
            assertNotNull(outcome);
            if (!outcome.isAccepted()) {
                assertNotNull(outcome.getRejectionReason());
            }
        }
    }

    private static String randomString(Random random, int len) {
        StringBuilder sb = new StringBuilder(len);
        for (int i = 0; i < len; i++) {
            sb.append(ALPHABET.charAt(random.nextInt(ALPHABET.length())));
        }
        return sb.toString();
    }
}

package org.searlelab.jchronologer.api;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;

class ChronologerLibraryOptionsBuilderTest {

    @Test
    void builderAppliesAllCustomFields() {
        ChronologerLibraryOptions options = ChronologerLibraryOptions.builder()
                .chronologerModelResource("models/chronologer-custom.pt")
                .chronologerPreprocessingResource("models/chronologer-custom.json")
                .cartographerModelResource("models/cartographer-custom.pt")
                .cartographerPreprocessingResource("models/cartographer-custom.json")
                .electricianModelResource("models/electrician-custom.pt")
                .electricianPreprocessingResource("models/electrician-custom.json")
                .batchSize(321)
                .cartographerBatchSize(123)
                .inferenceThreads(4)
                .massMatchEpsilon(1e-4)
                .minimumReportedIntensity(0.25f)
                .verboseLogging(true)
                .build();

        assertEquals("models/chronologer-custom.pt", options.getChronologerModelResource());
        assertEquals("models/chronologer-custom.json", options.getChronologerPreprocessingResource());
        assertEquals("models/cartographer-custom.pt", options.getCartographerModelResource());
        assertEquals("models/cartographer-custom.json", options.getCartographerPreprocessingResource());
        assertEquals("models/electrician-custom.pt", options.getElectricianModelResource());
        assertEquals("models/electrician-custom.json", options.getElectricianPreprocessingResource());
        assertEquals(321, options.getBatchSize());
        assertEquals(123, options.getCartographerBatchSize());
        assertEquals(4, options.getInferenceThreads());
        assertEquals(1e-4, options.getMassMatchEpsilon());
        assertEquals(0.25f, options.getMinimumReportedIntensity());
        assertTrue(options.isVerboseLogging());
    }

    @Test
    void builderRejectsInvalidStringFields() {
        assertThrows(IllegalArgumentException.class, () -> ChronologerLibraryOptions.builder()
                .chronologerModelResource(null)
                .build());
        assertThrows(IllegalArgumentException.class, () -> ChronologerLibraryOptions.builder()
                .chronologerModelResource(" ")
                .build());

        assertThrows(IllegalArgumentException.class, () -> ChronologerLibraryOptions.builder()
                .chronologerPreprocessingResource(null)
                .build());
        assertThrows(IllegalArgumentException.class, () -> ChronologerLibraryOptions.builder()
                .chronologerPreprocessingResource(" ")
                .build());

        assertThrows(IllegalArgumentException.class, () -> ChronologerLibraryOptions.builder()
                .cartographerModelResource(null)
                .build());
        assertThrows(IllegalArgumentException.class, () -> ChronologerLibraryOptions.builder()
                .cartographerModelResource(" ")
                .build());

        assertThrows(IllegalArgumentException.class, () -> ChronologerLibraryOptions.builder()
                .cartographerPreprocessingResource(null)
                .build());
        assertThrows(IllegalArgumentException.class, () -> ChronologerLibraryOptions.builder()
                .cartographerPreprocessingResource(" ")
                .build());

        assertThrows(IllegalArgumentException.class, () -> ChronologerLibraryOptions.builder()
                .electricianModelResource(null)
                .build());
        assertThrows(IllegalArgumentException.class, () -> ChronologerLibraryOptions.builder()
                .electricianModelResource(" ")
                .build());

        assertThrows(IllegalArgumentException.class, () -> ChronologerLibraryOptions.builder()
                .electricianPreprocessingResource(null)
                .build());
        assertThrows(IllegalArgumentException.class, () -> ChronologerLibraryOptions.builder()
                .electricianPreprocessingResource(" ")
                .build());
    }

    @Test
    void builderRejectsInvalidNumericFields() {
        assertThrows(IllegalArgumentException.class, () -> ChronologerLibraryOptions.builder()
                .batchSize(0)
                .build());
        assertThrows(IllegalArgumentException.class, () -> ChronologerLibraryOptions.builder()
                .cartographerBatchSize(0)
                .build());
        assertThrows(IllegalArgumentException.class, () -> ChronologerLibraryOptions.builder()
                .inferenceThreads(0)
                .build());
        assertThrows(IllegalArgumentException.class, () -> ChronologerLibraryOptions.builder()
                .massMatchEpsilon(0.0)
                .build());
        assertThrows(IllegalArgumentException.class, () -> ChronologerLibraryOptions.builder()
                .minimumReportedIntensity(-0.0001f)
                .build());
    }
}

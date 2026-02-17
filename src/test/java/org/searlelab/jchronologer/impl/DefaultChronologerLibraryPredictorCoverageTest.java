package org.searlelab.jchronologer.impl;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertInstanceOf;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import org.junit.jupiter.api.Test;
import org.searlelab.jchronologer.inference.CartographerSpectrumDecoder;

class DefaultChronologerLibraryPredictorCoverageTest {

    @Test
    void loadCartographerMetadataParsesCustomValues() throws Exception {
        Object metadata = invokeLoadMetadata("data/preprocessing/cartographer_metadata_custom.json");

        assertEquals(2, getIntField(metadata, "minPrecursorCharge"));
        assertEquals(5, getIntField(metadata, "maxPrecursorCharge"));
        assertEquals(222, getIntField(metadata, "outputWidth"));
    }

    @Test
    void loadCartographerMetadataUsesDefaultsWhenFieldsMissing() throws Exception {
        Object metadata = invokeLoadMetadata("data/preprocessing/cartographer_metadata_minimal.json");

        assertEquals(1, getIntField(metadata, "minPrecursorCharge"));
        assertEquals(6, getIntField(metadata, "maxPrecursorCharge"));
        assertEquals(CartographerSpectrumDecoder.VECTOR_LENGTH, getIntField(metadata, "outputWidth"));
    }

    @Test
    void loadCartographerMetadataRejectsMissingResource() throws Exception {
        Method method = loadMetadataMethod();
        InvocationTargetException error = assertThrows(
                InvocationTargetException.class,
                () -> method.invoke(null, "data/preprocessing/does_not_exist.json"));

        IllegalArgumentException cause = assertInstanceOf(IllegalArgumentException.class, error.getCause());
        assertTrue(cause.getMessage().contains("Missing cartographer metadata resource"));
    }

    @Test
    void loadCartographerMetadataWrapsMalformedJson() throws Exception {
        Method method = loadMetadataMethod();
        InvocationTargetException error = assertThrows(
                InvocationTargetException.class,
                () -> method.invoke(null, "data/preprocessing/cartographer_metadata_malformed.json"));

        IllegalStateException cause = assertInstanceOf(IllegalStateException.class, error.getCause());
        assertTrue(cause.getMessage().contains("Failed to parse cartographer metadata resource"));
    }

    private static Object invokeLoadMetadata(String resource) throws Exception {
        Method method = loadMetadataMethod();
        return method.invoke(null, resource);
    }

    private static Method loadMetadataMethod() throws NoSuchMethodException {
        Method method = DefaultChronologerLibraryPredictor.class.getDeclaredMethod(
                "loadCartographerMetadata",
                String.class);
        method.setAccessible(true);
        return method;
    }

    private static int getIntField(Object target, String fieldName) throws Exception {
        Field field = target.getClass().getDeclaredField(fieldName);
        field.setAccessible(true);
        return field.getInt(target);
    }
}

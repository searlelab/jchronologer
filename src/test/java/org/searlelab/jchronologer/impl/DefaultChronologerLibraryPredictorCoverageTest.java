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
        Object metadata = invokeLoadCartographerMetadata("data/preprocessing/cartographer_metadata_custom.json");

        assertEquals(2, getIntField(metadata, "minPrecursorCharge"));
        assertEquals(5, getIntField(metadata, "maxPrecursorCharge"));
        assertEquals(222, getIntField(metadata, "outputWidth"));
    }

    @Test
    void loadCartographerMetadataUsesDefaultsWhenFieldsMissing() throws Exception {
        Object metadata = invokeLoadCartographerMetadata("data/preprocessing/cartographer_metadata_minimal.json");

        assertEquals(1, getIntField(metadata, "minPrecursorCharge"));
        assertEquals(6, getIntField(metadata, "maxPrecursorCharge"));
        assertEquals(CartographerSpectrumDecoder.VECTOR_LENGTH, getIntField(metadata, "outputWidth"));
    }

    @Test
    void loadCartographerMetadataRejectsMissingResource() throws Exception {
        Method method = loadCartographerMetadataMethod();
        InvocationTargetException error = assertThrows(
                InvocationTargetException.class,
                () -> method.invoke(null, "data/preprocessing/does_not_exist.json"));

        IllegalArgumentException cause = assertInstanceOf(IllegalArgumentException.class, error.getCause());
        assertTrue(cause.getMessage().contains("Missing cartographer metadata resource"));
    }

    @Test
    void loadCartographerMetadataWrapsMalformedJson() throws Exception {
        Method method = loadCartographerMetadataMethod();
        InvocationTargetException error = assertThrows(
                InvocationTargetException.class,
                () -> method.invoke(null, "data/preprocessing/cartographer_metadata_malformed.json"));

        IllegalStateException cause = assertInstanceOf(IllegalStateException.class, error.getCause());
        assertTrue(cause.getMessage().contains("Failed to parse cartographer metadata resource"));
    }

    @Test
    void loadElectricianMetadataParsesCustomOutputWidth() throws Exception {
        Object metadata = invokeLoadElectricianMetadata("data/preprocessing/cartographer_metadata_custom.json");

        assertEquals(222, getIntField(metadata, "outputWidth"));
    }

    @Test
    void loadElectricianMetadataUsesDefaultWhenFieldsMissing() throws Exception {
        Object metadata = invokeLoadElectricianMetadata("data/preprocessing/cartographer_metadata_minimal.json");

        assertEquals(6, getIntField(metadata, "outputWidth"));
    }

    @Test
    void loadElectricianMetadataRejectsMissingResource() throws Exception {
        Method method = loadElectricianMetadataMethod();
        InvocationTargetException error = assertThrows(
                InvocationTargetException.class,
                () -> method.invoke(null, "data/preprocessing/does_not_exist.json"));

        IllegalArgumentException cause = assertInstanceOf(IllegalArgumentException.class, error.getCause());
        assertTrue(cause.getMessage().contains("Missing electrician metadata resource"));
    }

    @Test
    void loadElectricianMetadataWrapsMalformedJson() throws Exception {
        Method method = loadElectricianMetadataMethod();
        InvocationTargetException error = assertThrows(
                InvocationTargetException.class,
                () -> method.invoke(null, "data/preprocessing/cartographer_metadata_malformed.json"));

        IllegalStateException cause = assertInstanceOf(IllegalStateException.class, error.getCause());
        assertTrue(cause.getMessage().contains("Failed to parse electrician metadata resource"));
    }

    @Test
    void loadSculptorMetadataParsesBundledResource() throws Exception {
        Object metadata = invokeLoadSculptorMetadata("models/Sculptor_20260311095327.preprocessing.json");

        assertEquals(50, getIntField(metadata, "maxPeptideLength"));
        assertEquals(6, getIntField(metadata, "chargeStateCount"));
        assertEquals(1, getIntField(metadata, "outputWidth"));
        assertEquals(490.40022836690935, getDoubleField(metadata, "ccsMean"), 1e-12);
        assertEquals(126.1480333943971, getDoubleField(metadata, "ccsStd"), 1e-12);
    }

    @Test
    void loadSculptorMetadataRejectsMissingResource() throws Exception {
        Method method = loadSculptorMetadataMethod();
        InvocationTargetException error = assertThrows(
                InvocationTargetException.class,
                () -> method.invoke(null, "data/preprocessing/does_not_exist.json"));

        IllegalArgumentException cause = assertInstanceOf(IllegalArgumentException.class, error.getCause());
        assertTrue(cause.getMessage().contains("Missing sculptor metadata resource"));
    }

    @Test
    void loadSculptorMetadataRejectsInvalidShapeMetadata() throws Exception {
        Method method = loadSculptorMetadataMethod();
        InvocationTargetException error = assertThrows(
                InvocationTargetException.class,
                () -> method.invoke(null, "data/preprocessing/cartographer_metadata_minimal.json"));

        IllegalStateException cause = assertInstanceOf(IllegalStateException.class, error.getCause());
        assertTrue(cause.getMessage().contains("Invalid Sculptor max_peptide_len"));
    }

    private static Object invokeLoadCartographerMetadata(String resource) throws Exception {
        Method method = loadCartographerMetadataMethod();
        return method.invoke(null, resource);
    }

    private static Object invokeLoadElectricianMetadata(String resource) throws Exception {
        Method method = loadElectricianMetadataMethod();
        return method.invoke(null, resource);
    }

    private static Object invokeLoadSculptorMetadata(String resource) throws Exception {
        Method method = loadSculptorMetadataMethod();
        return method.invoke(null, resource);
    }

    private static Method loadCartographerMetadataMethod() throws NoSuchMethodException {
        Method method = DefaultChronologerLibraryPredictor.class.getDeclaredMethod(
                "loadCartographerMetadata",
                String.class);
        method.setAccessible(true);
        return method;
    }

    private static Method loadElectricianMetadataMethod() throws NoSuchMethodException {
        Method method = DefaultChronologerLibraryPredictor.class.getDeclaredMethod(
                "loadElectricianMetadata",
                String.class);
        method.setAccessible(true);
        return method;
    }

    private static Method loadSculptorMetadataMethod() throws NoSuchMethodException {
        Method method = DefaultChronologerLibraryPredictor.class.getDeclaredMethod(
                "loadSculptorMetadata",
                String.class);
        method.setAccessible(true);
        return method;
    }

    private static int getIntField(Object target, String fieldName) throws Exception {
        Field field = target.getClass().getDeclaredField(fieldName);
        field.setAccessible(true);
        return field.getInt(target);
    }

    private static double getDoubleField(Object target, String fieldName) throws Exception {
        Field field = target.getClass().getDeclaredField(fieldName);
        field.setAccessible(true);
        return field.getDouble(target);
    }
}

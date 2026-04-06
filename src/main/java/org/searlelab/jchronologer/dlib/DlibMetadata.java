package org.searlelab.jchronologer.dlib;

import java.util.LinkedHashMap;
import java.util.Map;

/**
 * Metadata defaults for generated DLIB files.
 */
public final class DlibMetadata {

    private final Map<String, String> values;

    private DlibMetadata(Map<String, String> values) {
        this.values = Map.copyOf(values);
    }

    public static DlibMetadata defaults() {
        LinkedHashMap<String, String> values = new LinkedHashMap<>();
        values.put("version", "0.1.18");
        values.put("EncyclopediaVersion", "Unknown");
        values.put("hasMatchedDecoySpectra", "false");
        values.put("isCompleteFASTALibrary", "false");
        return new DlibMetadata(values);
    }

    public Map<String, String> getValues() {
        return values;
    }
}

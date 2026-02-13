package org.searlelab.jchronologer.impl;

import org.searlelab.jchronologer.api.Chronologer;
import org.searlelab.jchronologer.api.ChronologerOptions;

public final class ChronologerFactory {

    private ChronologerFactory() {
    }

    public static Chronologer createDefault() {
        return new DefaultChronologer(ChronologerOptions.builder().build());
    }

    public static Chronologer create(ChronologerOptions options) {
        return new DefaultChronologer(options);
    }
}

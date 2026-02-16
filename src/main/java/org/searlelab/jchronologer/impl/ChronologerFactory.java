package org.searlelab.jchronologer.impl;

import org.searlelab.jchronologer.api.Chronologer;
import org.searlelab.jchronologer.api.ChronologerLibraryOptions;
import org.searlelab.jchronologer.api.ChronologerLibraryPredictor;
import org.searlelab.jchronologer.api.ChronologerOptions;

/**
 * Factory for creating Chronologer inference clients.
 *
 * <p>All instances are currently backed by {@link DefaultChronologer}.
 */
public final class ChronologerFactory {

    private ChronologerFactory() {
    }

    /**
     * Creates a Chronologer client using bundled default resources and default batch size.
     *
     * @return ready-to-use Chronologer client
     */
    public static Chronologer createDefault() {
        return new DefaultChronologer(ChronologerOptions.builder().build());
    }

    /**
     * Creates a Chronologer client with explicit runtime options.
     *
     * @param options model/preprocessing resource and batching configuration
     * @return ready-to-use Chronologer client
     */
    public static Chronologer create(ChronologerOptions options) {
        return new DefaultChronologer(options);
    }

    /**
     * Creates a tandem library predictor with default Chronologer + Cartographer resources.
     *
     * @return ready-to-use library predictor
     */
    public static ChronologerLibraryPredictor createLibraryPredictorDefault() {
        return new DefaultChronologerLibraryPredictor(ChronologerLibraryOptions.builder().build());
    }

    /**
     * Creates a tandem library predictor with explicit runtime options.
     *
     * @param options model/preprocessing resources and batching configuration
     * @return ready-to-use library predictor
     */
    public static ChronologerLibraryPredictor createLibraryPredictor(ChronologerLibraryOptions options) {
        return new DefaultChronologerLibraryPredictor(options);
    }
}

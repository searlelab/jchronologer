package org.searlelab.jchronologer.api;

import java.util.List;

/**
 * Tandem inference API that combines Chronologer RT and Cartographer MS2 predictions.
 */
public interface ChronologerLibraryPredictor extends AutoCloseable {

    /**
     * Optional lifecycle hook.
     */
    default void init() {
        // no-op by default
    }

    /**
     * Predicts complete library entries for each request.
     *
     * @param requests library prediction requests
     * @return predicted entries in request order, expanded by precursor conditions
     */
    List<ChronologerLibraryEntry> predict(List<LibraryPredictionRequest> requests);

    @Override
    void close();
}

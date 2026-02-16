package org.searlelab.jchronologer.api;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Input request for tandem library prediction.
 */
public final class LibraryPredictionRequest {

    private final String peptideSequence;
    private final List<PrecursorCondition> precursorConditions;

    public LibraryPredictionRequest(String peptideSequence, List<PrecursorCondition> precursorConditions) {
        this.peptideSequence = peptideSequence;
        this.precursorConditions = Collections.unmodifiableList(new ArrayList<>(precursorConditions));
    }

    public String getPeptideSequence() {
        return peptideSequence;
    }

    public List<PrecursorCondition> getPrecursorConditions() {
        return precursorConditions;
    }
}

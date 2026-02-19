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
    private final double precursorNce;
    private final double minimumChargeProbability;
    private final boolean automaticChargeSelection;

    /**
     * Creates a request with explicit precursor conditions.
     */
    public LibraryPredictionRequest(String peptideSequence, List<PrecursorCondition> precursorConditions) {
        if (peptideSequence == null || peptideSequence.isBlank()) {
            throw new IllegalArgumentException("Peptide sequence must be non-empty.");
        }
        if (precursorConditions == null) {
            throw new IllegalArgumentException("Precursor conditions must be non-null.");
        }
        this.peptideSequence = peptideSequence;
        this.precursorConditions = Collections.unmodifiableList(new ArrayList<>(precursorConditions));
        this.precursorNce = Double.NaN;
        this.minimumChargeProbability = Double.NaN;
        this.automaticChargeSelection = false;
    }

    /**
     * Creates a request that derives precursor charges from Electrician charge probabilities.
     */
    public LibraryPredictionRequest(String peptideSequence, double precursorNce, double minimumChargeProbability) {
        if (peptideSequence == null || peptideSequence.isBlank()) {
            throw new IllegalArgumentException("Peptide sequence must be non-empty.");
        }
        this.peptideSequence = peptideSequence;
        this.precursorConditions = List.of();
        this.precursorNce = precursorNce;
        this.minimumChargeProbability = minimumChargeProbability;
        this.automaticChargeSelection = true;
    }

    public String getPeptideSequence() {
        return peptideSequence;
    }

    public List<PrecursorCondition> getPrecursorConditions() {
        return precursorConditions;
    }

    /**
     * @return precursor NCE for automatic charge-selection requests
     */
    public double getPrecursorNce() {
        return precursorNce;
    }

    /**
     * @return minimum normalized charge probability in [0.0, 1.0] for automatic mode
     */
    public double getMinimumChargeProbability() {
        return minimumChargeProbability;
    }

    public boolean usesExplicitPrecursorConditions() {
        return !automaticChargeSelection;
    }

    public boolean usesAutomaticChargeSelection() {
        return automaticChargeSelection;
    }
}

package org.searlelab.jchronologer.api;

/**
 * Explicit precursor condition (charge and NCE) for one MS2 prediction.
 */
public final class PrecursorCondition {

    private final byte precursorCharge;
    private final double precursorNce;

    public PrecursorCondition(byte precursorCharge, double precursorNce) {
        this.precursorCharge = precursorCharge;
        this.precursorNce = precursorNce;
    }

    public byte getPrecursorCharge() {
        return precursorCharge;
    }

    public double getPrecursorNce() {
        return precursorNce;
    }
}

package org.searlelab.jchronologer.api;

/**
 * Explicit precursor condition (charge and NCE) for one MS2 prediction.
 *
 * <p>NCE (normalized collision energy) should be specified in standard instrument units,
 * typically between 20 and 50 (e.g. 27.0 for NCE 27). The value is automatically
 * normalized internally (divided by 100) before being passed to the Cartographer model.
 */
public final class PrecursorCondition {

    private final byte precursorCharge;
    private final double precursorNce;

    /**
     * @param precursorCharge precursor charge state (e.g. 2 for +2H)
     * @param precursorNce normalized collision energy in standard units (e.g. 27.0, not 0.27)
     */
    public PrecursorCondition(byte precursorCharge, double precursorNce) {
        this.precursorCharge = precursorCharge;
        this.precursorNce = precursorNce;
    }

    public byte getPrecursorCharge() {
        return precursorCharge;
    }

    /**
     * @return NCE in standard instrument units (e.g. 27.0)
     */
    public double getPrecursorNce() {
        return precursorNce;
    }
}

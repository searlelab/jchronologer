package org.searlelab.jchronologer.api;

import java.util.Arrays;
import java.util.Optional;
import org.searlelab.jchronologer.preprocessing.PeptideSequenceConverter;

/**
 * Predicted peptide library row combining RT and MS2 outputs.
 */
public final class ChronologerLibraryEntry {

    private final String unimodPeptideSequence;
    private final byte precursorCharge;
    private final double precursorNce;

    private final double precursorMz;
    private final float retentionTimeInSeconds;
    private final double[] massArray;
    private final float[] intensityArray;
    private final String[] ionTypeArray;
    private final Optional<Float> ccs;

    public ChronologerLibraryEntry(
            String unimodPeptideSequence,
            byte precursorCharge,
            double precursorNce,
            double precursorMz,
            float retentionTimeInSeconds,
            double[] massArray,
            float[] intensityArray,
            String[] ionTypeArray) {
        this(
                unimodPeptideSequence,
                precursorCharge,
                precursorNce,
                precursorMz,
                retentionTimeInSeconds,
                massArray,
                intensityArray,
                ionTypeArray,
                Optional.empty());
    }

    public ChronologerLibraryEntry(
            String unimodPeptideSequence,
            byte precursorCharge,
            double precursorNce,
            double precursorMz,
            float retentionTimeInSeconds,
            double[] massArray,
            float[] intensityArray,
            String[] ionTypeArray,
            Optional<Float> ccs) {
        this.unimodPeptideSequence = unimodPeptideSequence;
        this.precursorCharge = precursorCharge;
        this.precursorNce = precursorNce;
        this.precursorMz = precursorMz;
        this.retentionTimeInSeconds = retentionTimeInSeconds;
        this.massArray = Arrays.copyOf(massArray, massArray.length);
        this.intensityArray = Arrays.copyOf(intensityArray, intensityArray.length);
        this.ionTypeArray = Arrays.copyOf(ionTypeArray, ionTypeArray.length);
        this.ccs = ccs == null ? Optional.empty() : ccs;
    }

    public String getUnimodPeptideSequence() {
        return unimodPeptideSequence;
    }

    public byte getPrecursorCharge() {
        return precursorCharge;
    }

    public double getPrecursorNce() {
        return precursorNce;
    }

    public double getPrecursorMz() {
        return precursorMz;
    }

    public float getRetentionTimeInSeconds() {
        return retentionTimeInSeconds;
    }

    public double[] getMassArray() {
        return Arrays.copyOf(massArray, massArray.length);
    }

    public float[] getIntensityArray() {
        return Arrays.copyOf(intensityArray, intensityArray.length);
    }

    public String[] getIonTypeArray() {
        return Arrays.copyOf(ionTypeArray, ionTypeArray.length);
    }

    public Optional<Float> getCCS() {
        return ccs;
    }

    /**
     * Converts this entry back to legacy mass-encoded PeptideModSeq style.
     */
    public String getPeptideModSeq() {
        return PeptideSequenceConverter.unimodToMassEncoded(unimodPeptideSequence);
    }
}

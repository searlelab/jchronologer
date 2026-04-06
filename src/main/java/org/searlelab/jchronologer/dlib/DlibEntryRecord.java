package org.searlelab.jchronologer.dlib;

import java.util.Arrays;
import java.util.Optional;

/**
 * Immutable row model for the DLIB entries table.
 */
public final class DlibEntryRecord {

    private final double precursorMz;
    private final int precursorCharge;
    private final String peptideModSeq;
    private final String peptideSeq;
    private final int copies;
    private final float rtInSeconds;
    private final float score;
    private final double[] massArray;
    private final float[] intensityArray;
    private final Optional<Float> ionMobility;
    private final String sourceFile;

    public DlibEntryRecord(
            double precursorMz,
            int precursorCharge,
            String peptideModSeq,
            String peptideSeq,
            int copies,
            float rtInSeconds,
            float score,
            double[] massArray,
            float[] intensityArray,
            Optional<Float> ionMobility,
            String sourceFile) {
        this.precursorMz = precursorMz;
        this.precursorCharge = precursorCharge;
        this.peptideModSeq = peptideModSeq;
        this.peptideSeq = peptideSeq;
        this.copies = copies;
        this.rtInSeconds = rtInSeconds;
        this.score = score;
        this.massArray = Arrays.copyOf(massArray, massArray.length);
        this.intensityArray = Arrays.copyOf(intensityArray, intensityArray.length);
        this.ionMobility = ionMobility == null ? Optional.empty() : ionMobility;
        this.sourceFile = sourceFile;
    }

    public double getPrecursorMz() {
        return precursorMz;
    }

    public int getPrecursorCharge() {
        return precursorCharge;
    }

    public String getPeptideModSeq() {
        return peptideModSeq;
    }

    public String getPeptideSeq() {
        return peptideSeq;
    }

    public int getCopies() {
        return copies;
    }

    public float getRtInSeconds() {
        return rtInSeconds;
    }

    public float getScore() {
        return score;
    }

    public double[] getMassArray() {
        return Arrays.copyOf(massArray, massArray.length);
    }

    public float[] getIntensityArray() {
        return Arrays.copyOf(intensityArray, intensityArray.length);
    }

    public Optional<Float> getIonMobility() {
        return ionMobility;
    }

    public String getSourceFile() {
        return sourceFile;
    }
}

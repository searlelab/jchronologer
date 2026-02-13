package org.searlelab.jchronologer.api;

import java.util.Arrays;

public final class AcceptedPrediction {

    private final int rowIndex;
    private final String peptideModSeq;
    private final String patchedPeptideModSeq;
    private final String codedPeptideSeq;
    private final long[] tokenArray;
    private final float predHi;

    public AcceptedPrediction(
            int rowIndex,
            String peptideModSeq,
            String patchedPeptideModSeq,
            String codedPeptideSeq,
            long[] tokenArray,
            float predHi) {
        this.rowIndex = rowIndex;
        this.peptideModSeq = peptideModSeq;
        this.patchedPeptideModSeq = patchedPeptideModSeq;
        this.codedPeptideSeq = codedPeptideSeq;
        this.tokenArray = Arrays.copyOf(tokenArray, tokenArray.length);
        this.predHi = predHi;
    }

    public int getRowIndex() {
        return rowIndex;
    }

    public String getPeptideModSeq() {
        return peptideModSeq;
    }

    public String getPatchedPeptideModSeq() {
        return patchedPeptideModSeq;
    }

    public String getCodedPeptideSeq() {
        return codedPeptideSeq;
    }

    public long[] getTokenArray() {
        return Arrays.copyOf(tokenArray, tokenArray.length);
    }

    public float getPredHi() {
        return predHi;
    }
}

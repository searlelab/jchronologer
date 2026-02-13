package org.searlelab.jchronologer.api;

/**
 * Immutable rejection record for one input row that could not be scored by Chronologer.
 *
 * <p>Contains the normalized sequence attempted during preprocessing, a coarse rejection category,
 * and optional detailed error text for diagnostics.
 */
public final class RejectedPrediction {

    private final int rowIndex;
    private final String peptideModSeq;
    private final String patchedPeptideModSeq;
    private final RejectionReason rejectionReason;
    private final String errorDetail;

    public RejectedPrediction(
            int rowIndex,
            String peptideModSeq,
            String patchedPeptideModSeq,
            RejectionReason rejectionReason,
            String errorDetail) {
        this.rowIndex = rowIndex;
        this.peptideModSeq = peptideModSeq;
        this.patchedPeptideModSeq = patchedPeptideModSeq;
        this.rejectionReason = rejectionReason;
        this.errorDetail = errorDetail;
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

    public RejectionReason getRejectionReason() {
        return rejectionReason;
    }

    public String getErrorDetail() {
        return errorDetail;
    }
}

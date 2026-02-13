package org.searlelab.jchronologer.preprocessing;

import org.searlelab.jchronologer.api.RejectionReason;

public final class PreprocessingOutcome {

    private final boolean accepted;
    private final String patchedPeptideModSeq;
    private final String codedPeptideSeq;
    private final long[] tokenArray;
    private final RejectionReason rejectionReason;
    private final String errorDetail;

    private PreprocessingOutcome(
            boolean accepted,
            String patchedPeptideModSeq,
            String codedPeptideSeq,
            long[] tokenArray,
            RejectionReason rejectionReason,
            String errorDetail) {
        this.accepted = accepted;
        this.patchedPeptideModSeq = patchedPeptideModSeq;
        this.codedPeptideSeq = codedPeptideSeq;
        this.tokenArray = tokenArray;
        this.rejectionReason = rejectionReason;
        this.errorDetail = errorDetail;
    }

    public static PreprocessingOutcome accepted(String patched, String coded, long[] tokens) {
        return new PreprocessingOutcome(true, patched, coded, tokens, null, null);
    }

    public static PreprocessingOutcome rejected(String patched, RejectionReason reason, String errorDetail) {
        return new PreprocessingOutcome(false, patched, null, null, reason, errorDetail);
    }

    public boolean isAccepted() {
        return accepted;
    }

    public String getPatchedPeptideModSeq() {
        return patchedPeptideModSeq;
    }

    public String getCodedPeptideSeq() {
        return codedPeptideSeq;
    }

    public long[] getTokenArray() {
        return tokenArray;
    }

    public RejectionReason getRejectionReason() {
        return rejectionReason;
    }

    public String getErrorDetail() {
        return errorDetail;
    }
}

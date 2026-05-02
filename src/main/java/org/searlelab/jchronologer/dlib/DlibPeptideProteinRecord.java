package org.searlelab.jchronologer.dlib;

/**
 * Immutable row model for the DLIB peptide-to-protein mapping table.
 */
public final class DlibPeptideProteinRecord {

    private final String peptideSeq;
    private final boolean isDecoy;
    private final String proteinAccession;

    public DlibPeptideProteinRecord(String peptideSeq, boolean isDecoy, String proteinAccession) {
        this.peptideSeq = peptideSeq;
        this.isDecoy = isDecoy;
        this.proteinAccession = proteinAccession;
    }

    public String peptideSeq() {
        return peptideSeq;
    }

    public boolean isDecoy() {
        return isDecoy;
    }

    public String proteinAccession() {
        return proteinAccession;
    }
}

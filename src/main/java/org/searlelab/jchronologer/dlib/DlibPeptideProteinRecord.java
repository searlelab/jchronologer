package org.searlelab.jchronologer.dlib;

/**
 * Immutable row model for the DLIB peptide-to-protein mapping table.
 */
public record DlibPeptideProteinRecord(String peptideSeq, boolean isDecoy, String proteinAccession) {
}

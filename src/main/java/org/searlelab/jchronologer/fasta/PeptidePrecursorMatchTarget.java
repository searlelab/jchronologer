package org.searlelab.jchronologer.fasta;

import java.util.LinkedHashSet;
import java.util.Set;

/**
 * Trie payload used for peptide-to-accession mapping.
 */
public final class PeptidePrecursorMatchTarget {

    private final String peptideSeq;
    private final LinkedHashSet<String> proteinAccessions;

    public PeptidePrecursorMatchTarget(String peptideSeq) {
        this.peptideSeq = peptideSeq;
        this.proteinAccessions = new LinkedHashSet<>();
    }

    public String getPeptideSeq() {
        return peptideSeq;
    }

    public Set<String> getProteinAccessions() {
        return proteinAccessions;
    }

    public void addProteinAccession(String proteinAccession) {
        proteinAccessions.add(proteinAccession);
    }
}

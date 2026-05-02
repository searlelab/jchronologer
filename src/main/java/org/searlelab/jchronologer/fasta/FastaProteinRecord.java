package org.searlelab.jchronologer.fasta;

/**
 * Immutable FASTA protein record.
 */
public final class FastaProteinRecord {

    private final String accession;
    private final String sequence;

    public FastaProteinRecord(String accession, String sequence) {
        this.accession = accession;
        this.sequence = sequence;
    }

    public String accession() {
        return accession;
    }

    public String sequence() {
        return sequence;
    }
}

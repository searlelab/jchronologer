package org.searlelab.jchronologer.fasta;

/**
 * Immutable FASTA protein record.
 */
public record FastaProteinRecord(String accession, String sequence) {
}

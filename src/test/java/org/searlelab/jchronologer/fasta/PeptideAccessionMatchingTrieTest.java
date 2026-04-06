package org.searlelab.jchronologer.fasta;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.util.List;
import java.util.Map;
import java.util.Set;
import org.junit.jupiter.api.Test;

class PeptideAccessionMatchingTrieTest {

    @Test
    void matchesExactSubstringsAcrossProteins() {
        PeptideAccessionMatchingTrie trie = new PeptideAccessionMatchingTrie(List.of("PEPTIDE", "TASEFDSAIAQDK", "MISS"));

        Map<String, PeptidePrecursorMatchTarget> matches = trie.match(List.of(
                new FastaProteinRecord("P1", "XXPEPTIDEXX"),
                new FastaProteinRecord("P2", "MTASEFDSAIAQDKYY"),
                new FastaProteinRecord("P3", "PEPTIDE")));

        assertEquals(Set.of("P1", "P3"), matches.get("PEPTIDE").getProteinAccessions());
        assertEquals(Set.of("P2"), matches.get("TASEFDSAIAQDK").getProteinAccessions());
        assertEquals(Set.of(), matches.get("MISS").getProteinAccessions());
    }
}

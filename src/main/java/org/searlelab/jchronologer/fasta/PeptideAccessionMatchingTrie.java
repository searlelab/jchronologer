package org.searlelab.jchronologer.fasta;

import java.util.Collection;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 * Reverse trie for exact peptide-to-protein substring matching.
 */
public final class PeptideAccessionMatchingTrie {

    private final TrieNode head;
    private final Map<String, PeptidePrecursorMatchTarget> targetByPeptide;

    public PeptideAccessionMatchingTrie(Collection<String> peptideSequences) {
        this.head = new TrieNode();
        this.targetByPeptide = new LinkedHashMap<>();
        for (String peptideSequence : peptideSequences) {
            PeptidePrecursorMatchTarget target = new PeptidePrecursorMatchTarget(peptideSequence);
            targetByPeptide.put(peptideSequence, target);
            TrieNode node = head;
            for (int i = peptideSequence.length() - 1; i >= 0; i--) {
                node = node.getOrCreate(peptideSequence.charAt(i));
            }
            node.addTarget(target);
        }
    }

    public Map<String, PeptidePrecursorMatchTarget> match(Iterable<FastaProteinRecord> proteins) {
        for (FastaProteinRecord protein : proteins) {
            addProtein(protein);
        }
        return targetByPeptide;
    }

    private void addProtein(FastaProteinRecord protein) {
        char[] sequence = protein.sequence().toCharArray();
        for (int i = sequence.length - 1; i >= 0; i--) {
            TrieNode node = head.get(sequence[i]);
            if (node == null) {
                continue;
            }
            node.addMatches(protein.accession());
            for (int j = i - 1; j >= 0; j--) {
                node = node.get(sequence[j]);
                if (node == null) {
                    break;
                }
                node.addMatches(protein.accession());
            }
        }
    }

    private static final class TrieNode {
        private final Map<Character, TrieNode> children;
        private final java.util.List<PeptidePrecursorMatchTarget> targets;

        private TrieNode() {
            this.children = new LinkedHashMap<>();
            this.targets = new java.util.ArrayList<>();
        }

        private TrieNode get(char residue) {
            return children.get(residue);
        }

        private TrieNode getOrCreate(char residue) {
            return children.computeIfAbsent(residue, key -> new TrieNode());
        }

        private void addTarget(PeptidePrecursorMatchTarget target) {
            targets.add(target);
        }

        private void addMatches(String accession) {
            for (PeptidePrecursorMatchTarget target : targets) {
                target.addProteinAccession(accession);
            }
        }
    }
}

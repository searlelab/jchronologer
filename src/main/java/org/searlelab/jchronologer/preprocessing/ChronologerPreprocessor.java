package org.searlelab.jchronologer.preprocessing;

import java.util.Map;
import org.searlelab.jchronologer.api.RejectionReason;
import org.searlelab.jchronologer.preprocessing.PreprocessingMetadataLoader.CompiledPreprocessingMetadata;
import org.searlelab.jchronologer.preprocessing.PreprocessingMetadataLoader.CompiledPreprocessingMetadata.CompiledRegexRule;

public final class ChronologerPreprocessor {

    private static final int MIN_PEPTIDE_LEN = 6;
    private static final String E_CYCLO_PATCH = "E[-18.0105647]";
    private static final String C_CYCLO_PATCH = "C[+39.99491463]";

    private final CompiledPreprocessingMetadata metadata;

    public ChronologerPreprocessor(CompiledPreprocessingMetadata metadata) {
        this.metadata = metadata;
    }

    public PreprocessingOutcome preprocess(String peptideModSeq) {
        if (peptideModSeq == null || peptideModSeq.isBlank()) {
            return PreprocessingOutcome.rejected(peptideModSeq, RejectionReason.TOKENIZATION_ERROR,
                    "Input peptide is blank.");
        }

        String patched = patchModSequence(peptideModSeq);
        final String coded;
        try {
            coded = modseqToCodedseq(patched);
        } catch (RuntimeException e) {
            return PreprocessingOutcome.rejected(
                    patched,
                    RejectionReason.TOKENIZATION_ERROR,
                    e.getClass().getSimpleName() + ": " + e.getMessage());
        }

        if (coded == null) {
            return PreprocessingOutcome.rejected(patched, RejectionReason.UNSUPPORTED_MODIFICATION, null);
        }

        int peptideLen = coded.length() - 2;
        if (peptideLen < MIN_PEPTIDE_LEN || peptideLen > metadata.getMaxPeptideLen()) {
            return PreprocessingOutcome.rejected(patched, RejectionReason.LENGTH_OUT_OF_RANGE, null);
        }

        long[] tokenArray = new long[metadata.tokenArrayLength()];
        Map<String, Integer> aaToInt = metadata.getAaToInt();
        for (int i = 0; i < coded.length(); i++) {
            String residue = String.valueOf(coded.charAt(i));
            Integer value = aaToInt.get(residue);
            if (value == null) {
                return PreprocessingOutcome.rejected(
                        patched,
                        RejectionReason.INVALID_RESIDUE_TOKEN,
                        "Unknown residue token: " + residue);
            }
            tokenArray[i] = value;
        }

        return PreprocessingOutcome.accepted(patched, coded, tokenArray);
    }

    static String patchModSequence(String peptide) {
        String patched = peptide
                .replace("E[-18.0]", E_CYCLO_PATCH)
                .replace("C[-17.0]", C_CYCLO_PATCH);

        if (patched.length() >= 13 && patched.substring(1, 13).equals("[+42.010565]")) {
            patched = "[+42.010565]" + patched.charAt(0) + patched.substring(13);
        }
        return patched;
    }

    private String modseqToCodedseq(String peptideModSeq) {
        String sequence = peptideModSeq;
        for (CompiledRegexRule rule : metadata.getModRegexRules()) {
            sequence = rule.getPattern().matcher(sequence).replaceAll(rule.getToken());
        }
        if (sequence.isEmpty()) {
            throw new IllegalArgumentException("Peptide sequence is empty after applying modification rules.");
        }

        if (sequence.charAt(0) == 'd') {
            sequence = ")" + sequence;
        } else if (sequence.charAt(0) == 'e') {
            sequence = "(" + sequence;
        } else if (sequence.charAt(0) == '[') {
            int end = sequence.indexOf(']');
            if (end < 0) {
                throw new IllegalArgumentException("Missing closing bracket in N-term annotation.");
            }
            if (end < 7) {
                throw new IllegalArgumentException("N-term annotation is too short: " + sequence.substring(0, end + 1));
            }
            String key = sequence.substring(1, 7);
            String ntermToken = metadata.getNtermKeys().get(key);
            if (ntermToken == null) {
                throw new IllegalArgumentException("Unknown N-term key: " + key);
            }
            sequence = ntermToken + sequence.substring(end + 1);
        } else {
            sequence = "-" + sequence;
        }
        sequence = sequence + "_";

        return sequence.indexOf('[') >= 0 ? null : sequence;
    }
}

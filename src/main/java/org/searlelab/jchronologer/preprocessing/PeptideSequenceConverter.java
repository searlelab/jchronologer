package org.searlelab.jchronologer.preprocessing;

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;

/**
 * Converts between legacy mass-encoded peptide strings and terminal-aware UNIMOD strings.
 */
public final class PeptideSequenceConverter {

    private static final String TRIMETHYLATION="UNIMOD:37";
	private static final String DIMETHYLATION="UNIMOD:36";
	private static final String SUCCINYLATION="UNIMOD:64";
	private static final String DEAMIDATION="UNIMOD:7";
	private static final String HEXNAC="UNIMOD:43";
	private static final String PYROGLUE="UNIMOD:27";
	private static final String PYROGLU="UNIMOD:28";
	private static final String METHYLATION="UNIMOD:34";
	private static final String ACETYLATION="UNIMOD:1";
	private static final String GLYGLY="UNIMOD:121";
	private static final String PHOSPHO="UNIMOD:21";
	private static final String OXIDATION="UNIMOD:35";
	private static final String CARBAMIDOMETHYL="UNIMOD:4";
	private static final String TMT0PLEX="UNIMOD:739";
	private static final String TMT6PLEX="UNIMOD:737";
	private static final Map<String, Double> UNIMOD_MASS_BY_ID = createUnimodMassMap();
    private static final Set<String> NTERM_COMPATIBLE_UNIMODS = Set.of(
    		TMT6PLEX,
            ACETYLATION,
            PYROGLU,
            PYROGLUE,
            TMT0PLEX);

    private PeptideSequenceConverter() {
    }

    public static String normalizeToUnimod(String peptideSequence, double epsilon) {
        if (peptideSequence == null || peptideSequence.isBlank()) {
            throw new IllegalArgumentException("Peptide sequence must be non-empty.");
        }
        ParsedUnimodSequence parsed = looksLikeTerminalAwareUnimod(peptideSequence)
                ? parseUnimod(peptideSequence)
                : parseMassEncoded(peptideSequence, epsilon);
        return formatUnimod(parsed);
    }

    public static ParsedUnimodSequence parseNormalizedUnimod(String unimodSequence) {
        return parseUnimod(unimodSequence);
    }

    public static String unimodToMassEncoded(String unimodSequence) {
        ParsedUnimodSequence parsed = parseUnimod(unimodSequence);
        String residues = parsed.getResidues();
        List<String> ntermMods = new ArrayList<>(parsed.getPositionMods().get(0));
        List<List<String>> residueMods = new ArrayList<>(residues.length());
        for (int i = 0; i < residues.length(); i++) {
            residueMods.add(new ArrayList<>(parsed.getPositionMods().get(i + 1)));
        }
        foldNtermPyroCarbamidomethylCysteine(residues, ntermMods, residueMods);

        StringBuilder builder = new StringBuilder();

        if (!ntermMods.isEmpty()) {
            builder.append('[')
                    .append(formatSignedMass(sumUnimodMass(ntermMods)))
                    .append(']');
        }

        for (int i = 0; i < residues.length(); i++) {
            builder.append(residues.charAt(i));
            List<String> mods = residueMods.get(i);
            if (!mods.isEmpty()) {
                builder.append('[')
                        .append(formatSignedMass(sumUnimodMass(mods)))
                        .append(']');
            }
        }

        List<String> ctermMods = parsed.getPositionMods().get(residues.length() + 1);
        if (!ctermMods.isEmpty()) {
            builder.append('[')
                    .append(formatSignedMass(sumUnimodMass(ctermMods)))
                    .append(']');
        }

        return builder.toString();
    }

    public static double sumUnimodMass(List<String> unimodTokens) {
        double sum = 0.0;
        for (String token : unimodTokens) {
            Double mass = UNIMOD_MASS_BY_ID.get(token);
            if (mass == null) {
                throw new IllegalArgumentException("Unsupported UNIMOD token: " + token);
            }
            sum += mass;
        }
        return sum;
    }

    /**
     * Folds first-residue pyroglu modifications onto the N-terminus.
     *
     * <p>Some datasets report pyro-glutamate as a modification on the first residue
     * (e.g. {@code []-Q[UNIMOD:28]PEPTIDE-[]}) rather than on the N-terminus
     * ({@code [UNIMOD:28]-QPEPTIDE-[]}). This method normalizes the first form
     * to the second, which is required for correct Cartographer tokenization.
     */
    public static void foldFirstResiduePyrogluToNterm(
            String residues,
            List<String> ntermMods,
            List<List<String>> residueMods) {
        if (residues.isEmpty() || !ntermMods.isEmpty()) {
            return;
        }
        char firstResidue = residues.charAt(0);
        List<String> firstMods = residueMods.get(0);
        if (firstResidue == 'Q' && firstMods.contains(PYROGLU)) {
            removeOne(firstMods, PYROGLU);
            ntermMods.add(PYROGLU);
        } else if (firstResidue == 'E' && firstMods.contains(PYROGLUE)) {
            removeOne(firstMods, PYROGLUE);
            ntermMods.add(PYROGLUE);
        }
    }

    private static void foldNtermPyroCarbamidomethylCysteine(
            String residues,
            List<String> ntermMods,
            List<List<String>> residueMods) {
        if (residues.isEmpty() || residues.charAt(0) != 'C') {
            return;
        }

        List<String> firstResidueMods = residueMods.get(0);
        boolean removedNterm = removeOne(ntermMods, PYROGLU);
        if (removedNterm && firstResidueMods.contains(CARBAMIDOMETHYL)) {
            firstResidueMods.add(PYROGLU);
            return;
        }
        if (removedNterm) {
            ntermMods.add(PYROGLU);
        }
    }

    private static boolean removeOne(List<String> tokens, String token) {
        int idx = tokens.indexOf(token);
        if (idx < 0) {
            return false;
        }
        tokens.remove(idx);
        return true;
    }

    private static boolean looksLikeTerminalAwareUnimod(String sequence) {
        return sequence.startsWith("[") && sequence.contains("-") && sequence.lastIndexOf("-") > 0;
    }

    private static ParsedUnimodSequence parseUnimod(String sequence) {
        int idx = 0;

        TokenParse nterm = parseBracketTokenBlock(sequence, idx, true);
        idx = nterm.nextIndex;
        if (idx >= sequence.length() || sequence.charAt(idx) != '-') {
            throw new IllegalArgumentException("UNIMOD sequence is missing '-' after N-term block: " + sequence);
        }
        idx++;

        StringBuilder residues = new StringBuilder();
        List<List<String>> residueMods = new ArrayList<>();
        while (idx < sequence.length() && sequence.charAt(idx) != '-') {
            char residue = sequence.charAt(idx);
            if (residue < 'A' || residue > 'Z') {
                throw new IllegalArgumentException("Invalid residue '" + residue + "' in sequence: " + sequence);
            }
            residues.append(residue);
            idx++;

            if (idx < sequence.length() && sequence.charAt(idx) == '[') {
                TokenParse residueTokenParse = parseBracketTokenBlock(sequence, idx, false);
                residueMods.add(residueTokenParse.tokens);
                idx = residueTokenParse.nextIndex;
            } else {
                residueMods.add(List.of());
            }
        }

        if (idx >= sequence.length() || sequence.charAt(idx) != '-') {
            throw new IllegalArgumentException("UNIMOD sequence is missing '-' before C-term block: " + sequence);
        }
        idx++;

        TokenParse cterm = parseBracketTokenBlock(sequence, idx, true);
        idx = cterm.nextIndex;
        if (idx != sequence.length()) {
            throw new IllegalArgumentException("Unexpected trailing characters in sequence: " + sequence);
        }

        List<List<String>> positionMods = new ArrayList<>(residueMods.size() + 2);
        positionMods.add(nterm.tokens);
        positionMods.addAll(residueMods);
        positionMods.add(cterm.tokens);
        return new ParsedUnimodSequence(residues.toString(), positionMods);
    }

    private static ParsedUnimodSequence parseMassEncoded(String sequence, double epsilon) {
        int idx = 0;
        List<String> ntermMods = new ArrayList<>();

        while (idx < sequence.length() && sequence.charAt(idx) == '[') {
            SingleToken token = parseSingleBracketToken(sequence, idx);
            ntermMods.addAll(resolveMassToUnimod(token.numericMass(sequence), epsilon, sequence));
            idx = token.nextIndex;
        }

        StringBuilder residues = new StringBuilder();
        List<List<String>> residueMods = new ArrayList<>();
        while (idx < sequence.length() && isResidue(sequence.charAt(idx))) {
            residues.append(sequence.charAt(idx));
            idx++;

            List<String> currentMods = new ArrayList<>();
            while (idx < sequence.length() && sequence.charAt(idx) == '[') {
                SingleToken token = parseSingleBracketToken(sequence, idx);
                List<String> mapped = resolveMassToUnimod(token.numericMass(sequence), epsilon, sequence);

                if (residues.length() == 1
                        && ntermMods.isEmpty()
                        && mapped.size() == 2
                        && mapped.get(0).equals(mapped.get(1))
                        && NTERM_COMPATIBLE_UNIMODS.contains(mapped.get(0))) {
                    ntermMods.add(mapped.get(0));
                    currentMods.add(mapped.get(1));
                } else {
                    currentMods.addAll(mapped);
                }
                idx = token.nextIndex;
            }
            residueMods.add(canonicalizeTokens(currentMods));
        }

        if (residues.length() == 0) {
            throw new IllegalArgumentException("Mass-encoded sequence has no residues: " + sequence);
        }

        List<String> ctermMods = new ArrayList<>();
        while (idx < sequence.length() && sequence.charAt(idx) == '[') {
            SingleToken token = parseSingleBracketToken(sequence, idx);
            ctermMods.addAll(resolveMassToUnimod(token.numericMass(sequence), epsilon, sequence));
            idx = token.nextIndex;
        }

        if (idx != sequence.length()) {
            throw new IllegalArgumentException("Invalid mass-encoded sequence (unexpected token near index " + idx + "): " + sequence);
        }

        List<List<String>> positionMods = new ArrayList<>(residueMods.size() + 2);
        positionMods.add(canonicalizeTokens(ntermMods));
        positionMods.addAll(residueMods);
        positionMods.add(canonicalizeTokens(ctermMods));
        return new ParsedUnimodSequence(residues.toString(), positionMods);
    }

    private static String formatUnimod(ParsedUnimodSequence parsed) {
        StringBuilder builder = new StringBuilder();
        appendTerminalBlock(builder, parsed.getPositionMods().get(0));
        builder.append('-');

        String residues = parsed.getResidues();
        for (int i = 0; i < residues.length(); i++) {
            builder.append(residues.charAt(i));
            appendInlineMods(builder, parsed.getPositionMods().get(i + 1));
        }

        builder.append('-');
        appendTerminalBlock(builder, parsed.getPositionMods().get(residues.length() + 1));
        return builder.toString();
    }

    private static void appendTerminalBlock(StringBuilder builder, List<String> tokens) {
        if (tokens.isEmpty()) {
            builder.append("[]");
            return;
        }
        for (String token : tokens) {
            builder.append('[').append(token).append(']');
        }
    }

    private static void appendInlineMods(StringBuilder builder, List<String> tokens) {
        for (String token : tokens) {
            builder.append('[').append(token).append(']');
        }
    }

    private static TokenParse parseBracketTokenBlock(String sequence, int startIndex, boolean allowEmpty) {
        int idx = startIndex;
        if (idx >= sequence.length() || sequence.charAt(idx) != '[') {
            throw new IllegalArgumentException("Expected bracket block at index " + startIndex + " in sequence: " + sequence);
        }

        List<String> tokens = new ArrayList<>();
        while (idx < sequence.length() && sequence.charAt(idx) == '[') {
            SingleToken single = parseSingleBracketToken(sequence, idx);
            if (!single.token.isEmpty()) {
                tokens.add(normalizeUnimodToken(single.token, sequence));
            }
            idx = single.nextIndex;
            if (!allowEmpty && (idx >= sequence.length() || sequence.charAt(idx) != '[')) {
                break;
            }
        }
        return new TokenParse(canonicalizeTokens(tokens), idx);
    }

    private static SingleToken parseSingleBracketToken(String sequence, int startIndex) {
        int close = sequence.indexOf(']', startIndex);
        if (close < 0) {
            throw new IllegalArgumentException("Missing closing bracket in sequence: " + sequence);
        }
        String token = sequence.substring(startIndex + 1, close).trim();
        return new SingleToken(token, close + 1);
    }

    private static String normalizeUnimodToken(String token, String sequence) {
        if (!token.startsWith("UNIMOD:")) {
            throw new IllegalArgumentException("Unsupported modification token '" + token + "' in sequence: " + sequence);
        }
        if (!UNIMOD_MASS_BY_ID.containsKey(token)) {
            throw new IllegalArgumentException("Unknown UNIMOD token '" + token + "' in sequence: " + sequence);
        }
        return token;
    }

    private static List<String> resolveMassToUnimod(double mass, double epsilon, String sequence) {
        List<String> singleMatches = new ArrayList<>();
        for (Map.Entry<String, Double> entry : UNIMOD_MASS_BY_ID.entrySet()) {
            if (Math.abs(entry.getValue() - mass) <= epsilon) {
                singleMatches.add(entry.getKey());
            }
        }

        if (singleMatches.size() == 1) {
            return List.of(singleMatches.get(0));
        }
        if (singleMatches.size() > 1) {
            throw new IllegalArgumentException(
                    "Ambiguous mass " + mass + " in sequence " + sequence + " mapped to multiple UNIMOD ids: " + singleMatches);
        }

        List<Map.Entry<String, Double>> mods = new ArrayList<>(UNIMOD_MASS_BY_ID.entrySet());
        Set<String> pairMatches = new LinkedHashSet<>();
        for (int i = 0; i < mods.size(); i++) {
            for (int j = i; j < mods.size(); j++) {
                double sum = mods.get(i).getValue() + mods.get(j).getValue();
                if (Math.abs(sum - mass) <= epsilon) {
                    List<String> pair = new ArrayList<>();
                    pair.add(mods.get(i).getKey());
                    pair.add(mods.get(j).getKey());
                    Collections.sort(pair);
                    pairMatches.add(String.join("|", pair));
                }
            }
        }

        if (pairMatches.size() == 1) {
            String[] parts = pairMatches.iterator().next().split("\\|");
            return List.of(parts[0], parts[1]);
        }
        if (!pairMatches.isEmpty()) {
            throw new IllegalArgumentException(
                    "Ambiguous compound mass " + mass + " in sequence " + sequence + " mapped to multiple UNIMOD pairs: " + pairMatches);
        }

        throw new IllegalArgumentException(
                "Mass " + mass + " in sequence " + sequence + " does not map to known UNIMOD ids within epsilon=" + epsilon);
    }

    private static List<String> canonicalizeTokens(List<String> tokens) {
        if (tokens.isEmpty()) {
            return List.of();
        }
        List<String> copy = new ArrayList<>(tokens);
        Collections.sort(copy);
        return List.copyOf(copy);
    }

    private static boolean isResidue(char c) {
        return c >= 'A' && c <= 'Z';
    }

    private static String formatSignedMass(double mass) {
        double rounded = Math.abs(mass) < 5e-7 ? 0.0 : mass;
        return String.format(Locale.US, "%+.6f", rounded);
    }

    private static Map<String, Double> createUnimodMassMap() {
        Map<String, Double> values = new LinkedHashMap<>();
        values.put(TMT6PLEX, 229.162932);
        values.put(OXIDATION, 15.994915);
        values.put(CARBAMIDOMETHYL, 57.021464);
        values.put(PHOSPHO, 79.966331);
        values.put(GLYGLY, 114.042927);
        values.put(ACETYLATION, 42.010565);
        values.put(METHYLATION, 14.015650);
        values.put(PYROGLU, -17.026549);
        values.put(HEXNAC, 203.079373);
        values.put(DEAMIDATION, 0.984016);
        values.put(PYROGLUE, -18.010565);
        values.put(SUCCINYLATION, 100.016044);
        values.put(DIMETHYLATION, 28.0313);
        values.put(TRIMETHYLATION, 42.04695);
        values.put(TMT0PLEX, 224.152478);
        return Collections.unmodifiableMap(values);
    }

    public static final class ParsedUnimodSequence {
        private final String residues;
        private final List<List<String>> positionMods;

        private ParsedUnimodSequence(String residues, List<List<String>> positionMods) {
            this.residues = residues;
            this.positionMods = Collections.unmodifiableList(new ArrayList<>(positionMods));
        }

        public String getResidues() {
            return residues;
        }

        public List<List<String>> getPositionMods() {
            return positionMods;
        }
    }

    private static final class SingleToken {
        private final String token;
        private final int nextIndex;

        private SingleToken(String token, int nextIndex) {
            this.token = token;
            this.nextIndex = nextIndex;
        }

        private double numericMass(String sequence) {
            String normalizedToken = token.replaceAll("\\s+", "");
            try {
                return Double.parseDouble(normalizedToken);
            } catch (NumberFormatException e) {
                throw new IllegalArgumentException(
                        "Expected numeric mass token but found '" + token + "' in sequence: " + sequence,
                        e);
            }
        }
    }

    private static final class TokenParse {
        private final List<String> tokens;
        private final int nextIndex;

        private TokenParse(List<String> tokens, int nextIndex) {
            this.tokens = tokens;
            this.nextIndex = nextIndex;
        }
    }
}

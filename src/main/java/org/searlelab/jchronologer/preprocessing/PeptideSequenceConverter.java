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

    private static final String UNIMOD_737 = "UNIMOD:737";
    private static final String UNIMOD_28 = "UNIMOD:28";
    private static final String UNIMOD_4 = "UNIMOD:4";
    private static final Map<String, Double> UNIMOD_MASS_BY_ID = createUnimodMassMap();
    private static final Set<String> NTERM_COMPATIBLE_UNIMODS = Set.of(
            UNIMOD_737,
            "UNIMOD:1",
            UNIMOD_28,
            "UNIMOD:27");

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

    private static void foldNtermPyroCarbamidomethylCysteine(
            String residues,
            List<String> ntermMods,
            List<List<String>> residueMods) {
        if (residues.isEmpty() || residues.charAt(0) != 'C') {
            return;
        }

        List<String> firstResidueMods = residueMods.get(0);
        boolean removedNterm = removeOne(ntermMods, UNIMOD_28);
        if (removedNterm && firstResidueMods.contains(UNIMOD_4)) {
            firstResidueMods.add(UNIMOD_28);
            return;
        }
        if (removedNterm) {
            ntermMods.add(UNIMOD_28);
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

        if (residues.isEmpty()) {
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
        values.put(UNIMOD_737, 229.162932);
        values.put("UNIMOD:35", 15.994915);
        values.put("UNIMOD:4", 57.021464);
        values.put("UNIMOD:21", 79.966331);
        values.put("UNIMOD:121", 114.042927);
        values.put("UNIMOD:1", 42.010565);
        values.put("UNIMOD:34", 14.015650);
        values.put("UNIMOD:28", -17.026549);
        values.put("UNIMOD:43", 203.079373);
        values.put("UNIMOD:7", 0.984016);
        values.put("UNIMOD:27", -18.010565);
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

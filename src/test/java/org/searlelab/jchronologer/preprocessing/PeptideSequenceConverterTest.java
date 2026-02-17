package org.searlelab.jchronologer.preprocessing;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.List;
import org.junit.jupiter.api.Test;

class PeptideSequenceConverterTest {

    @Test
    void normalizeRejectsNullOrBlankInput() {
        assertThrows(IllegalArgumentException.class, () -> PeptideSequenceConverter.normalizeToUnimod(null, 1e-5));
        assertThrows(IllegalArgumentException.class, () -> PeptideSequenceConverter.normalizeToUnimod("   ", 1e-5));
    }

    @Test
    void normalizeUnimodCanonicalizesTokenOrdering() {
        String normalized = PeptideSequenceConverter.normalizeToUnimod(
                "[UNIMOD:737][UNIMOD:1]-AK[UNIMOD:34][UNIMOD:121]-[]",
                1e-5);

        assertEquals("[UNIMOD:1][UNIMOD:737]-AK[UNIMOD:121][UNIMOD:34]-[]", normalized);
    }

    @Test
    void normalizeMassEncodedHandlesCompoundFirstResidueTmt() {
        String normalized = PeptideSequenceConverter.normalizeToUnimod(
                "K[458.325864]PGLAITFAK[229.162932]",
                1e-5);

        assertEquals("[UNIMOD:737]-K[UNIMOD:737]PGLAITFAK[UNIMOD:737]-[]", normalized);
        assertEquals("[+229.162932]K[+229.162932]PGLAITFAK[+229.162932]",
                PeptideSequenceConverter.unimodToMassEncoded(normalized));
    }

    @Test
    void normalizeMassEncodedHandlesNewUnimodLookups() {
        String normalized = PeptideSequenceConverter.normalizeToUnimod(
                "[+100.016044]K[+28.0313]R[+42.04695]K[+224.152478]",
                1e-5);

        assertEquals("[UNIMOD:64]-K[UNIMOD:36]R[UNIMOD:37]K[UNIMOD:739]-[]", normalized);
        assertEquals(
                "[+100.016044]K[+28.031300]R[+42.046950]K[+224.152478]",
                PeptideSequenceConverter.unimodToMassEncoded(normalized));
    }

    @Test
    void normalizeMassEncodedMapsLegacyTmtAliasMassToTmt0() {
        String normalized = PeptideSequenceConverter.normalizeToUnimod(
                "[+224.152478]PEPTIDEK[+224.152478]",
                1e-5);

        assertEquals("[UNIMOD:739]-PEPTIDEK[UNIMOD:739]-[]", normalized);
        assertEquals("[+224.152478]PEPTIDEK[+224.152478]",
                PeptideSequenceConverter.unimodToMassEncoded(normalized));
    }

    @Test
    void normalizeMassEncodedMapsProtonatedLegacyTmtAliasMassToUnimod737() {
        IllegalArgumentException error = assertThrows(
                IllegalArgumentException.class,
                () -> PeptideSequenceConverter.normalizeToUnimod(
                        "[+225.155833]PEPTIDEK[+225.155833]",
                        1e-5));
        assertTrue(error.getMessage().contains("does not map to known UNIMOD ids"));
    }

    @Test
    void normalizeMassEncodedMapsProtonatedLegacyTmtAliasWithWhitespaceToUnimod737() {
        IllegalArgumentException error = assertThrows(
                IllegalArgumentException.class,
                () -> PeptideSequenceConverter.normalizeToUnimod(
                        "[+\t225.155833]PEPTIDEK[+ 225.155833]",
                        1e-5));
        assertTrue(error.getMessage().contains("does not map to known UNIMOD ids"));
    }

    @Test
    void unsupportedMassFailsWithClearError() {
        IllegalArgumentException error = assertThrows(
                IllegalArgumentException.class,
                () -> PeptideSequenceConverter.normalizeToUnimod("PEP[+1.234567]TIDE", 1e-5));

        assertTrue(error.getMessage().contains("does not map to known UNIMOD ids"));
    }

    @Test
    void parseNormalizedUnimodKeepsTerminalAndResidueModCounts() {
        PeptideSequenceConverter.ParsedUnimodSequence parsed = PeptideSequenceConverter.parseNormalizedUnimod(
                "[UNIMOD:737]-K[UNIMOD:737]PEPTIDE-[]");

        assertEquals("KPEPTIDE", parsed.getResidues());
        assertEquals(parsed.getResidues().length() + 2, parsed.getPositionMods().size());
        assertEquals("UNIMOD:737", parsed.getPositionMods().get(0).get(0));
        assertEquals("UNIMOD:737", parsed.getPositionMods().get(1).get(0));
    }

    @Test
    void unimodToMassEncodedFoldsNtermPyroCarbamidomethylCysteineForChronologerCompatibility() {
        String massEncoded = PeptideSequenceConverter.unimodToMassEncoded(
                "[UNIMOD:28]-C[UNIMOD:4]AILTVLTAQDVSIFPNVHSDDSQVK-[]");

        assertEquals("C[+39.994915]AILTVLTAQDVSIFPNVHSDDSQVK", massEncoded);
    }

    @Test
    void unimodToMassEncodedRestoresNtermPyroWhenFirstResidueIsNotCarbamidomethylated() {
        String massEncoded = PeptideSequenceConverter.unimodToMassEncoded("[UNIMOD:28]-CPEPTIDE-[]");
        assertEquals("[-17.026549]CPEPTIDE", massEncoded);
    }

    @Test
    void unimodToMassEncodedIncludesCtermModMass() {
        String massEncoded = PeptideSequenceConverter.unimodToMassEncoded("[]-PEPTIDE-[UNIMOD:1]");
        assertEquals("PEPTIDE[+42.010565]", massEncoded);
    }

    @Test
    void normalizeMassEncodedSupportsCtermMassMods() {
        String normalized = PeptideSequenceConverter.normalizeToUnimod("PEPTIDE[+42.010565]", 1e-5);
        assertEquals("[]-PEPTIDE[UNIMOD:1]-[]", normalized);
    }

    @Test
    void sumUnimodMassRejectsUnsupportedTokens() {
        IllegalArgumentException error = assertThrows(
                IllegalArgumentException.class,
                () -> PeptideSequenceConverter.sumUnimodMass(List.of("UNIMOD:999")));
        assertTrue(error.getMessage().contains("Unsupported UNIMOD token"));
    }

    @Test
    void foldFirstResiduePyrogluHandlesQAndEAndNoOpCases() {
        List<String> ntermModsQ = new ArrayList<>();
        List<List<String>> residueModsQ = new ArrayList<>();
        residueModsQ.add(new ArrayList<>(List.of("UNIMOD:28")));
        PeptideSequenceConverter.foldFirstResiduePyrogluToNterm("QPEPTIDE", ntermModsQ, residueModsQ);
        assertEquals(List.of("UNIMOD:28"), ntermModsQ);
        assertEquals(List.of(), residueModsQ.get(0));

        List<String> ntermModsE = new ArrayList<>();
        List<List<String>> residueModsE = new ArrayList<>();
        residueModsE.add(new ArrayList<>(List.of("UNIMOD:27")));
        PeptideSequenceConverter.foldFirstResiduePyrogluToNterm("EPEPTIDE", ntermModsE, residueModsE);
        assertEquals(List.of("UNIMOD:27"), ntermModsE);
        assertEquals(List.of(), residueModsE.get(0));

        List<String> ntermExisting = new ArrayList<>(List.of("UNIMOD:1"));
        List<List<String>> residueModsExisting = new ArrayList<>();
        residueModsExisting.add(new ArrayList<>(List.of("UNIMOD:28")));
        PeptideSequenceConverter.foldFirstResiduePyrogluToNterm("QPEPTIDE", ntermExisting, residueModsExisting);
        assertEquals(List.of("UNIMOD:1"), ntermExisting);
        assertEquals(List.of("UNIMOD:28"), residueModsExisting.get(0));
    }

    @Test
    void parseUnimodRejectsMalformedSequences() {
        assertThrows(IllegalArgumentException.class, () -> PeptideSequenceConverter.parseNormalizedUnimod("[]PEPTIDE-[]"));
        assertThrows(IllegalArgumentException.class, () -> PeptideSequenceConverter.parseNormalizedUnimod("[]-PEPTIDE[]"));
        assertThrows(IllegalArgumentException.class, () -> PeptideSequenceConverter.parseNormalizedUnimod("[]-PEPTIDE-[]x"));
        assertThrows(IllegalArgumentException.class, () -> PeptideSequenceConverter.parseNormalizedUnimod("[]-PEP*IDE-[]"));
        assertThrows(IllegalArgumentException.class, () -> PeptideSequenceConverter.parseNormalizedUnimod("PEPTIDE-[]"));
        assertThrows(IllegalArgumentException.class, () -> PeptideSequenceConverter.parseNormalizedUnimod("[]-PEPTIDE-[UNIMOD:1"));
    }

    @Test
    void parseUnimodRejectsUnsupportedAndUnknownTokens() {
        IllegalArgumentException unsupported = assertThrows(
                IllegalArgumentException.class,
                () -> PeptideSequenceConverter.parseNormalizedUnimod("[ABC]-PEPTIDE-[]"));
        assertTrue(unsupported.getMessage().contains("Unsupported modification token"));

        IllegalArgumentException unknown = assertThrows(
                IllegalArgumentException.class,
                () -> PeptideSequenceConverter.parseNormalizedUnimod("[UNIMOD:999]-PEPTIDE-[]"));
        assertTrue(unknown.getMessage().contains("Unknown UNIMOD token"));
    }

    @Test
    void parseMassEncodedRejectsMalformedSequences() {
        assertThrows(IllegalArgumentException.class, () -> PeptideSequenceConverter.normalizeToUnimod("[+42.010565]", 1e-5));
        assertThrows(IllegalArgumentException.class, () -> PeptideSequenceConverter.normalizeToUnimod("PEPTIDE@", 1e-5));
        IllegalArgumentException invalidToken = assertThrows(
                IllegalArgumentException.class,
                () -> PeptideSequenceConverter.normalizeToUnimod("PEP[UNIMOD:1]TIDE", 1e-5));
        assertTrue(invalidToken.getMessage().contains("Expected numeric mass token"));
    }

    @Test
    void normalizeMassEncodedRejectsAmbiguousMappings() {
        IllegalArgumentException ambiguousSingle = assertThrows(
                IllegalArgumentException.class,
                () -> PeptideSequenceConverter.normalizeToUnimod("PEP[+42.020000]TIDE", 0.05));
        assertTrue(ambiguousSingle.getMessage().contains("Ambiguous mass"));

        IllegalArgumentException ambiguousPair = assertThrows(
                IllegalArgumentException.class,
                () -> PeptideSequenceConverter.normalizeToUnimod("PEP[+56.062600]TIDE", 1e-5));
        assertTrue(ambiguousPair.getMessage().contains("Ambiguous compound mass"));
    }

    @Test
    void formatSignedMassRoundsNearZeroToPositiveZero() throws Exception {
        Method method = PeptideSequenceConverter.class.getDeclaredMethod("formatSignedMass", double.class);
        method.setAccessible(true);
        try {
            String value = (String) method.invoke(null, 1e-8);
            assertEquals("+0.000000", value);
        } catch (InvocationTargetException e) {
            throw new AssertionError("formatSignedMass invocation failed", e.getCause());
        }
    }
}

package org.searlelab.jchronologer.preprocessing;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;

class PeptideSequenceConverterTest {

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
}

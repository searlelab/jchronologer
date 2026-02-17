package org.searlelab.jchronologer.inference;

import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;
import org.searlelab.jchronologer.preprocessing.PeptideSequenceConverter;
import org.searlelab.jchronologer.preprocessing.PeptideSequenceConverter.ParsedUnimodSequence;

class PeptideMassCalculatorTest {

    @Test
    void precursorAndFragmentMzArePositiveForModifiedPeptide() {
        ParsedUnimodSequence peptide = PeptideSequenceConverter.parseNormalizedUnimod(
                "[UNIMOD:1]-ACDM[UNIMOD:35]EK-[]");

        double precursorMz = PeptideMassCalculator.calculatePrecursorMz(peptide, 2);
        double b4 = PeptideMassCalculator.calculateFragmentMz(peptide, 4, 1, false);
        double y3 = PeptideMassCalculator.calculateFragmentMz(peptide, 3, 1, true);

        assertTrue(precursorMz > 0.0);
        assertTrue(b4 > 0.0);
        assertTrue(y3 > 0.0);
    }

    @Test
    void calculatePrecursorMzRejectsNonPositiveCharge() {
        ParsedUnimodSequence peptide = PeptideSequenceConverter.parseNormalizedUnimod("[]-PEPTIDE-[]");
        assertThrows(
                IllegalArgumentException.class,
                () -> PeptideMassCalculator.calculatePrecursorMz(peptide, 0));
    }

    @Test
    void calculateFragmentMzValidatesArguments() {
        ParsedUnimodSequence peptide = PeptideSequenceConverter.parseNormalizedUnimod("[]-PEPTIDE-[]");

        assertThrows(
                IllegalArgumentException.class,
                () -> PeptideMassCalculator.calculateFragmentMz(peptide, 0, 1, true));
        assertThrows(
                IllegalArgumentException.class,
                () -> PeptideMassCalculator.calculateFragmentMz(peptide, 7, 1, false));
        assertThrows(
                IllegalArgumentException.class,
                () -> PeptideMassCalculator.calculateFragmentMz(peptide, 3, 0, false));
    }

    @Test
    void calculateMassRejectsUnsupportedResidues() {
        ParsedUnimodSequence peptide = PeptideSequenceConverter.parseNormalizedUnimod("[]-ABCD-[]");
        assertThrows(
                IllegalArgumentException.class,
                () -> PeptideMassCalculator.calculatePrecursorMz(peptide, 2));
    }

    @Test
    void isFragmentPossibleValidatesBoundaries() {
        assertFalse(PeptideMassCalculator.isFragmentPossible(1, 2, 1, 1));
        assertFalse(PeptideMassCalculator.isFragmentPossible(8, 2, 0, 1));
        assertFalse(PeptideMassCalculator.isFragmentPossible(8, 2, 8, 1));
        assertFalse(PeptideMassCalculator.isFragmentPossible(8, 2, 3, 0));
        assertFalse(PeptideMassCalculator.isFragmentPossible(8, 2, 3, 4));
        assertFalse(PeptideMassCalculator.isFragmentPossible(8, 2, 3, 3));
        assertTrue(PeptideMassCalculator.isFragmentPossible(8, 3, 3, 3));
        assertTrue(PeptideMassCalculator.isFragmentPossible(8, 2, 3, 2));
    }
}

package org.searlelab.jchronologer.inference;

import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import org.searlelab.jchronologer.preprocessing.PeptideSequenceConverter;
import org.searlelab.jchronologer.preprocessing.PeptideSequenceConverter.ParsedUnimodSequence;

/**
 * Mass calculations for precursor and b/y fragment ions.
 */
public final class PeptideMassCalculator {

    private static final double PROTON = 1.00727646688;
    private static final double H2O = 18.0105647;
    private static final Map<Character, Double> RESIDUE_MASS = createResidueMassMap();

    private PeptideMassCalculator() {
    }

    public static double calculatePrecursorMz(ParsedUnimodSequence peptide, int precursorCharge) {
        if (precursorCharge <= 0) {
            throw new IllegalArgumentException("Precursor charge must be positive.");
        }
        double neutralMass = calculatePeptideNeutralMass(peptide);
        return (neutralMass + PROTON * precursorCharge) / precursorCharge;
    }

    public static double calculateFragmentMz(
            ParsedUnimodSequence peptide,
            int ionNumber,
            int fragmentCharge,
            boolean yIon) {
        if (ionNumber <= 0) {
            throw new IllegalArgumentException("Fragment ion number must be positive.");
        }
        if (fragmentCharge <= 0) {
            throw new IllegalArgumentException("Fragment charge must be positive.");
        }

        String residues = peptide.getResidues();
        int peptideLength = residues.length();
        if (ionNumber >= peptideLength) {
            throw new IllegalArgumentException("Fragment ion number must be less than peptide length.");
        }

        double[] residueMassWithMods = residueMassWithMods(peptide);
        double ntermDelta = PeptideSequenceConverter.sumUnimodMass(peptide.getPositionMods().get(0));
        double totalNeutralMass = calculatePeptideNeutralMass(peptide);

        double prefixMass = ntermDelta;
        for (int i = 0; i < ionNumber; i++) {
            prefixMass += residueMassWithMods[i];
        }

        double ionNeutralMass;
        if (yIon) {
            int leftPrefixLength = peptideLength - ionNumber;
            double leftPrefixMass = ntermDelta;
            for (int i = 0; i < leftPrefixLength; i++) {
                leftPrefixMass += residueMassWithMods[i];
            }
            ionNeutralMass = totalNeutralMass - leftPrefixMass;
        } else {
            ionNeutralMass = prefixMass;
        }

        return (ionNeutralMass + PROTON * fragmentCharge) / fragmentCharge;
    }

    public static boolean isFragmentPossible(int peptideLength, int precursorCharge, int ionNumber, int fragmentCharge) {
        if (peptideLength <= 1) {
            return false;
        }
        if (ionNumber <= 0 || ionNumber >= peptideLength) {
            return false;
        }
        if (fragmentCharge <= 0 || fragmentCharge > 3) {
            return false;
        }
        return fragmentCharge <= precursorCharge;
    }

    private static double calculatePeptideNeutralMass(ParsedUnimodSequence peptide) {
        String residues = peptide.getResidues();
        List<List<String>> mods = peptide.getPositionMods();

        double mass = H2O;
        mass += PeptideSequenceConverter.sumUnimodMass(mods.get(0));
        mass += PeptideSequenceConverter.sumUnimodMass(mods.get(residues.length() + 1));

        for (int i = 0; i < residues.length(); i++) {
            char residue = residues.charAt(i);
            Double residueMass = RESIDUE_MASS.get(residue);
            if (residueMass == null) {
                throw new IllegalArgumentException("Unsupported residue for mass calculation: " + residue);
            }
            mass += residueMass;
            mass += PeptideSequenceConverter.sumUnimodMass(mods.get(i + 1));
        }
        return mass;
    }

    private static double[] residueMassWithMods(ParsedUnimodSequence peptide) {
        String residues = peptide.getResidues();
        List<List<String>> mods = peptide.getPositionMods();

        double[] values = new double[residues.length()];
        for (int i = 0; i < residues.length(); i++) {
            char residue = residues.charAt(i);
            Double residueMass = RESIDUE_MASS.get(residue);
            if (residueMass == null) {
                throw new IllegalArgumentException("Unsupported residue for mass calculation: " + residue);
            }
            values[i] = residueMass + PeptideSequenceConverter.sumUnimodMass(mods.get(i + 1));
        }
        return values;
    }

    private static Map<Character, Double> createResidueMassMap() {
        Map<Character, Double> values = new LinkedHashMap<>();
        values.put('A', 71.037113805);
        values.put('C', 103.009184505);
        values.put('D', 115.026943065);
        values.put('E', 129.042593135);
        values.put('F', 147.068413945);
        values.put('G', 57.021463735);
        values.put('H', 137.058911875);
        values.put('I', 113.084064015);
        values.put('K', 128.094963050);
        values.put('L', 113.084064015);
        values.put('M', 131.040484645);
        values.put('N', 114.042927470);
        values.put('P', 97.052763875);
        values.put('Q', 128.058577540);
        values.put('R', 156.101111050);
        values.put('S', 87.032028435);
        values.put('T', 101.047678505);
        values.put('V', 99.068413945);
        values.put('W', 186.079312980);
        values.put('Y', 163.063328575);
        return Collections.unmodifiableMap(values);
    }
}

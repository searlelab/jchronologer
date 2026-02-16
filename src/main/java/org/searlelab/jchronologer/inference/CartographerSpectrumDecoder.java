package org.searlelab.jchronologer.inference;

import java.util.ArrayList;
import java.util.List;
import org.searlelab.jchronologer.preprocessing.PeptideSequenceConverter.ParsedUnimodSequence;

/**
 * Converts 174-bin Cartographer outputs into filtered ion arrays.
 */
public final class CartographerSpectrumDecoder {

    public static final int VECTOR_LENGTH = 174;
    private static final int MAX_FRAGMENT_ION_NUMBER = 29;

    private CartographerSpectrumDecoder() {
    }

    public static DecodedSpectrum decode(
            ParsedUnimodSequence peptide,
            byte precursorCharge,
            float[] intensities,
            float minimumReportedIntensity) {
        if (intensities.length != VECTOR_LENGTH) {
            throw new IllegalArgumentException(
                    "Expected Cartographer vector length " + VECTOR_LENGTH + " but received " + intensities.length);
        }

        int peptideLength = peptide.getResidues().length();
        List<Double> massValues = new ArrayList<>();
        List<Float> intensityValues = new ArrayList<>();
        List<String> ionTypes = new ArrayList<>();

        for (int ionNumber = 1; ionNumber <= MAX_FRAGMENT_ION_NUMBER; ionNumber++) {
            int baseIndex = (ionNumber - 1) * 6;
            appendIfPresent(peptide, precursorCharge, ionNumber, true, 1, intensities[baseIndex], minimumReportedIntensity,
                    peptideLength, massValues, intensityValues, ionTypes);
            appendIfPresent(peptide, precursorCharge, ionNumber, true, 2, intensities[baseIndex + 1], minimumReportedIntensity,
                    peptideLength, massValues, intensityValues, ionTypes);
            appendIfPresent(peptide, precursorCharge, ionNumber, true, 3, intensities[baseIndex + 2], minimumReportedIntensity,
                    peptideLength, massValues, intensityValues, ionTypes);
            appendIfPresent(peptide, precursorCharge, ionNumber, false, 1, intensities[baseIndex + 3], minimumReportedIntensity,
                    peptideLength, massValues, intensityValues, ionTypes);
            appendIfPresent(peptide, precursorCharge, ionNumber, false, 2, intensities[baseIndex + 4], minimumReportedIntensity,
                    peptideLength, massValues, intensityValues, ionTypes);
            appendIfPresent(peptide, precursorCharge, ionNumber, false, 3, intensities[baseIndex + 5], minimumReportedIntensity,
                    peptideLength, massValues, intensityValues, ionTypes);
        }

        return new DecodedSpectrum(
                toDoubleArray(massValues),
                toFloatArray(intensityValues),
                ionTypes.toArray(new String[0]));
    }

    private static void appendIfPresent(
            ParsedUnimodSequence peptide,
            byte precursorCharge,
            int ionNumber,
            boolean yIon,
            int fragmentCharge,
            float intensity,
            float minimumReportedIntensity,
            int peptideLength,
            List<Double> massValues,
            List<Float> intensityValues,
            List<String> ionTypes) {
        if (!Float.isFinite(intensity) || intensity < 0.0f || intensity < minimumReportedIntensity) {
            return;
        }
        if (!PeptideMassCalculator.isFragmentPossible(peptideLength, precursorCharge, ionNumber, fragmentCharge)) {
            return;
        }

        double mz = PeptideMassCalculator.calculateFragmentMz(peptide, ionNumber, fragmentCharge, yIon);
        massValues.add(mz);
        intensityValues.add(intensity);
        ionTypes.add(fragmentCharge + "+" + (yIon ? "y" : "b") + ionNumber);
    }

    private static double[] toDoubleArray(List<Double> values) {
        double[] out = new double[values.size()];
        for (int i = 0; i < values.size(); i++) {
            out[i] = values.get(i);
        }
        return out;
    }

    private static float[] toFloatArray(List<Float> values) {
        float[] out = new float[values.size()];
        for (int i = 0; i < values.size(); i++) {
            out[i] = values.get(i);
        }
        return out;
    }

    public static final class DecodedSpectrum {
        private final double[] massArray;
        private final float[] intensityArray;
        private final String[] ionTypeArray;

        private DecodedSpectrum(double[] massArray, float[] intensityArray, String[] ionTypeArray) {
            this.massArray = massArray;
            this.intensityArray = intensityArray;
            this.ionTypeArray = ionTypeArray;
        }

        public double[] getMassArray() {
            return massArray;
        }

        public float[] getIntensityArray() {
            return intensityArray;
        }

        public String[] getIonTypeArray() {
            return ionTypeArray;
        }
    }
}

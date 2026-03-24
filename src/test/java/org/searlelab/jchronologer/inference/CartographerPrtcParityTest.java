package org.searlelab.jchronologer.inference;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.junit.jupiter.api.Test;
import org.searlelab.jchronologer.api.ChronologerLibraryEntry;
import org.searlelab.jchronologer.api.ChronologerLibraryPredictor;
import org.searlelab.jchronologer.api.ChronologerOptions;
import org.searlelab.jchronologer.api.LibraryPredictionRequest;
import org.searlelab.jchronologer.api.PrecursorCondition;
import org.searlelab.jchronologer.impl.ChronologerFactory;
import org.searlelab.jchronologer.preprocessing.ChronologerPreprocessor;
import org.searlelab.jchronologer.preprocessing.PeptideSequenceConverter;
import org.searlelab.jchronologer.preprocessing.PeptideSequenceConverter.ParsedUnimodSequence;
import org.searlelab.jchronologer.preprocessing.PreprocessingMetadataLoader;
import org.searlelab.jchronologer.preprocessing.PreprocessingOutcome;

/**
 * Compares Cartographer predictions for 15 PRTC peptides against Python reference
 * vectors (epoch 12 from full_prtc_log.tsv).
 *
 * <p>These peptides are all unmodified, predicted at charge 2 and NCE 33.
 */
class CartographerPrtcParityTest {

	private static final String REFERENCE_RESOURCE = "data/golden/prtc_reference_vectors.tsv";
	private static final int VECTOR_LENGTH = CartographerSpectrumDecoder.VECTOR_LENGTH;
	private static final Pattern ION_TYPE_PATTERN = Pattern.compile("([1-3])\\+([yb])(\\d+)");

	/** PRTC peptides in the same order as the Python trainer. */
	private static final String[] PRTC_PEPTIDES = {
		"SSAAPPPPPR", "GISNEGQNASIK", "HVLTSIGEK", "DIPVPKPK",
		"IGDYAGIK", "TASEFDSAIAQDK", "SAAGAFGPELSR", "ELGQSGVDTYLQTK",
		"GLILVGGYGTR", "GILFVGSGVSGGEEGAR", "SFANQPLEVVYSK",
		"LTILEELR", "NGFILDGFPR", "ELASGLSFPVGFK", "LSSEAPALFQFDLK",
	};

	/** Model-level NCE (already normalized, for direct CartographerBatchPredictor calls). */
	private static final float PRTC_NCE_NORMALIZED = 0.33f;
	/** Standard NCE units (for PrecursorCondition / full-pipeline calls). */
	private static final double PRTC_NCE = 33.0;
	private static final byte PRTC_CHARGE = 2;

	/**
	 * Verifies that the raw 174-element Cartographer output vectors from Java/DJL
	 * match the Python reference within a cosine similarity threshold.
	 */
	@Test
	void rawVectorsShouldMatchPythonReference() throws IOException {
		Map<String, float[]> referenceByPeptide = loadReferenceVectors();
		float[][] predictions = predictVectors(PRTC_PEPTIDES);
		assertEquals(PRTC_PEPTIDES.length, predictions.length);

		for (int i = 0; i < PRTC_PEPTIDES.length; i++) {
			String peptide = PRTC_PEPTIDES[i];
			float[] predicted = predictions[i];
			float[] reference = referenceByPeptide.get(peptide);
			ParsedUnimodSequence parsed = PeptideSequenceConverter.parseNormalizedUnimod("[]-" + peptide + "-[]");
			List<String> validIonOrder = validIonTypes(parsed.getResidues().length(), PRTC_CHARGE);

			assertTrue(reference != null, "Missing reference for " + peptide);
			assertEquals(VECTOR_LENGTH, predicted.length, "Wrong output width for " + peptide);

			Map<String, Float> referenceByIon = decodeByIon(parsed, reference, 0.0f);
			Map<String, Float> predictedByIon = decodeByIon(parsed, predicted, 0.0f);
			float[] referenceNormalized = normalizeDecodedIntensities(referenceByIon, validIonOrder);
			float[] predictedNormalized = normalizeDecodedIntensities(predictedByIon, validIonOrder);

			double cosine = cosineSimilarity(predictedNormalized, referenceNormalized);
			int sigPred = countAbove(predictedNormalized, 0.01f);
			int sigRef = countAbove(referenceNormalized, 0.01f);
			assertTrue(cosine > 0.80,
					peptide + ": cosine similarity " + String.format("%.4f", cosine)
							+ " is below 0.80 threshold (decoded + masked + renormalized).");
			assertTrue(sigPred >= 5, peptide + ": only " + sigPred + " significant decoded intensities in prediction.");
			assertTrue(sigRef >= 5, peptide + ": reference unexpectedly has only " + sigRef + " significant decoded intensities.");

			int[] topReference = topIndices(referenceNormalized, 12);
			int[] topPredicted = topIndices(predictedNormalized, 12);
			int topOverlap = overlapCount(topReference, topPredicted);
			assertTrue(topOverlap >= 7,
					peptide + ": top-12 overlap too small (" + topOverlap + "), expected >= 7.");

			int[] referenceTop3 = topIndices(referenceNormalized, 3);
			assertTrue(contains(referenceTop3, topPredicted[0]),
					peptide + ": strongest predicted peak is not in reference top-3.");
		}
	}

	/**
	 * Verifies that each PRTC peptide produces a spectrum with multiple fragment ions,
	 * not just a single ion like y1.
	 */
	@Test
	void decodedSpectraShouldHaveMultipleIons() {
		float[][] predictions = predictVectors(PRTC_PEPTIDES);
		for (int i = 0; i < PRTC_PEPTIDES.length; i++) {
			String peptide = PRTC_PEPTIDES[i];
			float[] vector = predictions[i];

			ParsedUnimodSequence parsed = PeptideSequenceConverter.parseNormalizedUnimod(
					"[]-" + peptide + "-[]");
			CartographerSpectrumDecoder.DecodedSpectrum decoded =
					CartographerSpectrumDecoder.decode(parsed, PRTC_CHARGE, vector, 0.01f);

			assertIonArraysLookValid(decoded.getMassArray(),
					decoded.getIntensityArray(),
					decoded.getIonTypeArray(),
					parsed.getResidues().length(),
					PRTC_CHARGE,
					5);

			int highIntensityCount = countAbove(decoded.getIntensityArray(), 0.10f);
			assertTrue(highIntensityCount >= 2,
					peptide + ": expected >= 2 high-intensity ions, found " + highIntensityCount);
		}
	}

	/**
	 * Verifies that the tokenization of PRTC peptides produces the expected coded
	 * sequences and token arrays matching the Python preprocessing.
	 */
	@Test
	void prtcTokenizationMatchesPythonConvention() {
		ChronologerPreprocessor preprocessor = new ChronologerPreprocessor(
				PreprocessingMetadataLoader.loadFromClasspath(
						ChronologerOptions.DEFAULT_CARTOGRAPHER_PREPROCESSING_RESOURCE));

		for (String peptide : PRTC_PEPTIDES) {
			PreprocessingOutcome outcome = preprocessor.preprocess(peptide);
			assertTrue(outcome.isAccepted(),
					"Failed to tokenize " + peptide + ": " + outcome.getRejectionReason());

			// Coded sequence should be "-" + peptide + "_"
			String expectedCoded = "-" + peptide + "_";
			assertEquals(expectedCoded, outcome.getCodedPeptideSeq(),
					"Wrong coded sequence for " + peptide);

			// Token array should have length 33 (max_peptide_len=31 + 2)
			assertEquals(33, outcome.getTokenArray().length,
					"Wrong token array length for " + peptide);

			// First token should be '-' (45), last non-padding should be '_' (51)
			assertEquals(45, outcome.getTokenArray()[0],
					"First token should be '-' (45) for " + peptide);
			assertEquals(51, outcome.getTokenArray()[peptide.length() + 1],
					"Terminal token should be '_' (51) for " + peptide);

			// Remaining should be padding zeros
			for (int j = peptide.length() + 2; j < 33; j++) {
				assertEquals(0, outcome.getTokenArray()[j],
						"Position " + j + " should be padding for " + peptide);
			}
		}
	}

	/**
	 * Full-pipeline test: runs PRTC peptides through DefaultChronologerLibraryPredictor
	 * (the same path a real user would call) and verifies rich fragmentation.
	 */
	@Test
	void fullPipelinePrtcPeptidesShouldHaveRichSpectra() {
		List<LibraryPredictionRequest> requests = new ArrayList<>();
		for (String peptide : PRTC_PEPTIDES) {
			requests.add(new LibraryPredictionRequest(
					peptide,
					List.of(new PrecursorCondition(PRTC_CHARGE, PRTC_NCE))));
		}

		try (ChronologerLibraryPredictor predictor = ChronologerFactory.createLibraryPredictorDefault()) {
			List<ChronologerLibraryEntry> entries = predictor.predict(requests);
			assertEquals(PRTC_PEPTIDES.length, entries.size());

			for (int i = 0; i < entries.size(); i++) {
				ChronologerLibraryEntry entry = entries.get(i);
				String peptide = PRTC_PEPTIDES[i];
				ParsedUnimodSequence parsed = PeptideSequenceConverter.parseNormalizedUnimod(
						entry.getUnimodPeptideSequence());

				assertIonArraysLookValid(
						entry.getMassArray(),
						entry.getIntensityArray(),
						entry.getIonTypeArray(),
						parsed.getResidues().length(),
						entry.getPrecursorCharge(),
						5);
				assertTrue(entry.getPrecursorMz() > 0.0, peptide + ": precursor m/z must be positive.");
				assertTrue(entry.getRetentionTimeInSeconds() > 0.0f, peptide + ": RT must be positive.");
			}
		}
	}

	/**
	 * Verifies that specific positions in the raw TASEFDSAIAQDK vector match
	 * Prosit ordering (ion-interleaved, stride 6), not channel-first (stride 29).
	 *
	 * <p>Prosit index formula: {@code (ionNumber - 1) * 6 + channelOffset}
	 * where channelOffset: y+1=0, y+2=1, y+3=2, b+1=3, b+2=4, b+3=5
	 */
	@Test
	void decodedHighIntensityIonsShouldPreserveReferenceOrderingForTasef() throws IOException {
		Map<String, float[]> refs = loadReferenceVectors();
		float[] referenceVector = refs.get("TASEFDSAIAQDK");
		assertTrue(referenceVector != null, "Missing reference for TASEFDSAIAQDK");

		float[] predictedVector = predictVectors(new String[] {"TASEFDSAIAQDK"})[0];
		ParsedUnimodSequence parsed = PeptideSequenceConverter.parseNormalizedUnimod(
				"[]-TASEFDSAIAQDK-[]");

		Map<String, Float> referenceByIon = decodeByIon(parsed, referenceVector, 0.01f);
		Map<String, Float> predictedByIon = decodeByIon(parsed, predictedVector, 0.01f);

		List<String> referenceTop6 = topIonTypes(referenceByIon, 6);
		List<String> predictedTop6 = topIonTypes(predictedByIon, 6);
		assertEquals(6, referenceTop6.size());
		assertEquals(6, predictedTop6.size());

		// Stable reference ordering sanity checks.
		assertEquals("1+y4", referenceTop6.get(0));
		assertTrue(referenceTop6.contains("1+b2"));

		// High-intensity ions should remain stable after model refreshes.
		assertTrue(predictedTop6.contains(referenceTop6.get(0)),
				"Predicted top-6 should contain the strongest reference ion " + referenceTop6.get(0));
		assertTrue(predictedTop6.contains(referenceTop6.get(1)),
				"Predicted top-6 should contain the second strongest reference ion " + referenceTop6.get(1));
		int overlap = overlapCount(referenceTop6, predictedTop6);
		assertTrue(overlap >= 4, "Top-6 decoded ion overlap is too small: " + overlap + " (expected >= 4)");
	}

	@Test
	void electricianPrtcPeptidesShouldStronglyFavorChargeTwo() {
		float[][] distributions = predictChargeDistributions(PRTC_PEPTIDES);
		assertEquals(PRTC_PEPTIDES.length, distributions.length);

		for (int i = 0; i < PRTC_PEPTIDES.length; i++) {
			String peptide = PRTC_PEPTIDES[i];
			float[] normalized = normalizeDistribution(distributions[i]);
			assertEquals(6, normalized.length, "Electrician output must have 6 charge-state probabilities.");

			float plusTwo = normalized[PRTC_CHARGE - 1];
			float maxProbability = -1.0f;
			int maxCharge = -1;
			for (int charge = 1; charge <= normalized.length; charge++) {
				float probability = normalized[charge - 1];
				if (probability > maxProbability) {
					maxProbability = probability;
					maxCharge = charge;
				}
			}

			assertEquals(PRTC_CHARGE, maxCharge, peptide + ": max Electrician probability is not at charge +2.");
			assertTrue(plusTwo >= 0.90f, peptide + ": expected P(+2) >= 0.9, found " + plusTwo);
		}
	}

	// ── Helpers ──────────────────────────────────────────────────────────

	private static float[][] predictVectors(String[] peptides) {
		ChronologerPreprocessor preprocessor = new ChronologerPreprocessor(
				PreprocessingMetadataLoader.loadFromClasspath(
						ChronologerOptions.DEFAULT_CARTOGRAPHER_PREPROCESSING_RESOURCE));
		try (CartographerBatchPredictor predictor = new CartographerBatchPredictor(
				ChronologerOptions.DEFAULT_CARTOGRAPHER_MODEL_RESOURCE)) {
			long[][] tokenBatch = new long[peptides.length][];
			float[][] chargeBatch = new float[peptides.length][6];
			float[][] nceBatch = new float[peptides.length][1];

			for (int i = 0; i < peptides.length; i++) {
				PreprocessingOutcome outcome = preprocessor.preprocess(peptides[i]);
				assertTrue(outcome.isAccepted(),
						"Failed to tokenize peptide " + peptides[i] + ": " + outcome.getRejectionReason());
				tokenBatch[i] = outcome.getTokenArray();
				chargeBatch[i][PRTC_CHARGE - 1] = 1.0f;
				nceBatch[i][0] = PRTC_NCE_NORMALIZED;
			}
			return predictor.predict(tokenBatch, chargeBatch, nceBatch);
		}
	}

	private static float[][] predictChargeDistributions(String[] peptides) {
		ChronologerPreprocessor preprocessor = new ChronologerPreprocessor(
				PreprocessingMetadataLoader.loadFromClasspath(
						ChronologerOptions.DEFAULT_ELECTRICIAN_PREPROCESSING_RESOURCE));
		try (ElectricianBatchPredictor predictor = new ElectricianBatchPredictor(
				ChronologerOptions.DEFAULT_ELECTRICIAN_MODEL_RESOURCE)) {
			long[][] tokenBatch = new long[peptides.length][];
			for (int i = 0; i < peptides.length; i++) {
				PreprocessingOutcome outcome = preprocessor.preprocess(peptides[i]);
				assertTrue(outcome.isAccepted(),
						"Failed to tokenize peptide " + peptides[i] + ": " + outcome.getRejectionReason());
				tokenBatch[i] = outcome.getTokenArray();
			}
			return predictor.predict(tokenBatch);
		}
	}

	private static float[] normalizeDistribution(float[] rawDistribution) {
		float[] normalized = new float[rawDistribution.length];
		double sum = 0.0;
		for (int i = 0; i < rawDistribution.length; i++) {
			float value = rawDistribution[i];
			if (Float.isFinite(value) && value > 0.0f) {
				normalized[i] = value;
				sum += value;
			}
		}
		if (sum <= 0.0) {
			return normalized;
		}
		for (int i = 0; i < normalized.length; i++) {
			normalized[i] = (float) (normalized[i] / sum);
		}
		return normalized;
	}

	private static void assertIonArraysLookValid(
			double[] masses,
			float[] intensities,
			String[] ionTypes,
			int peptideLength,
			byte precursorCharge,
			int minimumIonCount) {
		assertEquals(masses.length, intensities.length);
		assertEquals(masses.length, ionTypes.length);
		assertTrue(ionTypes.length >= minimumIonCount,
				"Expected at least " + minimumIonCount + " ions, found " + ionTypes.length);

		for (int i = 0; i < ionTypes.length; i++) {
			assertTrue(masses[i] > 0.0, "Expected positive m/z for ion " + ionTypes[i]);
			assertTrue(intensities[i] >= 0.01f, "Expected filtered intensity for ion " + ionTypes[i]);

			Matcher matcher = ION_TYPE_PATTERN.matcher(ionTypes[i]);
			assertTrue(matcher.matches(), "Unexpected ion label format: " + ionTypes[i]);
			int fragmentCharge = Integer.parseInt(matcher.group(1));
			int ionNumber = Integer.parseInt(matcher.group(3));
			assertTrue(fragmentCharge <= precursorCharge, "Fragment charge exceeds precursor in " + ionTypes[i]);
			assertTrue(ionNumber > 0 && ionNumber < peptideLength, "Ion number out of range: " + ionTypes[i]);
		}
	}

	private static int[] topIndices(float[] vector, int count) {
		Integer[] indices = new Integer[vector.length];
		for (int i = 0; i < vector.length; i++) {
			indices[i] = i;
		}
		Arrays.sort(indices, (left, right) -> Float.compare(sortableIntensity(vector[right]), sortableIntensity(vector[left])));
		int size = Math.min(count, indices.length);
		int[] top = new int[size];
		for (int i = 0; i < size; i++) {
			top[i] = indices[i];
		}
		return top;
	}

	private static boolean contains(int[] values, int query) {
		for (int value : values) {
			if (value == query) {
				return true;
			}
		}
		return false;
	}

	private static int overlapCount(int[] left, int[] right) {
		Set<Integer> rightSet = Arrays.stream(right).boxed().collect(java.util.stream.Collectors.toSet());
		int overlap = 0;
		for (int value : left) {
			if (rightSet.contains(value)) {
				overlap++;
			}
		}
		return overlap;
	}

	private static int overlapCount(List<String> left, List<String> right) {
		Set<String> rightSet = Set.copyOf(right);
		int overlap = 0;
		for (String value : left) {
			if (rightSet.contains(value)) {
				overlap++;
			}
		}
		return overlap;
	}

	private static Map<String, Float> decodeByIon(ParsedUnimodSequence parsed, float[] vector, float minimumReportedIntensity) {
		CartographerSpectrumDecoder.DecodedSpectrum decoded =
				CartographerSpectrumDecoder.decode(parsed, PRTC_CHARGE, vector, minimumReportedIntensity);
		Map<String, Float> decodedByIon = new LinkedHashMap<>();
		for (int i = 0; i < decoded.getIonTypeArray().length; i++) {
			decodedByIon.put(decoded.getIonTypeArray()[i], decoded.getIntensityArray()[i]);
		}
		return decodedByIon;
	}

	private static List<String> topIonTypes(Map<String, Float> byIon, int count) {
		return byIon.entrySet().stream()
				.sorted((left, right) -> Float.compare(right.getValue(), left.getValue()))
				.limit(count)
				.map(Map.Entry::getKey)
				.toList();
	}

	private static List<String> validIonTypes(int peptideLength, byte precursorCharge) {
		List<String> ions = new ArrayList<>();
		for (int ionNumber = 1; ionNumber <= 29; ionNumber++) {
			for (int fragmentCharge = 1; fragmentCharge <= 3; fragmentCharge++) {
				if (PeptideMassCalculator.isFragmentPossible(
						peptideLength,
						precursorCharge,
						ionNumber,
						fragmentCharge)) {
					ions.add(fragmentCharge + "+y" + ionNumber);
					ions.add(fragmentCharge + "+b" + ionNumber);
				}
			}
		}
		return ions;
	}

	private static float[] normalizeDecodedIntensities(Map<String, Float> byIon, List<String> ionOrder) {
		float[] values = new float[ionOrder.size()];
		double sum = 0.0;
		for (int i = 0; i < ionOrder.size(); i++) {
			float intensity = byIon.getOrDefault(ionOrder.get(i), 0.0f);
			if (Float.isFinite(intensity) && intensity > 0.0f) {
				values[i] = intensity;
				sum += intensity;
			}
		}
		if (sum <= 0.0) {
			return values;
		}
		for (int i = 0; i < values.length; i++) {
			values[i] = (float) (values[i] / sum);
		}
		return values;
	}

	private static float sortableIntensity(float value) {
		return Float.isFinite(value) ? value : -Float.MAX_VALUE;
	}

	private static Map<String, float[]> loadReferenceVectors() throws IOException {
		Map<String, float[]> result = new LinkedHashMap<>();
		ClassLoader loader = Thread.currentThread().getContextClassLoader();
		try (InputStream stream = loader.getResourceAsStream(REFERENCE_RESOURCE)) {
			if (stream == null) {
				fail("Missing test resource: " + REFERENCE_RESOURCE);
			}
			BufferedReader reader = new BufferedReader(new InputStreamReader(stream));
			reader.readLine(); // skip header
			String line;
			while ((line = reader.readLine()) != null) {
				String[] parts = line.split("\t");
				// parts[0] = epoch, parts[1] = peptide, parts[2..175] = i0..i173
				String peptide = parts[1];
				float[] vector = new float[VECTOR_LENGTH];
				for (int j = 0; j < VECTOR_LENGTH; j++) {
					vector[j] = Float.parseFloat(parts[j + 2]);
				}
				result.put(peptide, vector);
			}
		}
		return result;
	}

	private static int countAbove(float[] v, float threshold) {
		int count = 0;
		for (float val : v) {
			if (val > threshold) {
				count++;
			}
		}
		return count;
	}

	private static double cosineSimilarity(float[] a, float[] b) {
		if (a.length != b.length) {
			throw new IllegalArgumentException("Vector length mismatch");
		}
		double dot = 0.0, normA = 0.0, normB = 0.0;
		for (int i = 0; i < a.length; i++) {
			dot += (double) a[i] * b[i];
			normA += (double) a[i] * a[i];
			normB += (double) b[i] * b[i];
		}
		double denom = Math.sqrt(normA) * Math.sqrt(normB);
		return denom == 0.0 ? 0.0 : dot / denom;
	}
}

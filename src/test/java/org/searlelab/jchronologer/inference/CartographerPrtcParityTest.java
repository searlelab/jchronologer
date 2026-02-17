package org.searlelab.jchronologer.inference;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import org.junit.jupiter.api.Test;
import org.searlelab.jchronologer.api.ChronologerLibraryEntry;
import org.searlelab.jchronologer.api.ChronologerLibraryOptions;
import org.searlelab.jchronologer.api.ChronologerLibraryPredictor;
import org.searlelab.jchronologer.api.LibraryPredictionRequest;
import org.searlelab.jchronologer.api.PrecursorCondition;
import org.searlelab.jchronologer.impl.ChronologerFactory;
import org.searlelab.jchronologer.preprocessing.ChronologerPreprocessor;
import org.searlelab.jchronologer.preprocessing.PreprocessingMetadataLoader;
import org.searlelab.jchronologer.preprocessing.PreprocessingOutcome;
import org.searlelab.jchronologer.preprocessing.PeptideSequenceConverter;
import org.searlelab.jchronologer.preprocessing.PeptideSequenceConverter.ParsedUnimodSequence;

/**
 * Compares Cartographer predictions for 15 PRTC peptides against Python reference
 * vectors (epoch 12 from full_prtc_log.tsv).
 *
 * <p>These peptides are all unmodified, predicted at charge 2 and NCE 33.
 */
class CartographerPrtcParityTest {

	private static final String REFERENCE_RESOURCE = "data/golden/prtc_reference_vectors.tsv";
	private static final int VECTOR_LENGTH = CartographerSpectrumDecoder.VECTOR_LENGTH;

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

		ChronologerPreprocessor preprocessor = new ChronologerPreprocessor(
				PreprocessingMetadataLoader.loadFromClasspath(
						ChronologerLibraryOptions.DEFAULT_CARTOGRAPHER_PREPROCESSING_RESOURCE));

		try (CartographerBatchPredictor predictor = new CartographerBatchPredictor(
				ChronologerLibraryOptions.DEFAULT_CARTOGRAPHER_MODEL_RESOURCE)) {

			long[][] tokenBatch = new long[PRTC_PEPTIDES.length][];
			float[][] chargeBatch = new float[PRTC_PEPTIDES.length][6];
			float[][] nceBatch = new float[PRTC_PEPTIDES.length][1];

			for (int i = 0; i < PRTC_PEPTIDES.length; i++) {
				PreprocessingOutcome outcome = preprocessor.preprocess(PRTC_PEPTIDES[i]);
				assertTrue(outcome.isAccepted(),
						"Failed to tokenize PRTC peptide " + PRTC_PEPTIDES[i] + ": " + outcome.getRejectionReason());
				tokenBatch[i] = outcome.getTokenArray();
				chargeBatch[i][PRTC_CHARGE - 1] = 1.0f;
				nceBatch[i][0] = PRTC_NCE_NORMALIZED;
			}

			float[][] predictions = predictor.predict(tokenBatch, chargeBatch, nceBatch);
			assertEquals(PRTC_PEPTIDES.length, predictions.length);

			System.out.println("=== PRTC Parity Diagnostics ===");
			for (int i = 0; i < PRTC_PEPTIDES.length; i++) {
				String peptide = PRTC_PEPTIDES[i];
				float[] predicted = predictions[i];
				float[] reference = referenceByPeptide.get(peptide);

				assertTrue(reference != null, "Missing reference for " + peptide);
				assertEquals(VECTOR_LENGTH, predicted.length,
						"Wrong output width for " + peptide);

				double cosine = cosineSimilarity(predicted, reference);
				int sigPred = countAbove(predicted, 0.01f);
				int sigRef = countAbove(reference, 0.01f);
				System.out.printf("  %s: cosine=%.6f, sig_pred=%d, sig_ref=%d%n",
						peptide, cosine, sigPred, sigRef);

				assertTrue(cosine > 0.90,
						peptide + ": cosine similarity " + String.format("%.4f", cosine)
								+ " is below 0.90 threshold (sig_pred=" + sigPred
								+ ", sig_ref=" + sigRef + ")");
			}
		}
	}

	/**
	 * Verifies that each PRTC peptide produces a spectrum with multiple fragment ions,
	 * not just a single ion like y1.
	 */
	@Test
	void decodedSpectraShouldHaveMultipleIons() {
		ChronologerPreprocessor preprocessor = new ChronologerPreprocessor(
				PreprocessingMetadataLoader.loadFromClasspath(
						ChronologerLibraryOptions.DEFAULT_CARTOGRAPHER_PREPROCESSING_RESOURCE));

		try (CartographerBatchPredictor predictor = new CartographerBatchPredictor(
				ChronologerLibraryOptions.DEFAULT_CARTOGRAPHER_MODEL_RESOURCE)) {

			long[][] tokenBatch = new long[PRTC_PEPTIDES.length][];
			float[][] chargeBatch = new float[PRTC_PEPTIDES.length][6];
			float[][] nceBatch = new float[PRTC_PEPTIDES.length][1];

			for (int i = 0; i < PRTC_PEPTIDES.length; i++) {
				PreprocessingOutcome outcome = preprocessor.preprocess(PRTC_PEPTIDES[i]);
				assertTrue(outcome.isAccepted());
				tokenBatch[i] = outcome.getTokenArray();
				chargeBatch[i][PRTC_CHARGE - 1] = 1.0f;
				nceBatch[i][0] = PRTC_NCE_NORMALIZED;
			}

			float[][] predictions = predictor.predict(tokenBatch, chargeBatch, nceBatch);

			for (int i = 0; i < PRTC_PEPTIDES.length; i++) {
				String peptide = PRTC_PEPTIDES[i];
				float[] vector = predictions[i];

				// Count how many positions in the raw vector have significant intensity
				int significantCount = 0;
				for (float v : vector) {
					if (v > 0.01f) {
						significantCount++;
					}
				}
				assertTrue(significantCount >= 5,
						peptide + ": only " + significantCount
								+ " significant intensities in raw vector (expected >= 5)");

				// Decode and verify multiple ions survive
				ParsedUnimodSequence parsed = PeptideSequenceConverter.parseNormalizedUnimod(
						"[]-" + peptide + "-[]");
				CartographerSpectrumDecoder.DecodedSpectrum decoded =
						CartographerSpectrumDecoder.decode(parsed, PRTC_CHARGE, vector, 0.01f);

				assertTrue(decoded.getIonTypeArray().length >= 5,
						peptide + ": decoded spectrum has only "
								+ decoded.getIonTypeArray().length + " ions (expected >= 5)");
			}
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
						ChronologerLibraryOptions.DEFAULT_CARTOGRAPHER_PREPROCESSING_RESOURCE));

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

				assertTrue(entry.getIonTypeArray().length >= 5,
						peptide + " (full pipeline): only " + entry.getIonTypeArray().length
								+ " ions decoded (expected >= 5)");
				assertTrue(entry.getMassArray().length >= 5,
						peptide + " (full pipeline): only " + entry.getMassArray().length
								+ " m/z values (expected >= 5)");
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
	void vectorLayoutShouldFollowPrositOrdering() throws IOException {
		Map<String, float[]> refs = loadReferenceVectors();
		float[] vector = refs.get("TASEFDSAIAQDK");
		assertTrue(vector != null, "Missing reference for TASEFDSAIAQDK");

		// Verify known positions from Python Prosit-ordering analysis.
		// index = (ionNumber - 1) * 6 + channelOffset
		assertEquals(1.0f,      vector[18], 0.001f, "i18: y4+1  = (4-1)*6+0 = 18");
		assertEquals(0.966279f, vector[9],  0.001f, "i9:  b2+1  = (2-1)*6+3 = 9");
		assertEquals(0.733529f, vector[48], 0.001f, "i48: y9+1  = (9-1)*6+0 = 48");
		assertEquals(0.674392f, vector[42], 0.001f, "i42: y8+1  = (8-1)*6+0 = 42");
		assertEquals(0.158662f, vector[0],  0.001f, "i0:  y1+1  = (1-1)*6+0 = 0");
		assertEquals(0.233861f, vector[15], 0.001f, "i15: b3+1  = (3-1)*6+3 = 15");
		assertEquals(0.004905f, vector[58], 0.001f, "i58: b10+2 = (10-1)*6+4 = 58");
	}

	/**
	 * Verifies that decoding the TASEFDSAIAQDK reference vector produces the
	 * correct ion type assignments. The expected ions are derived from the
	 * Python training output organized in Prosit ordering.
	 *
	 * <p>TASEFDSAIAQDK: 13 residues, charge 2, max ion number 12.
	 */
	@Test
	void decodedIonsShouldMatchExpectedAssignments() throws IOException {
		Map<String, float[]> refs = loadReferenceVectors();
		float[] vector = refs.get("TASEFDSAIAQDK");
		assertTrue(vector != null, "Missing reference for TASEFDSAIAQDK");

		ParsedUnimodSequence parsed = PeptideSequenceConverter.parseNormalizedUnimod(
				"[]-TASEFDSAIAQDK-[]");

		// Use threshold low enough to capture all 18 expected ions (min is b10+2 = 0.004905)
		CartographerSpectrumDecoder.DecodedSpectrum decoded =
				CartographerSpectrumDecoder.decode(parsed, (byte) 2, vector, 0.004f);

		Map<String, Float> decodedByIon = new LinkedHashMap<>();
		for (int i = 0; i < decoded.getIonTypeArray().length; i++) {
			decodedByIon.put(decoded.getIonTypeArray()[i], decoded.getIntensityArray()[i]);
		}

		// Print decoded ions for diagnostic purposes
		System.out.println("=== TASEFDSAIAQDK Decoded Ion Assignments ===");
		for (Map.Entry<String, Float> entry : decodedByIon.entrySet()) {
			System.out.printf("  %s: %.6f%n", entry.getKey(), entry.getValue());
		}

		// Expected: 18 unmasked ions with intensity >= 0.004 from Python reference
		assertIonIntensity(decodedByIon, "1+y4",  1.0f);
		assertIonIntensity(decodedByIon, "1+b2",  0.966279f);
		assertIonIntensity(decodedByIon, "1+y9",  0.733529f);
		assertIonIntensity(decodedByIon, "1+y8",  0.674392f);
		assertIonIntensity(decodedByIon, "1+y7",  0.557423f);
		assertIonIntensity(decodedByIon, "1+y11", 0.456221f);
		assertIonIntensity(decodedByIon, "1+y5",  0.444093f);
		assertIonIntensity(decodedByIon, "1+y2",  0.32067f);
		assertIonIntensity(decodedByIon, "1+y3",  0.277611f);
		assertIonIntensity(decodedByIon, "1+y6",  0.24649f);
		assertIonIntensity(decodedByIon, "1+b3",  0.233861f);
		assertIonIntensity(decodedByIon, "1+y10", 0.208115f);
		assertIonIntensity(decodedByIon, "1+y1",  0.158662f);
		assertIonIntensity(decodedByIon, "1+b4",  0.092071f);
		assertIonIntensity(decodedByIon, "1+b5",  0.031726f);
		assertIonIntensity(decodedByIon, "1+b10", 0.013218f);
		assertIonIntensity(decodedByIon, "1+b6",  0.0081f);
		assertIonIntensity(decodedByIon, "2+b10", 0.004905f);
	}

	// ── Helpers ──────────────────────────────────────────────────────────

	private static void assertIonIntensity(Map<String, Float> decodedByIon, String ionType, float expected) {
		Float actual = decodedByIon.get(ionType);
		assertTrue(actual != null, "Missing expected ion: " + ionType
				+ " (decoded ions: " + decodedByIon.keySet() + ")");
		assertEquals(expected, actual, 0.002f, "Wrong intensity for " + ionType);
	}

	private static Map<String, float[]> loadReferenceVectors() throws IOException {
		Map<String, float[]> result = new LinkedHashMap<>();
		ClassLoader loader = Thread.currentThread().getContextClassLoader();
		try (InputStream stream = loader.getResourceAsStream(REFERENCE_RESOURCE)) {
			if (stream == null) {
				fail("Missing test resource: " + REFERENCE_RESOURCE);
			}
			BufferedReader reader = new BufferedReader(new InputStreamReader(stream));
			String header = reader.readLine(); // skip header
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

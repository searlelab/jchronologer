package org.searlelab.jchronologer.impl;

import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import org.searlelab.jchronologer.api.AcceptedPrediction;
import org.searlelab.jchronologer.api.Chronologer;
import org.searlelab.jchronologer.api.ChronologerLibraryEntry;
import org.searlelab.jchronologer.api.ChronologerLibraryOptions;
import org.searlelab.jchronologer.api.ChronologerLibraryPredictor;
import org.searlelab.jchronologer.api.ChronologerOptions;
import org.searlelab.jchronologer.api.LibraryPredictionRequest;
import org.searlelab.jchronologer.api.PrecursorCondition;
import org.searlelab.jchronologer.api.PredictionResult;
import org.searlelab.jchronologer.api.RejectedPrediction;
import org.searlelab.jchronologer.inference.CartographerBatchPredictor;
import org.searlelab.jchronologer.inference.CartographerSpectrumDecoder;
import org.searlelab.jchronologer.inference.CartographerSpectrumDecoder.DecodedSpectrum;
import org.searlelab.jchronologer.inference.PeptideMassCalculator;
import org.searlelab.jchronologer.preprocessing.ChronologerPreprocessor;
import org.searlelab.jchronologer.preprocessing.PeptideSequenceConverter;
import org.searlelab.jchronologer.preprocessing.PeptideSequenceConverter.ParsedUnimodSequence;
import org.searlelab.jchronologer.preprocessing.PreprocessingMetadataLoader;
import org.searlelab.jchronologer.preprocessing.PreprocessingOutcome;

/**
 * Default tandem predictor combining Chronologer RT and Cartographer MS2 inference.
 */
public final class DefaultChronologerLibraryPredictor implements ChronologerLibraryPredictor {

    private static final ObjectMapper OBJECT_MAPPER = new ObjectMapper();

    /** Minimum accepted NCE value (standard units, e.g. 20 for 20%). */
    private static final double MIN_NCE = 10.0;
    /** Maximum accepted NCE value (standard units, e.g. 50 for 50%). */
    private static final double MAX_NCE = 60.0;
    /** Divisor to convert standard NCE to model-normalized NCE (NCE / 100). */
    private static final double NCE_NORMALIZATION_DIVISOR = 100.0;

    private final ChronologerLibraryOptions options;
    private final Chronologer chronologer;
    private final ChronologerPreprocessor tokenPreprocessor;
    private final CartographerBatchPredictor cartographerPredictor;
    private final int minPrecursorCharge;
    private final int maxPrecursorCharge;
    private final int cartographerOutputWidth;
    private volatile boolean closed;

    public DefaultChronologerLibraryPredictor(ChronologerLibraryOptions options) {
        this.options = options;
        ChronologerOptions chronologerOptions = ChronologerOptions.builder()
                .modelResource(options.getChronologerModelResource())
                .preprocessingResource(options.getChronologerPreprocessingResource())
                .batchSize(options.getBatchSize())
                .inferenceThreads(options.getInferenceThreads())
                .verboseLogging(options.isVerboseLogging())
                .build();

        this.chronologer = ChronologerFactory.create(chronologerOptions);
        this.tokenPreprocessor = new ChronologerPreprocessor(
                PreprocessingMetadataLoader.loadFromClasspath(options.getCartographerPreprocessingResource()));
        this.cartographerPredictor = new CartographerBatchPredictor(options.getCartographerModelResource());

        CartographerMetadata metadata = loadCartographerMetadata(options.getCartographerPreprocessingResource());
        this.minPrecursorCharge = metadata.minPrecursorCharge;
        this.maxPrecursorCharge = metadata.maxPrecursorCharge;
        this.cartographerOutputWidth = metadata.outputWidth;
    }

    @Override
    public void init() {
        if (closed) {
            throw new IllegalStateException("ChronologerLibraryPredictor has been closed.");
        }
    }

    @Override
    public List<ChronologerLibraryEntry> predict(List<LibraryPredictionRequest> requests) {
        init();
        if (requests == null || requests.isEmpty()) {
            return List.of();
        }

        List<NormalizedRequest> normalizedRequests = new ArrayList<>(requests.size());
        LinkedHashMap<String, String> massByUnimod = new LinkedHashMap<>();

        for (int i = 0; i < requests.size(); i++) {
            LibraryPredictionRequest request = requests.get(i);
            if (request == null) {
                throw new IllegalArgumentException("Request at index " + i + " is null.");
            }
            if (request.getPrecursorConditions().isEmpty()) {
                throw new IllegalArgumentException("Request at index " + i + " has no precursor conditions.");
            }

            String normalizedUnimod = PeptideSequenceConverter.normalizeToUnimod(
                    request.getPeptideSequence(),
                    options.getMassMatchEpsilon());
            String massEncoded = PeptideSequenceConverter.unimodToMassEncoded(normalizedUnimod);

            List<PrecursorCondition> validatedConditions = new ArrayList<>(request.getPrecursorConditions().size());
            for (PrecursorCondition condition : request.getPrecursorConditions()) {
                validateCondition(condition);
                validatedConditions.add(condition);
            }

            normalizedRequests.add(new NormalizedRequest(normalizedUnimod, validatedConditions));
            massByUnimod.putIfAbsent(normalizedUnimod, massEncoded);
        }

        LinkedHashMap<String, long[]> tokenByUnimod = new LinkedHashMap<>();
        List<String> massEncodedUnique = new ArrayList<>(massByUnimod.size());
        for (Map.Entry<String, String> entry : massByUnimod.entrySet()) {
            String unimod = entry.getKey();
            String massEncoded = entry.getValue();

            // Fold first-residue pyroglu to nterm for Cartographer tokenization.
            // This produces e.g. [-17.026549]QPEPTIDE instead of Q[-17.026549]PEPTIDE,
            // matching the encoding used during Cartographer training.
            ParsedUnimodSequence parsed = PeptideSequenceConverter.parseNormalizedUnimod(unimod);
            List<String> ntermMods = new ArrayList<>(parsed.getPositionMods().get(0));
            List<List<String>> residueMods = new ArrayList<>(parsed.getResidues().length());
            for (int ri = 0; ri < parsed.getResidues().length(); ri++) {
                residueMods.add(new ArrayList<>(parsed.getPositionMods().get(ri + 1)));
            }
            PeptideSequenceConverter.foldFirstResiduePyrogluToNterm(
                    parsed.getResidues(), ntermMods, residueMods);
            // Rebuild mass-encoded with the fold applied
            StringBuilder cartographerMassEncoded = new StringBuilder();
            if (!ntermMods.isEmpty()) {
                cartographerMassEncoded.append('[')
                        .append(String.format(java.util.Locale.US, "%+.6f",
                                PeptideSequenceConverter.sumUnimodMass(ntermMods)))
                        .append(']');
            }
            for (int ri = 0; ri < parsed.getResidues().length(); ri++) {
                cartographerMassEncoded.append(parsed.getResidues().charAt(ri));
                List<String> mods = residueMods.get(ri);
                if (!mods.isEmpty()) {
                    cartographerMassEncoded.append('[')
                            .append(String.format(java.util.Locale.US, "%+.6f",
                                    PeptideSequenceConverter.sumUnimodMass(mods)))
                            .append(']');
                }
            }

            PreprocessingOutcome outcome = tokenPreprocessor.preprocess(cartographerMassEncoded.toString());
            if (!outcome.isAccepted()) {
                throw new IllegalArgumentException(
                        "Failed to tokenize peptide " + massEncoded + " (normalized " + unimod + "): "
                                + outcome.getRejectionReason());
            }
            tokenByUnimod.put(unimod, outcome.getTokenArray());
            massEncodedUnique.add(massEncoded);
        }

        PredictionResult rtResult = chronologer.predict(massEncodedUnique);
        if (!rtResult.getRejected().isEmpty()) {
            RejectedPrediction first = rtResult.getRejected().get(0);
            throw new IllegalArgumentException(
                    "RT prediction rejected peptide before library inference: "
                            + first.getPeptideModSeq()
                            + " reason="
                            + first.getRejectionReason());
        }

        Map<String, Float> predHiByMass = new LinkedHashMap<>();
        for (AcceptedPrediction accepted : rtResult.getAccepted()) {
            predHiByMass.put(accepted.getPeptideModSeq(), accepted.getPredHi());
        }

        LinkedHashMap<String, ParsedUnimodSequence> parsedByUnimod = new LinkedHashMap<>();
        List<PredictionJob> jobs = new ArrayList<>();
        for (NormalizedRequest normalizedRequest : normalizedRequests) {
            String unimod = normalizedRequest.unimodPeptide;
            String massEncoded = massByUnimod.get(unimod);
            Float predHi = predHiByMass.get(massEncoded);
            if (predHi == null) {
                throw new IllegalStateException("Missing RT prediction for peptide: " + massEncoded);
            }

            ParsedUnimodSequence parsed = parsedByUnimod.computeIfAbsent(
                    unimod,
                    PeptideSequenceConverter::parseNormalizedUnimod);
            long[] tokens = tokenByUnimod.get(unimod);

            for (PrecursorCondition condition : normalizedRequest.conditions) {
                jobs.add(new PredictionJob(
                        unimod,
                        parsed,
                        tokens,
                        condition.getPrecursorCharge(),
                        condition.getPrecursorNce(),
                        predHi * 60.0f));
            }
        }

        List<ChronologerLibraryEntry> results = new ArrayList<>(jobs.size());
        for (int start = 0; start < jobs.size(); start += options.getCartographerBatchSize()) {
            int end = Math.min(start + options.getCartographerBatchSize(), jobs.size());
            List<PredictionJob> batch = jobs.subList(start, end);

            long[][] tokenBatch = new long[batch.size()][];
            float[][] chargeBatch = new float[batch.size()][maxPrecursorCharge - minPrecursorCharge + 1];
            float[][] nceBatch = new float[batch.size()][1];

            for (int i = 0; i < batch.size(); i++) {
                PredictionJob job = batch.get(i);
                tokenBatch[i] = job.tokenArray;
                chargeBatch[i][job.charge - minPrecursorCharge] = 1.0f;
                nceBatch[i][0] = (float) (job.nce / NCE_NORMALIZATION_DIVISOR);
            }

            float[][] predictions = cartographerPredictor.predict(tokenBatch, chargeBatch, nceBatch);
            for (int i = 0; i < batch.size(); i++) {
                PredictionJob job = batch.get(i);
                float[] vector = predictions[i];
                if (vector.length != cartographerOutputWidth) {
                    throw new IllegalStateException(
                            "Unexpected Cartographer output width: " + vector.length + " expected=" + cartographerOutputWidth);
                }

                DecodedSpectrum decoded = CartographerSpectrumDecoder.decode(
                        job.parsedPeptide,
                        job.charge,
                        vector,
                        options.getMinimumReportedIntensity());

                double precursorMz = PeptideMassCalculator.calculatePrecursorMz(job.parsedPeptide, job.charge);
                results.add(new ChronologerLibraryEntry(
                        job.unimodPeptide,
                        job.charge,
                        job.nce,
                        precursorMz,
                        job.retentionTimeSeconds,
                        decoded.getMassArray(),
                        decoded.getIntensityArray(),
                        decoded.getIonTypeArray()));
            }
        }

        return results;
    }

    private void validateCondition(PrecursorCondition condition) {
        if (condition == null) {
            throw new IllegalArgumentException("Precursor condition must be non-null.");
        }
        int charge = condition.getPrecursorCharge();
        if (charge < minPrecursorCharge || charge > maxPrecursorCharge) {
            throw new IllegalArgumentException(
                    "Precursor charge " + charge + " is outside supported range "
                            + minPrecursorCharge + "-" + maxPrecursorCharge);
        }
        double nce = condition.getPrecursorNce();
        if (!Double.isFinite(nce)) {
            throw new IllegalArgumentException("NCE must be a finite number.");
        }
        if (nce < MIN_NCE || nce > MAX_NCE) {
            throw new IllegalArgumentException(
                    "NCE " + nce + " is outside the supported range " + MIN_NCE + "-" + MAX_NCE
                            + ". NCE should be in standard units (e.g. 27.0 for 27%).");
        }
    }

    @Override
    public void close() {
        if (closed) {
            return;
        }
        closed = true;
        cartographerPredictor.close();
        chronologer.close();
    }

    private static CartographerMetadata loadCartographerMetadata(String resource) {
        ClassLoader loader = Thread.currentThread().getContextClassLoader();
        try (InputStream input = loader.getResourceAsStream(resource)) {
            if (input == null) {
                throw new IllegalArgumentException("Missing cartographer metadata resource: " + resource);
            }
            JsonNode root = OBJECT_MAPPER.readTree(input);
            int minCharge = root.path("min_precursor_charge").asInt(1);
            int maxCharge = root.path("max_precursor_charge").asInt(6);
            JsonNode outputShape = root.path("output_shape");
            int outputWidth = outputShape.size() > 1
                    ? outputShape.get(1).asInt(CartographerSpectrumDecoder.VECTOR_LENGTH)
                    : CartographerSpectrumDecoder.VECTOR_LENGTH;
            return new CartographerMetadata(minCharge, maxCharge, outputWidth);
        } catch (IOException e) {
            throw new IllegalStateException("Failed to parse cartographer metadata resource: " + resource, e);
        }
    }

    private static final class NormalizedRequest {
        private final String unimodPeptide;
        private final List<PrecursorCondition> conditions;

        private NormalizedRequest(String unimodPeptide, List<PrecursorCondition> conditions) {
            this.unimodPeptide = unimodPeptide;
            this.conditions = conditions;
        }
    }

    private static final class PredictionJob {
        private final String unimodPeptide;
        private final ParsedUnimodSequence parsedPeptide;
        private final long[] tokenArray;
        private final byte charge;
        private final double nce;
        private final float retentionTimeSeconds;

        private PredictionJob(
                String unimodPeptide,
                ParsedUnimodSequence parsedPeptide,
                long[] tokenArray,
                byte charge,
                double nce,
                float retentionTimeSeconds) {
            this.unimodPeptide = unimodPeptide;
            this.parsedPeptide = parsedPeptide;
            this.tokenArray = tokenArray;
            this.charge = charge;
            this.nce = nce;
            this.retentionTimeSeconds = retentionTimeSeconds;
        }
    }

    private static final class CartographerMetadata {
        private final int minPrecursorCharge;
        private final int maxPrecursorCharge;
        private final int outputWidth;

        private CartographerMetadata(int minPrecursorCharge, int maxPrecursorCharge, int outputWidth) {
            this.minPrecursorCharge = minPrecursorCharge;
            this.maxPrecursorCharge = maxPrecursorCharge;
            this.outputWidth = outputWidth;
        }
    }
}

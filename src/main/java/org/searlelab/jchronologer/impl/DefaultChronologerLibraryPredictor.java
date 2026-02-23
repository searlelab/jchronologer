package org.searlelab.jchronologer.impl;

import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;
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
import org.searlelab.jchronologer.inference.ElectricianBatchPredictor;
import org.searlelab.jchronologer.inference.PeptideMassCalculator;
import org.searlelab.jchronologer.preprocessing.ChronologerPreprocessor;
import org.searlelab.jchronologer.preprocessing.PeptideSequenceConverter;
import org.searlelab.jchronologer.preprocessing.PeptideSequenceConverter.ParsedUnimodSequence;
import org.searlelab.jchronologer.preprocessing.PreprocessingMetadataLoader;
import org.searlelab.jchronologer.preprocessing.PreprocessingOutcome;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Default tandem predictor combining Chronologer RT and Cartographer MS2 inference.
 */
public final class DefaultChronologerLibraryPredictor implements ChronologerLibraryPredictor {

    private static final ObjectMapper OBJECT_MAPPER = new ObjectMapper();
    private static final Logger LOGGER = LoggerFactory.getLogger(DefaultChronologerLibraryPredictor.class);

    /** Minimum accepted NCE value (standard units, e.g. 20 for 20%). */
    private static final double MIN_NCE = 10.0;
    /** Maximum accepted NCE value (standard units, e.g. 50 for 50%). */
    private static final double MAX_NCE = 60.0;
    /** Divisor to convert standard NCE to model-normalized NCE (NCE / 100). */
    private static final double NCE_NORMALIZATION_DIVISOR = 100.0;
    /** Number of Electrician outputs corresponding to z=1..6. */
    private static final int ELECTRICIAN_CHARGE_STATE_COUNT = 6;

    private final ChronologerLibraryOptions options;
    private final Chronologer chronologer;
    private final ChronologerPreprocessor tokenPreprocessor;
    private final CartographerBatchPredictor cartographerPredictor;
    private final ElectricianBatchPredictor electricianPredictor;
    private final int minPrecursorCharge;
    private final int maxPrecursorCharge;
    private final int cartographerOutputWidth;
    private final int electricianOutputWidth;
    private final ExecutorService cartographerExecutor;
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
        this.electricianPredictor = new ElectricianBatchPredictor(options.getElectricianModelResource());

        CartographerMetadata cartographerMetadata = loadCartographerMetadata(options.getCartographerPreprocessingResource());
        this.minPrecursorCharge = cartographerMetadata.minPrecursorCharge;
        this.maxPrecursorCharge = cartographerMetadata.maxPrecursorCharge;
        this.cartographerOutputWidth = cartographerMetadata.outputWidth;

        ElectricianMetadata electricianMetadata = loadElectricianMetadata(options.getElectricianPreprocessingResource());
        this.electricianOutputWidth = electricianMetadata.outputWidth;
        if (this.electricianOutputWidth != ELECTRICIAN_CHARGE_STATE_COUNT) {
            throw new IllegalStateException(
                    "Unexpected Electrician metadata output width: "
                            + this.electricianOutputWidth
                            + " expected="
                            + ELECTRICIAN_CHARGE_STATE_COUNT);
        }
        this.cartographerExecutor = options.getInferenceThreads() > 1
                ? newCartographerExecutor(options.getInferenceThreads())
                : null;
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

        long totalStartNanos = System.nanoTime();
        long normalizeStartNanos = System.nanoTime();
        List<NormalizedRequest> normalizedRequests = new ArrayList<>(requests.size());
        LinkedHashMap<String, String> massByUnimod = new LinkedHashMap<>();

        for (int i = 0; i < requests.size(); i++) {
            LibraryPredictionRequest request = requests.get(i);
            if (request == null) {
                throw new IllegalArgumentException("Request at index " + i + " is null.");
            }

            String normalizedUnimod = PeptideSequenceConverter.normalizeToUnimod(
                    request.getPeptideSequence(),
                    options.getMassMatchEpsilon());
            String massEncoded = PeptideSequenceConverter.unimodToMassEncoded(normalizedUnimod);

            if (request.usesAutomaticChargeSelection()) {
                validateAutomaticRequest(request.getPrecursorNce(), request.getMinimumChargeProbability());
                normalizedRequests.add(NormalizedRequest.automatic(
                        normalizedUnimod,
                        request.getPrecursorNce(),
                        request.getMinimumChargeProbability()));
            } else {
                if (request.getPrecursorConditions().isEmpty()) {
                    throw new IllegalArgumentException("Request at index " + i + " has no precursor conditions.");
                }

                List<PrecursorCondition> validatedConditions = new ArrayList<>(request.getPrecursorConditions().size());
                for (PrecursorCondition condition : request.getPrecursorConditions()) {
                    validateCondition(condition);
                    validatedConditions.add(condition);
                }
                normalizedRequests.add(NormalizedRequest.explicit(normalizedUnimod, validatedConditions));
            }

            massByUnimod.putIfAbsent(normalizedUnimod, massEncoded);
        }
        long normalizeNanos = System.nanoTime() - normalizeStartNanos;

        long tokenizationStartNanos = System.nanoTime();
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
            // Rebuild mass-encoded with the fold applied.
            StringBuilder cartographerMassEncoded = new StringBuilder();
            if (!ntermMods.isEmpty()) {
                cartographerMassEncoded.append('[')
                        .append(String.format(Locale.US, "%+.6f", PeptideSequenceConverter.sumUnimodMass(ntermMods)))
                        .append(']');
            }
            for (int ri = 0; ri < parsed.getResidues().length(); ri++) {
                cartographerMassEncoded.append(parsed.getResidues().charAt(ri));
                List<String> mods = residueMods.get(ri);
                if (!mods.isEmpty()) {
                    cartographerMassEncoded.append('[')
                            .append(String.format(Locale.US, "%+.6f", PeptideSequenceConverter.sumUnimodMass(mods)))
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
        long tokenizationNanos = System.nanoTime() - tokenizationStartNanos;

        long rtStartNanos = System.nanoTime();
        PredictionResult rtResult = chronologer.predict(massEncodedUnique);
        long rtNanos = System.nanoTime() - rtStartNanos;
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

        long chargeStartNanos = System.nanoTime();
        Map<String, float[]> chargeDistributionByUnimod = inferChargeDistributions(normalizedRequests, tokenByUnimod);
        long chargeNanos = System.nanoTime() - chargeStartNanos;

        long jobExpansionStartNanos = System.nanoTime();
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

            if (normalizedRequest.automaticChargeSelection) {
                float[] distribution = chargeDistributionByUnimod.get(unimod);
                if (distribution == null) {
                    throw new IllegalStateException("Missing Electrician distribution for peptide: " + unimod);
                }

                for (int charge = minPrecursorCharge; charge <= maxPrecursorCharge; charge++) {
                    int chargeIndex = charge - 1;
                    if (chargeIndex < 0 || chargeIndex >= distribution.length) {
                        continue;
                    }
                    if (distribution[chargeIndex] >= normalizedRequest.minimumChargeProbability) {
                        jobs.add(new PredictionJob(
                                unimod,
                                parsed,
                                tokens,
                                (byte) charge,
                                normalizedRequest.precursorNce,
                                predHi * 60.0f));
                    }
                }
            } else {
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
        }
        long jobExpansionNanos = System.nanoTime() - jobExpansionStartNanos;

        int cartographerBatchCount = jobs.isEmpty()
                ? 0
                : (int) Math.ceil(jobs.size() / (double) options.getCartographerBatchSize());
        BatchPredictionResult[] orderedPredictionBatches = new BatchPredictionResult[cartographerBatchCount];
        List<Future<BatchPredictionResult>> futures = new ArrayList<>(cartographerBatchCount);

        if (!jobs.isEmpty()) {
            int batchIndex = 0;
            for (int start = 0; start < jobs.size(); start += options.getCartographerBatchSize()) {
                int end = Math.min(start + options.getCartographerBatchSize(), jobs.size());
                final int startFinal = start;
                final int endFinal = end;
                final int batchIndexFinal = batchIndex;
                if (cartographerExecutor == null || cartographerBatchCount == 1) {
                    orderedPredictionBatches[batchIndexFinal] = predictCartographerBatch(jobs, startFinal, endFinal, batchIndexFinal);
                } else {
                    futures.add(cartographerExecutor.submit(
                            () -> predictCartographerBatch(jobs, startFinal, endFinal, batchIndexFinal)));
                }
                batchIndex++;
            }
        }

        if (!futures.isEmpty()) {
            try {
                for (Future<BatchPredictionResult> future : futures) {
                    BatchPredictionResult batchResult = future.get();
                    orderedPredictionBatches[batchResult.batchIndex] = batchResult;
                }
            } catch (InterruptedException e) {
                Thread.currentThread().interrupt();
                futures.forEach(future -> future.cancel(true));
                throw new IllegalStateException("Cartographer inference was interrupted.", e);
            } catch (ExecutionException e) {
                futures.forEach(future -> future.cancel(true));
                Throwable cause = e.getCause() == null ? e : e.getCause();
                throw new IllegalStateException("Failed to run Cartographer inference.", cause);
            }
        }

        long cartographerInferenceNanos = 0L;
        List<ChronologerLibraryEntry> results = new ArrayList<>(jobs.size());
        long decodeStartNanos = System.nanoTime();
        for (int batchIndex = 0; batchIndex < orderedPredictionBatches.length; batchIndex++) {
            BatchPredictionResult batchResult = orderedPredictionBatches[batchIndex];
            if (batchResult == null) {
                throw new IllegalStateException("Missing Cartographer batch result for batch index " + batchIndex);
            }
            cartographerInferenceNanos += batchResult.inferenceNanos;
            for (int i = 0; i < batchResult.batchSize; i++) {
                PredictionJob job = jobs.get(batchResult.startJobIndex + i);
                float[] vector = batchResult.predictions[i];
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
        long decodeNanos = System.nanoTime() - decodeStartNanos;

        LOGGER.info(
                "Library predict completed: requests={}, uniquePeptides={}, jobs={}, cartographerBatches={}, normalizeMs={}, tokenizationMs={}, rtMs={}, chargeMs={}, jobExpandMs={}, cartographerInferenceMs={}, decodeMs={}, totalMs={}",
                requests.size(),
                massByUnimod.size(),
                jobs.size(),
                cartographerBatchCount,
                elapsedMillis(normalizeNanos),
                elapsedMillis(tokenizationNanos),
                elapsedMillis(rtNanos),
                elapsedMillis(chargeNanos),
                elapsedMillis(jobExpansionNanos),
                elapsedMillis(cartographerInferenceNanos),
                elapsedMillis(decodeNanos),
                elapsedMillis(System.nanoTime() - totalStartNanos));

        return results;
    }

    private BatchPredictionResult predictCartographerBatch(
            List<PredictionJob> jobs,
            int start,
            int end,
            int batchIndex) {
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

        long inferenceStartNanos = System.nanoTime();
        float[][] predictions = cartographerPredictor.predict(tokenBatch, chargeBatch, nceBatch);
        long inferenceNanos = System.nanoTime() - inferenceStartNanos;
        return new BatchPredictionResult(batchIndex, start, batch.size(), predictions, inferenceNanos);
    }

    private static ExecutorService newCartographerExecutor(int threads) {
        AtomicInteger threadCounter = new AtomicInteger(1);
        return Executors.newFixedThreadPool(threads, runnable -> {
            Thread thread = new Thread(runnable);
            thread.setName("cartographer-inference-" + threadCounter.getAndIncrement());
            thread.setDaemon(true);
            return thread;
        });
    }

    private static void shutdownCartographerExecutor(ExecutorService executor) {
        if (executor == null) {
            return;
        }
        executor.shutdown();
        try {
            if (!executor.awaitTermination(30, TimeUnit.SECONDS)) {
                executor.shutdownNow();
            }
        } catch (InterruptedException e) {
            executor.shutdownNow();
            Thread.currentThread().interrupt();
        }
    }

    private static long elapsedMillis(long nanos) {
        return nanos / 1_000_000L;
    }

    private Map<String, float[]> inferChargeDistributions(
            List<NormalizedRequest> normalizedRequests,
            Map<String, long[]> tokenByUnimod) {
        LinkedHashMap<String, long[]> tokensByAutoUnimod = new LinkedHashMap<>();
        for (NormalizedRequest request : normalizedRequests) {
            if (request.automaticChargeSelection) {
                tokensByAutoUnimod.putIfAbsent(request.unimodPeptide, tokenByUnimod.get(request.unimodPeptide));
            }
        }
        if (tokensByAutoUnimod.isEmpty()) {
            return Map.of();
        }

        long[][] tokenBatch = new long[tokensByAutoUnimod.size()][];
        List<String> unimodBatch = new ArrayList<>(tokensByAutoUnimod.size());
        int index = 0;
        for (Map.Entry<String, long[]> entry : tokensByAutoUnimod.entrySet()) {
            unimodBatch.add(entry.getKey());
            tokenBatch[index++] = entry.getValue();
        }

        float[][] distributions = electricianPredictor.predict(tokenBatch);
        if (distributions.length != unimodBatch.size()) {
            throw new IllegalStateException(
                    "Unexpected Electrician batch size: " + distributions.length + " expected=" + unimodBatch.size());
        }

        LinkedHashMap<String, float[]> normalizedByUnimod = new LinkedHashMap<>();
        for (int i = 0; i < unimodBatch.size(); i++) {
            float[] distribution = distributions[i];
            if (distribution.length != electricianOutputWidth) {
                throw new IllegalStateException(
                        "Unexpected Electrician output width: " + distribution.length + " expected=" + electricianOutputWidth);
            }
            normalizedByUnimod.put(unimodBatch.get(i), normalizeChargeDistribution(distribution));
        }
        return normalizedByUnimod;
    }

    private float[] normalizeChargeDistribution(float[] distribution) {
        float[] normalized = new float[distribution.length];
        double sum = 0.0;
        for (int i = 0; i < distribution.length; i++) {
            float value = distribution[i];
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
        validateNce(condition.getPrecursorNce());
    }

    private void validateAutomaticRequest(double nce, double minimumChargeProbability) {
        validateNce(nce);
        if (!Double.isFinite(minimumChargeProbability)) {
            throw new IllegalArgumentException("Minimum charge probability must be a finite number.");
        }
        if (minimumChargeProbability < 0.0 || minimumChargeProbability > 1.0) {
            throw new IllegalArgumentException(
                    "Minimum charge probability "
                            + minimumChargeProbability
                            + " is outside the supported range 0.0-1.0.");
        }
    }

    private void validateNce(double nce) {
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
        shutdownCartographerExecutor(cartographerExecutor);
        electricianPredictor.close();
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

    private static ElectricianMetadata loadElectricianMetadata(String resource) {
        ClassLoader loader = Thread.currentThread().getContextClassLoader();
        try (InputStream input = loader.getResourceAsStream(resource)) {
            if (input == null) {
                throw new IllegalArgumentException("Missing electrician metadata resource: " + resource);
            }
            JsonNode root = OBJECT_MAPPER.readTree(input);
            JsonNode outputShape = root.path("output_shape");
            int outputWidth = outputShape.size() > 1
                    ? outputShape.get(1).asInt(ELECTRICIAN_CHARGE_STATE_COUNT)
                    : ELECTRICIAN_CHARGE_STATE_COUNT;
            return new ElectricianMetadata(outputWidth);
        } catch (IOException e) {
            throw new IllegalStateException("Failed to parse electrician metadata resource: " + resource, e);
        }
    }

    private static final class NormalizedRequest {
        private final String unimodPeptide;
        private final List<PrecursorCondition> conditions;
        private final boolean automaticChargeSelection;
        private final double precursorNce;
        private final double minimumChargeProbability;

        private NormalizedRequest(
                String unimodPeptide,
                List<PrecursorCondition> conditions,
                boolean automaticChargeSelection,
                double precursorNce,
                double minimumChargeProbability) {
            this.unimodPeptide = unimodPeptide;
            this.conditions = conditions;
            this.automaticChargeSelection = automaticChargeSelection;
            this.precursorNce = precursorNce;
            this.minimumChargeProbability = minimumChargeProbability;
        }

        private static NormalizedRequest explicit(String unimodPeptide, List<PrecursorCondition> conditions) {
            return new NormalizedRequest(unimodPeptide, conditions, false, Double.NaN, Double.NaN);
        }

        private static NormalizedRequest automatic(
                String unimodPeptide,
                double precursorNce,
                double minimumChargeProbability) {
            return new NormalizedRequest(unimodPeptide, List.of(), true, precursorNce, minimumChargeProbability);
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

    private static final class BatchPredictionResult {
        private final int batchIndex;
        private final int startJobIndex;
        private final int batchSize;
        private final float[][] predictions;
        private final long inferenceNanos;

        private BatchPredictionResult(
                int batchIndex,
                int startJobIndex,
                int batchSize,
                float[][] predictions,
                long inferenceNanos) {
            this.batchIndex = batchIndex;
            this.startJobIndex = startJobIndex;
            this.batchSize = batchSize;
            this.predictions = predictions;
            this.inferenceNanos = inferenceNanos;
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

    private static final class ElectricianMetadata {
        private final int outputWidth;

        private ElectricianMetadata(int outputWidth) {
            this.outputWidth = outputWidth;
        }
    }
}

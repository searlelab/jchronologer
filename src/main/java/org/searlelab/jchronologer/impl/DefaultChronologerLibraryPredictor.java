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
import java.util.Optional;
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
    private final CartographerBatchPredictor sculptorPredictor;
    private final SculptorTokenizer sculptorTokenizer;
    private final int sculptorChargeStateCount;
    private final int sculptorOutputWidth;
    private final float sculptorCCSMean;
    private final float sculptorCCSStd;
    private final boolean ccsPredictionEnabled;
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
                .sculptorModelResource(options.getSculptorModelResource())
                .sculptorPreprocessingResource(options.getSculptorPreprocessingResource())
                .ccsPredictionEnabled(options.isCCSPredictionEnabled())
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
        this.ccsPredictionEnabled = options.isCCSPredictionEnabled();
        if (ccsPredictionEnabled) {
            SculptorMetadata sculptorMetadata = loadSculptorMetadata(options.getSculptorPreprocessingResource());
            this.sculptorTokenizer = new SculptorTokenizer(sculptorMetadata);
            this.sculptorChargeStateCount = sculptorMetadata.chargeStateCount;
            this.sculptorOutputWidth = sculptorMetadata.outputWidth;
            this.sculptorCCSMean = (float) sculptorMetadata.ccsMean;
            this.sculptorCCSStd = (float) sculptorMetadata.ccsStd;
            if (this.sculptorOutputWidth != 1) {
                throw new IllegalStateException(
                        "Unexpected Sculptor metadata output width: " + this.sculptorOutputWidth + " expected=1");
            }
            this.sculptorPredictor = new CartographerBatchPredictor(options.getSculptorModelResource());
        } else {
            this.sculptorTokenizer = null;
            this.sculptorChargeStateCount = 0;
            this.sculptorOutputWidth = 0;
            this.sculptorCCSMean = 0.0f;
            this.sculptorCCSStd = 1.0f;
            this.sculptorPredictor = null;
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
                                predHi * 60.0f,
                                Optional.of(distribution[chargeIndex])));
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
                            predHi * 60.0f,
                            Optional.empty()));
                }
            }
        }
        long jobExpansionNanos = System.nanoTime() - jobExpansionStartNanos;
        long ccsStartNanos = System.nanoTime();
        Map<String, Float> ccsByUnimodAndCharge = inferCCSByJob(jobs);
        long ccsNanos = System.nanoTime() - ccsStartNanos;

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
                Optional<Float> ccs = Optional.ofNullable(
                        ccsByUnimodAndCharge.get(ccsKey(job.unimodPeptide, job.charge)));
                results.add(new ChronologerLibraryEntry(
                        job.unimodPeptide,
                        job.charge,
                        job.nce,
                        precursorMz,
                        job.retentionTimeSeconds,
                        decoded.getMassArray(),
                        decoded.getIntensityArray(),
                        decoded.getIonTypeArray(),
                        ccs,
                        job.chargeProbability));
            }
        }
        long decodeNanos = System.nanoTime() - decodeStartNanos;

        LOGGER.info(
                "Library predict completed: requests={}, uniquePeptides={}, jobs={}, cartographerBatches={}, normalizeMs={}, tokenizationMs={}, rtMs={}, chargeMs={}, jobExpandMs={}, ccsMs={}, cartographerInferenceMs={}, decodeMs={}, totalMs={}",
                requests.size(),
                massByUnimod.size(),
                jobs.size(),
                cartographerBatchCount,
                elapsedMillis(normalizeNanos),
                elapsedMillis(tokenizationNanos),
                elapsedMillis(rtNanos),
                elapsedMillis(chargeNanos),
                elapsedMillis(jobExpansionNanos),
                elapsedMillis(ccsNanos),
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

    private Map<String, Float> inferCCSByJob(List<PredictionJob> jobs) {
        if (!ccsPredictionEnabled || jobs.isEmpty()) {
            return Map.of();
        }

        LinkedHashMap<String, CCSCandidate> candidatesByKey = new LinkedHashMap<>();
        for (PredictionJob job : jobs) {
            String key = ccsKey(job.unimodPeptide, job.charge);
            if (candidatesByKey.containsKey(key)) {
                continue;
            }
            if (job.charge < 1 || job.charge > sculptorChargeStateCount) {
                continue;
            }
            long[] tokenArray = sculptorTokenizer.tokenize(job.parsedPeptide);
            if (tokenArray == null) {
                continue;
            }
            candidatesByKey.put(key, new CCSCandidate(key, tokenArray, job.charge));
        }

        if (candidatesByKey.isEmpty()) {
            return Map.of();
        }

        long[][] tokenBatch = new long[candidatesByKey.size()][];
        float[][] chargeBatch = new float[candidatesByKey.size()][sculptorChargeStateCount];
        List<String> keys = new ArrayList<>(candidatesByKey.size());
        int index = 0;
        for (CCSCandidate candidate : candidatesByKey.values()) {
            keys.add(candidate.key);
            tokenBatch[index] = candidate.tokenArray;
            chargeBatch[index][candidate.charge - 1] = 1.0f;
            index++;
        }

        float[][] predictions = sculptorPredictor.predict(tokenBatch, chargeBatch);
        if (predictions.length != keys.size()) {
            throw new IllegalStateException(
                    "Unexpected Sculptor batch size: " + predictions.length + " expected=" + keys.size());
        }

        LinkedHashMap<String, Float> ccsByKey = new LinkedHashMap<>();
        for (int i = 0; i < predictions.length; i++) {
            float[] output = predictions[i];
            if (output.length != sculptorOutputWidth) {
                throw new IllegalStateException(
                        "Unexpected Sculptor output width: " + output.length + " expected=" + sculptorOutputWidth);
            }
            // Sculptor outputs normalized CCS (ccs_norm); convert back to physical CCS units.
            float ccs = (output[0] * sculptorCCSStd) + sculptorCCSMean;
            if (Float.isFinite(ccs)) {
                ccsByKey.put(keys.get(i), ccs);
            }
        }
        return ccsByKey;
    }

    private static String ccsKey(String unimodPeptide, byte charge) {
        return unimodPeptide + "|" + charge;
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
        if (sculptorPredictor != null) {
            sculptorPredictor.close();
        }
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

    private static SculptorMetadata loadSculptorMetadata(String resource) {
        ClassLoader loader = Thread.currentThread().getContextClassLoader();
        try (InputStream input = loader.getResourceAsStream(resource)) {
            if (input == null) {
                throw new IllegalArgumentException("Missing sculptor metadata resource: " + resource);
            }
            JsonNode root = OBJECT_MAPPER.readTree(input);

            int maxPeptideLength = root.path("max_peptide_len").asInt(0);
            if (maxPeptideLength <= 0) {
                throw new IllegalStateException("Invalid Sculptor max_peptide_len in resource: " + resource);
            }

            JsonNode aaToIntNode = root.path("aa_to_int");
            if (!aaToIntNode.isObject() || aaToIntNode.size() == 0) {
                throw new IllegalStateException("Invalid Sculptor aa_to_int in resource: " + resource);
            }
            LinkedHashMap<Character, Integer> tokenByResidue = new LinkedHashMap<>();
            aaToIntNode.fields().forEachRemaining(entry -> {
                String token = entry.getKey();
                if (token.length() != 1) {
                    return;
                }
                int value = entry.getValue().asInt(0);
                if (value <= 0) {
                    throw new IllegalStateException(
                            "Invalid Sculptor aa_to_int value for token '" + token + "' in resource: " + resource);
                }
                tokenByResidue.put(token.charAt(0), value);
            });
            if (tokenByResidue.isEmpty()) {
                throw new IllegalStateException("Sculptor aa_to_int has no single-character tokens in resource: " + resource);
            }

            int chargeStateCount = root.path("charge_dist_len").asInt(0);
            if (chargeStateCount <= 0) {
                JsonNode inputShapes = root.path("input_shapes");
                if (inputShapes.size() > 1 && inputShapes.get(1).size() > 1) {
                    chargeStateCount = inputShapes.get(1).get(1).asInt(0);
                }
            }
            if (chargeStateCount <= 0) {
                throw new IllegalStateException("Invalid Sculptor charge_dist_len in resource: " + resource);
            }

            JsonNode outputShape = root.path("output_shape");
            int outputWidth = outputShape.size() > 1 ? outputShape.get(1).asInt(1) : 1;
            double ccsMean = root.path("train_ccs_mean").asDouble(0.0);
            double ccsStd = root.path("train_ccs_std").asDouble(1.0);
            if (!Double.isFinite(ccsMean)) {
                ccsMean = 0.0;
            }
            if (!Double.isFinite(ccsStd) || ccsStd <= 0.0) {
                ccsStd = 1.0;
            }

            LinkedHashMap<String, Character> ntermTokenByUnimod = new LinkedHashMap<>();
            JsonNode ntermMap = root.path("nterm_unimod_map");
            if (!ntermMap.isObject()) {
                throw new IllegalStateException("Invalid Sculptor nterm_unimod_map in resource: " + resource);
            }
            ntermMap.fields().forEachRemaining(entry -> {
                String unimod = entry.getKey();
                String token = entry.getValue().asText();
                if (token.length() != 1) {
                    throw new IllegalStateException(
                            "Invalid Sculptor nterm token for " + unimod + " in resource: " + resource);
                }
                ntermTokenByUnimod.put(unimod, token.charAt(0));
            });

            LinkedHashMap<String, Character> residueTokenByResidueAndUnimod = new LinkedHashMap<>();
            JsonNode residueMap = root.path("residue_unimod_map");
            if (!residueMap.isArray()) {
                throw new IllegalStateException("Invalid Sculptor residue_unimod_map in resource: " + resource);
            }
            for (int i = 0; i < residueMap.size(); i++) {
                JsonNode entry = residueMap.get(i);
                String residue = entry.path("residue").asText();
                String unimod = entry.path("unimod").asText();
                String token = entry.path("token").asText();
                if (residue.length() != 1 || unimod.isBlank() || token.length() != 1) {
                    throw new IllegalStateException(
                            "Invalid Sculptor residue_unimod_map entry at index " + i + " in resource: " + resource);
                }
                residueTokenByResidueAndUnimod.put(residue + "|" + unimod, token.charAt(0));
            }

            return new SculptorMetadata(
                    maxPeptideLength,
                    tokenByResidue,
                    ntermTokenByUnimod,
                    residueTokenByResidueAndUnimod,
                    chargeStateCount,
                    outputWidth,
                    ccsMean,
                    ccsStd);
        } catch (IOException e) {
            throw new IllegalStateException("Failed to parse sculptor metadata resource: " + resource, e);
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
        private final Optional<Float> chargeProbability;

        private PredictionJob(
                String unimodPeptide,
                ParsedUnimodSequence parsedPeptide,
                long[] tokenArray,
                byte charge,
                double nce,
                float retentionTimeSeconds,
                Optional<Float> chargeProbability) {
            this.unimodPeptide = unimodPeptide;
            this.parsedPeptide = parsedPeptide;
            this.tokenArray = tokenArray;
            this.charge = charge;
            this.nce = nce;
            this.retentionTimeSeconds = retentionTimeSeconds;
            this.chargeProbability = chargeProbability;
        }
    }

    private static final class CCSCandidate {
        private final String key;
        private final long[] tokenArray;
        private final byte charge;

        private CCSCandidate(String key, long[] tokenArray, byte charge) {
            this.key = key;
            this.tokenArray = tokenArray;
            this.charge = charge;
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

    private static final class SculptorTokenizer {
        private final Map<Character, Integer> tokenByResidue;
        private final Map<String, Character> ntermTokenByUnimod;
        private final Map<String, Character> residueTokenByResidueAndUnimod;
        private final int maxPeptideLength;
        private final int tokenArrayLength;

        private SculptorTokenizer(SculptorMetadata metadata) {
            this.tokenByResidue = metadata.tokenByResidue;
            this.ntermTokenByUnimod = metadata.ntermTokenByUnimod;
            this.residueTokenByResidueAndUnimod = metadata.residueTokenByResidueAndUnimod;
            this.maxPeptideLength = metadata.maxPeptideLength;
            this.tokenArrayLength = metadata.maxPeptideLength + 2;
        }

        private long[] tokenize(ParsedUnimodSequence parsedPeptide) {
            String residues = parsedPeptide.getResidues();
            if (residues.length() > maxPeptideLength) {
                return null;
            }

            List<String> ntermMods = parsedPeptide.getPositionMods().get(0);
            if (ntermMods.size() > 1) {
                return null;
            }

            List<String> ctermMods = parsedPeptide.getPositionMods().get(residues.length() + 1);
            if (!ctermMods.isEmpty()) {
                return null;
            }

            char ntermToken = '-';
            if (!ntermMods.isEmpty()) {
                Character mapped = ntermTokenByUnimod.get(ntermMods.get(0));
                if (mapped == null) {
                    return null;
                }
                ntermToken = mapped;
            }

            long[] tokenArray = new long[tokenArrayLength];
            Integer ntermValue = tokenByResidue.get(ntermToken);
            if (ntermValue == null) {
                return null;
            }
            tokenArray[0] = ntermValue;

            for (int i = 0; i < residues.length(); i++) {
                List<String> residueMods = parsedPeptide.getPositionMods().get(i + 1);
                char residueToken;
                if (residueMods.isEmpty()) {
                    residueToken = residues.charAt(i);
                } else if (residueMods.size() == 1) {
                    Character mapped = residueTokenByResidueAndUnimod.get(
                            residues.charAt(i) + "|" + residueMods.get(0));
                    if (mapped == null) {
                        return null;
                    }
                    residueToken = mapped;
                } else {
                    return null;
                }

                Integer residueValue = tokenByResidue.get(residueToken);
                if (residueValue == null) {
                    return null;
                }
                tokenArray[i + 1] = residueValue;
            }

            Integer endTokenValue = tokenByResidue.get('_');
            if (endTokenValue == null) {
                return null;
            }
            tokenArray[residues.length() + 1] = endTokenValue;
            return tokenArray;
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

    private static final class SculptorMetadata {
        private final int maxPeptideLength;
        private final Map<Character, Integer> tokenByResidue;
        private final Map<String, Character> ntermTokenByUnimod;
        private final Map<String, Character> residueTokenByResidueAndUnimod;
        private final int chargeStateCount;
        private final int outputWidth;
        private final double ccsMean;
        private final double ccsStd;

        private SculptorMetadata(
                int maxPeptideLength,
                Map<Character, Integer> tokenByResidue,
                Map<String, Character> ntermTokenByUnimod,
                Map<String, Character> residueTokenByResidueAndUnimod,
                int chargeStateCount,
                int outputWidth,
                double ccsMean,
                double ccsStd) {
            this.maxPeptideLength = maxPeptideLength;
            this.tokenByResidue = tokenByResidue;
            this.ntermTokenByUnimod = ntermTokenByUnimod;
            this.residueTokenByResidueAndUnimod = residueTokenByResidueAndUnimod;
            this.chargeStateCount = chargeStateCount;
            this.outputWidth = outputWidth;
            this.ccsMean = ccsMean;
            this.ccsStd = ccsStd;
        }
    }
}

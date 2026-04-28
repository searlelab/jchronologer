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
import java.util.Objects;
import java.util.Optional;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;
import org.searlelab.jchronologer.api.ChronologerLibraryEntry;
import org.searlelab.jchronologer.api.ChronologerLibraryOptions;
import org.searlelab.jchronologer.api.ChronologerLibraryPredictor;
import org.searlelab.jchronologer.api.LibraryPredictionRequest;
import org.searlelab.jchronologer.api.PrecursorCondition;
import org.searlelab.jchronologer.inference.ElectricianBatchPredictor;
import org.searlelab.jchronologer.inference.PeptideMassCalculator;
import org.searlelab.jchronologer.inference.ScoutBatchPredictor;
import org.searlelab.jchronologer.inference.ScoutBatchPredictor.ScoutBatchPredictionResult;
import org.searlelab.jchronologer.inference.ScoutSpectrumDecoder;
import org.searlelab.jchronologer.preprocessing.ChronologerPreprocessor;
import org.searlelab.jchronologer.preprocessing.PeptideSequenceConverter;
import org.searlelab.jchronologer.preprocessing.PeptideSequenceConverter.ParsedUnimodSequence;
import org.searlelab.jchronologer.preprocessing.PreprocessingMetadataLoader;
import org.searlelab.jchronologer.preprocessing.PreprocessingOutcome;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public final class DefaultScoutLibraryPredictor implements ChronologerLibraryPredictor {

    private static final ObjectMapper OBJECT_MAPPER = new ObjectMapper();
    private static final Logger LOGGER = LoggerFactory.getLogger(DefaultScoutLibraryPredictor.class);

    private static final double MIN_NCE = 10.0;
    private static final double MAX_NCE = 60.0;
    private static final double NCE_NORMALIZATION_DIVISOR = 100.0;
    private static final float SCOUT_IRT_MEAN = 11.587185f;
    private static final float SCOUT_IRT_STD = 5.611839f;
    private static final float SCOUT_CCS_MEAN = 472.605222f;
    private static final float SCOUT_CCS_STD = 113.822638f;
    private static final double SCOUT_SEED_BASIC_RESIDUE_COEFFICIENT = 0.35;
    private static final double SCOUT_SEED_LENGTH_COEFFICIENT = 0.08;
    private static final int ELECTRICIAN_CHARGE_STATE_COUNT = 6;

    private final ChronologerLibraryOptions options;
    private final ChronologerPreprocessor tokenPreprocessor;
    private final ScoutBatchPredictor scoutPredictor;
    private final AutomaticChargeSelectionMode automaticChargeSelectionMode;
    private final ElectricianBatchPredictor electricianPredictor;
    private final int electricianOutputWidth;
    private final int minPrecursorCharge;
    private final int maxPrecursorCharge;
    private final int scoutMs2OutputWidth;
    private final int scoutIrtOutputWidth;
    private final int scoutCcsOutputWidth;
    private final int scoutChargeDistributionOutputWidth;
    private final ExecutorService scoutExecutor;
    private volatile boolean closed;

    public DefaultScoutLibraryPredictor(ChronologerLibraryOptions options) {
        this(options, AutomaticChargeSelectionMode.SCOUT_SEED_REUSE);
    }

    DefaultScoutLibraryPredictor(
            ChronologerLibraryOptions options,
            AutomaticChargeSelectionMode automaticChargeSelectionMode) {
        this.options = options;
        this.automaticChargeSelectionMode = automaticChargeSelectionMode;
        this.tokenPreprocessor = new ChronologerPreprocessor(
                PreprocessingMetadataLoader.loadFromClasspath(options.getCartographerPreprocessingResource()));
        this.scoutPredictor = new ScoutBatchPredictor(options.getScoutModelResource());

        ScoutMetadata scoutMetadata = loadScoutMetadata(options.getScoutPreprocessingResource());
        this.minPrecursorCharge = scoutMetadata.minPrecursorCharge;
        this.maxPrecursorCharge = scoutMetadata.maxPrecursorCharge;
        this.scoutMs2OutputWidth = scoutMetadata.ms2OutputWidth;
        this.scoutIrtOutputWidth = scoutMetadata.irtOutputWidth;
        this.scoutCcsOutputWidth = scoutMetadata.ccsOutputWidth;
        this.scoutChargeDistributionOutputWidth = scoutMetadata.chargeDistributionOutputWidth;
        if (automaticChargeSelectionMode == AutomaticChargeSelectionMode.ELECTRICIAN_SELECTION) {
            this.electricianPredictor = new ElectricianBatchPredictor(options.getElectricianModelResource());
            this.electricianOutputWidth = ELECTRICIAN_CHARGE_STATE_COUNT;
        } else {
            this.electricianPredictor = null;
            this.electricianOutputWidth = 0;
        }
        this.scoutExecutor = options.getInferenceThreads() > 1 ? newScoutExecutor(options.getInferenceThreads()) : null;
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
        LinkedHashMap<String, String> massByUnimod = new LinkedHashMap<>();
        List<NormalizedRequest> normalizedRequests = new ArrayList<>(requests.size());

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

        LinkedHashMap<String, long[]> tokenByUnimod = new LinkedHashMap<>();
        LinkedHashMap<String, ParsedUnimodSequence> parsedByUnimod = new LinkedHashMap<>();
        for (Map.Entry<String, String> entry : massByUnimod.entrySet()) {
            ParsedUnimodSequence parsed = PeptideSequenceConverter.parseNormalizedUnimod(entry.getKey());
            parsedByUnimod.put(entry.getKey(), parsed);
            PreprocessingOutcome outcome = tokenPreprocessor.preprocess(toMassEncodedTokenInput(parsed));
            if (!outcome.isAccepted()) {
                throw new IllegalArgumentException(
                        "Failed to tokenize peptide " + entry.getValue() + " (normalized " + entry.getKey() + "): "
                                + outcome.getRejectionReason());
            }
            tokenByUnimod.put(entry.getKey(), outcome.getTokenArray());
        }

        List<PredictionJob> finalJobs = new ArrayList<>();
        AutomaticChargeSelectionResult automaticSelection = selectAutomaticJobs(
                normalizedRequests,
                parsedByUnimod,
                tokenByUnimod,
                finalJobs);

        Map<AutomaticSeedKey, ScoutPredictionRecord> reusedSeedPredictions = automaticSelection.reusedSeedPredictions;
        List<PredictionJob> rerunJobs = automaticSelection.rerunJobs;
        int reusedSeedJobCount = automaticSelection.reusedSeedJobCount;

        Map<PredictionJobKey, ScoutPredictionRecord> rerunPredictions = predictJobs(rerunJobs);
        List<ChronologerLibraryEntry> results = new ArrayList<>(finalJobs.size());
        for (PredictionJob job : finalJobs) {
            ScoutPredictionRecord prediction;
            if (job.automaticSeedKey.isPresent() && job.reuseSeedPrediction) {
                prediction = reusedSeedPredictions.get(job.automaticSeedKey.get());
            } else {
                prediction = rerunPredictions.get(job.toPredictionJobKey());
            }
            if (prediction == null) {
                throw new IllegalStateException("Missing Scout prediction for job: " + job.toPredictionJobKey());
            }
            results.add(toLibraryEntry(job, prediction));
        }

        LOGGER.info(
                "Scout library predict completed: requests={}, uniquePeptides={}, jobs={}, seedJobs={}, reusedSeedJobs={}, rerunJobs={}, scoutBatches={}, totalMs={}",
                requests.size(),
                massByUnimod.size(),
                finalJobs.size(),
                automaticSelection.seedJobCount,
                reusedSeedJobCount,
                rerunJobs.size(),
                finalJobs.isEmpty() ? 0 : (int) Math.ceil(finalJobs.size() / (double) options.getBatchSize()),
                elapsedMillis(System.nanoTime() - totalStartNanos));
        if (!automaticSelection.missedSeedDiagnostics.isEmpty()) {
            LOGGER.info("Temporary Scout missed-seed diagnostics (count={}):\n{}",
                    automaticSelection.missedSeedDiagnostics.size(),
                    String.join(System.lineSeparator(), automaticSelection.missedSeedDiagnostics));
        }
        return results;
    }

    private AutomaticChargeSelectionResult selectAutomaticJobs(
            List<NormalizedRequest> normalizedRequests,
            LinkedHashMap<String, ParsedUnimodSequence> parsedByUnimod,
            LinkedHashMap<String, long[]> tokenByUnimod,
            List<PredictionJob> finalJobs) {
        if (automaticChargeSelectionMode == AutomaticChargeSelectionMode.ELECTRICIAN_SELECTION) {
            return selectAutomaticJobsWithElectrician(normalizedRequests, parsedByUnimod, tokenByUnimod, finalJobs);
        }
        return selectAutomaticJobsWithScoutSeed(normalizedRequests, parsedByUnimod, tokenByUnimod, finalJobs);
    }

    private AutomaticChargeSelectionResult selectAutomaticJobsWithScoutSeed(
            List<NormalizedRequest> normalizedRequests,
            LinkedHashMap<String, ParsedUnimodSequence> parsedByUnimod,
            LinkedHashMap<String, long[]> tokenByUnimod,
            List<PredictionJob> finalJobs) {
        LinkedHashMap<AutomaticSeedKey, SeedPredictionContext> automaticContexts = new LinkedHashMap<>();
        for (NormalizedRequest normalizedRequest : normalizedRequests) {
            ParsedUnimodSequence parsed = parsedByUnimod.get(normalizedRequest.unimodPeptide);
            long[] tokens = tokenByUnimod.get(normalizedRequest.unimodPeptide);
            if (normalizedRequest.automaticChargeSelection) {
                AutomaticSeedKey key = new AutomaticSeedKey(normalizedRequest.unimodPeptide, normalizedRequest.precursorNce);
                automaticContexts.computeIfAbsent(
                        key,
                        ignored -> new SeedPredictionContext(
                                key,
                                parsed,
                                tokens,
                                (byte) selectSeedCharge(parsed),
                                normalizedRequest.minimumChargeProbability));
            } else {
                for (PrecursorCondition condition : normalizedRequest.conditions) {
                    finalJobs.add(PredictionJob.explicit(
                            normalizedRequest.unimodPeptide,
                            parsed,
                            tokens,
                            condition.getPrecursorCharge(),
                            condition.getPrecursorNce()));
                }
            }
        }

        Map<AutomaticSeedKey, ScoutPredictionRecord> seedPredictionByKey = predictSeedPasses(automaticContexts);
        List<String> missedSeedDiagnostics = new ArrayList<>();
        for (NormalizedRequest normalizedRequest : normalizedRequests) {
            if (!normalizedRequest.automaticChargeSelection) {
                continue;
            }
            AutomaticSeedKey key = new AutomaticSeedKey(normalizedRequest.unimodPeptide, normalizedRequest.precursorNce);
            SeedPredictionContext context = automaticContexts.get(key);
            ScoutPredictionRecord seedPrediction = seedPredictionByKey.get(key);
            if (context == null || seedPrediction == null) {
                throw new IllegalStateException("Missing Scout seed prediction for peptide: " + normalizedRequest.unimodPeptide);
            }
            float[] normalizedDistribution = normalizeChargeDistribution(seedPrediction.chargeDist);
            boolean seedChargeSelected = false;
            for (int charge = minPrecursorCharge; charge <= maxPrecursorCharge; charge++) {
                int distributionIndex = charge - minPrecursorCharge;
                float probability = distributionIndex >= 0 && distributionIndex < normalizedDistribution.length
                        ? normalizedDistribution[distributionIndex]
                        : 0.0f;
                if (probability < normalizedRequest.minimumChargeProbability) {
                    continue;
                }
                if (charge == context.seedCharge) {
                    seedChargeSelected = true;
                }
                finalJobs.add(PredictionJob.automatic(
                        key,
                        normalizedRequest.unimodPeptide,
                        context.parsedPeptide,
                        context.tokenArray,
                        (byte) charge,
                        normalizedRequest.precursorNce,
                        probability,
                        charge == context.seedCharge));
            }
            if (!seedChargeSelected) {
                missedSeedDiagnostics.add(formatMissedSeedDiagnostic(
                        normalizedRequest.unimodPeptide,
                        context.seedCharge,
                        normalizedRequest.minimumChargeProbability,
                        normalizedDistribution));
            }
        }

        List<PredictionJob> rerunJobs = new ArrayList<>();
        int reusedSeedJobCount = 0;
        for (PredictionJob job : finalJobs) {
            if (job.automaticSeedKey.isPresent() && job.reuseSeedPrediction) {
                reusedSeedJobCount++;
                continue;
            }
            rerunJobs.add(job);
        }
        return new AutomaticChargeSelectionResult(
                automaticContexts.size(),
                reusedSeedJobCount,
                rerunJobs,
                seedPredictionByKey,
                missedSeedDiagnostics);
    }

    private AutomaticChargeSelectionResult selectAutomaticJobsWithElectrician(
            List<NormalizedRequest> normalizedRequests,
            LinkedHashMap<String, ParsedUnimodSequence> parsedByUnimod,
            LinkedHashMap<String, long[]> tokenByUnimod,
            List<PredictionJob> finalJobs) {
        Map<String, float[]> chargeDistributionByUnimod = inferElectricianChargeDistributions(normalizedRequests, tokenByUnimod);
        for (NormalizedRequest normalizedRequest : normalizedRequests) {
            ParsedUnimodSequence parsed = parsedByUnimod.get(normalizedRequest.unimodPeptide);
            long[] tokens = tokenByUnimod.get(normalizedRequest.unimodPeptide);
            if (!normalizedRequest.automaticChargeSelection) {
                for (PrecursorCondition condition : normalizedRequest.conditions) {
                    finalJobs.add(PredictionJob.explicit(
                            normalizedRequest.unimodPeptide,
                            parsed,
                            tokens,
                            condition.getPrecursorCharge(),
                            condition.getPrecursorNce()));
                }
                continue;
            }
            float[] distribution = chargeDistributionByUnimod.get(normalizedRequest.unimodPeptide);
            if (distribution == null) {
                throw new IllegalStateException("Missing Electrician distribution for peptide: " + normalizedRequest.unimodPeptide);
            }
            for (int charge = minPrecursorCharge; charge <= maxPrecursorCharge; charge++) {
                int chargeIndex = charge - 1;
                if (chargeIndex < 0 || chargeIndex >= distribution.length) {
                    continue;
                }
                if (distribution[chargeIndex] >= normalizedRequest.minimumChargeProbability) {
                    finalJobs.add(PredictionJob.electricianAutomatic(
                            normalizedRequest.unimodPeptide,
                            parsed,
                            tokens,
                            (byte) charge,
                            normalizedRequest.precursorNce,
                            distribution[chargeIndex]));
                }
            }
        }
        return new AutomaticChargeSelectionResult(
                0,
                0,
                new ArrayList<>(finalJobs),
                Map.of(),
                List.of());
    }

    private Map<AutomaticSeedKey, ScoutPredictionRecord> predictSeedPasses(
            LinkedHashMap<AutomaticSeedKey, SeedPredictionContext> contexts) {
        List<PredictionJob> seedJobs = new ArrayList<>(contexts.size());
        for (SeedPredictionContext context : contexts.values()) {
            seedJobs.add(PredictionJob.seed(
                    context.seedKey,
                    context.seedKey.unimodPeptide,
                    context.parsedPeptide,
                    context.tokenArray,
                    context.seedCharge,
                    context.seedKey.precursorNce));
        }
        Map<PredictionJobKey, ScoutPredictionRecord> jobPredictions = predictJobs(seedJobs);
        LinkedHashMap<AutomaticSeedKey, ScoutPredictionRecord> predictionsByKey = new LinkedHashMap<>();
        for (PredictionJob job : seedJobs) {
            ScoutPredictionRecord prediction = jobPredictions.get(job.toPredictionJobKey());
            if (prediction == null) {
                throw new IllegalStateException("Missing Scout seed prediction for key: " + job.automaticSeedKey.get());
            }
            predictionsByKey.put(job.automaticSeedKey.get(), prediction);
        }
        return predictionsByKey;
    }

    private Map<PredictionJobKey, ScoutPredictionRecord> predictJobs(List<PredictionJob> jobs) {
        if (jobs.isEmpty()) {
            return Map.of();
        }
        int batchCount = (int) Math.ceil(jobs.size() / (double) options.getBatchSize());
        BatchPredictionResult[] orderedBatches = new BatchPredictionResult[batchCount];
        List<Future<BatchPredictionResult>> futures = new ArrayList<>(batchCount);

        int batchIndex = 0;
        for (int start = 0; start < jobs.size(); start += options.getBatchSize()) {
            int end = Math.min(start + options.getBatchSize(), jobs.size());
            final int startFinal = start;
            final int endFinal = end;
            final int batchIndexFinal = batchIndex;
            if (scoutExecutor == null || batchCount == 1) {
                orderedBatches[batchIndexFinal] = predictScoutBatch(jobs, startFinal, endFinal, batchIndexFinal);
            } else {
                futures.add(scoutExecutor.submit(() -> predictScoutBatch(jobs, startFinal, endFinal, batchIndexFinal)));
            }
            batchIndex++;
        }

        if (!futures.isEmpty()) {
            try {
                for (Future<BatchPredictionResult> future : futures) {
                    BatchPredictionResult batchResult = future.get();
                    orderedBatches[batchResult.batchIndex] = batchResult;
                }
            } catch (InterruptedException e) {
                Thread.currentThread().interrupt();
                futures.forEach(future -> future.cancel(true));
                throw new IllegalStateException("Scout inference was interrupted.", e);
            } catch (ExecutionException e) {
                futures.forEach(future -> future.cancel(true));
                Throwable cause = e.getCause() == null ? e : e.getCause();
                throw new IllegalStateException("Failed to run Scout inference.", cause);
            }
        }

        LinkedHashMap<PredictionJobKey, ScoutPredictionRecord> predictionByJob = new LinkedHashMap<>();
        for (BatchPredictionResult batchResult : orderedBatches) {
            if (batchResult == null) {
                throw new IllegalStateException("Missing Scout batch result.");
            }
            for (int i = 0; i < batchResult.batchSize; i++) {
                PredictionJob job = jobs.get(batchResult.startJobIndex + i);
                float[] ms2 = batchResult.ms2[i];
                float[] irt = batchResult.irt[i];
                float[] ccs = batchResult.ccs[i];
                float[] chargeDist = batchResult.chargeDist[i];
                validatePredictionShapes(ms2, irt, ccs, chargeDist);
                predictionByJob.put(
                        job.toPredictionJobKey(),
                        new ScoutPredictionRecord(ms2, irt, ccs, chargeDist));
            }
        }
        return predictionByJob;
    }

    private BatchPredictionResult predictScoutBatch(List<PredictionJob> jobs, int start, int end, int batchIndex) {
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
        ScoutBatchPredictionResult predictions = scoutPredictor.predict(tokenBatch, chargeBatch, nceBatch);
        long inferenceNanos = System.nanoTime() - inferenceStartNanos;
        return new BatchPredictionResult(
                batchIndex,
                start,
                batch.size(),
                predictions.getMs2(),
                predictions.getIrt(),
                predictions.getCcs(),
                predictions.getChargeDist(),
                inferenceNanos);
    }

    private Map<String, float[]> inferElectricianChargeDistributions(
            List<NormalizedRequest> normalizedRequests,
            Map<String, long[]> tokenByUnimod) {
        if (electricianPredictor == null) {
            throw new IllegalStateException("Electrician predictor is not enabled for this Scout predictor.");
        }
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

    private ChronologerLibraryEntry toLibraryEntry(PredictionJob job, ScoutPredictionRecord prediction) {
        float denormalizedIrt = (prediction.irt[0] * SCOUT_IRT_STD) + SCOUT_IRT_MEAN;
        float retentionTimeSeconds = denormalizedIrt * 60.0f;
        float denormalizedCcs = (prediction.ccs[0] * SCOUT_CCS_STD) + SCOUT_CCS_MEAN;
        ScoutSpectrumDecoder.DecodedSpectrum decoded = ScoutSpectrumDecoder.decode(
                job.parsedPeptide,
                job.charge,
                prediction.ms2,
                options.getMinimumReportedIntensity());
        double precursorMz = PeptideMassCalculator.calculatePrecursorMz(job.parsedPeptide, job.charge);
        return new ChronologerLibraryEntry(
                job.unimodPeptide,
                job.charge,
                job.nce,
                precursorMz,
                retentionTimeSeconds,
                decoded.getMassArray(),
                decoded.getIntensityArray(),
                decoded.getIonTypeArray(),
                Float.isFinite(denormalizedCcs) ? Optional.of(denormalizedCcs) : Optional.empty(),
                job.chargeProbability);
    }

    private void validatePredictionShapes(float[] ms2, float[] irt, float[] ccs, float[] chargeDist) {
        if (ms2.length != scoutMs2OutputWidth) {
            throw new IllegalStateException(
                    "Unexpected Scout MS2 output width: " + ms2.length + " expected=" + scoutMs2OutputWidth);
        }
        if (irt.length != scoutIrtOutputWidth) {
            throw new IllegalStateException(
                    "Unexpected Scout IRT output width: " + irt.length + " expected=" + scoutIrtOutputWidth);
        }
        if (ccs.length != scoutCcsOutputWidth) {
            throw new IllegalStateException(
                    "Unexpected Scout CCS output width: " + ccs.length + " expected=" + scoutCcsOutputWidth);
        }
        if (chargeDist.length != scoutChargeDistributionOutputWidth) {
            throw new IllegalStateException(
                    "Unexpected Scout charge_dist output width: "
                            + chargeDist.length
                            + " expected="
                            + scoutChargeDistributionOutputWidth);
        }
    }

    private float[] normalizeChargeDistribution(float[] rawDistribution) {
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

    private String formatMissedSeedDiagnostic(
            String unimodPeptide,
            byte seedCharge,
            double minimumChargeProbability,
            float[] normalizedDistribution) {
        StringBuilder distribution = new StringBuilder();
        for (int i = 0; i < normalizedDistribution.length; i++) {
            if (i > 0) {
                distribution.append(", ");
            }
            distribution.append(i + minPrecursorCharge)
                    .append('=')
                    .append(String.format(Locale.US, "%.4f", normalizedDistribution[i]));
        }
        return String.format(
                Locale.US,
                "TEMP scout missed seed: peptide=%s seedCharge=%d minProbability=%.4f chargeDist=[%s]",
                unimodPeptide,
                seedCharge,
                minimumChargeProbability,
                distribution);
    }

    private int selectSeedCharge(ParsedUnimodSequence peptide) {
        int seedGuess = estimateLikelyCharge(peptide);
        return Math.max(minPrecursorCharge, Math.min(maxPrecursorCharge, seedGuess));
    }

    private int estimateLikelyCharge(ParsedUnimodSequence peptide) {
        String residues = peptide.getResidues();
        int basicResidueCount = 0;
        for (int i = 0; i < residues.length(); i++) {
            char residue = residues.charAt(i);
            if (residue == 'K' || residue == 'R') {
                basicResidueCount++;
            }
        }
        int length = residues.length();
        int heuristicCharge = (int) Math.round(
                1.0
                        + (SCOUT_SEED_BASIC_RESIDUE_COEFFICIENT * basicResidueCount)
                        + (SCOUT_SEED_LENGTH_COEFFICIENT * length));
        return Math.min(
                basicResidueCount + 1,
                Math.max(1, heuristicCharge));
    }

    private static String toMassEncodedTokenInput(ParsedUnimodSequence parsed) {
        List<String> ntermMods = new ArrayList<>(parsed.getPositionMods().get(0));
        List<List<String>> residueMods = new ArrayList<>(parsed.getResidues().length());
        for (int ri = 0; ri < parsed.getResidues().length(); ri++) {
            residueMods.add(new ArrayList<>(parsed.getPositionMods().get(ri + 1)));
        }
        PeptideSequenceConverter.foldFirstResiduePyrogluToNterm(parsed.getResidues(), ntermMods, residueMods);

        StringBuilder massEncoded = new StringBuilder();
        if (!ntermMods.isEmpty()) {
            massEncoded.append('[')
                    .append(String.format(Locale.US, "%+.6f", PeptideSequenceConverter.sumUnimodMass(ntermMods)))
                    .append(']');
        }
        for (int ri = 0; ri < parsed.getResidues().length(); ri++) {
            massEncoded.append(parsed.getResidues().charAt(ri));
            List<String> mods = residueMods.get(ri);
            if (!mods.isEmpty()) {
                massEncoded.append('[')
                        .append(String.format(Locale.US, "%+.6f", PeptideSequenceConverter.sumUnimodMass(mods)))
                        .append(']');
            }
        }
        return massEncoded.toString();
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
        shutdownScoutExecutor(scoutExecutor);
        if (electricianPredictor != null) {
            electricianPredictor.close();
        }
        scoutPredictor.close();
    }

    private static ExecutorService newScoutExecutor(int threads) {
        AtomicInteger threadCounter = new AtomicInteger(1);
        return Executors.newFixedThreadPool(threads, runnable -> {
            Thread thread = new Thread(runnable);
            thread.setName("scout-inference-" + threadCounter.getAndIncrement());
            thread.setDaemon(true);
            return thread;
        });
    }

    private static void shutdownScoutExecutor(ExecutorService executor) {
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

    private static ScoutMetadata loadScoutMetadata(String resource) {
        ClassLoader loader = Thread.currentThread().getContextClassLoader();
        try (InputStream input = loader.getResourceAsStream(resource)) {
            if (input == null) {
                throw new IllegalArgumentException("Missing Scout metadata resource: " + resource);
            }
            JsonNode root = OBJECT_MAPPER.readTree(input);
            int minCharge = root.path("min_precursor_charge").asInt(1);
            int maxCharge = root.path("max_precursor_charge").asInt(6);
            JsonNode outputShapes = root.path("output_shapes");
            if (!outputShapes.isArray() || outputShapes.size() < 4) {
                throw new IllegalStateException("Invalid Scout output_shapes in resource: " + resource);
            }
            int ms2Width = outputShapes.get(0).size() > 1
                    ? outputShapes.get(0).get(1).asInt(ScoutSpectrumDecoder.VECTOR_LENGTH)
                    : ScoutSpectrumDecoder.VECTOR_LENGTH;
            int irtWidth = outputShapes.get(1).size() > 1 ? outputShapes.get(1).get(1).asInt(1) : 1;
            int ccsWidth = outputShapes.get(2).size() > 1 ? outputShapes.get(2).get(1).asInt(1) : 1;
            int chargeDistributionWidth = outputShapes.get(3).size() > 1
                    ? outputShapes.get(3).get(1).asInt(maxCharge - minCharge + 1)
                    : maxCharge - minCharge + 1;
            int expectedChargeDistributionWidth = maxCharge - minCharge + 1;
            if (chargeDistributionWidth != expectedChargeDistributionWidth) {
                throw new IllegalStateException(
                        "Unexpected Scout charge_dist width: "
                                + chargeDistributionWidth
                                + " expected="
                                + expectedChargeDistributionWidth);
            }
            return new ScoutMetadata(minCharge, maxCharge, ms2Width, irtWidth, ccsWidth, chargeDistributionWidth);
        } catch (IOException e) {
            throw new IllegalStateException("Failed to parse Scout metadata resource: " + resource, e);
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

        private static NormalizedRequest automatic(String unimodPeptide, double precursorNce, double minimumChargeProbability) {
            return new NormalizedRequest(unimodPeptide, List.of(), true, precursorNce, minimumChargeProbability);
        }
    }

    enum AutomaticChargeSelectionMode {
        SCOUT_SEED_REUSE,
        ELECTRICIAN_SELECTION
    }

    private static final class PredictionJob {
        private final Optional<AutomaticSeedKey> automaticSeedKey;
        private final String unimodPeptide;
        private final ParsedUnimodSequence parsedPeptide;
        private final long[] tokenArray;
        private final byte charge;
        private final double nce;
        private final Optional<Float> chargeProbability;
        private final boolean reuseSeedPrediction;

        private PredictionJob(
                Optional<AutomaticSeedKey> automaticSeedKey,
                String unimodPeptide,
                ParsedUnimodSequence parsedPeptide,
                long[] tokenArray,
                byte charge,
                double nce,
                Optional<Float> chargeProbability,
                boolean reuseSeedPrediction) {
            this.automaticSeedKey = automaticSeedKey;
            this.unimodPeptide = unimodPeptide;
            this.parsedPeptide = parsedPeptide;
            this.tokenArray = tokenArray;
            this.charge = charge;
            this.nce = nce;
            this.chargeProbability = chargeProbability;
            this.reuseSeedPrediction = reuseSeedPrediction;
        }

        private static PredictionJob explicit(
                String unimodPeptide,
                ParsedUnimodSequence parsedPeptide,
                long[] tokenArray,
                byte charge,
                double nce) {
            return new PredictionJob(Optional.empty(), unimodPeptide, parsedPeptide, tokenArray, charge, nce, Optional.empty(), false);
        }

        private static PredictionJob automatic(
                AutomaticSeedKey automaticSeedKey,
                String unimodPeptide,
                ParsedUnimodSequence parsedPeptide,
                long[] tokenArray,
                byte charge,
                double nce,
                float chargeProbability,
                boolean reuseSeedPrediction) {
            return new PredictionJob(
                    Optional.of(automaticSeedKey),
                    unimodPeptide,
                    parsedPeptide,
                    tokenArray,
                    charge,
                    nce,
                    Optional.of(chargeProbability),
                    reuseSeedPrediction);
        }

        private static PredictionJob seed(
                AutomaticSeedKey automaticSeedKey,
                String unimodPeptide,
                ParsedUnimodSequence parsedPeptide,
                long[] tokenArray,
                byte charge,
                double nce) {
            return new PredictionJob(
                    Optional.of(automaticSeedKey),
                    unimodPeptide,
                    parsedPeptide,
                    tokenArray,
                    charge,
                    nce,
                    Optional.empty(),
                    false);
        }

        private static PredictionJob electricianAutomatic(
                String unimodPeptide,
                ParsedUnimodSequence parsedPeptide,
                long[] tokenArray,
                byte charge,
                double nce,
                float chargeProbability) {
            return new PredictionJob(
                    Optional.empty(),
                    unimodPeptide,
                    parsedPeptide,
                    tokenArray,
                    charge,
                    nce,
                    Optional.of(chargeProbability),
                    false);
        }

        private PredictionJobKey toPredictionJobKey() {
            return new PredictionJobKey(unimodPeptide, charge, nce);
        }
    }

    private static final class BatchPredictionResult {
        private final int batchIndex;
        private final int startJobIndex;
        private final int batchSize;
        private final float[][] ms2;
        private final float[][] irt;
        private final float[][] ccs;
        private final float[][] chargeDist;
        @SuppressWarnings("unused")
        private final long inferenceNanos;

        private BatchPredictionResult(
                int batchIndex,
                int startJobIndex,
                int batchSize,
                float[][] ms2,
                float[][] irt,
                float[][] ccs,
                float[][] chargeDist,
                long inferenceNanos) {
            this.batchIndex = batchIndex;
            this.startJobIndex = startJobIndex;
            this.batchSize = batchSize;
            this.ms2 = ms2;
            this.irt = irt;
            this.ccs = ccs;
            this.chargeDist = chargeDist;
            this.inferenceNanos = inferenceNanos;
        }
    }

    private static final class AutomaticChargeSelectionResult {
        private final int seedJobCount;
        private final int reusedSeedJobCount;
        private final List<PredictionJob> rerunJobs;
        private final Map<AutomaticSeedKey, ScoutPredictionRecord> reusedSeedPredictions;
        private final List<String> missedSeedDiagnostics;

        private AutomaticChargeSelectionResult(
                int seedJobCount,
                int reusedSeedJobCount,
                List<PredictionJob> rerunJobs,
                Map<AutomaticSeedKey, ScoutPredictionRecord> reusedSeedPredictions,
                List<String> missedSeedDiagnostics) {
            this.seedJobCount = seedJobCount;
            this.reusedSeedJobCount = reusedSeedJobCount;
            this.rerunJobs = rerunJobs;
            this.reusedSeedPredictions = reusedSeedPredictions;
            this.missedSeedDiagnostics = missedSeedDiagnostics;
        }
    }

    private static final class ScoutMetadata {
        private final int minPrecursorCharge;
        private final int maxPrecursorCharge;
        private final int ms2OutputWidth;
        private final int irtOutputWidth;
        private final int ccsOutputWidth;
        private final int chargeDistributionOutputWidth;

        private ScoutMetadata(
                int minPrecursorCharge,
                int maxPrecursorCharge,
                int ms2OutputWidth,
                int irtOutputWidth,
                int ccsOutputWidth,
                int chargeDistributionOutputWidth) {
            this.minPrecursorCharge = minPrecursorCharge;
            this.maxPrecursorCharge = maxPrecursorCharge;
            this.ms2OutputWidth = ms2OutputWidth;
            this.irtOutputWidth = irtOutputWidth;
            this.ccsOutputWidth = ccsOutputWidth;
            this.chargeDistributionOutputWidth = chargeDistributionOutputWidth;
        }
    }

    private static final class SeedPredictionContext {
        private final AutomaticSeedKey seedKey;
        private final ParsedUnimodSequence parsedPeptide;
        private final long[] tokenArray;
        private final byte seedCharge;
        @SuppressWarnings("unused")
        private final double minimumChargeProbability;

        private SeedPredictionContext(
                AutomaticSeedKey seedKey,
                ParsedUnimodSequence parsedPeptide,
                long[] tokenArray,
                byte seedCharge,
                double minimumChargeProbability) {
            this.seedKey = seedKey;
            this.parsedPeptide = parsedPeptide;
            this.tokenArray = tokenArray;
            this.seedCharge = seedCharge;
            this.minimumChargeProbability = minimumChargeProbability;
        }
    }

    private static final class ScoutPredictionRecord {
        private final float[] ms2;
        private final float[] irt;
        private final float[] ccs;
        private final float[] chargeDist;

        private ScoutPredictionRecord(float[] ms2, float[] irt, float[] ccs, float[] chargeDist) {
            this.ms2 = ms2;
            this.irt = irt;
            this.ccs = ccs;
            this.chargeDist = chargeDist;
        }
    }

    private static final class AutomaticSeedKey {
        private final String unimodPeptide;
        private final double precursorNce;

        private AutomaticSeedKey(String unimodPeptide, double precursorNce) {
            this.unimodPeptide = unimodPeptide;
            this.precursorNce = precursorNce;
        }

        @Override
        public boolean equals(Object object) {
            if (this == object) {
                return true;
            }
            if (!(object instanceof AutomaticSeedKey other)) {
                return false;
            }
            return unimodPeptide.equals(other.unimodPeptide)
                    && Double.doubleToLongBits(precursorNce) == Double.doubleToLongBits(other.precursorNce);
        }

        @Override
        public int hashCode() {
            return Objects.hash(unimodPeptide, Double.doubleToLongBits(precursorNce));
        }

        @Override
        public String toString() {
            return unimodPeptide + "@" + precursorNce;
        }
    }

    private static final class PredictionJobKey {
        private final String unimodPeptide;
        private final byte charge;
        private final double nce;

        private PredictionJobKey(String unimodPeptide, byte charge, double nce) {
            this.unimodPeptide = unimodPeptide;
            this.charge = charge;
            this.nce = nce;
        }

        @Override
        public boolean equals(Object object) {
            if (this == object) {
                return true;
            }
            if (!(object instanceof PredictionJobKey other)) {
                return false;
            }
            return charge == other.charge
                    && unimodPeptide.equals(other.unimodPeptide)
                    && Double.doubleToLongBits(nce) == Double.doubleToLongBits(other.nce);
        }

        @Override
        public int hashCode() {
            return Objects.hash(unimodPeptide, charge, Double.doubleToLongBits(nce));
        }

        @Override
        public String toString() {
            return unimodPeptide + "/" + charge + "/" + nce;
        }
    }
}

package org.searlelab.jchronologer;

import java.io.PrintStream;
import java.nio.file.Path;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.Statement;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Objects;
import java.util.Optional;
import java.util.TreeMap;
import java.util.function.BiConsumer;
import java.util.function.Function;
import org.searlelab.jchronologer.api.ChronologerLibraryEntry;
import org.searlelab.jchronologer.api.ChronologerLibraryOptions;
import org.searlelab.jchronologer.api.ChronologerLibraryPredictor;
import org.searlelab.jchronologer.api.LibraryPredictionRequest;
import org.searlelab.jchronologer.api.PrecursorCondition;
import org.searlelab.jchronologer.dlib.DlibCodec;
import org.searlelab.jchronologer.impl.ChronologerFactory;
import org.searlelab.jchronologer.preprocessing.PeptideSequenceConverter;

/**
 * CLI entry point for calibrating Cartographer NCE against an existing DLIB.
 */
public final class CalibrateNceMain {

    static final double DEFAULT_START_NCE = 33.0;
    static final double DEFAULT_PPM_TOLERANCE = 10.0;
    static final int DEFAULT_MAX_ENTRIES = 1000;
    static final int DEFAULT_GROUP_TARGET = 100;
    static final double DEFAULT_MZ_BIN_WIDTH = 100.0;
    static final double MIN_NCE = 10.0;
    static final double MAX_NCE = 60.0;
    static final double SEED_STEP = 3.0;
    static final double SEED_SPAN = 6.0;
    static final double NCE_TOLERANCE = 0.1;
    static final double MASS_MATCH_EPSILON = ChronologerLibraryOptions.DEFAULT_MASS_MATCH_EPSILON;
    private static final double SCORE_EPSILON = 1e-9;

    private CalibrateNceMain() {
    }

    public static void main(String[] args) {
        int code = run(args, System.out, System.err);
        if (code != 0) {
            System.exit(code);
        }
    }

    static int run(String[] args, PrintStream out, PrintStream err) {
        CliArgs cliArgs;
        try {
            cliArgs = parseArgs(args);
        } catch (IllegalArgumentException e) {
            err.println(e.getMessage());
            printUsage(err);
            return 2;
        }

        if (cliArgs.help()) {
            printUsage(out);
            return 0;
        }

        try {
            CalibrationOutcome outcome = calibrate(
                    cliArgs,
                    (phase, evaluation) -> printEvaluationStatus(out, phase, evaluation, evaluation.rows().size()));
            printFinalReport(out, outcome.searchResult(), outcome.selectedRows());
            return 0;
        } catch (Exception e) {
            err.println("DLIB NCE calibration failed: " + e.getMessage());
            return 1;
        }
    }

    static CalibrationOutcome calibrate(CliArgs cliArgs) throws Exception {
        return calibrate(cliArgs, (phase, evaluation) -> {
        });
    }

    static CalibrationOutcome calibrate(CliArgs cliArgs, BiConsumer<String, Evaluation> evaluationConsumer) throws Exception {
        List<CalibrationRow> rows = readCalibrationRows(cliArgs.input(), cliArgs.mzBinWidth());
        if (rows.isEmpty()) {
            throw new IllegalArgumentException("No valid calibration rows were found in DLIB: " + cliArgs.input());
        }

        List<CalibrationRow> selectedRows = selectRows(rows, cliArgs.maxEntries(), cliArgs.groupTarget());
        if (selectedRows.isEmpty()) {
            throw new IllegalArgumentException("No calibration rows were selected from DLIB: " + cliArgs.input());
        }

        ChronologerLibraryOptions options = ChronologerLibraryOptions.builder()
                .batchSize(cliArgs.batchSize())
                .cartographerBatchSize(cliArgs.batchSize())
                .inferenceThreads(cliArgs.inferenceThreads())
                .verboseLogging(cliArgs.verbose())
                .build();

        try (ChronologerLibraryPredictor predictor = ChronologerFactory.createLibraryPredictor(options)) {
            predictor.init();
            Calibrator calibrator = new Calibrator(selectedRows, cliArgs.ppmTolerance(), predictor);
            SearchResult searchResult = searchBestNce(
                    cliArgs.startNce(),
                    calibrator::evaluate,
                    evaluationConsumer);
            return new CalibrationOutcome(rows, selectedRows, searchResult);
        }
    }

    static List<CalibrationRow> readCalibrationRows(Path input, double mzBinWidth) throws Exception {
        Objects.requireNonNull(input, "input");
        if (mzBinWidth <= 0.0) {
            throw new IllegalArgumentException("m/z bin width must be positive.");
        }

        try {
            Class.forName("org.sqlite.JDBC");
        } catch (ClassNotFoundException e) {
            throw new IllegalStateException("SQLite JDBC driver is not available.", e);
        }

        Map<String, CalibrationRow> byPeptideAndCharge = new LinkedHashMap<>();
        try (Connection connection = DriverManager.getConnection("jdbc:sqlite:" + input.toAbsolutePath());
                Statement statement = connection.createStatement();
                ResultSet resultSet = statement.executeQuery(
                        "select PeptideModSeq, PrecursorCharge, PrecursorMz, MassArray, IntensityArray "
                                + "from entries")) {
            while (resultSet.next()) {
                String peptideModSeq = resultSet.getString(1);
                int charge = resultSet.getInt(2);
                double precursorMz = resultSet.getDouble(3);
                byte[] massBlob = resultSet.getBytes(4);
                byte[] intensityBlob = resultSet.getBytes(5);

                Optional<CalibrationRow> candidate = decodeCalibrationRow(
                        peptideModSeq,
                        charge,
                        precursorMz,
                        massBlob,
                        intensityBlob,
                        mzBinWidth);
                if (candidate.isEmpty()) {
                    continue;
                }

                CalibrationRow row = candidate.orElseThrow();
                String key = calibrationKey(row.unimodPeptideSequence(), row.precursorCharge());
                CalibrationRow existing = byPeptideAndCharge.get(key);
                if (existing == null || row.summedObservedIntensity() > existing.summedObservedIntensity()) {
                    byPeptideAndCharge.put(key, row);
                }
            }
        }

        return byPeptideAndCharge.values().stream()
                .sorted(Comparator.comparingDouble(CalibrationRow::summedObservedIntensity).reversed()
                        .thenComparing(CalibrationRow::unimodPeptideSequence)
                        .thenComparingInt(row -> row.precursorCharge()))
                .toList();
    }

    private static Optional<CalibrationRow> decodeCalibrationRow(
            String peptideModSeq,
            int charge,
            double precursorMz,
            byte[] massBlob,
            byte[] intensityBlob,
            double mzBinWidth) {
        if (peptideModSeq == null || peptideModSeq.isBlank() || charge <= 0 || !Double.isFinite(precursorMz)) {
            return Optional.empty();
        }
        if (massBlob == null || intensityBlob == null || massBlob.length == 0 || intensityBlob.length == 0) {
            return Optional.empty();
        }

        try {
            String unimodPeptide = PeptideSequenceConverter.normalizeToUnimod(peptideModSeq, MASS_MATCH_EPSILON);
            double[] observedMasses = decodeMassArray(DlibCodec.decompress(massBlob));
            float[] observedIntensities = decodeIntensityArray(DlibCodec.decompress(intensityBlob));
            if (observedMasses.length == 0 || observedMasses.length != observedIntensities.length) {
                return Optional.empty();
            }

            double summedIntensity = 0.0;
            for (float intensity : observedIntensities) {
                if (Float.isFinite(intensity) && intensity > 0.0f) {
                    summedIntensity += intensity;
                }
            }
            if (!(summedIntensity > 0.0)) {
                return Optional.empty();
            }

            return Optional.of(new CalibrationRow(
                    unimodPeptide,
                    (byte) charge,
                    precursorMz,
                    observedMasses,
                    observedIntensities,
                    summedIntensity,
                    mzBinStart(precursorMz, mzBinWidth),
                    formatMzBinLabel(mzBinStart(precursorMz, mzBinWidth), mzBinWidth)));
        } catch (RuntimeException e) {
            return Optional.empty();
        }
    }

    static List<CalibrationRow> selectRows(List<CalibrationRow> rows, int maxEntries, int groupTarget) {
        if (rows == null || rows.isEmpty()) {
            return List.of();
        }
        if (maxEntries <= 0) {
            throw new IllegalArgumentException("Max entries must be positive.");
        }
        if (groupTarget <= 0) {
            throw new IllegalArgumentException("Group target must be positive.");
        }

        List<CalibrationRow> rankedRows = rows.stream()
                .sorted(Comparator.comparingDouble(CalibrationRow::summedObservedIntensity).reversed()
                        .thenComparing(CalibrationRow::unimodPeptideSequence)
                        .thenComparingInt(row -> row.precursorCharge()))
                .toList();

        Map<Byte, List<CalibrationRow>> byCharge = groupBy(rankedRows, CalibrationRow::precursorCharge);
        LinkedHashSet<CalibrationRow> chargeCoverage = new LinkedHashSet<>();
        for (Map.Entry<Byte, List<CalibrationRow>> entry : byCharge.entrySet()) {
            List<CalibrationRow> candidates = entry.getValue();
            for (int i = 0; i < Math.min(groupTarget, candidates.size()); i++) {
                chargeCoverage.add(candidates.get(i));
            }
        }

        LinkedHashSet<CalibrationRow> selected = new LinkedHashSet<>();
        selected.addAll(chargeCoverage);

        Map<Double, List<CalibrationRow>> byMzBin = groupBy(rankedRows, CalibrationRow::mzBinStart);
        List<CalibrationRow> mzCoverageCandidates = new ArrayList<>();
        for (Map.Entry<Double, List<CalibrationRow>> entry : byMzBin.entrySet()) {
            int added = 0;
            for (CalibrationRow row : entry.getValue()) {
                if (selected.contains(row)) {
                    continue;
                }
                mzCoverageCandidates.add(row);
                added++;
                if (added >= groupTarget) {
                    break;
                }
            }
        }
        mzCoverageCandidates.sort(Comparator.comparingDouble(CalibrationRow::summedObservedIntensity).reversed()
                .thenComparing(CalibrationRow::unimodPeptideSequence)
                .thenComparingInt(row -> row.precursorCharge()));
        for (CalibrationRow row : mzCoverageCandidates) {
            selected.add(row);
        }

        if (selected.size() < maxEntries) {
            for (CalibrationRow row : rankedRows) {
                if (selected.size() >= maxEntries) {
                    break;
                }
                selected.add(row);
            }
        }

        return List.copyOf(selected);
    }

    static double scorePredictedIonCosine(
            double[] predictedMasses,
            float[] predictedIntensities,
            double[] observedMasses,
            float[] observedIntensities,
            double ppmTolerance) {
        if (predictedMasses.length != predictedIntensities.length) {
            throw new IllegalArgumentException("Predicted masses and intensities must have the same length.");
        }
        if (observedMasses.length != observedIntensities.length) {
            throw new IllegalArgumentException("Observed masses and intensities must have the same length.");
        }

        double dot = 0.0;
        double predictedNorm = 0.0;
        double observedNorm = 0.0;
        for (int i = 0; i < predictedMasses.length; i++) {
            float predictedIntensity = predictedIntensities[i];
            double predictedMz = predictedMasses[i];
            if (!Float.isFinite(predictedIntensity) || predictedIntensity <= 0.0f || !Double.isFinite(predictedMz)) {
                continue;
            }
            double observedIntensity = nearestObservedIntensity(predictedMz, observedMasses, observedIntensities, ppmTolerance);
            dot += predictedIntensity * observedIntensity;
            predictedNorm += predictedIntensity * predictedIntensity;
            observedNorm += observedIntensity * observedIntensity;
        }

        if (predictedNorm <= 0.0 || observedNorm <= 0.0) {
            return 0.0;
        }
        return dot / (Math.sqrt(predictedNorm) * Math.sqrt(observedNorm));
    }

    static double[] decodeMassArray(byte[] bytes) {
        return chooseMorePlausibleMasses(
                decodeDoubles(bytes, ByteOrder.LITTLE_ENDIAN),
                decodeDoubles(bytes, ByteOrder.BIG_ENDIAN));
    }

    static float[] decodeIntensityArray(byte[] bytes) {
        return chooseMorePlausibleIntensities(
                decodeFloats(bytes, ByteOrder.LITTLE_ENDIAN),
                decodeFloats(bytes, ByteOrder.BIG_ENDIAN));
    }

    private static double[] decodeDoubles(byte[] bytes, ByteOrder byteOrder) {
        ByteBuffer buffer = ByteBuffer.wrap(bytes).order(byteOrder);
        double[] values = new double[bytes.length / Double.BYTES];
        for (int i = 0; i < values.length; i++) {
            values[i] = buffer.getDouble();
        }
        return values;
    }

    private static float[] decodeFloats(byte[] bytes, ByteOrder byteOrder) {
        ByteBuffer buffer = ByteBuffer.wrap(bytes).order(byteOrder);
        float[] values = new float[bytes.length / Float.BYTES];
        for (int i = 0; i < values.length; i++) {
            values[i] = buffer.getFloat();
        }
        return values;
    }

    private static double[] chooseMorePlausibleMasses(double[] first, double[] second) {
        return massPlausibilityScore(first) >= massPlausibilityScore(second) ? first : second;
    }

    private static float[] chooseMorePlausibleIntensities(float[] first, float[] second) {
        return intensityPlausibilityScore(first) >= intensityPlausibilityScore(second) ? first : second;
    }

    private static int massPlausibilityScore(double[] values) {
        int score = 0;
        double previous = Double.NEGATIVE_INFINITY;
        for (double value : values) {
            if (Double.isFinite(value)) {
                score += 2;
            }
            if (value > 0.0 && value < 4000.0) {
                score += 4;
            }
            if (value >= previous) {
                score += 1;
            }
            previous = value;
        }
        return score;
    }

    private static int intensityPlausibilityScore(float[] values) {
        int score = 0;
        for (float value : values) {
            if (Float.isFinite(value)) {
                score += 2;
            }
            if (value >= 0.0f) {
                score += 3;
            }
            if (value <= 10_000_000.0f) {
                score += 1;
            }
        }
        return score;
    }

    private static double nearestObservedIntensity(
            double predictedMz,
            double[] observedMasses,
            float[] observedIntensities,
            double ppmTolerance) {
        double bestAbsDelta = Double.POSITIVE_INFINITY;
        double bestIntensity = 0.0;
        for (int i = 0; i < observedMasses.length; i++) {
            double observedMz = observedMasses[i];
            float observedIntensity = observedIntensities[i];
            if (!Double.isFinite(observedMz) || !Float.isFinite(observedIntensity) || observedIntensity <= 0.0f) {
                continue;
            }

            double absDelta = Math.abs(observedMz - predictedMz);
            double ppm = (absDelta / predictedMz) * 1_000_000.0;
            if (ppm <= ppmTolerance && absDelta < bestAbsDelta) {
                bestAbsDelta = absDelta;
                bestIntensity = observedIntensity;
            }
        }
        return bestIntensity;
    }

    static SearchResult searchBestNce(
            double startNce,
            Function<Double, Evaluation> evaluator) {
        return searchBestNce(startNce, evaluator, (phase, evaluation) -> {
        });
    }

    static SearchResult searchBestNce(
            double startNce,
            Function<Double, Evaluation> evaluator,
            BiConsumer<String, Evaluation> evaluationConsumer) {
        TreeMap<Double, Evaluation> evaluations = new TreeMap<>();
        Function<Double, Evaluation> cachedEvaluator = nce -> evaluations.computeIfAbsent(nce, evaluator::apply);

        List<Double> seedPoints = new ArrayList<>();
        seedPoints.add(clampNce(startNce - SEED_SPAN));
        seedPoints.add(clampNce(startNce - SEED_STEP));
        seedPoints.add(clampNce(startNce));
        seedPoints.add(clampNce(startNce + SEED_STEP));
        seedPoints.add(clampNce(startNce + SEED_SPAN));

        for (double point : new LinkedHashSet<>(seedPoints)) {
            cachedEvaluator.apply(point);
            evaluationConsumer.accept("seed", contextualizeEvaluation(evaluations, point));
        }

        while (true) {
            Evaluation best = currentBestEvaluation(evaluations);
            List<Double> ordered = new ArrayList<>(evaluations.keySet());
            int bestIndex = ordered.indexOf(best.nce());
            if (bestIndex == 0 && best.nce() > MIN_NCE) {
                double next = clampNce(best.nce() - SEED_STEP);
                if (evaluations.containsKey(next)) {
                    break;
                }
                cachedEvaluator.apply(next);
                evaluationConsumer.accept("expand", contextualizeEvaluation(evaluations, next));
                continue;
            }
            if (bestIndex == ordered.size() - 1 && best.nce() < MAX_NCE) {
                double next = clampNce(best.nce() + SEED_STEP);
                if (evaluations.containsKey(next)) {
                    break;
                }
                cachedEvaluator.apply(next);
                evaluationConsumer.accept("expand", contextualizeEvaluation(evaluations, next));
                continue;
            }
            break;
        }

        Evaluation incumbent = currentBestEvaluation(evaluations);

        while (true) {
            List<Double> ordered = new ArrayList<>(evaluations.keySet());
            int bestIndex = ordered.indexOf(incumbent.nce());
            if (bestIndex < 0) {
                throw new IllegalStateException("Best NCE is missing from cached evaluations.");
            }

            if (bestIndex == 0) {
                double right = ordered.get(1);
                if (right - incumbent.nce() <= NCE_TOLERANCE) {
                    break;
                }
                double mid = roundNce((incumbent.nce() + right) / 2.0);
                if (evaluations.containsKey(mid)) {
                    break;
                }
                cachedEvaluator.apply(mid);
                Evaluation evaluation = contextualizeEvaluation(evaluations, mid);
                evaluationConsumer.accept("refine", evaluation);
                if (compareEvaluationsAscending(evaluation, incumbent) > 0) {
                    incumbent = evaluation;
                }
                continue;
            }

            if (bestIndex == ordered.size() - 1) {
                double left = ordered.get(bestIndex - 1);
                if (incumbent.nce() - left <= NCE_TOLERANCE) {
                    break;
                }
                double mid = roundNce((incumbent.nce() + left) / 2.0);
                if (evaluations.containsKey(mid)) {
                    break;
                }
                cachedEvaluator.apply(mid);
                Evaluation evaluation = contextualizeEvaluation(evaluations, mid);
                evaluationConsumer.accept("refine", evaluation);
                if (compareEvaluationsAscending(evaluation, incumbent) > 0) {
                    incumbent = evaluation;
                }
                continue;
            }

            double left = ordered.get(bestIndex - 1);
            double right = ordered.get(bestIndex + 1);
            if (incumbent.nce() - left <= NCE_TOLERANCE && right - incumbent.nce() <= NCE_TOLERANCE) {
                break;
            }

            double leftMid = roundNce((left + incumbent.nce()) / 2.0);
            double rightMid = roundNce((incumbent.nce() + right) / 2.0);
            if (!evaluations.containsKey(leftMid)) {
                cachedEvaluator.apply(leftMid);
                Evaluation leftMidEvaluation = contextualizeEvaluation(evaluations, leftMid);
                evaluationConsumer.accept("refine", leftMidEvaluation);
            }
            if (!evaluations.containsKey(rightMid)) {
                cachedEvaluator.apply(rightMid);
                Evaluation rightMidEvaluation = contextualizeEvaluation(evaluations, rightMid);
                evaluationConsumer.accept("refine", rightMidEvaluation);
            }

            Evaluation leftMidEvaluation = contextualizeEvaluation(evaluations, leftMid);
            Evaluation rightMidEvaluation = contextualizeEvaluation(evaluations, rightMid);

            Evaluation nextIncumbent = incumbent;
            if (compareEvaluationsAscending(leftMidEvaluation, nextIncumbent) > 0) {
                nextIncumbent = leftMidEvaluation;
            }
            if (compareEvaluationsAscending(rightMidEvaluation, nextIncumbent) > 0) {
                nextIncumbent = rightMidEvaluation;
            }
            incumbent = nextIncumbent;
        }

        List<Evaluation> orderedEvaluations = new ArrayList<>(evaluations.values());
        orderedEvaluations.sort(Comparator.comparingDouble(Evaluation::nce));
        DistributionReport distributionReport = buildDistributionReport(orderedEvaluations, incumbent);
        return new SearchResult(incumbent, orderedEvaluations, distributionReport);
    }

    private static Evaluation currentBestEvaluation(TreeMap<Double, Evaluation> evaluations) {
        return evaluations.values().stream()
                .map(evaluation -> contextualizeEvaluation(evaluations, evaluation.nce()))
                .max(CalibrateNceMain::compareEvaluationsAscending)
                .orElseThrow();
    }

    private static Evaluation contextualizeEvaluation(TreeMap<Double, Evaluation> evaluations, double nce) {
        Evaluation base = evaluations.get(nce);
        if (base == null) {
            throw new IllegalStateException("Missing cached evaluation for NCE " + nce);
        }
        Evaluation incumbent = evaluations.values().stream()
                .max(CalibrateNceMain::compareEvaluationsAscending)
                .orElse(base);
        Double left = evaluations.lowerKey(nce);
        Double right = evaluations.higherKey(nce);
        return new Evaluation(
                base.nce(),
                base.meanCosine(),
                base.perRowScores(),
                base.rows(),
                base.evaluationIndex(),
                incumbent.nce(),
                Optional.ofNullable(left),
                Optional.ofNullable(right));
    }

    private static int compareEvaluationsDescending(Evaluation left, Evaluation right) {
        return -compareEvaluationsAscending(left, right);
    }

    static int compareEvaluationsAscending(Evaluation left, Evaluation right) {
        int scoreComparison = Double.compare(left.meanCosine(), right.meanCosine());
        if (Math.abs(left.meanCosine() - right.meanCosine()) <= SCORE_EPSILON) {
            scoreComparison = 0;
        }
        if (scoreComparison != 0) {
            return scoreComparison;
        }
        int nceComparison = Double.compare(right.nce(), left.nce());
        if (nceComparison != 0) {
            return nceComparison;
        }
        return Integer.compare(right.evaluationIndex(), left.evaluationIndex());
    }

    static DistributionReport buildDistributionReport(List<Evaluation> evaluations, Evaluation bestEvaluation) {
        if (evaluations.isEmpty()) {
            throw new IllegalArgumentException("Evaluations must be non-empty.");
        }

        int rowCount = bestEvaluation.perRowScores().size();
        List<Double> bestNces = new ArrayList<>(rowCount);
        Map<Byte, List<Double>> byCharge = new TreeMap<>();
        Map<String, List<Double>> byMzBin = new TreeMap<>();

        for (int i = 0; i < rowCount; i++) {
            Evaluation bestForRow = null;
            for (Evaluation evaluation : evaluations) {
                if (bestForRow == null) {
                    bestForRow = evaluation;
                    continue;
                }
                double score = evaluation.perRowScores().get(i);
                double bestScore = bestForRow.perRowScores().get(i);
                if (score > bestScore + SCORE_EPSILON
                        || (Math.abs(score - bestScore) <= SCORE_EPSILON
                                && (evaluation.nce() < bestForRow.nce()
                                        || (Double.compare(evaluation.nce(), bestForRow.nce()) == 0
                                                && evaluation.evaluationIndex() < bestForRow.evaluationIndex())))) {
                    bestForRow = evaluation;
                }
            }

            double bestNce = bestForRow.nce();
            bestNces.add(bestNce);

            CalibrationRow row = bestEvaluation.rows().get(i);
            byCharge.computeIfAbsent(row.precursorCharge(), ignored -> new ArrayList<>()).add(bestNce);
            byMzBin.computeIfAbsent(row.mzBinLabel(), ignored -> new ArrayList<>()).add(bestNce);
        }

        return new DistributionReport(
                summarize(bestNces),
                summarizeGroups(byCharge),
                summarizeGroups(byMzBin),
                buildHistogram(bestNces));
    }

    static SummaryStats summarize(Collection<Double> values) {
        if (values.isEmpty()) {
            throw new IllegalArgumentException("Values must be non-empty.");
        }
        List<Double> sorted = values.stream().sorted().toList();
        return new SummaryStats(
                sorted.size(),
                percentile(sorted, 0.05),
                percentile(sorted, 0.25),
                percentile(sorted, 0.50),
                percentile(sorted, 0.75),
                percentile(sorted, 0.95));
    }

    private static <T> Map<T, SummaryStats> summarizeGroups(Map<T, List<Double>> groupedValues) {
        Map<T, SummaryStats> summary = new LinkedHashMap<>();
        for (Map.Entry<T, List<Double>> entry : groupedValues.entrySet()) {
            summary.put(entry.getKey(), summarize(entry.getValue()));
        }
        return summary;
    }

    static List<HistogramBin> buildHistogram(Collection<Double> values) {
        if (values.isEmpty()) {
            return List.of();
        }
        double min = values.stream().mapToDouble(Double::doubleValue).min().orElseThrow();
        double max = values.stream().mapToDouble(Double::doubleValue).max().orElseThrow();
        int start = (int) Math.floor(min);
        int end = (int) Math.floor(max);
        Map<Integer, Integer> counts = new LinkedHashMap<>();
        for (int bin = start; bin <= end; bin++) {
            counts.put(bin, 0);
        }
        for (double value : values) {
            int bin = (int) Math.floor(value);
            counts.computeIfPresent(bin, (ignored, current) -> current + 1);
        }
        List<HistogramBin> bins = new ArrayList<>();
        int maxCount = counts.values().stream().mapToInt(Integer::intValue).max().orElse(0);
        for (Map.Entry<Integer, Integer> entry : counts.entrySet()) {
            int barLength = maxCount == 0 ? 0 : (int) Math.round((entry.getValue() / (double) maxCount) * 40.0);
            bins.add(new HistogramBin(entry.getKey(), entry.getKey() + 1, entry.getValue(), "#".repeat(barLength)));
        }
        return bins;
    }

    static double percentile(List<Double> sortedValues, double percentile) {
        if (sortedValues.isEmpty()) {
            throw new IllegalArgumentException("Sorted values must be non-empty.");
        }
        if (percentile <= 0.0) {
            return sortedValues.get(0);
        }
        if (percentile >= 1.0) {
            return sortedValues.get(sortedValues.size() - 1);
        }
        double position = (sortedValues.size() - 1) * percentile;
        int lower = (int) Math.floor(position);
        int upper = (int) Math.ceil(position);
        if (lower == upper) {
            return sortedValues.get(lower);
        }
        double fraction = position - lower;
        return sortedValues.get(lower) + ((sortedValues.get(upper) - sortedValues.get(lower)) * fraction);
    }

    private static void printEvaluationStatus(PrintStream out, String phase, Evaluation evaluation, int selectedRowCount) {
        String bracket;
        if (evaluation.leftBracket().isPresent() && evaluation.rightBracket().isPresent()) {
            bracket = String.format(
                    Locale.US,
                    "[%.4f, %.4f]",
                    evaluation.leftBracket().orElseThrow(),
                    evaluation.rightBracket().orElseThrow());
        } else if (evaluation.leftBracket().isPresent()) {
            bracket = String.format(Locale.US, "[%.4f, upper-bound]", evaluation.leftBracket().orElseThrow());
        } else if (evaluation.rightBracket().isPresent()) {
            bracket = String.format(Locale.US, "[lower-bound, %.4f]", evaluation.rightBracket().orElseThrow());
        } else {
            bracket = "unbracketed";
        }
        out.printf(
                Locale.US,
                "phase=%s nce=%.4f entries=%d mean_cosine=%.6f incumbent_best=%.4f bracket=%s%n",
                phase,
                evaluation.nce(),
                selectedRowCount,
                evaluation.meanCosine(),
                evaluation.bestKnownNce(),
                bracket);
    }

    private static void printFinalReport(PrintStream out, SearchResult searchResult, List<CalibrationRow> rows) {
        out.println();
        out.println("Global calibration summary");
        out.printf(
                Locale.US,
                "best_nce=%.4f mean_cosine=%.6f entries=%d evaluations=%d%n",
                searchResult.bestEvaluation().nce(),
                searchResult.bestEvaluation().meanCosine(),
                rows.size(),
                searchResult.evaluations().size());

        out.println();
        out.println("Overall best-sampled NCE distribution");
        printSummaryStats(out, searchResult.distributionReport().overall());

        out.println();
        out.println("Best-sampled NCE distribution by charge");
        for (Map.Entry<Byte, SummaryStats> entry : searchResult.distributionReport().byCharge().entrySet()) {
            out.printf(Locale.US, "charge=%d ", entry.getKey());
            printSummaryStats(out, entry.getValue());
        }

        out.println();
        out.println("Best-sampled NCE distribution by precursor m/z bin");
        for (Map.Entry<String, SummaryStats> entry : searchResult.distributionReport().byMzBin().entrySet()) {
            out.printf(Locale.US, "mz_bin=%s ", entry.getKey());
            printSummaryStats(out, entry.getValue());
        }

        out.println();
        out.println("Histogram of best-sampled NCE values");
        for (HistogramBin bin : searchResult.distributionReport().histogram()) {
            out.printf(
                    Locale.US,
                    "%.1f-%.1f count=%d %s%n",
                    (double) bin.startInclusive(),
                    (double) bin.endExclusive(),
                    bin.count(),
                    bin.bar());
        }
    }

    private static void printSummaryStats(PrintStream out, SummaryStats summaryStats) {
        out.printf(
                Locale.US,
                "count=%d p05=%.4f p25=%.4f p50=%.4f p75=%.4f p95=%.4f%n",
                summaryStats.count(),
                summaryStats.p05(),
                summaryStats.p25(),
                summaryStats.p50(),
                summaryStats.p75(),
                summaryStats.p95());
    }

    static CliArgs parseArgs(String[] args) {
        Path input = null;
        double startNce = DEFAULT_START_NCE;
        double ppmTolerance = DEFAULT_PPM_TOLERANCE;
        int maxEntries = DEFAULT_MAX_ENTRIES;
        int groupTarget = DEFAULT_GROUP_TARGET;
        double mzBinWidth = DEFAULT_MZ_BIN_WIDTH;
        int batchSize = ChronologerLibraryOptions.builder().build().getBatchSize();
        int inferenceThreads = ChronologerLibraryOptions.builder().build().getInferenceThreads();
        boolean verbose = false;
        boolean help = false;

        for (int i = 0; i < args.length; i++) {
            String arg = args[i];
            switch (arg) {
                case "--help":
                case "-h":
                    help = true;
                    break;
                case "--start_nce":
                case "--start-nce":
                    startNce = parseDoubleOption(requireValue(args, ++i, arg), arg);
                    break;
                case "--ppm_tolerance":
                case "--ppm-tolerance":
                    ppmTolerance = parseDoubleOption(requireValue(args, ++i, arg), arg);
                    break;
                case "--max_entries":
                case "--max-entries":
                    maxEntries = parseIntOption(requireValue(args, ++i, arg), arg);
                    break;
                case "--group_target":
                case "--group-target":
                    groupTarget = parseIntOption(requireValue(args, ++i, arg), arg);
                    break;
                case "--mz_bin_width":
                case "--mz-bin-width":
                    mzBinWidth = parseDoubleOption(requireValue(args, ++i, arg), arg);
                    break;
                case "--batch_size":
                case "--batch-size":
                    batchSize = parseIntOption(requireValue(args, ++i, arg), arg);
                    break;
                case "--inference_threads":
                case "--inference-threads":
                    inferenceThreads = parseIntOption(requireValue(args, ++i, arg), arg);
                    break;
                case "--verbose":
                    verbose = true;
                    break;
                default:
                    if (arg.startsWith("--")) {
                        throw new IllegalArgumentException("Unknown option: " + arg);
                    }
                    if (input != null) {
                        throw new IllegalArgumentException("Expected a single input DLIB path.");
                    }
                    input = Path.of(arg);
                    break;
            }
        }

        if (help) {
            return new CliArgs(null, Double.NaN, Double.NaN, 0, 0, Double.NaN, 0, 0, false, true);
        }
        if (input == null) {
            throw new IllegalArgumentException("Expected input DLIB path.");
        }
        if (!Double.isFinite(startNce) || startNce < MIN_NCE || startNce > MAX_NCE) {
            throw new IllegalArgumentException("Start NCE must be a finite number in the range 10.0-60.0.");
        }
        if (!Double.isFinite(ppmTolerance) || ppmTolerance <= 0.0) {
            throw new IllegalArgumentException("ppm tolerance must be a finite positive number.");
        }
        if (maxEntries <= 0) {
            throw new IllegalArgumentException("Max entries must be positive.");
        }
        if (groupTarget <= 0) {
            throw new IllegalArgumentException("Group target must be positive.");
        }
        if (!Double.isFinite(mzBinWidth) || mzBinWidth <= 0.0) {
            throw new IllegalArgumentException("m/z bin width must be a finite positive number.");
        }
        if (batchSize <= 0) {
            throw new IllegalArgumentException("Batch size must be positive.");
        }
        if (inferenceThreads <= 0) {
            throw new IllegalArgumentException("Inference threads must be positive.");
        }

        return new CliArgs(
                input,
                startNce,
                ppmTolerance,
                maxEntries,
                groupTarget,
                mzBinWidth,
                batchSize,
                inferenceThreads,
                verbose,
                false);
    }

    private static String requireValue(String[] args, int index, String flag) {
        if (index >= args.length) {
            throw new IllegalArgumentException("Missing value for option: " + flag);
        }
        String value = args[index];
        if (value.startsWith("--")) {
            throw new IllegalArgumentException("Missing value for option: " + flag);
        }
        return value;
    }

    private static double parseDoubleOption(String value, String flag) {
        try {
            return Double.parseDouble(value);
        } catch (NumberFormatException e) {
            throw new IllegalArgumentException("Invalid numeric value for option: " + flag + " (" + value + ")");
        }
    }

    private static int parseIntOption(String value, String flag) {
        try {
            return Integer.parseInt(value);
        } catch (NumberFormatException e) {
            throw new IllegalArgumentException("Invalid integer value for option: " + flag + " (" + value + ")");
        }
    }

    private static void printUsage(PrintStream out) {
        out.println("Usage: calibrate-nce [options] <input.dlib>");
        out.println();
        out.println("Calibrate the best Cartographer NCE against an existing DLIB.");
        out.println();
        out.println("Options:");
        out.println("  --start_nce <value>          Starting NCE probe center (default: 33.0)");
        out.println("  --ppm_tolerance <value>      Fragment ppm matching tolerance (default: 10.0)");
        out.println("  --max_entries <value>        Target selected DLIB rows; may exceed for coverage (default: 1000)");
        out.println("  --group_target <value>       Per-charge and per-m/z-bin coverage target (default: 100)");
        out.println("  --mz_bin_width <value>       Precursor m/z bin width (default: 100.0)");
        out.println("  --batch_size <value>         Prediction batch size");
        out.println("  --inference_threads <value>  Number of inference threads");
        out.println("  --verbose                    Enable verbose predictor logging");
        out.println("  --help, -h                   Show this help message");
    }

    private static <T> Map<T, List<CalibrationRow>> groupBy(List<CalibrationRow> rows, Function<CalibrationRow, T> keyExtractor) {
        Map<T, List<CalibrationRow>> grouped = new LinkedHashMap<>();
        for (CalibrationRow row : rows) {
            grouped.computeIfAbsent(keyExtractor.apply(row), ignored -> new ArrayList<>()).add(row);
        }
        return grouped;
    }

    static double mzBinStart(double precursorMz, double mzBinWidth) {
        return Math.floor(precursorMz / mzBinWidth) * mzBinWidth;
    }

    static String formatMzBinLabel(double start, double width) {
        double end = start + width - 0.0001;
        return String.format(Locale.US, "%.1f-%.4f", start, end);
    }

    private static String calibrationKey(String unimodPeptideSequence, byte precursorCharge) {
        return unimodPeptideSequence + "|" + precursorCharge;
    }

    private static double clampNce(double nce) {
        return roundNce(Math.max(MIN_NCE, Math.min(MAX_NCE, nce)));
    }

    private static double roundNce(double nce) {
        return Math.round(nce * 10_000.0) / 10_000.0;
    }

    record CliArgs(
            Path input,
            double startNce,
            double ppmTolerance,
            int maxEntries,
            int groupTarget,
            double mzBinWidth,
            int batchSize,
            int inferenceThreads,
            boolean verbose,
            boolean help) {
    }

    record CalibrationRow(
            String unimodPeptideSequence,
            byte precursorCharge,
            double precursorMz,
            double[] observedMasses,
            float[] observedIntensities,
            double summedObservedIntensity,
            double mzBinStart,
            String mzBinLabel) {
    }

    record Evaluation(
            double nce,
            double meanCosine,
            List<Double> perRowScores,
            List<CalibrationRow> rows,
            int evaluationIndex,
            double bestKnownNce,
            Optional<Double> leftBracket,
            Optional<Double> rightBracket) {
    }

    record SearchResult(Evaluation bestEvaluation, List<Evaluation> evaluations, DistributionReport distributionReport) {
    }

    record CalibrationOutcome(List<CalibrationRow> allRows, List<CalibrationRow> selectedRows, SearchResult searchResult) {
    }

    record SummaryStats(int count, double p05, double p25, double p50, double p75, double p95) {
    }

    record HistogramBin(int startInclusive, int endExclusive, int count, String bar) {
    }

    record DistributionReport(
            SummaryStats overall,
            Map<Byte, SummaryStats> byCharge,
            Map<String, SummaryStats> byMzBin,
            List<HistogramBin> histogram) {
    }

    static final class Calibrator {
        private final List<CalibrationRow> rows;
        private final double ppmTolerance;
        private final ChronologerLibraryPredictor predictor;
        private int evaluationCounter;
        private final Map<Double, Evaluation> cache = new HashMap<>();

        Calibrator(List<CalibrationRow> rows, double ppmTolerance, ChronologerLibraryPredictor predictor) {
            this.rows = List.copyOf(rows);
            this.ppmTolerance = ppmTolerance;
            this.predictor = predictor;
        }

        Evaluation evaluate(double nce) {
            return cache.computeIfAbsent(nce, this::evaluateInternal);
        }

        private Evaluation evaluateInternal(double nce) {
            List<LibraryPredictionRequest> requests = new ArrayList<>(rows.size());
            for (CalibrationRow row : rows) {
                requests.add(new LibraryPredictionRequest(
                        row.unimodPeptideSequence(),
                        List.of(new PrecursorCondition(row.precursorCharge(), nce))));
            }

            List<ChronologerLibraryEntry> predictedEntries = predictor.predict(requests);
            Map<String, ChronologerLibraryEntry> predictionsByKey = new HashMap<>();
            for (ChronologerLibraryEntry entry : predictedEntries) {
                String normalizedPeptide = PeptideSequenceConverter.normalizeToUnimod(
                        entry.getUnimodPeptideSequence(),
                        MASS_MATCH_EPSILON);
                predictionsByKey.put(calibrationKey(normalizedPeptide, entry.getPrecursorCharge()), entry);
            }

            List<Double> perRowScores = new ArrayList<>(rows.size());
            double total = 0.0;
            for (CalibrationRow row : rows) {
                ChronologerLibraryEntry prediction = predictionsByKey.get(calibrationKey(
                        row.unimodPeptideSequence(),
                        row.precursorCharge()));
                double score = prediction == null
                        ? 0.0
                        : scorePredictedIonCosine(
                                prediction.getMassArray(),
                                prediction.getIntensityArray(),
                                row.observedMasses(),
                                row.observedIntensities(),
                                ppmTolerance);
                perRowScores.add(score);
                total += score;
            }

            evaluationCounter++;
            return new Evaluation(
                    nce,
                    total / rows.size(),
                    List.copyOf(perRowScores),
                    rows,
                    evaluationCounter,
                    nce,
                    Optional.empty(),
                    Optional.empty());
        }
    }
}

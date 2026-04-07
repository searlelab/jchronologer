package org.searlelab.jchronologer;

import java.io.IOException;
import java.io.PrintStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Batch wrapper for running NCE calibration across many DLIB files.
 */
public final class BatchCalibrateNceMain {

    private static final Path DEFAULT_INPUT_DIR = Path.of("/Users/searle.brian/Documents/astral/hela_test");
    private static final String DEFAULT_GLOB = "*.dlib";
    private static final double DEFAULT_START_OFFSET = 6.0;
    private static final Pattern NCE_PATTERN = Pattern.compile("_NCE(\\d+(?:\\.\\d+)?)_");

    private BatchCalibrateNceMain() {
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
            List<Path> dlibs = findDlibs(cliArgs.inputDir(), cliArgs.glob());
            if (dlibs.isEmpty()) {
                throw new IllegalArgumentException("No matching DLIB files found in " + cliArgs.inputDir());
            }

            List<BatchResult> results = new ArrayList<>();
            for (Path dlib : dlibs) {
                BatchResult result = calibrateOne(cliArgs, dlib, out);
                results.add(result);
            }
            printSummary(out, results);
            return 0;
        } catch (Exception e) {
            err.println("Batch DLIB NCE calibration failed: " + e.getMessage());
            return 1;
        }
    }

    private static BatchResult calibrateOne(CliArgs cliArgs, Path dlib, PrintStream out) throws Exception {
        String fileName = dlib.getFileName().toString();
        Double acquiredNce = parseAcquiredNce(fileName);
        double startNce = acquiredNce == null
                ? cliArgs.fallbackStartNce()
                : clampNce(acquiredNce + cliArgs.startOffset());

        out.printf(
                Locale.US,
                "== %s acquired_nce=%s start_nce=%.4f ==%n",
                fileName,
                acquiredNce == null ? "NA" : String.format(Locale.US, "%.4f", acquiredNce),
                startNce);

        CalibrateNceMain.CliArgs singleArgs = new CalibrateNceMain.CliArgs(
                dlib,
                startNce,
                cliArgs.ppmTolerance(),
                cliArgs.maxEntries(),
                cliArgs.groupTarget(),
                cliArgs.mzBinWidth(),
                cliArgs.batchSize(),
                cliArgs.inferenceThreads(),
                cliArgs.verbose(),
                false);

        CalibrateNceMain.CalibrationOutcome outcome = CalibrateNceMain.calibrate(singleArgs);
        for (CalibrateNceMain.Evaluation evaluation : outcome.searchResult().evaluations()) {
            printEvaluationStatus(out, evaluation, outcome.selectedRows().size());
        }
        out.printf(
                Locale.US,
                "completed %s best_nce=%.4f mean_cosine=%.6f selected_rows=%d%n%n",
                fileName,
                outcome.searchResult().bestEvaluation().nce(),
                outcome.searchResult().bestEvaluation().meanCosine(),
                outcome.selectedRows().size());
        out.printf(
                Locale.US,
                "charge_median_best_nce %s%n",
                formatChargeMedians(outcome.searchResult().distributionReport().byCharge()));
        out.printf(
                Locale.US,
                "mz_bin_median_best_nce %s%n%n",
                formatMzBinMedians(outcome.searchResult().distributionReport().byMzBin()));

        return new BatchResult(
                dlib,
                acquiredNce,
                startNce,
                outcome.searchResult().bestEvaluation().nce(),
                outcome.searchResult().bestEvaluation().meanCosine(),
                outcome.selectedRows().size(),
                outcome.allRows().size());
    }

    private static void printEvaluationStatus(PrintStream out, CalibrateNceMain.Evaluation evaluation, int selectedRows) {
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
                "nce=%.4f mean_cosine=%.6f incumbent_best=%.4f rows=%d bracket=%s%n",
                evaluation.nce(),
                evaluation.meanCosine(),
                evaluation.bestKnownNce(),
                selectedRows,
                bracket);
    }

    private static void printSummary(PrintStream out, List<BatchResult> results) {
        out.println("Batch summary");
        List<String> headers = List.of(
                "File",
                "AcquiredNCE",
                "StartNCE",
                "EstimatedNCE",
                "Delta",
                "MeanCosine",
                "SelectedRows",
                "ValidRows");
        int[] widths = headers.stream().mapToInt(String::length).toArray();
        List<List<String>> rows = new ArrayList<>();
        for (BatchResult result : results) {
            String acquired = result.acquiredNce() == null ? "NA" : formatDouble(result.acquiredNce());
            String delta = result.acquiredNce() == null ? "NA" : formatDouble(result.estimatedNce() - result.acquiredNce());
            List<String> row = List.of(
                    result.dlib().getFileName().toString(),
                    acquired,
                    formatDouble(result.startNce()),
                    formatDouble(result.estimatedNce()),
                    delta,
                    formatDouble(result.meanCosine()),
                    Integer.toString(result.selectedRows()),
                    Integer.toString(result.validRows()));
            rows.add(row);
            for (int i = 0; i < row.size(); i++) {
                widths[i] = Math.max(widths[i], row.get(i).length());
            }
        }

        printLine(out, headers, widths);
        printSeparator(out, widths);
        for (List<String> row : rows) {
            printLine(out, row, widths);
        }
    }

    private static void printLine(PrintStream out, List<String> values, int[] widths) {
        StringBuilder builder = new StringBuilder();
        for (int i = 0; i < values.size(); i++) {
            if (i > 0) {
                builder.append(" | ");
            }
            builder.append(padRight(values.get(i), widths[i]));
        }
        out.println(builder);
    }

    private static void printSeparator(PrintStream out, int[] widths) {
        StringBuilder builder = new StringBuilder();
        for (int i = 0; i < widths.length; i++) {
            if (i > 0) {
                builder.append("-+-");
            }
            builder.append("-".repeat(widths[i]));
        }
        out.println(builder);
    }

    private static String padRight(String value, int width) {
        if (value.length() >= width) {
            return value;
        }
        return value + " ".repeat(width - value.length());
    }

    private static String formatDouble(double value) {
        return String.format(Locale.US, "%.4f", value);
    }

    private static String formatChargeMedians(Map<Byte, CalibrateNceMain.SummaryStats> byCharge) {
        return joinMedians(byCharge, charge -> "z" + charge);
    }

    private static String formatMzBinMedians(Map<String, CalibrateNceMain.SummaryStats> byMzBin) {
        List<Map.Entry<String, CalibrateNceMain.SummaryStats>> entries = new ArrayList<>(byMzBin.entrySet());
        entries.sort(Comparator.comparingDouble(entry -> parseMzBinStart(entry.getKey())));
        List<String> parts = new ArrayList<>(entries.size());
        for (Map.Entry<String, CalibrateNceMain.SummaryStats> entry : entries) {
            parts.add(entry.getKey() + "=" + formatDouble(entry.getValue().p50()));
        }
        return String.join(" ", parts);
    }

    private static <T> String joinMedians(Map<T, CalibrateNceMain.SummaryStats> groups, java.util.function.Function<T, String> labeler) {
        List<String> parts = new ArrayList<>(groups.size());
        for (Map.Entry<T, CalibrateNceMain.SummaryStats> entry : groups.entrySet()) {
            parts.add(labeler.apply(entry.getKey()) + "=" + formatDouble(entry.getValue().p50()));
        }
        return String.join(" ", parts);
    }

    private static double parseMzBinStart(String label) {
        int separator = label.indexOf('-');
        if (separator < 0) {
            return Double.POSITIVE_INFINITY;
        }
        try {
            return Double.parseDouble(label.substring(0, separator));
        } catch (NumberFormatException e) {
            return Double.POSITIVE_INFINITY;
        }
    }

    static List<Path> findDlibs(Path inputDir, String glob) throws IOException {
        try (var stream = Files.list(inputDir)) {
            PathMatcher matcher = new PathMatcher(glob);
            return stream
                    .filter(Files::isRegularFile)
                    .filter(path -> matcher.matches(path.getFileName().toString()))
                    .sorted(Comparator.comparing(path -> path.getFileName().toString()))
                    .toList();
        }
    }

    static Double parseAcquiredNce(String fileName) {
        Matcher matcher = NCE_PATTERN.matcher(fileName);
        if (!matcher.find()) {
            return null;
        }
        return Double.parseDouble(matcher.group(1));
    }

    private static double clampNce(double nce) {
        return Math.max(CalibrateNceMain.MIN_NCE, Math.min(CalibrateNceMain.MAX_NCE, nce));
    }

    static CliArgs parseArgs(String[] args) {
        Path inputDir = DEFAULT_INPUT_DIR;
        String glob = DEFAULT_GLOB;
        double startOffset = DEFAULT_START_OFFSET;
        double fallbackStartNce = CalibrateNceMain.DEFAULT_START_NCE;
        double ppmTolerance = CalibrateNceMain.DEFAULT_PPM_TOLERANCE;
        int maxEntries = CalibrateNceMain.DEFAULT_MAX_ENTRIES;
        int groupTarget = CalibrateNceMain.DEFAULT_GROUP_TARGET;
        double mzBinWidth = CalibrateNceMain.DEFAULT_MZ_BIN_WIDTH;
        int batchSize = org.searlelab.jchronologer.api.ChronologerOptions.DEFAULT_BATCH_SIZE;
        int inferenceThreads = org.searlelab.jchronologer.api.ChronologerOptions.DEFAULT_INFERENCE_THREADS;
        boolean verbose = false;
        boolean help = false;

        for (int i = 0; i < args.length; i++) {
            String arg = args[i];
            switch (arg) {
                case "--help":
                case "-h":
                    help = true;
                    break;
                case "--glob":
                    glob = requireValue(args, ++i, arg);
                    break;
                case "--start_offset":
                case "--start-offset":
                    startOffset = parseDoubleOption(requireValue(args, ++i, arg), arg);
                    break;
                case "--fallback_start_nce":
                case "--fallback-start-nce":
                    fallbackStartNce = parseDoubleOption(requireValue(args, ++i, arg), arg);
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
                    inputDir = Path.of(arg);
                    break;
            }
        }

        if (help) {
            return new CliArgs(null, null, 0.0, 0.0, 0.0, 0, 0, 0.0, 0, 0, false, true);
        }
        if (inputDir == null || !Files.isDirectory(inputDir)) {
            throw new IllegalArgumentException("Input directory must exist: " + inputDir);
        }
        if (glob == null || glob.isBlank()) {
            throw new IllegalArgumentException("Glob must be non-empty.");
        }
        if (!Double.isFinite(startOffset)) {
            throw new IllegalArgumentException("Start offset must be finite.");
        }
        if (!Double.isFinite(fallbackStartNce)
                || fallbackStartNce < CalibrateNceMain.MIN_NCE
                || fallbackStartNce > CalibrateNceMain.MAX_NCE) {
            throw new IllegalArgumentException("Fallback start NCE must be in the range 10.0-60.0.");
        }
        if (!Double.isFinite(ppmTolerance) || ppmTolerance <= 0.0) {
            throw new IllegalArgumentException("ppm tolerance must be positive.");
        }
        if (maxEntries <= 0 || groupTarget <= 0 || !Double.isFinite(mzBinWidth) || mzBinWidth <= 0.0
                || batchSize <= 0 || inferenceThreads <= 0) {
            throw new IllegalArgumentException("Batch calibration numeric options must be positive.");
        }

        return new CliArgs(
                inputDir,
                glob,
                startOffset,
                fallbackStartNce,
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
        out.println("Usage: batch-calibrate-nce [options] [input-dir]");
        out.println();
        out.println("Run CalibrateNceMain sequentially across all matching DLIBs in a directory.");
        out.println("If a file name contains _NCE<value>_, the default start NCE is acquired+6.0.");
        out.println();
        out.println("Options:");
        out.println("  --glob <pattern>                File-name glob within the directory (default: *.dlib)");
        out.println("  --start_offset <value>          Added to filename NCE for the per-file start NCE (default: 6.0)");
        out.println("  --fallback_start_nce <value>    Used when the file name has no NCE token (default: 33.0)");
        out.println("  --ppm_tolerance <value>         Fragment ppm matching tolerance");
        out.println("  --max_entries <value>           Target selected rows; may exceed for coverage");
        out.println("  --group_target <value>          Per-charge and per-m/z-bin coverage target");
        out.println("  --mz_bin_width <value>          Precursor m/z bin width");
        out.println("  --batch_size <value>            Prediction batch size");
        out.println("  --inference_threads <value>     Number of inference threads");
        out.println("  --verbose                       Enable verbose predictor logging");
        out.println("  --help, -h                      Show this help message");
        out.println();
        out.printf(Locale.US, "Default input dir: %s%n", DEFAULT_INPUT_DIR);
    }

    record CliArgs(
            Path inputDir,
            String glob,
            double startOffset,
            double fallbackStartNce,
            double ppmTolerance,
            int maxEntries,
            int groupTarget,
            double mzBinWidth,
            int batchSize,
            int inferenceThreads,
            boolean verbose,
            boolean help) {
    }

    record BatchResult(
            Path dlib,
            Double acquiredNce,
            double startNce,
            double estimatedNce,
            double meanCosine,
            int selectedRows,
            int validRows) {
    }

    static final class PathMatcher {
        private final Pattern pattern;

        PathMatcher(String glob) {
            this.pattern = Pattern.compile(globToRegex(glob));
        }

        boolean matches(String fileName) {
            return pattern.matcher(fileName).matches();
        }

        private static String globToRegex(String glob) {
            StringBuilder regex = new StringBuilder("^");
            for (int i = 0; i < glob.length(); i++) {
                char c = glob.charAt(i);
                switch (c) {
                    case '*':
                        regex.append(".*");
                        break;
                    case '?':
                        regex.append('.');
                        break;
                    case '.':
                        regex.append("\\.");
                        break;
                    default:
                        regex.append(Pattern.quote(Character.toString(c)));
                        break;
                }
            }
            regex.append('$');
            return regex.toString();
        }
    }
}

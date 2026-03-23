package org.searlelab.jchronologer;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Locale;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.searlelab.jchronologer.api.ChronologerLibraryEntry;
import org.searlelab.jchronologer.api.ChronologerLibraryOptions;
import org.searlelab.jchronologer.api.ChronologerLibraryPredictor;
import org.searlelab.jchronologer.api.LibraryPredictionRequest;
import org.searlelab.jchronologer.impl.ChronologerFactory;
import org.searlelab.jchronologer.preprocessing.PeptideSequenceConverter;

/**
 * CLI entrypoint for quick library-accuracy inspection with one row per peptide+charge.
 */
public final class LibraryReportMain {

    private static final double DEFAULT_NCE = 33.0;
    private static final double DEFAULT_MINIMUM_CHARGE_PROBABILITY = 0.05;
    private static final double MASS_MATCH_EPSILON = ChronologerLibraryOptions.DEFAULT_MASS_MATCH_EPSILON;
    private static final int TOP_ION_COUNT = 3;
    private static final int MIN_NCE = 10;
    private static final int MAX_NCE = 60;

    private static final String ERROR_TOKEN = "ERROR";
    private static final String NO_CHARGE_TOKEN = "NO_CHARGE";

    private static final Pattern RAW_ION_PATTERN = Pattern.compile("(\\d+)\\+([yb])(\\d+)");

    private static final List<String> DEFAULT_PEPTIDES = List.of(
            "[]-HC[UNIMOD:4]VDPAVIAAIISR-[]",
            "[]-TLLISSLSPALPAEHLEDR-[]",
            "[]-TPIGSFLGSLS-[]",
            "[]-ANAEKTSGSNVKIVKVKKE-[]",
            "[UNIMOD:737]-K[UNIMOD:737]PGLAITFAK[UNIMOD:737]-[]");

    private static final List<String> HEADERS = List.of(
            "Peptide sequence (formatted with mods)",
            "Retention time (%ACN)",
            "Precursor charge",
            "CCS (by charge)",
            "Top 3 ions (by charge)");

    private LibraryReportMain() {
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

        if (cliArgs.help) {
            printUsage(out);
            return 0;
        }

        try {
            List<ReportRow> rows = buildRows(cliArgs);
            writeTable(out, rows);
            return 0;
        } catch (Exception e) {
            err.println("Library report failed: " + e.getMessage());
            return 1;
        }
    }

    private static List<ReportRow> buildRows(CliArgs cliArgs) {
        List<String> requestedPeptides = cliArgs.peptides.isEmpty() ? DEFAULT_PEPTIDES : cliArgs.peptides;

        List<ReportRow> rows = new ArrayList<>();
        List<String> validPeptides = new ArrayList<>();
        for (String peptide : requestedPeptides) {
            String input = peptide == null ? "" : peptide.trim();
            if (input.isEmpty()) {
                rows.add(ReportRow.error("<blank>", "Peptide sequence must be non-empty."));
                continue;
            }
            try {
                String canonical = PeptideSequenceConverter.normalizeToUnimod(input, MASS_MATCH_EPSILON);
                if (!input.equals(canonical)) {
                    throw new IllegalArgumentException(
                            "Sequence must be canonical terminal-aware UNIMOD (expected " + canonical + ")");
                }
                // Re-parse explicitly so strict terminal-aware UNIMOD input is required.
                PeptideSequenceConverter.parseNormalizedUnimod(input);
                validPeptides.add(input);
            } catch (RuntimeException e) {
                rows.add(ReportRow.error(input, e.getMessage()));
            }
        }

        if (validPeptides.isEmpty()) {
            return rows;
        }

        try (ChronologerLibraryPredictor predictor = ChronologerFactory.createLibraryPredictorDefault()) {
            for (String peptide : validPeptides) {
                try {
                    LibraryPredictionRequest request = new LibraryPredictionRequest(
                            peptide,
                            cliArgs.nce,
                            cliArgs.minimumChargeProbability);
                    List<ChronologerLibraryEntry> entries = predictor.predict(List.of(request));
                    if (entries.isEmpty()) {
                        rows.add(ReportRow.noCharge(peptide));
                        continue;
                    }

                    entries.stream()
                            .sorted(Comparator.comparingInt(entry -> entry.getPrecursorCharge()))
                            .forEach(entry -> rows.add(ReportRow.fromEntry(entry)));
                } catch (RuntimeException e) {
                    rows.add(ReportRow.error(peptide, e.getMessage()));
                }
            }
        }

        return rows;
    }

    private static void writeTable(PrintStream out, List<ReportRow> rows) {
        int[] widths = new int[HEADERS.size()];
        for (int i = 0; i < HEADERS.size(); i++) {
            widths[i] = HEADERS.get(i).length();
        }
        for (ReportRow row : rows) {
            List<String> columns = row.columns();
            for (int i = 0; i < columns.size(); i++) {
                widths[i] = Math.max(widths[i], columns.get(i).length());
            }
        }

        out.println(formatLine(HEADERS, widths));
        out.println(formatSeparator(widths));
        for (ReportRow row : rows) {
            out.println(formatLine(row.columns(), widths));
        }
    }

    private static String formatLine(List<String> values, int[] widths) {
        StringBuilder builder = new StringBuilder();
        for (int i = 0; i < values.size(); i++) {
            if (i > 0) {
                builder.append(" | ");
            }
            builder.append(padRight(values.get(i), widths[i]));
        }
        return builder.toString();
    }

    private static String formatSeparator(int[] widths) {
        StringBuilder builder = new StringBuilder();
        for (int i = 0; i < widths.length; i++) {
            if (i > 0) {
                builder.append("-+-");
            }
            builder.append("-".repeat(widths[i]));
        }
        return builder.toString();
    }

    private static String padRight(String value, int width) {
        if (value.length() >= width) {
            return value;
        }
        return value + " ".repeat(width - value.length());
    }

    private static CliArgs parseArgs(String[] args) {
        double nce = DEFAULT_NCE;
        double minimumChargeProbability = DEFAULT_MINIMUM_CHARGE_PROBABILITY;
        boolean help = false;
        List<String> peptides = new ArrayList<>();

        for (int i = 0; i < args.length; i++) {
            String arg = args[i];
            switch (arg) {
                case "--help":
                case "-h":
                    help = true;
                    break;
                case "--nce":
                    nce = parseDoubleOption(requireValue(args, ++i, arg), arg);
                    break;
                case "--minimum_charge_probability":
                case "--minimum-charge-probability":
                    minimumChargeProbability = parseDoubleOption(requireValue(args, ++i, arg), arg);
                    break;
                default:
                    if (arg.startsWith("--")) {
                        throw new IllegalArgumentException("Unknown option: " + arg);
                    }
                    peptides.add(arg);
                    break;
            }
        }

        if (help) {
            return new CliArgs(Double.NaN, Double.NaN, List.of(), true);
        }

        if (!Double.isFinite(nce) || nce < MIN_NCE || nce > MAX_NCE) {
            throw new IllegalArgumentException(
                    "NCE must be a finite number in the range " + MIN_NCE + "-" + MAX_NCE + ".");
        }
        if (!Double.isFinite(minimumChargeProbability)
                || minimumChargeProbability < 0.0
                || minimumChargeProbability > 1.0) {
            throw new IllegalArgumentException(
                    "Minimum charge probability must be a finite number in the range 0.0-1.0.");
        }

        return new CliArgs(nce, minimumChargeProbability, List.copyOf(peptides), false);
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

    private static void printUsage(PrintStream stream) {
        stream.println("Usage: library-report [options] [peptide1 peptide2 ...]");
        stream.println();
        stream.println("Generate a quick console report from DefaultChronologerLibraryPredictor.");
        stream.println("Each peptide must be a canonical terminal-aware UNIMOD sequence.");
        stream.println("If no peptides are provided, a built-in default peptide set is used.");
        stream.println();
        stream.println("Options:");
        stream.println("  --nce <value>                              Automatic-mode NCE (default: 33.0)");
        stream.println("  --minimum_charge_probability <value>       Minimum charge probability (default: 0.05)");
        stream.println("  --minimum-charge-probability <value>       Alias for --minimum_charge_probability");
        stream.println("  --help, -h                                 Show this help");
    }

    private static String formatTopIons(ChronologerLibraryEntry entry, int topN) {
        float[] intensities = entry.getIntensityArray();
        String[] ionTypes = entry.getIonTypeArray();
        List<Integer> indices = new ArrayList<>(intensities.length);
        for (int i = 0; i < intensities.length; i++) {
            indices.add(i);
        }

        indices.sort((left, right) -> {
            int intensityCompare = Float.compare(intensities[right], intensities[left]);
            if (intensityCompare != 0) {
                return intensityCompare;
            }
            int ionTypeCompare = ionTypes[left].compareTo(ionTypes[right]);
            if (ionTypeCompare != 0) {
                return ionTypeCompare;
            }
            return Integer.compare(left, right);
        });

        int limit = Math.min(topN, indices.size());
        if (limit == 0) {
            return "N/A";
        }

        List<String> formatted = new ArrayList<>(limit);
        for (int i = 0; i < limit; i++) {
            formatted.add(formatIonTypeForDisplay(ionTypes[indices.get(i)]));
        }
        return String.join(", ", formatted);
    }

    static String formatIonTypeForDisplay(String rawIonType) {
        Matcher matcher = RAW_ION_PATTERN.matcher(rawIonType);
        if (!matcher.matches()) {
            return rawIonType;
        }

        int fragmentCharge = Integer.parseInt(matcher.group(1));
        StringBuilder builder = new StringBuilder();
        builder.append(matcher.group(2)).append(matcher.group(3));
        if (fragmentCharge > 1) {
            builder.append("+".repeat(fragmentCharge));
        }
        return builder.toString();
    }

    private static String summarizeErrorMessage(String message) {
        String raw = message == null || message.isBlank() ? "Unknown error" : message.trim();
        if (raw.length() <= 120) {
            return raw;
        }
        return raw.substring(0, 117) + "...";
    }

    private record CliArgs(double nce, double minimumChargeProbability, List<String> peptides, boolean help) {
    }

    private static final class ReportRow {
        private final String peptide;
        private final String retentionTimeAcn;
        private final String precursorCharge;
        private final String ccs;
        private final String topIons;

        private ReportRow(
                String peptide,
                String retentionTimeAcn,
                String precursorCharge,
                String ccs,
                String topIons) {
            this.peptide = peptide;
            this.retentionTimeAcn = retentionTimeAcn;
            this.precursorCharge = precursorCharge;
            this.ccs = ccs;
            this.topIons = topIons;
        }

        private static ReportRow fromEntry(ChronologerLibraryEntry entry) {
            String acn = String.format(Locale.US, "%.4f", entry.getRetentionTimeInSeconds() / 60.0f);
            String charge = Byte.toString(entry.getPrecursorCharge());
            String ccsValue = entry.getCCS().isPresent()
                    ? String.format(Locale.US, "%.4f", entry.getCCS().get())
                    : "N/A";
            String topIonsValue = formatTopIons(entry, TOP_ION_COUNT);
            return new ReportRow(entry.getUnimodPeptideSequence(), acn, charge, ccsValue, topIonsValue);
        }

        private static ReportRow noCharge(String peptide) {
            return new ReportRow(peptide, NO_CHARGE_TOKEN, NO_CHARGE_TOKEN, NO_CHARGE_TOKEN, NO_CHARGE_TOKEN);
        }

        private static ReportRow error(String peptide, String message) {
            String errorText = ERROR_TOKEN + ": " + summarizeErrorMessage(message);
            return new ReportRow(peptide, ERROR_TOKEN, ERROR_TOKEN, ERROR_TOKEN, errorText);
        }

        private List<String> columns() {
            return List.of(peptide, retentionTimeAcn, precursorCharge, ccs, topIons);
        }
    }
}

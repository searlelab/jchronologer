package org.searlelab.jchronologer.cli;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.searlelab.jchronologer.api.AcceptedPrediction;
import org.searlelab.jchronologer.api.Chronologer;
import org.searlelab.jchronologer.api.ChronologerOptions;
import org.searlelab.jchronologer.api.PredictionResult;
import org.searlelab.jchronologer.api.RejectedPrediction;
import org.searlelab.jchronologer.impl.ChronologerFactory;
import org.searlelab.jchronologer.util.TsvTable;

/**
 * Command-style CLI entry point with explicit {@code predict} subcommand semantics.
 *
 * <p>Compared with {@link org.searlelab.jchronologer.Main}, this variant requires explicit input
 * and output options and can optionally retain rejected rows with diagnostics.
 */
public final class Main {

    private Main() {
    }

    public static void main(String[] args) {
        int code = run(args);
        if (code != 0) {
            System.exit(code);
        }
    }

    static int run(String[] args) {
        if (args.length == 0 || isHelp(args[0])) {
            printUsage();
            return args.length == 0 ? 1 : 0;
        }

        if (!"predict".equals(args[0])) {
            System.err.println("Unknown command: " + args[0]);
            printUsage();
            return 2;
        }

        PredictArgs predictArgs;
        try {
            predictArgs = parsePredictArgs(args);
        } catch (IllegalArgumentException e) {
            System.err.println(e.getMessage());
            printUsage();
            return 2;
        }

        try {
            runPredict(predictArgs);
            return 0;
        } catch (Exception e) {
            System.err.println("Prediction failed: " + e.getMessage());
            return 1;
        }
    }

    private static void runPredict(PredictArgs predictArgs) throws IOException {
        TsvTable inputTable = TsvTable.read(predictArgs.input);
        int peptideColumnIndex = inputTable.columnIndex(predictArgs.peptideColumn);
        if (peptideColumnIndex < 0) {
            throw new IllegalArgumentException("Input TSV is missing peptide column: " + predictArgs.peptideColumn);
        }

        List<String> peptides = new ArrayList<>();
        for (String[] row : inputTable.getRows()) {
            peptides.add(row[peptideColumnIndex]);
        }

        ChronologerOptions options = ChronologerOptions.builder()
                .batchSize(predictArgs.batchSize)
                .build();

        PredictionResult predictionResult;
        try (Chronologer chronologer = ChronologerFactory.create(options)) {
            predictionResult = chronologer.predict(peptides);
        }

        Map<Integer, AcceptedPrediction> acceptedByRow = new HashMap<>();
        for (AcceptedPrediction accepted : predictionResult.getAccepted()) {
            acceptedByRow.put(accepted.getRowIndex(), accepted);
        }

        Map<Integer, RejectedPrediction> rejectedByRow = new HashMap<>();
        for (RejectedPrediction rejected : predictionResult.getRejected()) {
            rejectedByRow.put(rejected.getRowIndex(), rejected);
        }

        List<String> outputHeaders = new ArrayList<>(inputTable.getHeaders());
        outputHeaders.add("Pred_HI");
        if (predictArgs.keepRejected) {
            outputHeaders.add("ChronologerStatus");
            outputHeaders.add("RejectionReason");
            outputHeaders.add("ErrorDetail");
        }

        List<String[]> outputRows = new ArrayList<>();
        List<String[]> rows = inputTable.getRows();
        int inputWidth = inputTable.getHeaders().size();

        for (int i = 0; i < rows.size(); i++) {
            String[] inputRow = rows.get(i);
            AcceptedPrediction accepted = acceptedByRow.get(i);
            if (accepted == null && !predictArgs.keepRejected) {
                continue;
            }

            int extra = predictArgs.keepRejected ? 4 : 1;
            String[] outputRow = new String[inputWidth + extra];
            System.arraycopy(inputRow, 0, outputRow, 0, Math.min(inputRow.length, inputWidth));

            if (accepted != null) {
                outputRow[inputWidth] = Float.toString(accepted.getPredHi());
                if (predictArgs.keepRejected) {
                    outputRow[inputWidth + 1] = "ACCEPTED";
                    outputRow[inputWidth + 2] = "";
                    outputRow[inputWidth + 3] = "";
                }
            } else {
                outputRow[inputWidth] = "";
                if (predictArgs.keepRejected) {
                    RejectedPrediction rejected = rejectedByRow.get(i);
                    outputRow[inputWidth + 1] = "REJECTED";
                    outputRow[inputWidth + 2] =
                            rejected == null ? "" : rejected.getRejectionReason().name();
                    outputRow[inputWidth + 3] =
                            rejected == null || rejected.getErrorDetail() == null
                                    ? ""
                                    : rejected.getErrorDetail();
                }
            }
            outputRows.add(outputRow);
        }

        TsvTable.write(predictArgs.output, outputHeaders, outputRows);
    }

    private static PredictArgs parsePredictArgs(String[] args) {
        Path input = null;
        Path output = null;
        String peptideColumn = "PeptideModSeq";
        int batchSize = ChronologerOptions.DEFAULT_BATCH_SIZE;
        boolean keepRejected = false;

        for (int i = 1; i < args.length; i++) {
            String arg = args[i];
            switch (arg) {
                case "--input":
                    input = Path.of(requireValue(args, ++i, arg));
                    break;
                case "--output":
                    output = Path.of(requireValue(args, ++i, arg));
                    break;
                case "--peptide-column":
                    peptideColumn = requireValue(args, ++i, arg);
                    break;
                case "--batch-size":
                    batchSize = parseIntegerOption(requireValue(args, ++i, arg), arg);
                    break;
                case "--keep-rejected":
                    keepRejected = true;
                    break;
                case "--help":
                case "-h":
                    throw new IllegalArgumentException("Use `jchronologer --help` to print usage.");
                default:
                    throw new IllegalArgumentException("Unknown option: " + arg);
            }
        }

        if (input == null) {
            throw new IllegalArgumentException("Missing required option: --input");
        }
        if (output == null) {
            throw new IllegalArgumentException("Missing required option: --output");
        }
        if (batchSize <= 0) {
            throw new IllegalArgumentException("Batch size must be positive.");
        }

        return new PredictArgs(input, output, peptideColumn, batchSize, keepRejected);
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

    private static int parseIntegerOption(String value, String flag) {
        try {
            return Integer.parseInt(value);
        } catch (NumberFormatException e) {
            throw new IllegalArgumentException("Invalid integer value for option: " + flag + " (" + value + ")");
        }
    }

    private static boolean isHelp(String arg) {
        return "--help".equals(arg) || "-h".equals(arg);
    }

    private static void printUsage() {
        System.out.println("Usage:");
        System.out.println("  jchronologer predict --input <input.tsv> --output <output.tsv> [options]");
        System.out.println();
        System.out.println("Options:");
        System.out.println("  --peptide-column <name>   Input column containing peptide mod sequences (default: PeptideModSeq)");
        System.out.println("  --batch-size <n>          Inference batch size (default: 2048)");
        System.out.println("  --keep-rejected           Keep rejected rows and add rejection diagnostics");
    }

    /**
     * Parsed arguments for {@code predict} command execution.
     */
    private record PredictArgs(
            Path input,
            Path output,
            String peptideColumn,
            int batchSize,
            boolean keepRejected) {
    }
}

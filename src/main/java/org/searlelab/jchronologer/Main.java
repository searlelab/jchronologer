package org.searlelab.jchronologer;

import java.io.IOException;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.searlelab.jchronologer.api.AcceptedPrediction;
import org.searlelab.jchronologer.api.Chronologer;
import org.searlelab.jchronologer.api.ChronologerOptions;
import org.searlelab.jchronologer.api.PredictionResult;
import org.searlelab.jchronologer.impl.ChronologerFactory;
import org.searlelab.jchronologer.util.TsvTable;

public final class Main {

    private static final String DEFAULT_PEPTIDE_COLUMN = "PeptideModSeq";

    private Main() {
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
            runPrediction(cliArgs, out);
            return 0;
        } catch (Exception e) {
            err.println("Prediction failed: " + e.getMessage());
            return 1;
        }
    }

    private static void runPrediction(CliArgs cliArgs, PrintStream out) throws IOException {
        InputTable inputTable = readInputTable(cliArgs.input, cliArgs.peptideColumn);
        List<String> peptides = new ArrayList<>(inputTable.rows.size());
        for (String[] row : inputTable.rows) {
            peptides.add(row[inputTable.peptideColumnIndex]);
        }

        ChronologerOptions options = ChronologerOptions.builder()
                .batchSize(cliArgs.batchSize)
                .build();

        PredictionResult predictionResult;
        try (Chronologer chronologer = ChronologerFactory.create(options)) {
            predictionResult = chronologer.predict(peptides);
        }

        Map<Integer, AcceptedPrediction> acceptedByRow = new HashMap<>();
        for (AcceptedPrediction accepted : predictionResult.getAccepted()) {
            acceptedByRow.put(accepted.getRowIndex(), accepted);
        }

        List<String> headers = new ArrayList<>(inputTable.headers);
        headers.add("Pred_HI");

        List<String[]> outputRows = new ArrayList<>();
        int width = inputTable.headers.size();
        for (int rowIndex = 0; rowIndex < inputTable.rows.size(); rowIndex++) {
            AcceptedPrediction accepted = acceptedByRow.get(rowIndex);
            if (accepted == null) {
                continue;
            }
            String[] source = inputTable.rows.get(rowIndex);
            String[] row = new String[width + 1];
            System.arraycopy(source, 0, row, 0, Math.min(source.length, width));
            row[width] = Float.toString(accepted.getPredHi());
            outputRows.add(row);
        }

        if (cliArgs.output == null) {
            writeToStdout(out, headers, outputRows);
        } else {
            TsvTable.write(cliArgs.output, headers, outputRows);
        }
    }

    private static InputTable readInputTable(Path input, String peptideColumn) throws IOException {
        List<String> lines = Files.readAllLines(input, StandardCharsets.UTF_8);
        String firstNonBlank = firstNonBlank(lines);
        if (firstNonBlank == null) {
            throw new IllegalArgumentException("Input file is empty: " + input);
        }

        if (looksLikeTsv(firstNonBlank, peptideColumn)) {
            TsvTable table = TsvTable.read(input);
            int peptideColumnIndex = table.columnIndex(peptideColumn);
            if (peptideColumnIndex < 0) {
                throw new IllegalArgumentException("Input TSV is missing peptide column: " + peptideColumn);
            }
            return new InputTable(table.getHeaders(), table.getRows(), peptideColumnIndex);
        }

        List<String[]> rows = new ArrayList<>();
        for (String line : lines) {
            String peptide = line.trim();
            if (!peptide.isEmpty()) {
                rows.add(new String[] {peptide});
            }
        }
        return new InputTable(List.of(peptideColumn), rows, 0);
    }

    private static String firstNonBlank(List<String> lines) {
        for (String line : lines) {
            if (!line.isBlank()) {
                return line;
            }
        }
        return null;
    }

    private static boolean looksLikeTsv(String firstNonBlankLine, String peptideColumn) {
        return firstNonBlankLine.contains("\t") || firstNonBlankLine.trim().equals(peptideColumn);
    }

    private static void writeToStdout(PrintStream out, List<String> headers, List<String[]> rows) {
        out.println(String.join("\t", headers));
        for (String[] row : rows) {
            out.println(String.join("\t", row));
        }
    }

    private static CliArgs parseArgs(String[] args) {
        if (args.length == 0) {
            throw new IllegalArgumentException("Missing required input file.");
        }

        Path input = null;
        Path output = null;
        String peptideColumn = DEFAULT_PEPTIDE_COLUMN;
        int batchSize = ChronologerOptions.DEFAULT_BATCH_SIZE;
        boolean help = false;

        List<String> positional = new ArrayList<>();
        for (int i = 0; i < args.length; i++) {
            String arg = args[i];
            switch (arg) {
                case "--help":
                case "-h":
                    help = true;
                    break;
                case "--batch_size":
                case "--batch-size":
                    batchSize = parseIntegerOption(requireValue(args, ++i, arg), arg);
                    break;
                case "--peptide_column":
                case "--peptide-column":
                    peptideColumn = requireValue(args, ++i, arg);
                    break;
                default:
                    if (arg.startsWith("--")) {
                        throw new IllegalArgumentException("Unknown option: " + arg);
                    }
                    positional.add(arg);
                    break;
            }
        }

        if (help) {
            return new CliArgs(null, null, peptideColumn, batchSize, true);
        }

        if (positional.isEmpty()) {
            throw new IllegalArgumentException("Missing required input file.");
        }
        if (positional.size() > 2) {
            throw new IllegalArgumentException("Too many positional arguments.");
        }
        if (batchSize <= 0) {
            throw new IllegalArgumentException("Batch size must be positive.");
        }

        input = Path.of(positional.get(0));
        if (positional.size() == 2) {
            output = Path.of(positional.get(1));
        }
        return new CliArgs(input, output, peptideColumn, batchSize, false);
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

    private static void printUsage(PrintStream stream) {
        stream.println("Usage: jchronologer [options] <input.tsv|input.txt> [output.tsv]");
        stream.println();
        stream.println("Predict RTs using Chronologer.");
        stream.println("Input may be:");
        stream.println("  1) TSV with a peptide column (default: PeptideModSeq), or");
        stream.println("  2) Plain text with one peptide per line.");
        stream.println("Incompatible peptides are dropped.");
        stream.println();
        stream.println("Options:");
        stream.println("  --batch_size <n>          Inference batch size (default: 2048)");
        stream.println("  --peptide_column <name>   Peptide column name for TSV input (default: "+DEFAULT_PEPTIDE_COLUMN+")");
        stream.println("  --help, -h                Show this help");
    }

    private record CliArgs(
            Path input,
            Path output,
            String peptideColumn,
            int batchSize,
            boolean help) {
    }

    private record InputTable(
            List<String> headers,
            List<String[]> rows,
            int peptideColumnIndex) {
    }
}

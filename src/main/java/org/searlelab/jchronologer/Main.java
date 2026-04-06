package org.searlelab.jchronologer;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import org.searlelab.jchronologer.api.ChronologerLibraryEntry;
import org.searlelab.jchronologer.api.ChronologerLibraryOptions;
import org.searlelab.jchronologer.api.ChronologerLibraryPredictor;
import org.searlelab.jchronologer.api.ChronologerOptions;
import org.searlelab.jchronologer.api.LibraryPredictionRequest;
import org.searlelab.jchronologer.dlib.DlibDatabase;
import org.searlelab.jchronologer.dlib.DlibEntryRecord;
import org.searlelab.jchronologer.dlib.DlibMapper;
import org.searlelab.jchronologer.dlib.DlibMetadata;
import org.searlelab.jchronologer.dlib.DlibPeptideProteinRecord;
import org.searlelab.jchronologer.fasta.FastaProteinRecord;
import org.searlelab.jchronologer.fasta.PeptideAccessionMatchingTrie;
import org.searlelab.jchronologer.fasta.PeptidePrecursorMatchTarget;
import org.searlelab.jchronologer.fasta.SimpleFastaReader;
import org.searlelab.jchronologer.impl.ChronologerFactory;
import org.searlelab.jchronologer.preprocessing.PeptideSequenceConverter;
import org.searlelab.jchronologer.util.TsvTable;

/**
 * DLIB generator CLI entry point for jchronologer.
 */
public final class Main {

    private static final String DEFAULT_PEPTIDE_COLUMN = "PeptideModSeq";
    private static final double DEFAULT_NCE = 33.0;
    private static final double DEFAULT_MIN_CHARGE_PROBABILITY = 0.01;
    private static final double UNIMOD_MASS_MATCH_EPSILON = 1e-5;
    private static final int MIN_NCE = 10;
    private static final int MAX_NCE = 60;
    private static final String SUPPORTED_RESIDUES = "ACDEFGHIKLMNPQRSTVWY";
    private static final int MIN_LIBRARY_PEPTIDE_LENGTH = 7;
    private static final int MAX_LIBRARY_PEPTIDE_LENGTH = 31;

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
            runLibraryGeneration(cliArgs, err);
            return 0;
        } catch (Exception e) {
            err.println("DLIB generation failed: " + e.getMessage());
            return 1;
        }
    }

    private static void runLibraryGeneration(CliArgs cliArgs, PrintStream err) throws IOException, SQLException {
        List<String> peptides = readInputPeptides(cliArgs.input, cliArgs.peptideColumn);
        List<FastaProteinRecord> proteins = SimpleFastaReader.read(cliArgs.fasta);
        if (peptides.isEmpty()) {
            throw new IllegalArgumentException("Input file is empty: " + cliArgs.input);
        }

        ChronologerLibraryOptions options = ChronologerLibraryOptions.builder()
                .batchSize(cliArgs.batchSize)
                .cartographerBatchSize(cliArgs.batchSize)
                .verboseLogging(cliArgs.verbose)
                .build();

        Summary summary = new Summary(peptides.size());
        boolean success = false;
        try (ChronologerLibraryPredictor predictor = ChronologerFactory.createLibraryPredictor(options);
                DlibDatabase database = new DlibDatabase(cliArgs.output, DlibMetadata.defaults())) {
            for (int start = 0; start < peptides.size(); start += cliArgs.batchSize) {
                int end = Math.min(start + cliArgs.batchSize, peptides.size());
                processBatch(peptides.subList(start, end), proteins, predictor, database, cliArgs, summary);
            }
            if (summary.entriesWritten == 0) {
                throw new IllegalArgumentException("No library entries were written.");
            }
            database.createIndices();
            success = true;
        } finally {
            if (!success || summary.entriesWritten == 0) {
                Files.deleteIfExists(cliArgs.output);
            }
        }

        err.println("Read " + summary.peptidesRead + " peptides from input.");
        err.println("Processed " + summary.batchesProcessed + " peptide batches.");
        err.println("Generated " + summary.predictedEntries + " predicted precursor entries.");
        err.println("Wrote " + summary.entriesWritten + " DLIB entries to " + cliArgs.output + ".");
        err.println("Dropped " + summary.invalidPeptides + " peptides for invalid input.");
        err.println("Dropped " + summary.noChargePeptides + " peptides with no qualifying charge.");
        err.println("Dropped " + summary.unmatchedPeptides + " peptides and "
                + summary.unmatchedEntries + " entries with no FASTA accession match.");
    }

    private static void processBatch(
            List<String> rawPeptides,
            List<FastaProteinRecord> proteins,
            ChronologerLibraryPredictor predictor,
            DlibDatabase database,
            CliArgs cliArgs,
            Summary summary) throws SQLException {
        List<String> normalizedPeptides = new ArrayList<>(rawPeptides.size());
        Map<String, Integer> occurrencesByUnimod = new LinkedHashMap<>();
        for (String peptide : rawPeptides) {
            String input = peptide == null ? "" : peptide.trim();
            if (input.isEmpty()) {
                summary.invalidPeptides++;
                continue;
            }
            try {
                String normalized = PeptideSequenceConverter.normalizeToUnimod(input, UNIMOD_MASS_MATCH_EPSILON);
                if (!isSupportedForLibraryPrediction(normalized)) {
                    summary.invalidPeptides++;
                    continue;
                }
                normalizedPeptides.add(normalized);
                occurrencesByUnimod.merge(normalized, 1, Integer::sum);
            } catch (RuntimeException e) {
                summary.invalidPeptides++;
            }
        }

        summary.batchesProcessed++;
        if (normalizedPeptides.isEmpty()) {
            return;
        }

        List<LibraryPredictionRequest> requests = new ArrayList<>(normalizedPeptides.size());
        for (String normalizedPeptide : normalizedPeptides) {
            requests.add(new LibraryPredictionRequest(
                    normalizedPeptide,
                    cliArgs.nce,
                    cliArgs.minimumChargeProbability));
        }

        List<ChronologerLibraryEntry> predictions = predictor.predict(requests);
        summary.predictedEntries += predictions.size();

        LinkedHashSet<String> peptidesWithPredictions = new LinkedHashSet<>();
        for (ChronologerLibraryEntry entry : predictions) {
            peptidesWithPredictions.add(entry.getUnimodPeptideSequence());
        }
        for (Map.Entry<String, Integer> occurrence : occurrencesByUnimod.entrySet()) {
            if (!peptidesWithPredictions.contains(occurrence.getKey())) {
                summary.noChargePeptides += occurrence.getValue();
            }
        }

        LinkedHashSet<String> peptideSequences = new LinkedHashSet<>();
        Map<String, Integer> unmatchedOccurrencesByPeptideSeq = new LinkedHashMap<>();
        for (String normalizedPeptide : occurrencesByUnimod.keySet()) {
            String peptideSeq = PeptideSequenceConverter.parseNormalizedUnimod(normalizedPeptide).getResidues();
            peptideSequences.add(peptideSeq);
            unmatchedOccurrencesByPeptideSeq.merge(peptideSeq, occurrencesByUnimod.get(normalizedPeptide), Integer::sum);
        }

        Map<String, PeptidePrecursorMatchTarget> matches = new PeptideAccessionMatchingTrie(peptideSequences).match(proteins);
        List<DlibEntryRecord> entriesToWrite = new ArrayList<>();
        LinkedHashSet<DlibPeptideProteinRecord> peptideProteinRows = new LinkedHashSet<>();
        LinkedHashSet<String> matchedPeptideSeqs = new LinkedHashSet<>();
        for (ChronologerLibraryEntry prediction : predictions) {
            DlibEntryRecord dlibEntry = DlibMapper.toDlibEntry(prediction);
            PeptidePrecursorMatchTarget matchTarget = matches.get(dlibEntry.getPeptideSeq());
            if (matchTarget == null || matchTarget.getProteinAccessions().isEmpty()) {
                summary.unmatchedEntries++;
                continue;
            }
            matchedPeptideSeqs.add(dlibEntry.getPeptideSeq());
            entriesToWrite.add(dlibEntry);
            for (String accession : matchTarget.getProteinAccessions()) {
                peptideProteinRows.add(DlibMapper.toPeptideProtein(dlibEntry.getPeptideSeq(), accession));
            }
        }

        for (String peptideSeq : peptideSequences) {
            if (!matchedPeptideSeqs.contains(peptideSeq)) {
                summary.unmatchedPeptides += unmatchedOccurrencesByPeptideSeq.getOrDefault(peptideSeq, 0);
            }
        }

        if (!entriesToWrite.isEmpty()) {
            database.writeBatch(entriesToWrite, new ArrayList<>(peptideProteinRows));
            summary.entriesWritten += entriesToWrite.size();
        }
    }

    private static List<String> readInputPeptides(Path input, String peptideColumn) throws IOException {
        try (BufferedReader reader = Files.newBufferedReader(input, StandardCharsets.UTF_8)) {
            String line;
            String firstNonBlank = null;
            while ((line = reader.readLine()) != null) {
                if (!line.isBlank()) {
                    firstNonBlank = line;
                    break;
                }
            }
            if (firstNonBlank == null) {
                return List.of();
            }

            if (looksLikeTsv(firstNonBlank, peptideColumn)) {
                TsvTable table = TsvTable.read(input);
                int peptideColumnIndex = table.columnIndex(peptideColumn);
                if (peptideColumnIndex < 0) {
                    throw new IllegalArgumentException("Input TSV is missing peptide column: " + peptideColumn);
                }
                List<String> peptides = new ArrayList<>(table.getRows().size());
                for (String[] row : table.getRows()) {
                    peptides.add(row[peptideColumnIndex]);
                }
                return peptides;
            }

            List<String> peptides = new ArrayList<>();
            peptides.add(firstNonBlank.trim());
            while ((line = reader.readLine()) != null) {
                String peptide = line.trim();
                if (!peptide.isEmpty()) {
                    peptides.add(peptide);
                }
            }
            return peptides;
        }
    }

    private static boolean looksLikeTsv(String firstNonBlankLine, String peptideColumn) {
        return firstNonBlankLine.contains("\t") || firstNonBlankLine.trim().equals(peptideColumn);
    }

    private static boolean isSupportedForLibraryPrediction(String normalizedPeptide) {
        String residues = PeptideSequenceConverter.parseNormalizedUnimod(normalizedPeptide).getResidues();
        if (residues.length() < MIN_LIBRARY_PEPTIDE_LENGTH || residues.length() > MAX_LIBRARY_PEPTIDE_LENGTH) {
            return false;
        }
        for (int i = 0; i < residues.length(); i++) {
            if (SUPPORTED_RESIDUES.indexOf(residues.charAt(i)) < 0) {
                return false;
            }
        }
        return true;
    }

    private static CliArgs parseArgs(String[] args) {
        if (args.length == 0) {
            throw new IllegalArgumentException("Missing required input peptide file.");
        }

        String peptideColumn = DEFAULT_PEPTIDE_COLUMN;
        int batchSize = ChronologerOptions.DEFAULT_BATCH_SIZE;
        double nce = DEFAULT_NCE;
        double minimumChargeProbability = DEFAULT_MIN_CHARGE_PROBABILITY;
        boolean verbose = false;
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
                case "--nce":
                    nce = parseDoubleOption(requireValue(args, ++i, arg), arg);
                    break;
                case "--min_charge_probability":
                case "--min-charge-probability":
                    minimumChargeProbability = parseDoubleOption(requireValue(args, ++i, arg), arg);
                    break;
                case "--peptide_column":
                case "--peptide-column":
                    peptideColumn = requireValue(args, ++i, arg);
                    break;
                case "--verbose":
                    verbose = true;
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
            return new CliArgs(null, null, null, peptideColumn, batchSize, nce, minimumChargeProbability, verbose, true);
        }
        if (positional.size() != 3) {
            throw new IllegalArgumentException(
                    positional.size() < 3
                            ? "Expected input peptides, FASTA, and output DLIB."
                            : "Too many positional arguments.");
        }
        if (batchSize <= 0) {
            throw new IllegalArgumentException("Batch size must be positive.");
        }
        if (!Double.isFinite(nce) || nce < MIN_NCE || nce > MAX_NCE) {
            throw new IllegalArgumentException("NCE must be a finite number in the range " + MIN_NCE + "-" + MAX_NCE + ".");
        }
        if (!Double.isFinite(minimumChargeProbability)
                || minimumChargeProbability < 0.0
                || minimumChargeProbability > 1.0) {
            throw new IllegalArgumentException(
                    "Minimum charge probability must be a finite number in the range 0.0-1.0.");
        }

        Path input = Path.of(positional.get(0));
        Path fasta = Path.of(positional.get(1));
        Path output = Path.of(positional.get(2));
        if (!Files.isRegularFile(input)) {
            throw new IllegalArgumentException("Input peptide file does not exist: " + input);
        }
        if (!Files.isRegularFile(fasta)) {
            throw new IllegalArgumentException("Protein FASTA file does not exist: " + fasta);
        }
        return new CliArgs(input, fasta, output, peptideColumn, batchSize, nce, minimumChargeProbability, verbose, false);
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

    private static double parseDoubleOption(String value, String flag) {
        try {
            return Double.parseDouble(value);
        } catch (NumberFormatException e) {
            throw new IllegalArgumentException("Invalid numeric value for option: " + flag + " (" + value + ")");
        }
    }

    private static void printUsage(PrintStream stream) {
        stream.println("Usage: jchronologer [options] <input.tsv|input.txt> <proteins.fasta> <output.dlib>");
        stream.println();
        stream.println("Generate EncyclopeDIA DLIB libraries from peptide lists.");
        stream.println("Input peptides may be mass-encoded or terminal-aware UNIMOD sequences.");
        stream.println("Protein FASTA is required for peptide-to-protein mapping.");
        stream.println();
        stream.println("Options:");
        stream.println("  --nce <value>                           Precursor NCE (default: 33.0)");
        stream.println("  --min_charge_probability <value>        Minimum charge probability (default: 0.01)");
        stream.println("  --min-charge-probability <value>        Alias for --min_charge_probability");
        stream.println("  --batch_size <n>                        Shared prediction/write batch size (default: 2048)");
        stream.println("  --peptide_column <name>                 Peptide column for TSV input (default: " + DEFAULT_PEPTIDE_COLUMN + ")");
        stream.println("  --verbose                               Enable detailed predictor diagnostics");
        stream.println("  --help, -h                              Show this help");
    }

    private record CliArgs(
            Path input,
            Path fasta,
            Path output,
            String peptideColumn,
            int batchSize,
            double nce,
            double minimumChargeProbability,
            boolean verbose,
            boolean help) {
    }

    private static final class Summary {
        private final int peptidesRead;
        private int batchesProcessed;
        private int predictedEntries;
        private int entriesWritten;
        private int invalidPeptides;
        private int noChargePeptides;
        private int unmatchedPeptides;
        private int unmatchedEntries;

        private Summary(int peptidesRead) {
            this.peptidesRead = peptidesRead;
        }
    }
}

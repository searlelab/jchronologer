package org.searlelab.jchronologer;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.ByteArrayOutputStream;
import java.io.InputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.Statement;
import org.junit.jupiter.api.Test;
import org.searlelab.jchronologer.api.AcceptedPrediction;
import org.searlelab.jchronologer.api.Chronologer;
import org.searlelab.jchronologer.api.ChronologerLibraryEntry;
import org.searlelab.jchronologer.api.ChronologerLibraryPredictor;
import org.searlelab.jchronologer.api.LibraryPredictionRequest;
import org.searlelab.jchronologer.api.PredictionResult;
import org.searlelab.jchronologer.impl.ChronologerFactory;
import org.searlelab.jchronologer.preprocessing.PeptideSequenceConverter;
import org.searlelab.jchronologer.util.TsvTable;

class UnimodPeptidesOutSmokeTest {

    private static final String UNIMOD_PEPTIDES_OUT_RESOURCE = "data/demo/unimod_peptides_out.txt";
    private static final double UNIMOD_MASS_MATCH_EPSILON = 1e-5;
    private static final int SMOKE_ROW_LIMIT = 128;
    private static final String SUPPORTED_RESIDUES = "ACDEFGHIKLMNPQRSTVWY";
    private static final int MIN_LIBRARY_PEPTIDE_LENGTH = 7;
    private static final int MAX_LIBRARY_PEPTIDE_LENGTH = 31;

    @Test
    void mainHandlesUnimodInputConsistentlyWithMassCodedChronologerPath() throws Exception {
        Path input = copyResourceToTemp(UNIMOD_PEPTIDES_OUT_RESOURCE);
        TsvTable inputTable = TsvTable.read(input);
        TsvTable smokeTable = new TsvTable(
                inputTable.getHeaders(),
                new ArrayList<>(inputTable.getRows().subList(0, Math.min(SMOKE_ROW_LIMIT, inputTable.getRows().size()))));
        Files.deleteIfExists(input);
        input = Files.createTempFile("jchronologer-unimod-peptides-out-smoke", ".tsv");
        TsvTable.write(input, smokeTable.getHeaders(), smokeTable.getRows());
        inputTable = smokeTable;
        int peptideColumn = inputTable.columnIndex("PeptideModSeq");
        assertTrue(peptideColumn >= 0, "Missing PeptideModSeq column in unimod_peptides_out resource.");

        List<String> converted = new ArrayList<>(inputTable.getRows().size());
        for (String[] row : inputTable.getRows()) {
            converted.add(toChronologerInputSequence(row[peptideColumn]));
        }

        PredictionResult baseline;
        try (Chronologer chronologer = ChronologerFactory.createDefault()) {
            baseline = chronologer.predict(converted);
        }
        assertEquals(inputTable.getRows().size(), baseline.getAcceptedCount(), "Expected all unimod rows to be accepted.");
        assertEquals(0, baseline.getRejectedCount(), "Expected no rejected unimod rows.");

        LinkedHashSet<String> residues = new LinkedHashSet<>();
        for (String[] row : inputTable.getRows()) {
            residues.add(PeptideSequenceConverter.parseNormalizedUnimod(
                    PeptideSequenceConverter.normalizeToUnimod(row[peptideColumn], UNIMOD_MASS_MATCH_EPSILON)).getResidues());
        }
        Path fasta = Files.createTempFile("jchronologer-unimod-main-fasta", ".fasta");
        Files.writeString(fasta, buildFasta(residues), StandardCharsets.UTF_8);

        Path output = Files.createTempFile("jchronologer-unimod-main-out", ".dlib");
        ByteArrayOutputStream stdoutBytes = new ByteArrayOutputStream();
        ByteArrayOutputStream stderrBytes = new ByteArrayOutputStream();
        PrintStream stdout = new PrintStream(stdoutBytes, true, StandardCharsets.UTF_8);
        PrintStream stderr = new PrintStream(stderrBytes, true, StandardCharsets.UTF_8);
        int code = Main.run(new String[] {input.toString(), fasta.toString(), output.toString()}, stdout, stderr);
        assertEquals(0, code, "CLI failed: " + stderrBytes.toString(StandardCharsets.UTF_8));

        LinkedHashSet<String> expectedPeptides = new LinkedHashSet<>();
        try (ChronologerLibraryPredictor predictor = ChronologerFactory.createLibraryPredictorDefault()) {
            Map<Integer, AcceptedPrediction> acceptedByRow = baseline.getAcceptedByRowIndex();
            for (int rowIndex = 0; rowIndex < inputTable.getRows().size(); rowIndex++) {
                AcceptedPrediction accepted = acceptedByRow.get(rowIndex);
                assertNotNull(accepted, "Missing accepted prediction for row " + rowIndex);
                String peptide = inputTable.getRows().get(rowIndex)[peptideColumn];
                String normalized = PeptideSequenceConverter.normalizeToUnimod(peptide, UNIMOD_MASS_MATCH_EPSILON);
                if (!isSupportedForLibraryPrediction(normalized)) {
                    continue;
                }
                List<ChronologerLibraryEntry> entries = predictor.predict(List.of(
                        new LibraryPredictionRequest(normalized, 33.0, 0.01)));
                if (!entries.isEmpty()) {
                    expectedPeptides.add(entries.get(0).getPeptideModSeq());
                }
            }
        }

        LinkedHashSet<String> observedPeptides = new LinkedHashSet<>();
        int writtenRows = 0;
        try (Connection connection = DriverManager.getConnection("jdbc:sqlite:" + output.toAbsolutePath());
                Statement statement = connection.createStatement();
                ResultSet resultSet = statement.executeQuery("select distinct PeptideModSeq from entries order by PeptideModSeq")) {
            while (resultSet.next()) {
                observedPeptides.add(resultSet.getString(1));
                writtenRows++;
            }
        }

        assertTrue(writtenRows > 0, "Expected at least one DLIB entry for UNIMOD input.");
        assertEquals(expectedPeptides, observedPeptides);
    }

    private static String toChronologerInputSequence(String peptide) {
        try {
            String unimod = PeptideSequenceConverter.normalizeToUnimod(peptide, UNIMOD_MASS_MATCH_EPSILON);
            return PeptideSequenceConverter.unimodToMassEncoded(unimod);
        } catch (RuntimeException e) {
            return peptide;
        }
    }

    private static Path copyResourceToTemp(String resource) throws Exception {
        try (InputStream stream = Thread.currentThread().getContextClassLoader().getResourceAsStream(resource)) {
            assertNotNull(stream, "Missing test resource: " + resource);
            Path tempFile = Files.createTempFile("jchronologer-unimod-peptides-out", ".tsv");
            Files.copy(stream, tempFile, StandardCopyOption.REPLACE_EXISTING);
            return tempFile;
        }
    }

    private static String buildFasta(Iterable<String> peptideSequences) {
        StringBuilder builder = new StringBuilder();
        int index = 1;
        for (String peptideSequence : peptideSequences) {
            builder.append(">P").append(index++).append('\n');
            builder.append("M").append(peptideSequence).append("K").append('\n');
        }
        return builder.toString();
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
}

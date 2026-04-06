package org.searlelab.jchronologer.fasta;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

/**
 * Minimal FASTA parser for protein accession matching.
 */
public final class SimpleFastaReader {

    private SimpleFastaReader() {
    }

    public static List<FastaProteinRecord> read(Path path) throws IOException {
        List<FastaProteinRecord> proteins = new ArrayList<>();
        try (BufferedReader reader = Files.newBufferedReader(path, StandardCharsets.UTF_8)) {
            String header = null;
            StringBuilder sequence = new StringBuilder();
            String line;
            while ((line = reader.readLine()) != null) {
                String trimmed = line.trim();
                if (trimmed.isEmpty()) {
                    continue;
                }
                if (trimmed.startsWith(">")) {
                    if (header != null) {
                        proteins.add(toRecord(header, sequence.toString(), path));
                    }
                    header = trimmed.substring(1).trim();
                    sequence.setLength(0);
                    continue;
                }
                if (header == null) {
                    throw new IOException("Malformed FASTA: sequence encountered before header in " + path);
                }
                sequence.append(trimmed.toUpperCase());
            }
            if (header != null) {
                proteins.add(toRecord(header, sequence.toString(), path));
            }
        }
        if (proteins.isEmpty()) {
            throw new IOException("FASTA file is empty: " + path);
        }
        return proteins;
    }

    private static FastaProteinRecord toRecord(String header, String sequence, Path path) throws IOException {
        if (header.isBlank()) {
            throw new IOException("Malformed FASTA: blank header in " + path);
        }
        if (sequence.isBlank()) {
            throw new IOException("Malformed FASTA: missing sequence for header " + header + " in " + path);
        }
        String[] parts = header.split("\\s+", 2);
        String accession = parts[0].trim();
        if (accession.isEmpty()) {
            throw new IOException("Malformed FASTA: missing accession in header " + header + " in " + path);
        }
        return new FastaProteinRecord(accession, sequence);
    }
}

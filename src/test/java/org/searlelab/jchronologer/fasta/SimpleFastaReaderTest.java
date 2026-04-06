package org.searlelab.jchronologer.fasta;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import org.junit.jupiter.api.Test;

class SimpleFastaReaderTest {

    @Test
    void readsUppercasesSequencesAndUsesFirstHeaderTokenAsAccession() throws Exception {
        Path fasta = Files.createTempFile("simple-fasta-reader", ".fasta");
        Files.writeString(
                fasta,
                "\n>sp|P1|desc protein one\nacde\nfg\n\n>P2 second protein\nmnpq\n",
                StandardCharsets.UTF_8);

        List<FastaProteinRecord> proteins = SimpleFastaReader.read(fasta);

        assertEquals(2, proteins.size());
        assertEquals("sp|P1|desc", proteins.get(0).accession());
        assertEquals("ACDEFG", proteins.get(0).sequence());
        assertEquals("P2", proteins.get(1).accession());
        assertEquals("MNPQ", proteins.get(1).sequence());
    }

    @Test
    void rejectsSequenceBeforeHeader() throws Exception {
        Path fasta = Files.createTempFile("simple-fasta-reader-invalid", ".fasta");
        Files.writeString(fasta, "ACDE\n>P1\nFG\n", StandardCharsets.UTF_8);

        IOException error = assertThrows(IOException.class, () -> SimpleFastaReader.read(fasta));
        assertEquals(true, error.getMessage().contains("sequence encountered before header"));
    }

    @Test
    void rejectsBlankHeaderAndMissingSequence() throws Exception {
        Path blankHeader = Files.createTempFile("simple-fasta-reader-blank", ".fasta");
        Files.writeString(blankHeader, ">\nACDE\n", StandardCharsets.UTF_8);
        IOException blankHeaderError = assertThrows(IOException.class, () -> SimpleFastaReader.read(blankHeader));
        assertEquals(true, blankHeaderError.getMessage().contains("blank header"));

        Path missingSequence = Files.createTempFile("simple-fasta-reader-missing-seq", ".fasta");
        Files.writeString(missingSequence, ">P1\n>P2\nACDE\n", StandardCharsets.UTF_8);
        IOException missingSequenceError = assertThrows(IOException.class, () -> SimpleFastaReader.read(missingSequence));
        assertEquals(true, missingSequenceError.getMessage().contains("missing sequence"));
    }

    @Test
    void rejectsEmptyFile() throws Exception {
        Path fasta = Files.createTempFile("simple-fasta-reader-empty", ".fasta");
        IOException error = assertThrows(IOException.class, () -> SimpleFastaReader.read(fasta));
        assertEquals(true, error.getMessage().contains("FASTA file is empty"));
    }
}

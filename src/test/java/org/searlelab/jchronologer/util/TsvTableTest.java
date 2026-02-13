package org.searlelab.jchronologer.util;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import org.junit.jupiter.api.Test;

class TsvTableTest {

    @Test
    void readThrowsOnEmptyFile() throws IOException {
        Path file = Files.createTempFile("tsv-empty", ".tsv");
        IllegalArgumentException error = assertThrows(IllegalArgumentException.class, () -> TsvTable.read(file));
        assertTrue(error.getMessage().contains("TSV file is empty"));
    }

    @Test
    void readPadsRowsShorterThanHeaderAndKeepsLongRows() throws IOException {
        Path file = Files.createTempFile("tsv-pad", ".tsv");
        Files.writeString(
                file,
                "A\tB\tC\n1\t2\n3\t4\t5\t6\n",
                StandardCharsets.UTF_8);

        TsvTable table = TsvTable.read(file);
        assertEquals(List.of("A", "B", "C"), table.getHeaders());
        assertEquals(2, table.getRows().size());
        assertEquals(3, table.getRows().get(0).length);
        assertEquals("", table.getRows().get(0)[2]);
        assertEquals(4, table.getRows().get(1).length);
        assertEquals("6", table.getRows().get(1)[3]);
        assertEquals(1, table.columnIndex("B"));
        assertEquals(-1, table.columnIndex("missing"));
    }

    @Test
    void headersAndRowsAreUnmodifiableViews() {
        List<String[]> rows = new ArrayList<>();
        rows.add(new String[] {"x"});
        TsvTable table = new TsvTable(List.of("A"), rows);
        assertThrows(UnsupportedOperationException.class, () -> table.getHeaders().add("B"));
        assertThrows(UnsupportedOperationException.class, () -> table.getRows().add(new String[] {"y"}));
    }

    @Test
    void writeCreatesParentDirectoriesAndRoundTrips() throws IOException {
        Path root = Files.createTempDirectory("tsv-write");
        Path output = root.resolve("nested").resolve("out.tsv");
        List<String> headers = List.of("PeptideModSeq", "Pred_HI");
        List<String[]> rows = new ArrayList<>();
        rows.add(new String[] {"VATVSLPR", "8.3945856"});

        TsvTable.write(output, headers, rows);
        assertTrue(Files.exists(output));

        TsvTable readBack = TsvTable.read(output);
        assertEquals(headers, readBack.getHeaders());
        assertEquals(1, readBack.getRows().size());
        assertEquals("VATVSLPR", readBack.getRows().get(0)[0]);
        assertEquals("8.3945856", readBack.getRows().get(0)[1]);
    }
}

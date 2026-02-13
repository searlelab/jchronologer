package org.searlelab.jchronologer.cli;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.util.List;
import org.junit.jupiter.api.Test;
import org.searlelab.jchronologer.util.TsvTable;

class CliIntegrationTest {

    @Test
    void defaultModeDropsRejectedRows() throws IOException {
        Path input = copyResourceToTemp("data/golden/parity_cases.tsv", "-input.tsv");
        Path output = Files.createTempFile("jchronologer-cli-default", ".tsv");

        int code = Main.run(new String[] {
                "predict",
                "--input", input.toString(),
                "--output", output.toString(),
                "--peptide-column", "PeptideModSeq"
        });

        assertEquals(0, code);
        TsvTable table = TsvTable.read(output);
        assertTrue(table.getHeaders().contains("Pred_HI"));
        assertTrue(!table.getHeaders().contains("ChronologerStatus"));
        assertEquals(20, table.getRows().size());
    }

    @Test
    void keepRejectedModePreservesAllRowsWithDiagnostics() throws IOException {
        Path input = copyResourceToTemp("data/golden/parity_cases.tsv", "-input.tsv");
        Path output = Files.createTempFile("jchronologer-cli-keep", ".tsv");

        int code = Main.run(new String[] {
                "predict",
                "--input", input.toString(),
                "--output", output.toString(),
                "--peptide-column", "PeptideModSeq",
                "--keep-rejected"
        });

        assertEquals(0, code);
        TsvTable table = TsvTable.read(output);
        List<String> headers = table.getHeaders();
        int statusColumn = headers.indexOf("ChronologerStatus");
        int reasonColumn = headers.indexOf("RejectionReason");

        assertTrue(statusColumn >= 0);
        assertTrue(reasonColumn >= 0);
        assertEquals(24, table.getRows().size());

        int accepted = 0;
        int rejected = 0;
        for (String[] row : table.getRows()) {
            if ("ACCEPTED".equals(row[statusColumn])) {
                accepted++;
            } else if ("REJECTED".equals(row[statusColumn])) {
                rejected++;
                assertTrue(!row[reasonColumn].isBlank());
            }
        }

        assertEquals(20, accepted);
        assertEquals(4, rejected);
    }

    private static Path copyResourceToTemp(String resource, String suffix) throws IOException {
        try (InputStream stream = Thread.currentThread().getContextClassLoader().getResourceAsStream(resource)) {
            if (stream == null) {
                throw new IllegalStateException("Missing resource: " + resource);
            }
            Path tempFile = Files.createTempFile("jchronologer-cli", suffix);
            Files.copy(stream, tempFile, StandardCopyOption.REPLACE_EXISTING);
            return tempFile;
        }
    }
}

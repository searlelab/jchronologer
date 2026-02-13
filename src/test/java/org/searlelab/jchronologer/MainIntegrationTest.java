package org.searlelab.jchronologer;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import org.junit.jupiter.api.Test;
import org.searlelab.jchronologer.util.TsvTable;

class MainIntegrationTest {

    @Test
    void plainTextInputWritesPredictionsToStdoutWhenOutputIsOmitted() throws IOException {
        Path input = Files.createTempFile("jchronologer-main-plain", ".txt");
        Files.writeString(
                input,
                "VATVSLPR\n[42.010565]ACDEFGHIK\n",
                StandardCharsets.UTF_8);

        ByteArrayOutputStream stdoutBytes = new ByteArrayOutputStream();
        ByteArrayOutputStream stderrBytes = new ByteArrayOutputStream();
        PrintStream stdout = new PrintStream(stdoutBytes, true, StandardCharsets.UTF_8);
        PrintStream stderr = new PrintStream(stderrBytes, true, StandardCharsets.UTF_8);

        int code = Main.run(new String[] {input.toString()}, stdout, stderr);

        assertEquals(0, code);
        String output = stdoutBytes.toString(StandardCharsets.UTF_8);
        String[] lines = output.split("\\R");
        assertEquals("PeptideModSeq\tPred_HI", lines[0]);
        assertEquals(2, lines.length);
        assertTrue(lines[1].startsWith("VATVSLPR\t"));
        assertTrue(!output.contains("[42.010565]ACDEFGHIK"));
    }

    @Test
    void tsvInputWritesPredictionsToOutputFileAndDropsRejectedRows() throws IOException {
        Path input = copyResourceToTemp("data/golden/parity_cases.tsv", "-input.tsv");
        Path output = Files.createTempFile("jchronologer-main-tsv", ".tsv");

        int code = Main.run(new String[] {input.toString(), output.toString()}, System.out, System.err);

        assertEquals(0, code);
        TsvTable table = TsvTable.read(output);
        assertTrue(table.getHeaders().contains("Pred_HI"));
        assertEquals(20, table.getRows().size());
    }

    private static Path copyResourceToTemp(String resource, String suffix) throws IOException {
        try (InputStream stream = Thread.currentThread().getContextClassLoader().getResourceAsStream(resource)) {
            if (stream == null) {
                throw new IllegalStateException("Missing resource: " + resource);
            }
            Path tempFile = Files.createTempFile("jchronologer-main", suffix);
            Files.copy(stream, tempFile, StandardCopyOption.REPLACE_EXISTING);
            return tempFile;
        }
    }
}

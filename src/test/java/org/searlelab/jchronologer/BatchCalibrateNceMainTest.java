package org.searlelab.jchronologer;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import org.junit.jupiter.api.Test;

class BatchCalibrateNceMainTest {

    @Test
    void parsesAcquiredNceFromFileName() {
        assertEquals(
                33.0,
                BatchCalibrateNceMain.parseAcquiredNce("20260403_no10_tne_hela_NCE33_OrbiDDA_01.dia.dlib"),
                1e-9);
        assertEquals(
                27.5,
                BatchCalibrateNceMain.parseAcquiredNce("sample_NCE27.5_AstDDA_01.dlib"),
                1e-9);
        assertNull(BatchCalibrateNceMain.parseAcquiredNce("sample_without_nce.dlib"));
    }

    @Test
    void findsOnlyMatchingDlibs() throws Exception {
        Path dir = Files.createTempDirectory("batch-calibrate-nce");
        Files.writeString(dir.resolve("a.dlib"), "");
        Files.writeString(dir.resolve("b.txt"), "");
        Files.writeString(dir.resolve("c.dia.dlib"), "");

        List<Path> matches = BatchCalibrateNceMain.findDlibs(dir, "*.dlib");

        assertEquals(2, matches.size());
        assertTrue(matches.get(0).getFileName().toString().endsWith(".dlib"));
        assertTrue(matches.get(1).getFileName().toString().endsWith(".dlib"));
    }
}

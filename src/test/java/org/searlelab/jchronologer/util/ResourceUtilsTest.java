package org.searlelab.jchronologer.util;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import org.junit.jupiter.api.Test;

class ResourceUtilsTest {

    @Test
    void copyClasspathResourceToTempFileCopiesBytes() throws IOException {
        Path copied = ResourceUtils.copyClasspathResourceToTempFile("data/demo/demo_rt.txt", ".txt");
        assertTrue(Files.exists(copied));

        byte[] expected = Files.readAllBytes(Path.of("src/test/resources/data/demo/demo_rt.txt"));
        byte[] actual = Files.readAllBytes(copied);
        assertArrayEquals(expected, actual);
    }

    @Test
    void copyClasspathResourceToTempFileThrowsForMissingResource() {
        IllegalArgumentException error = assertThrows(
                IllegalArgumentException.class,
                () -> ResourceUtils.copyClasspathResourceToTempFile("does/not/exist.txt", ".txt"));
        assertTrue(error.getMessage().contains("Missing classpath resource"));
    }
}

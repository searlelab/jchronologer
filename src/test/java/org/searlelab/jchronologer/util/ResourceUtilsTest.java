package org.searlelab.jchronologer.util;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.IOException;
import java.io.InputStream;
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

    @Test
    void copyClasspathResourceToTempFileWrapsReadIoExceptions() {
        ClassLoader original = Thread.currentThread().getContextClassLoader();
        ClassLoader failing = new ClassLoader(original) {
            @Override
            public InputStream getResourceAsStream(String name) {
                if (!"broken/resource.txt".equals(name)) {
                    return super.getResourceAsStream(name);
                }
                return new InputStream() {
                    @Override
                    public int read() throws IOException {
                        throw new IOException("simulated read failure");
                    }
                };
            }
        };

        Thread.currentThread().setContextClassLoader(failing);
        try {
            IllegalStateException error = assertThrows(
                    IllegalStateException.class,
                    () -> ResourceUtils.copyClasspathResourceToTempFile("broken/resource.txt", ".txt"));
            assertTrue(error.getMessage().contains("Failed to copy classpath resource"));
            assertTrue(error.getCause() instanceof IOException);
        } finally {
            Thread.currentThread().setContextClassLoader(original);
        }
    }
}

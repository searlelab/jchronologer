package org.searlelab.jchronologer.util;

import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;

/**
 * Utility helpers for working with classpath resources during runtime loading.
 */
public final class ResourceUtils {

    private ResourceUtils() {
    }

    /**
     * Copies a classpath resource to a temporary file and returns the resulting path.
     *
     * <p>This is used for DJL model loading paths that require filesystem-backed artifacts.
     *
     * @param resource classpath resource to copy
     * @param suffix suffix used when creating the temporary file
     * @return temporary file path containing resource bytes
     */
    public static Path copyClasspathResourceToTempFile(String resource, String suffix) {
        try (InputStream stream = Thread.currentThread().getContextClassLoader().getResourceAsStream(resource)) {
            if (stream == null) {
                throw new IllegalArgumentException("Missing classpath resource: " + resource);
            }
            Path tempFile = Files.createTempFile("jchronologer-", suffix);
            Files.copy(stream, tempFile, StandardCopyOption.REPLACE_EXISTING);
            tempFile.toFile().deleteOnExit();
            return tempFile;
        } catch (IOException e) {
            throw new IllegalStateException("Failed to copy classpath resource to temp file: " + resource, e);
        }
    }
}

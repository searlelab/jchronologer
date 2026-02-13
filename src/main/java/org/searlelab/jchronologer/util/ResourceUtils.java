package org.searlelab.jchronologer.util;

import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;

public final class ResourceUtils {

    private ResourceUtils() {
    }

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

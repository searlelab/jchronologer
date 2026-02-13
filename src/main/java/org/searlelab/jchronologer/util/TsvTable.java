package org.searlelab.jchronologer.util;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public final class TsvTable {

    private final List<String> headers;
    private final List<String[]> rows;

    public TsvTable(List<String> headers, List<String[]> rows) {
        this.headers = Collections.unmodifiableList(new ArrayList<>(headers));
        this.rows = Collections.unmodifiableList(new ArrayList<>(rows));
    }

    public List<String> getHeaders() {
        return headers;
    }

    public List<String[]> getRows() {
        return rows;
    }

    public int columnIndex(String columnName) {
        for (int i = 0; i < headers.size(); i++) {
            if (headers.get(i).equals(columnName)) {
                return i;
            }
        }
        return -1;
    }

    public static TsvTable read(Path path) throws IOException {
        try (BufferedReader reader = Files.newBufferedReader(path, StandardCharsets.UTF_8)) {
            String headerLine = reader.readLine();
            if (headerLine == null) {
                throw new IllegalArgumentException("TSV file is empty: " + path);
            }
            List<String> headers = List.of(headerLine.split("\\t", -1));
            List<String[]> rows = new ArrayList<>();
            String line;
            while ((line = reader.readLine()) != null) {
                String[] values = line.split("\\t", -1);
                if (values.length < headers.size()) {
                    String[] normalized = new String[headers.size()];
                    System.arraycopy(values, 0, normalized, 0, values.length);
                    for (int i = values.length; i < normalized.length; i++) {
                        normalized[i] = "";
                    }
                    rows.add(normalized);
                } else {
                    rows.add(values);
                }
            }
            return new TsvTable(headers, rows);
        }
    }

    public static void write(Path path, List<String> headers, List<String[]> rows) throws IOException {
        Path parent = path.getParent();
        if (parent != null) {
            Files.createDirectories(parent);
        }
        try (BufferedWriter writer = Files.newBufferedWriter(path, StandardCharsets.UTF_8)) {
            writer.write(String.join("\t", headers));
            writer.newLine();
            for (String[] row : rows) {
                writer.write(String.join("\t", row));
                writer.newLine();
            }
        }
    }
}

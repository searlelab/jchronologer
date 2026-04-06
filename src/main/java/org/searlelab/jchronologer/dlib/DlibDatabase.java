package org.searlelab.jchronologer.dlib;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.List;
import java.util.Map;

/**
 * Low-level SQLite DLIB writer.
 */
public final class DlibDatabase implements AutoCloseable {

    private final Connection connection;
    private final DlibBatchWriter batchWriter;

    public DlibDatabase(Path outputPath, DlibMetadata metadata) throws IOException, SQLException {
        try {
            Class.forName("org.sqlite.JDBC");
        } catch (ClassNotFoundException e) {
            throw new IllegalStateException("SQLite JDBC driver is not available.", e);
        }
        Files.deleteIfExists(outputPath);
        this.connection = DriverManager.getConnection("jdbc:sqlite:" + outputPath.toAbsolutePath());
        this.connection.setAutoCommit(false);
        this.batchWriter = new DlibBatchWriter();
        initializeSchema();
        insertMetadata(metadata.getValues());
        connection.commit();
    }

    private void initializeSchema() throws SQLException {
        try (Statement statement = connection.createStatement()) {
            statement.execute(DlibSchema.CREATE_METADATA_TABLE);
            statement.execute(DlibSchema.CREATE_ENTRIES_TABLE);
            statement.execute(DlibSchema.CREATE_PEPTIDE_TO_PROTEIN_TABLE);
        }
    }

    private void insertMetadata(Map<String, String> metadata) throws SQLException {
        try (PreparedStatement statement = connection.prepareStatement("INSERT INTO metadata (Key, Value) VALUES (?, ?)")) {
            for (Map.Entry<String, String> entry : metadata.entrySet()) {
                statement.setString(1, entry.getKey());
                statement.setString(2, entry.getValue());
                statement.addBatch();
            }
            statement.executeBatch();
        }
    }

    public void writeBatch(List<DlibEntryRecord> entries, List<DlibPeptideProteinRecord> peptideProteins) throws SQLException {
        batchWriter.writeBatch(connection, entries, peptideProteins);
        connection.commit();
    }

    public void createIndices() throws SQLException {
        try (Statement statement = connection.createStatement()) {
            for (String sql : DlibSchema.CREATE_INDICES) {
                statement.execute(sql);
            }
        }
        connection.commit();
    }

    public Connection getConnection() {
        return connection;
    }

    @Override
    public void close() throws SQLException {
        connection.close();
    }
}

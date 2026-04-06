package org.searlelab.jchronologer.dlib;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.nio.file.Files;
import java.nio.file.Path;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.Statement;
import java.util.List;
import java.util.Optional;
import org.junit.jupiter.api.Test;

class DlibDatabaseTest {

    @Test
    void writesSchemaMetadataAndRoundTripsEncodedArrays() throws Exception {
        Path output = Files.createTempFile("dlib-database-test", ".dlib");
        DlibEntryRecord entry = new DlibEntryRecord(
                500.25,
                2,
                "PEPTIDE",
                "PEPTIDE",
                1,
                120.5f,
                0.2f,
                new double[] {100.0, 200.0},
                new float[] {0.4f, 0.6f},
                Optional.of(300.5f),
                "jchronologer");
        DlibPeptideProteinRecord mapping = new DlibPeptideProteinRecord("PEPTIDE", false, "P12345");

        try (DlibDatabase database = new DlibDatabase(output, DlibMetadata.defaults())) {
            database.writeBatch(List.of(entry), List.of(mapping));
            database.createIndices();
        }

        try (Connection connection = DriverManager.getConnection("jdbc:sqlite:" + output.toAbsolutePath());
                Statement statement = connection.createStatement()) {
            ResultSet metadata = statement.executeQuery("select count(*) from metadata");
            assertTrue(metadata.next());
            assertEquals(4, metadata.getInt(1));

            ResultSet rows = statement.executeQuery(
                    "select MassEncodedLength, MassArray, IntensityEncodedLength, IntensityArray, IonMobility "
                            + "from entries");
            assertTrue(rows.next());
            byte[] massBytes = DlibCodec.decompress(rows.getBytes(2));
            byte[] intensityBytes = DlibCodec.decompress(rows.getBytes(4));
            assertEquals(rows.getInt(1), massBytes.length);
            assertEquals(rows.getInt(3), intensityBytes.length);
            assertArrayEquals(new double[] {100.0, 200.0}, DlibCodec.decodeDoubles(massBytes), 1e-9);
            assertArrayEquals(new float[] {0.4f, 0.6f}, DlibCodec.decodeFloats(intensityBytes), 1e-6f);
            assertEquals(300.5f, rows.getFloat(5), 1e-5f);

            ResultSet proteins = statement.executeQuery("select PeptideSeq, IsDecoy, ProteinAccession from peptidetoprotein");
            assertTrue(proteins.next());
            assertEquals("PEPTIDE", proteins.getString(1));
            assertEquals(false, proteins.getBoolean(2));
            assertEquals("P12345", proteins.getString(3));
        }
    }
}

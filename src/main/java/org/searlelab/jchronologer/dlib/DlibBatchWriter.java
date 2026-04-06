package org.searlelab.jchronologer.dlib;

import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.SQLException;
import java.sql.Types;
import java.util.List;

/**
 * Writes one flushed DLIB batch in a single transaction.
 */
final class DlibBatchWriter {

    private static final String INSERT_ENTRY =
            "INSERT INTO entries (PrecursorMz, PrecursorCharge, PeptideModSeq, PeptideSeq, Copies, RTInSeconds, "
                    + "Score, MassEncodedLength, MassArray, IntensityEncodedLength, IntensityArray, "
                    + "CorrelationEncodedLength, CorrelationArray, QuantifiedIonsArray, RTInSecondsStart, "
                    + "RTInSecondsStop, IonMobility, MedianChromatogramEncodedLength, MedianChromatogramArray, SourceFile) "
                    + "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)";

    private static final String INSERT_PEPTIDE_TO_PROTEIN =
            "INSERT OR IGNORE INTO peptidetoprotein (PeptideSeq, IsDecoy, ProteinAccession) VALUES (?,?,?)";

    void writeBatch(
            Connection connection,
            List<DlibEntryRecord> entries,
            List<DlibPeptideProteinRecord> peptideProteins) throws SQLException {
        try (PreparedStatement entryStatement = connection.prepareStatement(INSERT_ENTRY);
                PreparedStatement peptideProteinStatement = connection.prepareStatement(INSERT_PEPTIDE_TO_PROTEIN)) {
            for (DlibEntryRecord entry : entries) {
                byte[] massBytes = DlibCodec.encode(entry.getMassArray());
                byte[] intensityBytes = DlibCodec.encode(entry.getIntensityArray());
                byte[] quantifiedIonsBytes = DlibCodec.encode(DlibCodec.unitBooleanArray(entry.getMassArray().length));

                entryStatement.setDouble(1, entry.getPrecursorMz());
                entryStatement.setInt(2, entry.getPrecursorCharge());
                entryStatement.setString(3, entry.getPeptideModSeq());
                entryStatement.setString(4, entry.getPeptideSeq());
                entryStatement.setInt(5, entry.getCopies());
                entryStatement.setFloat(6, entry.getRtInSeconds());
                entryStatement.setFloat(7, entry.getScore());
                entryStatement.setInt(8, massBytes.length);
                entryStatement.setBytes(9, DlibCodec.compress(massBytes));
                entryStatement.setInt(10, intensityBytes.length);
                entryStatement.setBytes(11, DlibCodec.compress(intensityBytes));
                entryStatement.setNull(12, Types.INTEGER);
                entryStatement.setNull(13, Types.BLOB);
                entryStatement.setBytes(14, DlibCodec.compress(quantifiedIonsBytes));
                entryStatement.setNull(15, Types.FLOAT);
                entryStatement.setNull(16, Types.FLOAT);
                if (entry.getIonMobility().isPresent()) {
                    entryStatement.setFloat(17, entry.getIonMobility().get());
                } else {
                    entryStatement.setNull(17, Types.FLOAT);
                }
                entryStatement.setNull(18, Types.INTEGER);
                entryStatement.setNull(19, Types.BLOB);
                entryStatement.setString(20, entry.getSourceFile());
                entryStatement.addBatch();
            }
            entryStatement.executeBatch();

            for (DlibPeptideProteinRecord peptideProtein : peptideProteins) {
                peptideProteinStatement.setString(1, peptideProtein.peptideSeq());
                peptideProteinStatement.setBoolean(2, peptideProtein.isDecoy());
                peptideProteinStatement.setString(3, peptideProtein.proteinAccession());
                peptideProteinStatement.addBatch();
            }
            peptideProteinStatement.executeBatch();
        }
    }
}

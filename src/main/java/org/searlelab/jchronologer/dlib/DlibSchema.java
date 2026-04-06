package org.searlelab.jchronologer.dlib;

import java.util.List;

/**
 * SQL schema constants for generated DLIB files.
 */
public final class DlibSchema {

    public static final String CREATE_METADATA_TABLE =
            "CREATE TABLE IF NOT EXISTS metadata (Key string not null, Value string not null)";

    public static final String CREATE_ENTRIES_TABLE =
            "CREATE TABLE IF NOT EXISTS entries ( "
                    + "PrecursorMz double not null, "
                    + "PrecursorCharge int not null, "
                    + "PeptideModSeq string not null, "
                    + "PeptideSeq string not null, "
                    + "Copies int not null, "
                    + "RTInSeconds double not null, "
                    + "Score double not null, "
                    + "MassEncodedLength int not null, "
                    + "MassArray blob not null, "
                    + "IntensityEncodedLength int not null, "
                    + "IntensityArray blob not null, "
                    + "CorrelationEncodedLength int, "
                    + "CorrelationArray blob, "
                    + "QuantifiedIonsArray blob, "
                    + "RTInSecondsStart double, "
                    + "RTInSecondsStop double, "
                    + "IonMobility double, "
                    + "MedianChromatogramEncodedLength int, "
                    + "MedianChromatogramArray blob, "
                    + "SourceFile string not null"
                    + ")";

    public static final String CREATE_PEPTIDE_TO_PROTEIN_TABLE =
            "CREATE TABLE IF NOT EXISTS peptidetoprotein ("
                    + "PeptideSeq string not null,"
                    + "isDecoy boolean,"
                    + "ProteinAccession string not null"
                    + ")";

    public static final List<String> CREATE_INDICES = List.of(
            "create index if not exists 'PeptideModSeq_PrecursorCharge_SourceFile_Entries_index' "
                    + "on 'entries' ('PeptideModSeq' ASC, 'PrecursorCharge' ASC, 'SourceFile' ASC)",
            "create index if not exists 'PeptideSeq_Entries_index' on 'entries' ('PeptideSeq' ASC)",
            "create index if not exists 'PrecursorMz_Entries_index' on 'entries' ('PrecursorMz' ASC)",
            "create index if not exists 'ProteinAccession_PeptideToProtein_index' "
                    + "on 'peptidetoprotein' ('ProteinAccession' ASC)",
            "create index if not exists 'PeptideSeq_PeptideToProtein_index' "
                    + "on 'peptidetoprotein' ('PeptideSeq' ASC)");

    private DlibSchema() {
    }
}

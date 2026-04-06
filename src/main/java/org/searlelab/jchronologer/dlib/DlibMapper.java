package org.searlelab.jchronologer.dlib;

import java.util.Optional;
import org.searlelab.jchronologer.api.ChronologerLibraryEntry;
import org.searlelab.jchronologer.preprocessing.PeptideSequenceConverter;

/**
 * Maps predicted library entries to DLIB row models.
 */
public final class DlibMapper {

    private static final String SOURCE_FILE = "jchronologer";

    private DlibMapper() {
    }

    public static DlibEntryRecord toDlibEntry(ChronologerLibraryEntry entry) {
        float chargeProbability = entry.getChargeProbability()
                .orElseThrow(() -> new IllegalArgumentException(
                        "Automatic library entry is missing charge probability: "
                                + entry.getUnimodPeptideSequence()));
        String peptideSeq = PeptideSequenceConverter.parseNormalizedUnimod(entry.getUnimodPeptideSequence()).getResidues();
        return new DlibEntryRecord(
                entry.getPrecursorMz(),
                entry.getPrecursorCharge(),
                entry.getPeptideModSeq(),
                peptideSeq,
                1,
                entry.getRetentionTimeInSeconds(),
                1.0f - chargeProbability,
                entry.getMassArray(),
                entry.getIntensityArray(),
                entry.getCCS(),
                SOURCE_FILE);
    }

    public static DlibPeptideProteinRecord toPeptideProtein(String peptideSeq, String accession) {
        return new DlibPeptideProteinRecord(peptideSeq, false, accession);
    }
}

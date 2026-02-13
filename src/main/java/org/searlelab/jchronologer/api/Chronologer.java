package org.searlelab.jchronologer.api;

import java.util.List;

/**
 * Primary inference entry point for Chronologer peptide retention prediction.
 *
 * <p>Implementations preprocess each peptide, run batched model inference for accepted rows, and
 * return both accepted and rejected outcomes in a single {@link PredictionResult}.
 */
public interface Chronologer extends AutoCloseable {

    /**
     * Predicts hydrophobic index values ({@code Pred_HI}) for a batch of peptide modification
     * sequences.
     *
     * <p>Input ordering is preserved through {@code rowIndex} fields in accepted and rejected
     * records.
     *
     * @param peptideModSeqs peptide sequences in EncyclopeDIA-style modified format
     * @return accepted predictions and structured rejection diagnostics
     */
    PredictionResult predict(List<String> peptideModSeqs);

    @Override
    void close();
}

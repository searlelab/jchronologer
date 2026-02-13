package org.searlelab.jchronologer.api;

/**
 * High-level categories for why a peptide input could not be converted into a model prediction.
 */
public enum RejectionReason {
    UNSUPPORTED_MODIFICATION,
    NTERM_MAPPING_FAILED,
    INVALID_RESIDUE_TOKEN,
    LENGTH_OUT_OF_RANGE,
    TOKENIZATION_ERROR
}

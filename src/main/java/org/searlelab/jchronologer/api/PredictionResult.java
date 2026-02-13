package org.searlelab.jchronologer.api;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Immutable container for one Chronologer prediction batch result.
 *
 * <p>The accepted and rejected lists are aligned to input rows via {@code rowIndex}, allowing
 * callers to reconstruct original ordering or filter behavior as needed.
 */
public final class PredictionResult {

    private final List<AcceptedPrediction> accepted;
    private final List<RejectedPrediction> rejected;

    public PredictionResult(List<AcceptedPrediction> accepted, List<RejectedPrediction> rejected) {
        this.accepted = Collections.unmodifiableList(new ArrayList<>(accepted));
        this.rejected = Collections.unmodifiableList(new ArrayList<>(rejected));
    }

    public List<AcceptedPrediction> getAccepted() {
        return accepted;
    }

    public List<RejectedPrediction> getRejected() {
        return rejected;
    }

    public int getAcceptedCount() {
        return accepted.size();
    }

    public int getRejectedCount() {
        return rejected.size();
    }

    /**
     * Returns accepted predictions keyed by input row index.
     *
     * <p>If duplicate row indexes are present, the last value in iteration order is retained.
     */
    public Map<Integer, AcceptedPrediction> getAcceptedByRowIndex() {
        Map<Integer, AcceptedPrediction> byRow = new HashMap<>();
        for (AcceptedPrediction item : accepted) {
            byRow.put(item.getRowIndex(), item);
        }
        return Collections.unmodifiableMap(byRow);
    }

    /**
     * Returns rejected predictions keyed by input row index.
     *
     * <p>If duplicate row indexes are present, the last value in iteration order is retained.
     */
    public Map<Integer, RejectedPrediction> getRejectedByRowIndex() {
        Map<Integer, RejectedPrediction> byRow = new HashMap<>();
        for (RejectedPrediction item : rejected) {
            byRow.put(item.getRowIndex(), item);
        }
        return Collections.unmodifiableMap(byRow);
    }
}

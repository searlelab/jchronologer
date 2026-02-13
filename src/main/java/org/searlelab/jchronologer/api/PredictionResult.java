package org.searlelab.jchronologer.api;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

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
}

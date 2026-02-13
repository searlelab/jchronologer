package org.searlelab.jchronologer.api;

import java.util.List;

public interface Chronologer extends AutoCloseable {

    PredictionResult predict(List<String> peptideModSeqs);

    @Override
    void close();
}

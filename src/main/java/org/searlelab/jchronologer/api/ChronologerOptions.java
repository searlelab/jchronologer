package org.searlelab.jchronologer.api;

public final class ChronologerOptions {

    public static final String DEFAULT_MODEL_RESOURCE =
            "models/Chronologer_20220601193755.torchscript.pt";
    public static final String DEFAULT_PREPROCESSING_RESOURCE =
            "models/Chronologer_20220601193755.preprocessing.json";
    public static final int DEFAULT_BATCH_SIZE = 2048;

    private final String modelResource;
    private final String preprocessingResource;
    private final int batchSize;

    private ChronologerOptions(Builder builder) {
        this.modelResource = builder.modelResource;
        this.preprocessingResource = builder.preprocessingResource;
        this.batchSize = builder.batchSize;
    }

    public static Builder builder() {
        return new Builder();
    }

    public String getModelResource() {
        return modelResource;
    }

    public String getPreprocessingResource() {
        return preprocessingResource;
    }

    public int getBatchSize() {
        return batchSize;
    }

    public static final class Builder {
        private String modelResource = DEFAULT_MODEL_RESOURCE;
        private String preprocessingResource = DEFAULT_PREPROCESSING_RESOURCE;
        private int batchSize = DEFAULT_BATCH_SIZE;

        private Builder() {
        }

        public Builder modelResource(String modelResource) {
            this.modelResource = modelResource;
            return this;
        }

        public Builder preprocessingResource(String preprocessingResource) {
            this.preprocessingResource = preprocessingResource;
            return this;
        }

        public Builder batchSize(int batchSize) {
            this.batchSize = batchSize;
            return this;
        }

        public ChronologerOptions build() {
            if (modelResource == null || modelResource.isBlank()) {
                throw new IllegalArgumentException("Model resource must be non-empty.");
            }
            if (preprocessingResource == null || preprocessingResource.isBlank()) {
                throw new IllegalArgumentException("Preprocessing resource must be non-empty.");
            }
            if (batchSize <= 0) {
                throw new IllegalArgumentException("Batch size must be positive.");
            }
            return new ChronologerOptions(this);
        }
    }
}

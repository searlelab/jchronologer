package org.searlelab.jchronologer.api;

/**
 * Runtime configuration for creating a {@link Chronologer} instance.
 *
 * <p>Defaults point at bundled classpath resources for the reference model and preprocessing
 * metadata, a reference batch size, and a small inference thread pool so a caller can use
 * {@link Builder#build()} without overriding values.
 */
public final class ChronologerOptions {

    public static final String DEFAULT_MODEL_RESOURCE =
            "models/Chronologer_20220601193755.torchscript.pt";
    public static final String DEFAULT_PREPROCESSING_RESOURCE =
            "models/Chronologer_20220601193755.preprocessing.json";
    public static final int DEFAULT_BATCH_SIZE = 2048;
    public static final int DEFAULT_INFERENCE_THREADS =
            Math.max(1, Runtime.getRuntime().availableProcessors()-2);

    private final String modelResource;
    private final String preprocessingResource;
    private final int batchSize;
    private final int inferenceThreads;

    private ChronologerOptions(Builder builder) {
        this.modelResource = builder.modelResource;
        this.preprocessingResource = builder.preprocessingResource;
        this.batchSize = builder.batchSize;
        this.inferenceThreads = builder.inferenceThreads;
    }

    /**
     * Creates a new mutable builder initialized with default resource paths and batch size.
     *
     * @return options builder
     */
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

    public int getInferenceThreads() {
        return inferenceThreads;
    }

    /**
     * Builder for {@link ChronologerOptions}.
     *
     * <p>Validation is applied in {@link #build()} to enforce non-empty resources and positive
     * runtime sizing values.
     */
    public static final class Builder {
        private String modelResource = DEFAULT_MODEL_RESOURCE;
        private String preprocessingResource = DEFAULT_PREPROCESSING_RESOURCE;
        private int batchSize = DEFAULT_BATCH_SIZE;
        private int inferenceThreads = DEFAULT_INFERENCE_THREADS;

        private Builder() {
        }

        /**
         * Sets the classpath resource path for the TorchScript model artifact.
         *
         * @param modelResource model resource path
         * @return this builder
         */
        public Builder modelResource(String modelResource) {
            this.modelResource = modelResource;
            return this;
        }

        /**
         * Sets the classpath resource path for preprocessing metadata JSON.
         *
         * @param preprocessingResource preprocessing metadata resource path
         * @return this builder
         */
        public Builder preprocessingResource(String preprocessingResource) {
            this.preprocessingResource = preprocessingResource;
            return this;
        }

        /**
         * Sets the maximum number of accepted peptides per inference batch.
         *
         * @param batchSize batch size, must be positive
         * @return this builder
         */
        public Builder batchSize(int batchSize) {
            this.batchSize = batchSize;
            return this;
        }

        /**
         * Sets the maximum number of threads used to score inference batches in parallel.
         *
         * @param inferenceThreads inference threads, must be positive
         * @return this builder
         */
        public Builder inferenceThreads(int inferenceThreads) {
            this.inferenceThreads = inferenceThreads;
            return this;
        }

        /**
         * Builds immutable options after validating required fields.
         *
         * @return validated options instance
         */
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
            if (inferenceThreads <= 0) {
                throw new IllegalArgumentException("Inference threads must be positive.");
            }
            return new ChronologerOptions(this);
        }
    }
}

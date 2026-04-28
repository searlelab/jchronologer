package org.searlelab.jchronologer.api;

/**
 * Runtime configuration for creating a {@link Chronologer} instance.
 *
 * <p>Defaults point at bundled classpath resources for the reference model and preprocessing
 * metadata, a reference batch size, and a small inference thread pool so a caller can use
 * {@link Builder#build()} without overriding values.
 */
public final class ChronologerOptions {

	public static final String DEFAULT_MODEL_RESOURCE="models/Chronologer_20220601193755.torchscript.pt";
	public static final String DEFAULT_PREPROCESSING_RESOURCE="models/Chronologer_20220601193755.preprocessing.json";

	public static final String DEFAULT_SCULPTOR_MODEL_RESOURCE="models/Sculptor_20260311095327.torchscript.pt";
	public static final String DEFAULT_SCULPTOR_PREPROCESSING_RESOURCE="models/Sculptor_20260311095327.preprocessing.json";

	public static final String DEFAULT_CARTOGRAPHER_MODEL_RESOURCE="models/Cartographer_multistart_best_run01.torchscript.pt";
	public static final String DEFAULT_CARTOGRAPHER_PREPROCESSING_RESOURCE="models/Cartographer_multistart_best_run01.preprocessing.json";

	public static final String DEFAULT_SCOUT_MODEL_RESOURCE="models/Scout_20260427110352_run01.torchscript.pt";
	public static final String DEFAULT_SCOUT_PREPROCESSING_RESOURCE="models/Scout_20260427110352_run01.preprocessing.json";

	public static final String DEFAULT_ELECTRICIAN_MODEL_RESOURCE="models/Electrician_20260225110528.torchscript.pt";
	public static final String DEFAULT_ELECTRICIAN_PREPROCESSING_RESOURCE="models/Electrician_20260225110528.preprocessing.json";
    
    public static final boolean DEFAULT_CCS_PREDICTION_ENABLED = true;
    public static final int DEFAULT_BATCH_SIZE = 2048;
    public static final int DEFAULT_INFERENCE_THREADS =
            Math.max(1, Runtime.getRuntime().availableProcessors()-2);

    private final String modelResource;
    private final String preprocessingResource;
    private final String sculptorModelResource;
    private final String sculptorPreprocessingResource;
    private final boolean ccsPredictionEnabled;
    private final int batchSize;
    private final int inferenceThreads;
    private final boolean verboseLogging;

    private ChronologerOptions(Builder builder) {
        this.modelResource = builder.modelResource;
        this.preprocessingResource = builder.preprocessingResource;
        this.sculptorModelResource = builder.sculptorModelResource;
        this.sculptorPreprocessingResource = builder.sculptorPreprocessingResource;
        this.ccsPredictionEnabled = builder.ccsPredictionEnabled;
        this.batchSize = builder.batchSize;
        this.inferenceThreads = builder.inferenceThreads;
        this.verboseLogging = builder.verboseLogging;
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

    public String getSculptorModelResource() {
        return sculptorModelResource;
    }

    public String getSculptorPreprocessingResource() {
        return sculptorPreprocessingResource;
    }

    public boolean isCCSPredictionEnabled() {
        return ccsPredictionEnabled;
    }

    public int getBatchSize() {
        return batchSize;
    }

    public int getInferenceThreads() {
        return inferenceThreads;
    }

    public boolean isVerboseLogging() {
        return verboseLogging;
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
        private String sculptorModelResource = DEFAULT_SCULPTOR_MODEL_RESOURCE;
        private String sculptorPreprocessingResource = DEFAULT_SCULPTOR_PREPROCESSING_RESOURCE;
        private boolean ccsPredictionEnabled = DEFAULT_CCS_PREDICTION_ENABLED;
        private int batchSize = DEFAULT_BATCH_SIZE;
        private int inferenceThreads = DEFAULT_INFERENCE_THREADS;
        private boolean verboseLogging;

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

        public Builder sculptorModelResource(String sculptorModelResource) {
            this.sculptorModelResource = sculptorModelResource;
            return this;
        }

        public Builder sculptorPreprocessingResource(String sculptorPreprocessingResource) {
            this.sculptorPreprocessingResource = sculptorPreprocessingResource;
            return this;
        }

        public Builder ccsPredictionEnabled(boolean ccsPredictionEnabled) {
            this.ccsPredictionEnabled = ccsPredictionEnabled;
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
         * Enables verbose startup/inference logging for diagnostics.
         *
         * @param verboseLogging whether verbose logging is enabled
         * @return this builder
         */
        public Builder verboseLogging(boolean verboseLogging) {
            this.verboseLogging = verboseLogging;
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
            if (sculptorModelResource == null || sculptorModelResource.isBlank()) {
                throw new IllegalArgumentException("Sculptor model resource must be non-empty.");
            }
            if (sculptorPreprocessingResource == null || sculptorPreprocessingResource.isBlank()) {
                throw new IllegalArgumentException("Sculptor preprocessing resource must be non-empty.");
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

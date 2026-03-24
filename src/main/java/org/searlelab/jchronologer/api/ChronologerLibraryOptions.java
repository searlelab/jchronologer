package org.searlelab.jchronologer.api;

/**
 * Runtime configuration for tandem library prediction.
 */
public final class ChronologerLibraryOptions {

    public static final double DEFAULT_MASS_MATCH_EPSILON = 1e-5;
    public static final float DEFAULT_MINIMUM_REPORTED_INTENSITY = 0.01f;

    private final String chronologerModelResource;
    private final String chronologerPreprocessingResource;
    private final String cartographerModelResource;
    private final String cartographerPreprocessingResource;
    private final String electricianModelResource;
    private final String electricianPreprocessingResource;
    private final String sculptorModelResource;
    private final String sculptorPreprocessingResource;
    private final boolean ccsPredictionEnabled;
    private final int batchSize;
    private final int cartographerBatchSize;
    private final int inferenceThreads;
    private final double massMatchEpsilon;
    private final float minimumReportedIntensity;
    private final boolean verboseLogging;

    private ChronologerLibraryOptions(Builder builder) {
        this.chronologerModelResource = builder.chronologerModelResource;
        this.chronologerPreprocessingResource = builder.chronologerPreprocessingResource;
        this.cartographerModelResource = builder.cartographerModelResource;
        this.cartographerPreprocessingResource = builder.cartographerPreprocessingResource;
        this.electricianModelResource = builder.electricianModelResource;
        this.electricianPreprocessingResource = builder.electricianPreprocessingResource;
        this.sculptorModelResource = builder.sculptorModelResource;
        this.sculptorPreprocessingResource = builder.sculptorPreprocessingResource;
        this.ccsPredictionEnabled = builder.ccsPredictionEnabled;
        this.batchSize = builder.batchSize;
        this.cartographerBatchSize = builder.cartographerBatchSize;
        this.inferenceThreads = builder.inferenceThreads;
        this.massMatchEpsilon = builder.massMatchEpsilon;
        this.minimumReportedIntensity = builder.minimumReportedIntensity;
        this.verboseLogging = builder.verboseLogging;
    }

    public static Builder builder() {
        return new Builder();
    }

    public String getChronologerModelResource() {
        return chronologerModelResource;
    }

    public String getChronologerPreprocessingResource() {
        return chronologerPreprocessingResource;
    }

    public String getCartographerModelResource() {
        return cartographerModelResource;
    }

    public String getCartographerPreprocessingResource() {
        return cartographerPreprocessingResource;
    }

    public String getElectricianModelResource() {
        return electricianModelResource;
    }

    public String getElectricianPreprocessingResource() {
        return electricianPreprocessingResource;
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

    public int getCartographerBatchSize() {
        return cartographerBatchSize;
    }

    public int getInferenceThreads() {
        return inferenceThreads;
    }

    public double getMassMatchEpsilon() {
        return massMatchEpsilon;
    }

    public float getMinimumReportedIntensity() {
        return minimumReportedIntensity;
    }

    public boolean isVerboseLogging() {
        return verboseLogging;
    }

    public static final class Builder {
        private String chronologerModelResource = ChronologerOptions.DEFAULT_MODEL_RESOURCE;
        private String chronologerPreprocessingResource = ChronologerOptions.DEFAULT_PREPROCESSING_RESOURCE;
        private String cartographerModelResource = ChronologerOptions.DEFAULT_CARTOGRAPHER_MODEL_RESOURCE;
        private String cartographerPreprocessingResource = ChronologerOptions.DEFAULT_CARTOGRAPHER_PREPROCESSING_RESOURCE;
        private String electricianModelResource = ChronologerOptions.DEFAULT_ELECTRICIAN_MODEL_RESOURCE;
        private String electricianPreprocessingResource = ChronologerOptions.DEFAULT_ELECTRICIAN_PREPROCESSING_RESOURCE;
        private String sculptorModelResource = ChronologerOptions.DEFAULT_SCULPTOR_MODEL_RESOURCE;
        private String sculptorPreprocessingResource = ChronologerOptions.DEFAULT_SCULPTOR_PREPROCESSING_RESOURCE;
        private boolean ccsPredictionEnabled = ChronologerOptions.DEFAULT_CCS_PREDICTION_ENABLED;
        private int batchSize = ChronologerOptions.DEFAULT_BATCH_SIZE;
        private int cartographerBatchSize = ChronologerOptions.DEFAULT_BATCH_SIZE;
        private int inferenceThreads = ChronologerOptions.DEFAULT_INFERENCE_THREADS;
        private double massMatchEpsilon = DEFAULT_MASS_MATCH_EPSILON;
        private float minimumReportedIntensity = DEFAULT_MINIMUM_REPORTED_INTENSITY;
        private boolean verboseLogging;

        private Builder() {
        }

        public Builder chronologerModelResource(String chronologerModelResource) {
            this.chronologerModelResource = chronologerModelResource;
            return this;
        }

        public Builder chronologerPreprocessingResource(String chronologerPreprocessingResource) {
            this.chronologerPreprocessingResource = chronologerPreprocessingResource;
            return this;
        }

        public Builder cartographerModelResource(String cartographerModelResource) {
            this.cartographerModelResource = cartographerModelResource;
            return this;
        }

        public Builder cartographerPreprocessingResource(String cartographerPreprocessingResource) {
            this.cartographerPreprocessingResource = cartographerPreprocessingResource;
            return this;
        }

        public Builder electricianModelResource(String electricianModelResource) {
            this.electricianModelResource = electricianModelResource;
            return this;
        }

        public Builder electricianPreprocessingResource(String electricianPreprocessingResource) {
            this.electricianPreprocessingResource = electricianPreprocessingResource;
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

        public Builder batchSize(int batchSize) {
            this.batchSize = batchSize;
            return this;
        }

        public Builder cartographerBatchSize(int cartographerBatchSize) {
            this.cartographerBatchSize = cartographerBatchSize;
            return this;
        }

        public Builder inferenceThreads(int inferenceThreads) {
            this.inferenceThreads = inferenceThreads;
            return this;
        }

        public Builder massMatchEpsilon(double massMatchEpsilon) {
            this.massMatchEpsilon = massMatchEpsilon;
            return this;
        }

        public Builder minimumReportedIntensity(float minimumReportedIntensity) {
            this.minimumReportedIntensity = minimumReportedIntensity;
            return this;
        }

        public Builder verboseLogging(boolean verboseLogging) {
            this.verboseLogging = verboseLogging;
            return this;
        }

        public ChronologerLibraryOptions build() {
            if (chronologerModelResource == null || chronologerModelResource.isBlank()) {
                throw new IllegalArgumentException("Chronologer model resource must be non-empty.");
            }
            if (chronologerPreprocessingResource == null || chronologerPreprocessingResource.isBlank()) {
                throw new IllegalArgumentException("Chronologer preprocessing resource must be non-empty.");
            }
            if (cartographerModelResource == null || cartographerModelResource.isBlank()) {
                throw new IllegalArgumentException("Cartographer model resource must be non-empty.");
            }
            if (cartographerPreprocessingResource == null || cartographerPreprocessingResource.isBlank()) {
                throw new IllegalArgumentException("Cartographer preprocessing resource must be non-empty.");
            }
            if (electricianModelResource == null || electricianModelResource.isBlank()) {
                throw new IllegalArgumentException("Electrician model resource must be non-empty.");
            }
            if (electricianPreprocessingResource == null || electricianPreprocessingResource.isBlank()) {
                throw new IllegalArgumentException("Electrician preprocessing resource must be non-empty.");
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
            if (cartographerBatchSize <= 0) {
                throw new IllegalArgumentException("Cartographer batch size must be positive.");
            }
            if (inferenceThreads <= 0) {
                throw new IllegalArgumentException("Inference threads must be positive.");
            }
            if (massMatchEpsilon <= 0.0) {
                throw new IllegalArgumentException("Mass match epsilon must be positive.");
            }
            if (minimumReportedIntensity < 0.0f) {
                throw new IllegalArgumentException("Minimum reported intensity must be >= 0.");
            }
            return new ChronologerLibraryOptions(this);
        }
    }
}

package org.searlelab.jchronologer.inference;

import static org.junit.jupiter.api.Assertions.assertDoesNotThrow;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import ai.djl.Model;
import ai.djl.inference.Predictor;
import ai.djl.ndarray.NDList;
import ai.djl.translate.Translator;
import ai.djl.translate.TranslatorContext;
import java.lang.reflect.Field;
import java.util.Set;
import org.junit.jupiter.api.Test;

class BatchPredictorTest {

    private static final String DEFAULT_MODEL_RESOURCE =
            "models/Chronologer_20220601193755.torchscript.pt";

    @Test
    void oneArgConstructorLoadsBundledModel() {
        assertDoesNotThrow(() -> {
            try (BatchPredictor predictor = new BatchPredictor(DEFAULT_MODEL_RESOURCE)) {
                // no-op
            }
        });
    }

    @Test
    void oneArgConstructorRejectsMissingResource() {
        IllegalArgumentException error = assertThrows(
                IllegalArgumentException.class,
                () -> new BatchPredictor("does/not/exist.txt"));
        assertTrue(error.getMessage().contains("Missing classpath resource"));
    }

    @Test
    void constructorWrapsModelLoadFailures() {
        IllegalStateException error = assertThrows(
                IllegalStateException.class,
                () -> new BatchPredictor("data/demo/demo_peptides.txt", true));
        assertTrue(error.getMessage().contains("Failed to load Chronologer model from resource"));
    }

    @Test
    void constructorRejectsNonModelResource() {
        RuntimeException error = assertThrows(
                RuntimeException.class,
                () -> new BatchPredictor("data/demo/demo_peptides.txt"));
        assertTrue(!error.getMessage().isBlank());
    }

    @Test
    void predictWrapsTranslateException() throws Exception {
        try (BatchPredictor predictor = new BatchPredictor(DEFAULT_MODEL_RESOURCE)) {
            Field modelField = BatchPredictor.class.getDeclaredField("model");
            modelField.setAccessible(true);
            Model model = (Model) modelField.get(predictor);

            Predictor<NDList, NDList> failingPredictor = model.newPredictor(new ThrowingTranslator());

            Field threadLocalField = BatchPredictor.class.getDeclaredField("threadLocalPredictor");
            threadLocalField.setAccessible(true);
            @SuppressWarnings("unchecked")
            ThreadLocal<Predictor<NDList, NDList>> threadLocal =
                    (ThreadLocal<Predictor<NDList, NDList>>) threadLocalField.get(predictor);
            threadLocal.set(failingPredictor);

            Field predictorsField = BatchPredictor.class.getDeclaredField("predictors");
            predictorsField.setAccessible(true);
            @SuppressWarnings("unchecked")
            Set<Predictor<NDList, NDList>> predictors =
                    (Set<Predictor<NDList, NDList>>) predictorsField.get(predictor);
            predictors.add(failingPredictor);

            IllegalStateException error = assertThrows(
                    IllegalStateException.class,
                    () -> predictor.predict(new long[][] {{1L, 2L, 3L}}));
            assertTrue(error.getMessage().contains("Failed to run Chronologer inference."));
        }
    }

    private static final class ThrowingTranslator implements Translator<NDList, NDList> {
        @Override
        public NDList processInput(TranslatorContext ctx, NDList input) {
            return input;
        }

        @Override
        public NDList processOutput(TranslatorContext ctx, NDList list) {
            throw new IllegalStateException("simulated translator failure");
        }
    }
}

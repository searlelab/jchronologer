package org.searlelab.jchronologer.impl;

import static org.junit.jupiter.api.Assertions.assertDoesNotThrow;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertInstanceOf;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.AbstractExecutorService;
import java.util.concurrent.Callable;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.TimeoutException;
import org.junit.jupiter.api.Test;
import org.searlelab.jchronologer.api.Chronologer;
import org.searlelab.jchronologer.api.ChronologerOptions;
import org.searlelab.jchronologer.api.PredictionResult;
import org.searlelab.jchronologer.inference.BatchPredictor;

class DefaultChronologerCoverageTest {

    @Test
    void closeIsIdempotent() {
        Chronologer chronologer = ChronologerFactory.createDefault();
        chronologer.close();
        assertDoesNotThrow(chronologer::close);
    }

    @Test
    void predictWithAllRejectedRowsReturnsNoAcceptedPredictions() {
        ChronologerOptions options = ChronologerOptions.builder()
                .batchSize(4)
                .inferenceThreads(1)
                .build();

        try (Chronologer chronologer = ChronologerFactory.create(options)) {
            PredictionResult result = chronologer.predict(List.of(
                    "[42.010565]ACDEFGHIK",
                    "[42.010565]PEPTIDEK"));
            assertEquals(0, result.getAcceptedCount());
            assertEquals(2, result.getRejectedCount());
        }
    }

    @Test
    void scoreAcceptedDraftsWrapsExecutionException() throws Exception {
        ChronologerOptions options = ChronologerOptions.builder()
                .batchSize(1)
                .inferenceThreads(1)
                .build();

        try (DefaultChronologer chronologer = (DefaultChronologer) ChronologerFactory.create(options)) {
            Method method = scoreAcceptedDraftsMethod();
            BatchPredictor predictor = getBatchPredictor(chronologer);
            List<Object> acceptedDrafts = buildAcceptedDrafts(2);
            ExecutorService executor = new ExecutionFailingExecutor();

            InvocationTargetException error = assertThrows(
                    InvocationTargetException.class,
                    () -> method.invoke(chronologer, acceptedDrafts, predictor, executor));
            IllegalStateException cause = assertInstanceOf(IllegalStateException.class, error.getCause());
            assertTrue(cause.getMessage().contains("Failed to run Chronologer inference."));
            executor.shutdownNow();
        }
    }

    @Test
    void scoreAcceptedDraftsWrapsInterruptedExceptionAndReinterruptsThread() throws Exception {
        ChronologerOptions options = ChronologerOptions.builder()
                .batchSize(1)
                .inferenceThreads(1)
                .build();

        try (DefaultChronologer chronologer = (DefaultChronologer) ChronologerFactory.create(options)) {
            Method method = scoreAcceptedDraftsMethod();
            BatchPredictor predictor = getBatchPredictor(chronologer);
            List<Object> acceptedDrafts = buildAcceptedDrafts(2);
            ExecutorService executor = new InterruptedFailingExecutor();

            InvocationTargetException error = assertThrows(
                    InvocationTargetException.class,
                    () -> method.invoke(chronologer, acceptedDrafts, predictor, executor));
            IllegalStateException cause = assertInstanceOf(IllegalStateException.class, error.getCause());
            assertTrue(cause.getMessage().contains("Chronologer inference was interrupted."));
            assertTrue(Thread.currentThread().isInterrupted());

            // Clear interrupt status to avoid leaking thread state into other tests.
            Thread.interrupted();
            executor.shutdownNow();
        }
    }

    private static Method scoreAcceptedDraftsMethod() throws Exception {
        Method method = DefaultChronologer.class.getDeclaredMethod(
                "scoreAcceptedDrafts",
                List.class,
                BatchPredictor.class,
                ExecutorService.class);
        method.setAccessible(true);
        return method;
    }

    private static BatchPredictor getBatchPredictor(DefaultChronologer chronologer) throws Exception {
        Field field = DefaultChronologer.class.getDeclaredField("batchPredictor");
        field.setAccessible(true);
        return (BatchPredictor) field.get(chronologer);
    }

    private static List<Object> buildAcceptedDrafts(int count) throws Exception {
        Class<?> acceptedDraftClass =
                Class.forName("org.searlelab.jchronologer.impl.DefaultChronologer$AcceptedDraft");
        Constructor<?> constructor = acceptedDraftClass.getDeclaredConstructor(
                int.class,
                String.class,
                String.class,
                String.class,
                long[].class);
        constructor.setAccessible(true);

        List<Object> drafts = new ArrayList<>(count);
        for (int i = 0; i < count; i++) {
            drafts.add(constructor.newInstance(
                    i,
                    "VATVSLPR",
                    "VATVSLPR",
                    "_VATVSLPR_",
                    new long[] {1L, 2L, 3L}));
        }
        return drafts;
    }

    private static abstract class BaseFailingExecutor extends AbstractExecutorService {
        private boolean shutdown;

        @Override
        public void shutdown() {
            shutdown = true;
        }

        @Override
        public List<Runnable> shutdownNow() {
            shutdown = true;
            return List.of();
        }

        @Override
        public boolean isShutdown() {
            return shutdown;
        }

        @Override
        public boolean isTerminated() {
            return shutdown;
        }

        @Override
        public boolean awaitTermination(long timeout, TimeUnit unit) {
            return true;
        }

        @Override
        public void execute(Runnable command) {
            command.run();
        }
    }

    private static final class ExecutionFailingExecutor extends BaseFailingExecutor {
        @Override
        public <T> Future<T> submit(Callable<T> task) {
            return CompletableFuture.failedFuture(new RuntimeException("simulated execution failure"));
        }
    }

    private static final class InterruptedFailingExecutor extends BaseFailingExecutor {
        @Override
        public <T> Future<T> submit(Callable<T> task) {
            return new Future<>() {
                private boolean cancelled;

                @Override
                public boolean cancel(boolean mayInterruptIfRunning) {
                    cancelled = true;
                    return true;
                }

                @Override
                public boolean isCancelled() {
                    return cancelled;
                }

                @Override
                public boolean isDone() {
                    return cancelled;
                }

                @Override
                public T get() throws InterruptedException {
                    throw new InterruptedException("simulated interrupt");
                }

                @Override
                public T get(long timeout, TimeUnit unit) throws InterruptedException, TimeoutException {
                    throw new InterruptedException("simulated interrupt");
                }
            };
        }
    }
}

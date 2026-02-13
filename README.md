# jchronloger

Java 17 inference implementation of the Chronologer retention-time model using DJL (PyTorch engine).

## Scope

This project implements Chronologer inference only (no training) with preprocessing and tokenization behavior matched to the Python implementation.

## Build

Fast compile:

```bash
mvn -am -Dmaven.test.skip=true compile
```

Run all tests:

```bash
mvn -am test
```

Run one test class:

```bash
mvn -am -Dtest=GoldenInferenceParityTest test
```

Run the demo peptide smoke test (prints `PeptideModSeq` and `Pred_HI` to console):

```bash
mvn -am -Dtest=DemoPeptidesSmokeTest test -Dsurefire.useFile=false
```

Package thin + fat jars:

```bash
mvn -am package
```

Artifacts:
- Thin library jar: `target/jchronologer-0.1.0-SNAPSHOT.jar`
- Executable fat jar: `target/jchronologer-0.1.0-SNAPSHOT-all.jar`

## API Usage

```java
import java.util.List;
import org.searlelab.jchronologer.api.Chronologer;
import org.searlelab.jchronologer.api.PredictionResult;
import org.searlelab.jchronologer.impl.ChronologerFactory;

try (Chronologer chronologer = ChronologerFactory.createDefault()) {
    PredictionResult result = chronologer.predict(List.of("VATVSLPR", "ACDE[+123.456]FGHIK"));
    // result.getAccepted() contains Pred_HI values
    // result.getRejected() contains structured rejection reasons
}
```

## CLI Usage

```bash
java -jar target/jchronologer-0.1.0-SNAPSHOT-all.jar \
  predict \
  --input input.tsv \
  --output output.tsv \
  --peptide-column PeptideModSeq
```

By default rejected peptides are dropped (matching Python `Predict_RT.py`).
Use `--keep-rejected` to preserve all rows and include rejection diagnostics.

## Model Artifacts

Bundled classpath resources:
- `models/Chronologer_20220601193755.torchscript.pt`
- `models/Chronologer_20220601193755.preprocessing.json`

Golden parity fixture:
- `data/golden/chronologer_parity_cases.golden.json`

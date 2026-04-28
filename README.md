# jchronologer

`jchronologer` is a Java 17 inference implementation of the Chronologer peptide retention-time model using DJL (PyTorch engine).

The goal is parity with the Python Chronologer inference path (`Predict_RT.py`): same preprocessing conventions, same model artifacts, and deterministic acceptance/rejection behavior for peptide inputs.

## Scope

This repository provides inference only:
- peptide preprocessing and tokenization
- batched model execution (`Pred_HI` output)
- API and CLI entry points for TSV/plain-text workflows

Training and model export are intentionally out of scope for this Java project.

## High-Level Structure

Core packages:
- `org.searlelab.jchronologer.api`: public runtime API (`Chronologer`, options, accepted/rejected result models)
- `org.searlelab.jchronologer.impl`: default implementation and factory wiring
- `org.searlelab.jchronologer.preprocessing`: metadata loading, rule compilation, peptide patching/tokenization
- `org.searlelab.jchronologer.inference`: DJL model loading and batched prediction
- `org.searlelab.jchronologer.util`: TSV and classpath-resource utilities
- `org.searlelab.jchronologer.Main`: executable default CLI (plain-text or TSV input)

Prediction flow:
1. Load preprocessing metadata JSON and compile regex/N-term rules.
2. For each `PeptideModSeq` or Unimod encoding, tokenize to model vocabulary, and validate bounds.
3. Split accepted token arrays into batches.
4. Run TorchScript inference through DJL (parallelized across batch chunks) and map each score back to original row index.
5. Return `PredictionResult` containing accepted rows (`Pred_HI`) and rejected rows (`RejectionReason` + optional detail).

## Build And Test

Fast compile (skip tests):

```bash
mvn -am -Dmaven.test.skip=true compile
```

Run all tests:

```bash
mvn -am test
```

Run a single class (example: golden parity):

```bash
mvn -am -Dtest=GoldenInferenceParityTest test
```

Generate coverage:

```bash
mvn -am verify
```

Coverage output:
- XML: `target/site/jacoco/jacoco.xml`
- HTML: `target/site/jacoco/index.html`

Package default artifacts:

```bash
mvn -am package
```

Artifacts:
- Thin library jar: `target/jchronologer-1.0.0.jar`
- Javadocs jar: `target/jchronologer-1.0.0-javadoc.jar`

Optional fat jar (only when explicitly requested):

```bash
mvn -am -Pfat-jar package
```

Artifact:
- Executable fat jar: `target/jchronologer-1.0.0-all.jar`

## API Usage

```java
import java.util.List;
import org.searlelab.jchronologer.api.Chronologer;
import org.searlelab.jchronologer.api.ChronologerOptions;
import org.searlelab.jchronologer.api.PredictionResult;
import org.searlelab.jchronologer.impl.ChronologerFactory;

ChronologerOptions options = ChronologerOptions.builder()
        .batchSize(2048)
        .inferenceThreads(2)
        .verboseLogging(false)
        .build();

try (Chronologer chronologer = ChronologerFactory.create(options)) {
    // Optional explicit lifecycle hook.
    chronologer.init();

    PredictionResult first = chronologer.predict(List.of(
            "VATVSLPR",
            "[+42.010565]KGSPTPGFSTR",
            "ACDE[+123.456]FGHIK"));

    PredictionResult second = chronologer.predict(List.of(
            "LGEHNIDVLEGNEQFINAAK",
            "FLEQQNQVLQTK"));

    // Reuse one initialized Chronologer for multiple inference rounds.
    first.getAccepted();
    second.getRejected();
}
```

## CLI Usage

Default executable CLI (`org.searlelab.jchronologer.Main`, used by shaded jar):

```bash
java -jar target/jchronologer-1.0.0-all.jar input.tsv output.tsv
```

Notes:
- `input.tsv` can be:
  - a TSV containing a peptide column (default name `PeptideModSeq`), or
  - a plain text file with one peptide sequence per line.
- if `output.tsv` is omitted, predictions are written to stdout.
- default behavior drops rejected peptides (mirrors Python `Predict_RT.py` output semantics).
- options:
  - `--batch_size <n>`
  - `--peptide_column <name>`
  - `--verbose` (show detailed startup diagnostics)

## Input Compatibility

- expected peptide column: `PeptideModSeq` (configurable via CLI option)
- supported peptide length after tokenization: 6 to 50 amino acids
- modification formatting: EncyclopeDIA-style mass annotations (same style as Python Chronologer)
- unsupported/invalid peptides are surfaced as structured rejections (`PredictionResult`) and are dropped from CLI output

## Model Artifacts

Bundled classpath resources:
- `models/Chronologer_20220601193755.torchscript.pt`
- `models/Chronologer_20220601193755.preprocessing.json`

Golden parity fixtures:
- `src/test/resources/data/golden/chronologer_parity_cases.golden.json`
- `src/test/resources/data/golden/parity_cases.tsv`

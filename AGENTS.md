# AGENTS.md
Agent guidance for this repository.

## Read-first
- Before making changes, read `README.md` at least once to understand architecture, module boundaries, and existing patterns.

## Repo layout (high level)
- Maven project: parent (root `pom.xml`) 

## Build and test commands
### Fast compile, no native rebuild, no tests
mvn -am -Dmaven.test.skip=true compile

### Full test suite (core + gui), skip native rebuild
mvn -am test

### Recommended: run selected tests
mvn -am -Dtest=<TestClass> test

## What “good autonomy” looks like
- Reuse existing code paths and utilities, especially for file parsing, object models, windowing, and output writers.
- Do not speculate about code you have not opened. If a conclusion depends on specific behavior, open and read the relevant files first (source, tests, build config, scripts).
-  Minimize blast radius:
  - Change the fewest files possible.
  - Keep diffs small and easy to review.
  - Prefer additive changes over rewiring existing behavior.
- Make success checkable:
  - Add or update a runnable test, or a deterministic golden-output check using fixtures under `src/test/resources`.
  - Ensure outputs are deterministic (ordering, rounding/epsilon, fixed seeds if randomness exists).

## Coding standards
- Preserve formatting and conventions already present in nearby files.
- Preserve and maintain comments, update comments when behavior changes and use the same comment style as nearby code.
- Do not introduce new build steps, toolchains, or dependencies unless explicitly requested.
- Always use import statements at the top of the file and never use inline fully qualified names (FQNs).

## Final output expectation
Always attempt to compile the Java code before finishing. When done, report:
- Commands you ran (including focused tests and/or full tests)
- What changed (files/modules)
- How correctness was verified (test names, fixtures, golden checks)
- When you make a claim about behavior, include the evidence path: file names and the exact methods or sections you relied on.

# **Zinc**

_Zinc_ is a proof-of-concept implementation of a novel SNARK system engineered
by [Nethermind Research](https://www.nethermind.io/nethermind-research), based
on [Zinc: Succinct Arguments with Small Arithmetization
Overheads from IOPs of Proximity to the Integers](https://eprint.iacr.org/2025/316) by
Albert Garreta, Hendrik Waldner, Katerina Hristova and Luca Dall'Ava.

The primary goal of _Zinc_ is to address a significant performance bottleneck in modern zero-knowledge proof systems,
specifically the high cost of arithmetization.

⚠️ **Disclaimer:** This is a proof-of-concept prototype.
This implementation is provided for research and evaluation purposes. It has not undergone a formal security audit or
comprehensive code review and is NOT ready for production use. Use at your own risk.

## **Repository Structure**

The Zinc codebase is organized into several distinct modules, each responsible for a core component of the cryptographic
system.

| Path          | Description                                                                                                                                                                                                                                                         |
|:--------------|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| src/zinc/     | **Top-Level Protocol Logic.** Implements the main ZincProver and ZincVerifier structs. Orchestrates the Spartan protocol by coordinating the sumcheck and PCS layers.                                                                                               |
| src/zip/      | **The Zip Polynomial Commitment Scheme (PCS).** Contains the commitment (pcs/commit.rs), opening (pcs/open\_z.rs), and verification (pcs/verify\_z.rs) logic. Implements the underlying RAA linear codes (code\_raa.rs) and Merkle tree commitments (pcs/utils.rs). |
| src/ccs/      | **Customizable Constraint Systems.** Defines the structure for CCS over integers (ccs\_z.rs) and their projection to finite fields (ccs\_f.rs). Includes utilities for matrix operations (utils.rs).                                                                |
| src/sumcheck/ | **The Sumcheck Protocol.** Provides the core interactive proof (prover.rs, verifier.rs) for verifying polynomial identities, which serves as the engine for the Spartan protocol.                                                                                   |
| src/field/    | **Custom Arithmetic Primitives.** Defines the core data structures for integer (int.rs) and finite field (field.rs) arithmetic, including the crucial RandomField enum that manages the integer-to-field transition.                                                |
| src/poly/     | **Polynomial Utilities.** Contains definitions and operations for dense and sparse multilinear polynomial extensions (MLEs), separated for integer (poly\_z) and field (poly\_f) coefficients.                                                                      |
| src/traits/   | **Core Abstractions.** Defines the essential Rust traits like Field, Integer, and FieldMap that enable the library's generic and modular design.                                                                                                                    |
| benches/      | **Performance Benchmarks.** Contains the Criterion.rs benchmark suite for evaluating the performance of core components like sumcheck, Spartan, and the Zip PCS.                                                                                                    |

## **Getting Started**

Follow these instructions to build, test and benchmark the Zinc library.

### **Prerequisites**

Ensure you have the following installed:

- **Rust**: Install Rust using [rustup](https://rustup.rs/).
- **git**: Ensure you have Git installed to clone the repository.

### **Building the Library**

To build the Zinc library, run the following commands:

```bash
git clone https://github.com/NethermindEth/zinc
cd zinc
cargo build --release
````

#### **Building the Library with parallelization enabled**

For improved performance on multicore systems, you can enable parallelization by using the `parallel` feature flag,
which uses Rayon for parallel processing.

```bash
cargo build --release --features "parallel"
```

### **Testing the Library**

```bash
cargo test --release
```

### **Running Benchmarks**

The repository includes a suite of benchmarks using Criterion.rs. You can run the full suite or a lighter, core subset (
used in CI for pull requests).

- Full local run (all benches):

```bash
cargo bench
```

- Core subset only (fast):

```bash
CORE_BENCHMARKS=1 cargo bench --bench zip_benches
```

Notes:

- The environment variable `CORE_BENCHMARKS=1` switches zip_benches.rs into a reduced set of cases for faster iteration.
- We pass `--output-format bencher` to produce output compatible with our CI reporting.

#### CI Benchmarking

Our GitHub Actions workflow runs benchmarks automatically:

- On pushes to `main`: the full benchmark suite (`cargo bench -- --nocapture --output-format bencher`), and saves
  history via benchmark-action.
- On pull requests: a lighter suite focused on the Zip PCS (
  `CORE_BENCHMARKS=1 cargo bench --bench zip_benches -- --nocapture --output-format bencher`), and posts results to the
  PR.

Continuous benchmark results can be viewed at [GitHub pages](https://NethermindEth.github.io/zinc/dev/bench)

## **Example**

For usage example you can refer to the example [r1cs.rs](examples/simple_r1cs.rs) file.

## Usage

Import the library:

```toml
[dependencies]
zinc = { git = "https://github.com/NethermindEth/zinc.git" }
```

## License

The crates in this repository are licensed under the following licence.

* MIT license ([LICENSE](./LICENSE))

## Would like to contribute?

see [Contributing](./CONTRIBUTING.md).

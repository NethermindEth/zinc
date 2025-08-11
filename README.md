# Zinc

A proof-of-concept implementation of the Zinc SNARK engineered by [Nethermind](https://nethermind.io) based on the work
[Zinc: Succinct Arguments with Small Arithmetization Overheads from IOPs of Proximity to the Integers](https://eprint.iacr.org/2025/316) by Albert Garreta, Hendrik Waldner, Katerina Hristova and Luca Dall'Ava.

**DISCLAIMER:** This is a proof-of-concept prototype, and in particular has not received careful code review. This implementation is provided "as is" and NOT ready for production use. Use at your own risk.

## Building

The [rust-toolchain](https://github.com/NethermindEth/zinc/blob/main/rust-toolchain) file pins the version of the Rust toolchain, which the Zinc library builds with, to the specific version `nightly-2024-11-05`.

One can install the `nightly-2024-11-05` toolchain by invoking:
```bash
rustup install nightly-2024-11-05
```

After that, use `cargo`, the standard Rust build tool, to build the library:

```bash
git clone https://github.com/NethermindEth/zinc.git
cargo build --release
```

## Usage
Import the library:
```toml
[dependencies]
zinc = { git = "https://github.com/NethermindEth/zinc.git" }
```

## License
The crates in this repository are licensed under the following licence.

* MIT license ([LICENSE-MIT](LICENSE-MIT))

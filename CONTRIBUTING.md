# Contributor's guide

## Commit signing

Enable [commit signing](https://docs.github.com/en/authentication/managing-commit-signature-verification/signing-commits)

```sh
git config commit.gpgsign true
```

## Prerequisites

* [Rust](https://www.rust-lang.org/tools/install)
* [cargo deny](https://github.com/EmbarkStudios/cargo-deny)
* [typos](https://github.com/crate-ci/typos?tab=readme-ov-file#install)
* [cargo sort](https://github.com/DevinR528/cargo-sort)

## Code quality assurance

Install a pre-push git hook:

```sh
git config core.hooksPath .githooks
```

## Running the Rust Documentation Locally
After cloning the repository, follow the instructions below to run the documentation locally:

```sh
cargo doc
```

Docs for `TODO(template) template_crate`:

```sh
RUSTDOCFLAGS="--html-in-header katex-header.html" cargo doc --no-deps -p template_crate --open
```

## pre-commit Framework (Alternative)

As an alternative to the pre-push git hook, you can use the [pre-commit](https://pre-commit.com/) framework to manage and run hooks.
pre-commit is Python-based and allows you to define a set of hooks that will run automatically on certain git events, such as pre-commit or pre-push.

In addition to the hooks defined in `.githooks`, pre-commit configuration also cleansup the code by removing trailing whitespace and ensuring that all files end with a newline.

To set up pre-commit, follow these steps:

1. Install pre-commit:
   ```sh
   pip install pre-commit
   ```

2. Install the pre-commit hooks:
   ```sh
    pre-commit install --hook-type pre-push
   ```

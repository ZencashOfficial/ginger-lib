#!/bin/bash
set -euo pipefail

cd /ginger-lib

cargo check

cargo check --all-features --tests

cargo +nightly-2021-04-25 check --all-features --tests --benches

cd /ginger-lib/algebra

cargo check

cargo check --all-features --tests

cargo +nightly-2021-04-25 check --all-features --tests --benches

cargo check --features "parallel"

cargo check --features "fft"

cargo check --features "n_fold"

cargo check --features "llvm_asm"

cargo check --features "bls12_377"

cargo check --features "bls12_381"

cargo check --features "edwards_bls12"

cargo check --features "edwards_sw6"

cargo check --features "jubjub"

cargo check --features "sw6"

cargo check --features "mnt4_753"

cargo check --features "mnt6_298"

cargo check --features "mnt6_753"

cargo check --features "bn_382"

cargo check --features "tweedle"

cargo check --features "full"

cd /ginger-lib/primitives

cargo check

cargo check --all-features --tests

cargo +nightly-2021-04-25 check --all-features --tests --benches

cargo check --features "llvm_asm"

cargo check --features "commitment"

cargo check --features "merkle_tree"

cargo check --features "prf"

cargo check --features "signature"

cargo check --features "vrf"

cargo check --features "mnt4_753"

cargo check --features "mnt6_753"

cargo check --features "bn_382"

cargo check --features "tweedle"

cd /ginger-lib/proof-systems

cargo check

cargo check --all-features --tests --examples

cargo +nightly-2021-04-25 check --all-features --tests --benches --examples

cargo check --features "llvm_asm"

cargo check --features "gm17"

cargo check --features "groth16"

cargo check --features "darlin"

cd /ginger-lib/r1cs/core

cargo check

cd /ginger-lib/gadgets/crypto

cargo check

cargo check --all-features --tests

cargo check --features "llvm_asm"

cargo check --features "commitment"

cargo check --features "merkle_tree"

cargo check --features "prf"

cargo check --features "signature"

cargo check --features "vrf"

cargo check --features "nizk"

cargo check --features "mnt4_753"

cargo check --features "mnt6_753"

cargo check --features "bn_382"

cargo check --features "tweedle"


cd /ginger-lib/std

cargo check

cargo check --all-features --tests

cargo check --features "llvm_asm"

cargo check --features "full"

cargo check --features "bls12_377"

cargo check --features "bn_382"

cargo check --features "edwards_bls12"

cargo check --features "edwards_sw6"

cargo check --features "jubjub"

cargo check --features "mnt4_753"

cargo check --features "mnt6_753"

cargo check --features "tweedle"

cd /ginger-lib/

cargo test --all-features

RUSTFLAGS="-C target-feature=+bmi2,+adx --emit=asm" cargo +nightly-2021-04-25 test --all-features

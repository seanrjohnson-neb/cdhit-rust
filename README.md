# cdhit-rust

A Rust port of [CD-HIT](https://github.com/weizhongli/cdhit) — biological
sequence clustering — built to compile and run on **Windows** and **WebAssembly**,
which are difficult with the original C++/Perl implementation.

The port reproduces the CD-HIT clustering algorithm with **bit-for-bit output
fidelity**: for the supported programs and options, the `.clstr` file and the
representative FASTA are byte-identical to the reference C++ (`make openmp=no`).

## Status

All six clustering programs are implemented and verified byte-identical to the
C++ across a large option matrix and randomized differential fuzzing:

| Program | Description |
|---|---|
| `cd-hit` | Protein clustering |
| `cd-hit-est` | Nucleotide clustering (incl. reverse-complement strand) |
| `cd-hit-2d` | Protein: compare a second database against a reference |
| `cd-hit-est-2d` | Nucleotide 2D comparison |
| `cd-hit-454` | 454 pyrosequencing duplicate detection |
| `cd-hit-div` | Split a database into N segments |

Not yet ported (out of scope for now): the `cd-hit-auxtools` package
(`cd-hit-dup`, `cd-hit-lap`, `read-linker`), the Perl post-processing scripts,
and `psi-cd-hit`.

## Workspace layout

- `crates/cdhit-core` — the engine (a reusable, platform-agnostic library).
- `crates/cdhit-cli` — the `cdhit` command-line tool (native).
- `crates/cdhit-wasm` — WebAssembly bindings (`wasm-bindgen`).

## Building and running

### Native CLI

```sh
cargo build --release
./target/release/cdhit cd-hit     -i input.faa -o output -c 0.9 -n 5
./target/release/cdhit cd-hit-est -i reads.fna -o output -c 0.95 -n 10
./target/release/cdhit cd-hit-2d  -i db1.faa -i2 db2.faa -o novel -c 0.9 -n 5
./target/release/cdhit cd-hit-div -i input.faa -o part -div 4
```

Output is written to `<output>` (representatives) and `<output>.clstr`
(cluster membership). Gzipped input (`.gz`) is supported. The historical
program names also work if the binary is symlinked to them.

### WebAssembly

```sh
cargo build -p cdhit-wasm --target wasm32-unknown-unknown --release
wasm-bindgen --target web --out-dir pkg \
    target/wasm32-unknown-unknown/release/cdhit_wasm.wasm
```

```js
import init, { cluster, cluster_est } from "./pkg/cdhit_wasm.js";
await init();
const res = cluster(fastaText, "-c 0.9 -n 5");
console.log(res.num_clusters, res.clstr, res.rep_fasta);
```

### As a library

```toml
[dependencies]
cdhit-core = { path = "crates/cdhit-core" }
```

```rust
let out = cdhit_core::cluster_1d(fasta_bytes, &args, /* est = */ false)?;
```

## Fidelity and testing

Fidelity is enforced by golden tests (`crates/cdhit-core/tests/`) that diff the
Rust output against captured C++ reference output, plus a randomized
differential fuzz harness that generates protein/nucleotide inputs and asserts
Rust == C++ across many option combinations.

The port deliberately mirrors several C++ quirks required for exact output
(e.g. the word-count power-table integer overflow, the `-g 1` cutoff-update
behaviour, and float-promoted-to-double threshold comparisons).

## Portability notes / differences from upstream

- **Threading** (`-T`, optional `rayon` feature): OpenMP is replaced by rayon.
  - **2D clustering** parallelizes cleanly and is **output-neutral** (db1 is a
    fixed reference, so every db2 sequence is screened independently). Scales
    well (≈6× on 8 cores in tests).
  - **1D clustering** (`cd-hit`, `cd-hit-est`, `cd-hit-454`) overlaps serial
    representative-building of the current batch with parallel screening of the
    following sequences against the previous batch's representatives (CD-HIT's
    threaded strategy). Under `-T > 1` the cluster *assignments* are **not
    byte-identical to serial** — the greedy first/best-match is order-sensitive
    and cannot be reproduced exactly in parallel (CD-HIT's own threaded output
    likewise differs from its serial output). The clustering is of equivalent
    quality, and — unlike CD-HIT, whose OpenMP output is non-deterministic
    run-to-run — this implementation is **deterministic** (identical output for
    any `-T`). Speedup is modest and workload-dependent (best on large,
    high-redundancy datasets; for small inputs the serial representative-building
    and parallel overhead dominate, so it may not help).
  - **Use `-T 1` (the default) for output byte-identical to the reference C++
    serial algorithm.** WASM builds are single-threaded (`--no-default-features`).
- **`-B` disk swap**: dropped. The database is held in memory (this removes the
  non-portable temp-file path and is required for WASM). Very large inputs are
  bounded by available memory.
- **zlib**: replaced by the `flate2` crate (`gzip` feature, on by default for
  native, off for WASM).

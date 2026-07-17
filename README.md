# cdhit-rust

A Rust port of [CD-HIT](https://github.com/weizhongli/cdhit) — biological
sequence clustering — built to compile and run on **Windows** and **WebAssembly**,
which are difficult with the original C++/Perl implementation.

The port reproduces the CD-HIT clustering algorithm with **bit-for-bit output
fidelity**: for the supported programs and options, the `.clstr` file and the
representative FASTA are byte-identical to the reference C++ (`make openmp=no`).

## Status

All six clustering programs plus the three `cd-hit-auxtools` programs are
implemented and verified byte-identical to the C++ across an option matrix and
differential testing:

| Program | Description |
|---|---|
| `cd-hit` | Protein clustering |
| `cd-hit-est` | Nucleotide clustering (incl. reverse-complement strand) |
| `cd-hit-2d` | Protein: compare a second database against a reference |
| `cd-hit-est-2d` | Nucleotide 2D comparison |
| `cd-hit-454` | 454 pyrosequencing duplicate detection |
| `cd-hit-div` | Split a database into N segments |
| `cd-hit-dup` | Duplicate / near-duplicate read detection (single-end + paired-end, chimera filtering) |
| `cd-hit-lap` | Cluster reads that overlap end-to-end |
| `read-linker` | Join paired-end reads by 3'/5' overlap |

Many of the Perl `.clstr` post-processing scripts are also ported, as `cdhit
<name>` subcommands, each verified byte-identical to the reference Perl:

| Subcommand | Perl script | Description |
|---|---|---|
| `clstr_sort_by` | `clstr_sort_by.pl` | Sort clusters by size or representative length |
| `clstr_size_stat` | `clstr_size_stat.pl` | Distribution of cluster sizes |
| `clstr_size_histogram` | `clstr_size_histogram.pl` | Binned histogram of cluster sizes |
| `clstr2txt` | `clstr2txt.pl` | Tabular per-sequence view |
| `clstr_renumber` | `clstr_renumber.pl` | Renumber clusters and members |
| `clstr_select` | `clstr_select.pl` | Select clusters by size range |
| `clstr_cut` | `clstr_cut.pl` | Keep the top N members per cluster |
| `clstr_rep` | `clstr_rep.pl` | List representative of each cluster |
| `clstr2tree` | `clstr2tree.pl` | Newick-like tree |
| `cd-hit-clstr_2_blm8` | `cd-hit-clstr_2_blm8.pl` | Convert to BLAST tabular (m8) |
| `clstr_reduce` | `clstr_reduce.pl` | Sub-sample clusters by size segment |
| `clstr_rev` | `clstr_rev.pl` | Flatten a two-level hierarchical clustering |
| `clstr_merge` | `clstr_merge.pl` | Merge divided cluster files into a master |
| `clstr_merge_noorder` | `clstr_merge_noorder.pl` | Merge divided cluster files whose cluster order need not match |
| `clstr_reps_faa_rev` | `clstr_reps_faa_rev.pl` | Keep the top N sequences per cluster from a FASTA |
| `clstr_select_rep` | `clstr_select_rep.pl` | Print representatives of clusters within a size range |
| `clstr_sort_prot_by` | `clstr_sort_prot_by.pl` | Sort members within each cluster by length or id |
| `clstr_quality_eval_by_link` | `clstr_quality_eval_by_link.pl` | Sensitivity/specificity vs a benchmark (independent links) |
| `plot_len1` | `plot_len1.pl` | Tabular sequence/cluster counts by size × representative-length bins |
| `clstr_sql_tbl` | `clstr_sql_tbl.pl` | Build/extend a hierarchical cluster table (`id len cid rep`, two columns per level) |
| `clstr_sql_tbl_sort` | `clstr_sql_tbl_sort.pl` | Stable-sort a cluster SQL table by hierarchy columns |
| `make_multi_seq` | `make_multi_seq.pl` | Write one FASTA file per cluster above a size cutoff |
| `cd-hit-dup-PE-out` | `cd-hit-dup-PE-out.pl` | Export representative paired-end reads after `cd-hit-dup` |

Not ported: `psi-cd-hit` (orchestrates external BLAST/PSI-BLAST — not
WASM-relevant) and the `cd-hit-para.pl` / `cd-hit-2d-para.pl` grid wrappers
(superseded by the `-T` rayon parallelism). Some Perl scripts are intentionally
**not** ported because they cannot be made bit-for-bit reproducible or need
non-portable dependencies: `FET.pl` (external CPAN `Text::NSP` module + Perl
`Storable` + hash-ordered output), `clstr_quality_eval.pl` and `clstr2xml.pl`
(emit in Perl hash-iteration order), `clstr_list.pl` (Perl `Storable` binary
format), and `plot_2d.pl` (ImageMagick / non-deterministic drawing). See
[Remaining Perl scripts](#remaining-perl-scripts) for the full breakdown.

The `cd-hit-auxtools` programs work on FASTA/FASTQ, share the `cdhit-core`
crate, and — like the clustering programs — build for WebAssembly.

> **Note on `cd-hit-dup -f`/`-s` (chimera filtering):** the reference C++
> `DetectChimeric` contains out-of-bounds reads (undefined behaviour). The port
> reproduces the well-defined behaviour bit-for-bit and guards the UB, so it
> matches the reference on well-formed inputs; the de-duplication path is fully
> verified.

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

# auxtools (own flags; run with no args for help)
./target/release/cdhit cd-hit-dup  -i reads.fq -o uniq -e 1
./target/release/cdhit cd-hit-dup  -i R1.fq -i2 R2.fq -o uniq -o2 uniq.R2
./target/release/cdhit cd-hit-lap  -i reads.fa -o out -m 20
./target/release/cdhit read-linker -1 R1.fq -2 R2.fq -o contigs.fq -l 10 -e 1

# .clstr post-processing (read a file arg or stdin; write to stdout)
./target/release/cdhit clstr_sort_by len out.clstr > sorted.clstr
./target/release/cdhit clstr_size_stat out.clstr
./target/release/cdhit clstr2txt out.clstr > table.tsv
./target/release/cdhit clstr_select_rep 2 100 out.clstr > reps.txt
./target/release/cdhit plot_len1 out.clstr 1,2-5,6-up 1-100,101-200,201-up

# cd-hit-dup-PE-out uses getopts-style flags and writes two files
./target/release/cdhit cd-hit-dup-PE-out -i R1.fq -j R2.fq -c uniq.clstr \
    -o rep.R1.fq -p rep.R2.fq

# a few .clstr helpers write files rather than stdout (native only)
./target/release/cdhit make_multi_seq db.faa out.clstr multi-seq 20  # one FASTA/cluster
./target/release/cdhit clstr_sql_tbl db90.clstr tbl.txt              # create the table
./target/release/cdhit clstr_sql_tbl db60.clstr tbl.txt              # append the next level
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

## Reproduced upstream quirks and bugs

To stay byte-identical, the port faithfully reproduces a number of quirks —
including outright bugs and undefined behaviour — in the original C++/Perl.
These are matched deliberately; "fixing" them would break output parity.

**Clustering engine (`cdhit-common`)**
- Word-count power-table **integer overflow** is reproduced exactly.
- `-g 1` cutoff-update behaviour and float→double threshold promotion.

**cd-hit-auxtools**
- **`bioSequence` parser**: residues of *every record except the file's last*
  are upper-cased; the C++ omits `ToUpper()` on the final flush, so the last
  record keeps its original case.
- **`read-linker`**: the emitted `mismatch_no=` uses the overlap-check error
  count (which also counts `N`/`N` pairs), not the number of real substitutions
  used to enumerate contigs. Reads are also processed in fixed batches of
  10000, so a batch's tail is dropped when the two mates' files are of unequal
  length.
- **`cd-hit-dup -f`/`-s` (chimera filtering)**: the reference `DetectChimeric`
  reads out of bounds (e.g. it always inspects the top *two* candidate offsets
  even when only one exists) — undefined behaviour. The port reproduces the
  well-defined behaviour bit-for-bit and *guards* the OOB accesses, so it may
  diverge only on the degenerate cases the C++ leaves undefined.
- **`cd-hit-dup -u`/paired-end**: full-length reads are restored for output by
  re-reading the input and copying **by index** into the (possibly length-sorted)
  working list — a latent mismatch when sorting occurred. Reproduced as-is.

**Perl `.clstr` post-processing scripts**
- **`cd-hit-clstr_2_blm8`**: protein (`aa`) clusters and strand-only nucleotide
  lines have no alignment coordinates, so the script emits deterministic
  "garbage" — negative alignment lengths and bit scores, and a raw string
  (`97.50%` or `+`) in the `q_b` column — via Perl's string→number coercion.
  Matched exactly, including Perl's `%.15g` number formatting.
- **`clstr2tree`**: branch length `1-fr` is emitted with Perl's `%.15g`
  stringification (so `1-0.8` prints as `0.2`, not `0.199999…`).
- **`clstr_rep`**: only recognises protein (`aa`) representative lines; a
  non-`aa` representative is treated as a format error, as upstream.

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
- **`-B` disk swap** (and its `-tmp` companion): accepted for command-line
  compatibility but a no-op — exactly as in current upstream CD-HIT, which itself
  disabled disk-swap (its read paths are annotated "disable swap option"), so
  `-B 1` produces output identical to `-B 0`. The database is held in memory
  (this removes the non-portable temp-file path and is required for WASM). Very
  large inputs are bounded by available memory.
- **zlib**: replaced by the `flate2` crate (`gzip` feature, on by default for
  native, off for WASM).

## Remaining Perl scripts

The high-priority "best candidates" have all been ported (see the subcommand
table above), each verified byte-identical to the reference Perl:
`clstr_select_rep`, `clstr_sort_prot_by`, `clstr_merge_noorder`,
`cd-hit-dup-PE-out`, `clstr_quality_eval_by_link`, `plot_len1`, and
`clstr_sql_tbl_sort`.

Two more native-only helpers have since been ported (also byte-identical):
`make_multi_seq` (one FASTA per cluster) and `clstr_sql_tbl` (the table producer
that `clstr_sql_tbl_sort` consumes).

The following scripts remain unported. The first group is portable but each
carries a caveat; the second group is out of scope for the Windows/WASM story.

**Useful, but with caveats**

| Script | Feasibility | Notes |
|---|---|---|
| `cd-hit-div.pl` | High | Round-robin split without sorting; differs from the ported `cd-hit-div`. Mainly historical wrapper compatibility. |
| `clstr_quality_eval.pl` | Medium | Aggregate metrics are easy, but the pair-list sections iterate Perl hashes, so byte-identical output is not stable (only the link-based variant is ported). |
| `clstr2xml.pl` | Medium | Output order depends on Perl hash iteration; a deterministic XML writer would be better than byte-parity. |

**Out of scope**

| Script | Notes |
|---|---|
| `FET.pl` | Depends on `Text::NSP`, Perl `Storable`, and hash-ordered output. Fisher exact could be reimplemented, but not as byte-identical Perl compatibility. |
| `clstr_list.pl` / `clstr_list_sort.pl` | Revolve around Perl `Storable` binary structures — not portable or useful for WASM. |
| `plot_2d.pl` | Depends on ImageMagick and random drawing; output is not deterministic. |
| `cd-hit-para.pl` / `cd-hit-2d-para.pl` | Orchestrate distributed jobs via ssh/qsub; superseded by native `-T` rayon support. |
| `psi-cd-hit` | Requires BLAST/PSI-BLAST, `makeblastdb`, and local/qsub orchestration — a separate native-only project if desired. |
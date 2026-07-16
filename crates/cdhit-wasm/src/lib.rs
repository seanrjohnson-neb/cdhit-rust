//! WebAssembly bindings for the CD-HIT engine.
//!
//! Single-threaded, no filesystem, no gzip — the `cdhit-core` engine is fully
//! portable to `wasm32`. These bindings take an in-memory FASTA string and the
//! usual cd-hit argument list, and return the representative FASTA and `.clstr`
//! text as JS-visible fields.
//!
//! Example (JS):
//! ```js
//! import init, { cluster } from "./cdhit_wasm.js";
//! await init();
//! const res = cluster(fastaText, "-c 0.9 -n 5");   // protein
//! console.log(res.clstr, res.rep_fasta, res.num_clusters);
//! ```

use cdhit_core::{
    cd_hit_dup, cd_hit_lap, cluster_1d_program, cluster_2d, read_linker, DupParams, LapParams,
    LinkerParams, Program,
};
use wasm_bindgen::prelude::*;

/// Result of a clustering run, exposed to JavaScript.
#[wasm_bindgen(getter_with_clone)]
pub struct ClusterResult {
    /// Representative sequences in FASTA format (the `-o` output).
    pub rep_fasta: String,
    /// Cluster membership file (`.clstr`).
    pub clstr: String,
    /// Number of clusters.
    pub num_clusters: usize,
}

/// Split a whitespace-separated argument string into cd-hit flag/value tokens.
/// `-i`/`-o` are synthesised (the WASM API works on in-memory data) so callers
/// need only pass the algorithmic options.
fn build_args(extra: &str) -> Vec<String> {
    let mut args = vec![
        "-i".to_string(),
        "in".to_string(),
        "-o".to_string(),
        "out".to_string(),
    ];
    args.extend(extra.split_whitespace().map(|s| s.to_string()));
    args
}

fn run_1d(fasta: &str, opts: &str, program: Program) -> Result<ClusterResult, String> {
    let args = build_args(opts);
    let out = cluster_1d_program(fasta.as_bytes(), &args, program).map_err(|e| e.to_string())?;
    Ok(ClusterResult {
        rep_fasta: String::from_utf8_lossy(&out.rep_fasta).into_owned(),
        clstr: out.clstr,
        num_clusters: out.num_clusters,
    })
}

/// Cluster a protein FASTA (`cd-hit`). `opts` is a space-separated option
/// string, e.g. `"-c 0.9 -n 5"`.
#[wasm_bindgen]
pub fn cluster(fasta: &str, opts: &str) -> Result<ClusterResult, String> {
    run_1d(fasta, opts, Program::CdHit)
}

/// Cluster a nucleotide FASTA (`cd-hit-est`).
#[wasm_bindgen]
pub fn cluster_est(fasta: &str, opts: &str) -> Result<ClusterResult, String> {
    run_1d(fasta, opts, Program::CdHitEst)
}

/// 454 read duplicate detection (`cd-hit-454`).
#[wasm_bindgen]
pub fn cluster_454(fasta: &str, opts: &str) -> Result<ClusterResult, String> {
    run_1d(fasta, opts, Program::CdHit454)
}

/// 2D comparison (`cd-hit-2d` / `cd-hit-est-2d`): cluster `db2` against the
/// reference `db1`. Set `est = true` for nucleotides.
#[wasm_bindgen]
pub fn cluster_2d_wasm(
    db1: &str,
    db2: &str,
    opts: &str,
    est: bool,
) -> Result<ClusterResult, String> {
    let mut args = build_args(opts);
    args.push("-i2".to_string());
    args.push("in2".to_string());
    let out = cluster_2d(db1.as_bytes(), db2.as_bytes(), &args, est).map_err(|e| e.to_string())?;
    Ok(ClusterResult {
        rep_fasta: String::from_utf8_lossy(&out.rep_fasta).into_owned(),
        clstr: out.clstr,
        num_clusters: out.num_clusters,
    })
}

// ---------------------------------------------------------------------------
// cd-hit-auxtools bindings
// ---------------------------------------------------------------------------

/// Result of a `read-linker` run.
#[wasm_bindgen(getter_with_clone)]
pub struct LinkerResult {
    /// Linked contigs (FASTA or FASTQ).
    pub output: String,
    /// Progress / summary text (what the CLI prints to stdout).
    pub log: String,
}

/// Join paired-end reads by overlap (`read-linker`). `min` = minimum overlap
/// (`-l`), `error` = max mismatches (`-e`).
#[wasm_bindgen]
pub fn link_reads(first: &str, second: &str, min: i32, error: i32) -> LinkerResult {
    let out = read_linker(
        first.as_bytes(),
        second.as_bytes(),
        &LinkerParams {
            min,
            error,
            maxlen: 0,
        },
    );
    LinkerResult {
        output: String::from_utf8_lossy(&out.output).into_owned(),
        log: out.log,
    }
}

/// Result of a `cd-hit-lap` run.
#[wasm_bindgen(getter_with_clone)]
pub struct LapResult {
    /// Representative sequences (FASTA/FASTQ).
    pub rep: String,
    /// Cluster membership file (`.clstr`).
    pub clstr: String,
    /// Log text.
    pub log: String,
}

/// Cluster overlapping reads (`cd-hit-lap`). `minlen` = `-m`, `minper` = `-p`,
/// `deslen` = `-d`.
#[wasm_bindgen]
pub fn lap(fastaq: &str, minlen: i32, minper: f32, deslen: i32) -> LapResult {
    let out = cd_hit_lap(
        fastaq.as_bytes(),
        &LapParams {
            minlen,
            minper,
            deslen,
            seed: 0,
        },
    );
    LapResult {
        rep: String::from_utf8_lossy(&out.rep).into_owned(),
        clstr: String::from_utf8_lossy(&out.clstr).into_owned(),
        log: out.log,
    }
}

/// Result of a `cd-hit-dup` run (single-end).
#[wasm_bindgen(getter_with_clone)]
pub struct DupResult {
    /// Representative reads (FASTA/FASTQ).
    pub reps: String,
    /// Cluster membership file (`.clstr`).
    pub clstr: String,
    /// Chimeric-cluster file (`2.clstr`).
    pub clstr2: String,
    /// Log text.
    pub log: String,
}

/// Detect duplicate / near-duplicate reads (`cd-hit-dup`, single-end).
/// `errors` = `-e`, `match_length` = `-m`, `uselen` = `-u`, `nochimeric` = `-f`.
#[wasm_bindgen]
pub fn dedup(
    fastaq: &str,
    errors: i32,
    match_length: bool,
    uselen: i32,
    nochimeric: bool,
) -> Result<DupResult, String> {
    let out = cd_hit_dup(
        fastaq.as_bytes(),
        &DupParams {
            input_name: "in",
            input2: None,
            match_length,
            abundance: -1,
            deslen: 0,
            uselen,
            errors,
            errors2: errors as f32,
            nochimeric,
            shared: 30,
            abratio: 1.0,
            percent: 1.0,
        },
    )?;
    Ok(DupResult {
        reps: String::from_utf8_lossy(&out.reps_r1).into_owned(),
        clstr: String::from_utf8_lossy(&out.clstr).into_owned(),
        clstr2: String::from_utf8_lossy(&out.clstr2).into_owned(),
        log: out.log,
    })
}

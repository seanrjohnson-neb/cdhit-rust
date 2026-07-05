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

use cdhit_core::{cluster_1d_program, cluster_2d, Program};
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

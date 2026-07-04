//! Phase 9 milestone: nucleotide (EST) clustering must byte-match cd-hit-est,
//! including reverse-complement strand detection.

use cdhit_core::cluster_1d;

fn args(v: &[&str]) -> Vec<String> {
    v.iter().map(|s| s.to_string()).collect()
}

#[test]
fn dna_c09_n8_clstr_matches_cpp() {
    let fa = include_bytes!("data/dna.fasta");
    let expected = include_str!("data/dna_c0.9_n8.clstr");
    let out = cluster_1d(fa, &args(&["-i", "x", "-o", "y", "-c", "0.9", "-n", "8"]), true).unwrap();
    assert_eq!(out.clstr, expected);
    // The reverse-complement near-duplicate must be recorded on the minus strand.
    assert!(out.clstr.contains("-/"));
}

#[test]
fn dna_c09_n8_rep_fasta_matches_cpp() {
    let fa = include_bytes!("data/dna.fasta");
    let expected = include_bytes!("data/dna_c0.9_n8.fasta");
    let out = cluster_1d(fa, &args(&["-i", "x", "-o", "y", "-c", "0.9", "-n", "8"]), true).unwrap();
    assert_eq!(out.rep_fasta, expected.to_vec());
}

#[test]
fn dna_c09_n8_forward_only_matches_cpp() {
    let fa = include_bytes!("data/dna.fasta");
    let expected = include_str!("data/dna_c0.9_n8_r0.clstr");
    let out = cluster_1d(
        fa,
        &args(&["-i", "x", "-o", "y", "-c", "0.9", "-n", "8", "-r", "0"]),
        true,
    )
    .unwrap();
    assert_eq!(out.clstr, expected);
}

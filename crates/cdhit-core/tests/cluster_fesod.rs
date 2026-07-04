//! Phase 8 milestone: end-to-end protein clustering must byte-match the C++
//! `.clstr` and representative FASTA on the FeSOD fixture.

use cdhit_core::cluster_1d;

fn args(v: &[&str]) -> Vec<String> {
    v.iter().map(|s| s.to_string()).collect()
}

#[test]
fn fesod_c05_n3_clstr_matches_cpp() {
    let fa = include_bytes!("data/fesod.fasta");
    let expected = include_str!("data/fesod_c0.5_n3.clstr");
    let out = cluster_1d(
        fa,
        &args(&["-i", "x", "-o", "y", "-c", "0.5", "-n", "3"]),
        false,
    )
    .unwrap();
    assert_eq!(out.clstr, expected);
}

#[test]
fn fesod_c05_n3_rep_fasta_matches_cpp() {
    let fa = include_bytes!("data/fesod.fasta");
    let expected = include_bytes!("data/fesod_c0.5_n3.fasta");
    let out = cluster_1d(
        fa,
        &args(&["-i", "x", "-o", "y", "-c", "0.5", "-n", "3"]),
        false,
    )
    .unwrap();
    assert_eq!(out.rep_fasta, expected.to_vec());
}

#[test]
fn fesod_c09_n5_clstr_matches_cpp() {
    let fa = include_bytes!("data/fesod.fasta");
    let expected = include_str!("data/fesod_c0.9_n5.clstr");
    let out = cluster_1d(
        fa,
        &args(&["-i", "x", "-o", "y", "-c", "0.9", "-n", "5"]),
        false,
    )
    .unwrap();
    assert_eq!(out.clstr, expected);
}

#[test]
fn fesod_g1_best_mode_matches_cpp() {
    let fa = include_bytes!("data/fesod.fasta");
    let expected = include_str!("data/fesod_c0.5_n3_g1.clstr");
    let out = cluster_1d(
        fa,
        &args(&["-i", "x", "-o", "y", "-c", "0.5", "-n", "3", "-g", "1"]),
        false,
    )
    .unwrap();
    assert_eq!(out.clstr, expected);
}

#[test]
fn fesod_sc_sort_by_size_matches_cpp() {
    let fa = include_bytes!("data/fesod.fasta");
    let expected = include_str!("data/fesod_c0.5_n3_sc.clstr");
    let out = cluster_1d(
        fa,
        &args(&["-i", "x", "-o", "y", "-c", "0.5", "-n", "3", "-sc", "1"]),
        false,
    )
    .unwrap();
    assert_eq!(out.clstr, expected);
}

#[test]
fn fesod_local_identity_coverage_matches_cpp() {
    let fa = include_bytes!("data/fesod.fasta");
    let expected = include_str!("data/fesod_c0.7_G0_aS0.8.clstr");
    let out = cluster_1d(
        fa,
        &args(&["-i", "x", "-o", "y", "-c", "0.7", "-n", "5", "-G", "0", "-aS", "0.8"]),
        false,
    )
    .unwrap();
    assert_eq!(out.clstr, expected);
}

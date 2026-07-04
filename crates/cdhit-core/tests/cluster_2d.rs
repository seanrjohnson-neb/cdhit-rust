//! Phase 9 milestone: 2D clustering (`cd-hit-2d`) must byte-match the C++ for
//! both the `.clstr` and the representative FASTA of novel db2 sequences.

use cdhit_core::cluster_2d;

fn args(v: &[&str]) -> Vec<String> {
    v.iter().map(|s| s.to_string()).collect()
}

#[test]
fn fesod_2d_c06_n4_matches_cpp() {
    let db1 = include_bytes!("data/fesod_db1.fasta");
    let db2 = include_bytes!("data/fesod_db2.fasta");
    let clstr = include_str!("data/fesod_2d_c0.6_n4.clstr");
    let fasta = include_bytes!("data/fesod_2d_c0.6_n4.fasta");
    let out = cluster_2d(
        db1,
        db2,
        &args(&["-i", "a", "-i2", "b", "-o", "y", "-c", "0.6", "-n", "4"]),
        false,
    )
    .unwrap();
    assert_eq!(out.clstr, clstr);
    assert_eq!(out.rep_fasta, fasta.to_vec());
}

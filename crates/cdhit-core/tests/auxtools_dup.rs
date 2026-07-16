//! Golden tests for `cd-hit-dup`: assert the Rust port's representative FASTQ,
//! `.clstr` and `2.clstr` (chimeric) output are byte-identical to the reference
//! C++ (`cd-hit-auxtools/cd-hit-dup`). Fixtures contain exact duplicates,
//! 1-mismatch near-duplicates, and (for the chimera case) two high-abundance
//! parents plus a lower-abundance chimeric read.

use cdhit_core::{cd_hit_dup, DupParams};

const DUP_IN: &[u8] = include_bytes!("data/auxtools/dup_in.fq");
const CHIM_IN: &[u8] = include_bytes!("data/auxtools/dup_chim_in.fq");

fn params<'a>(name: &'a str) -> DupParams<'a> {
    DupParams {
        input_name: name,
        input2: None,
        match_length: true,
        abundance: -1,
        deslen: 0,
        uselen: 0,
        errors: 0,
        errors2: 0.0,
        nochimeric: false,
        shared: 30,
        abratio: 1.0,
        percent: 1.0,
    }
}

#[test]
fn dup_exact_matches_reference() {
    let out = cd_hit_dup(DUP_IN, &params("in")).unwrap();
    assert_eq!(out.reps_r1, include_bytes!("data/auxtools/dup_exact.rep"));
    assert_eq!(out.clstr, include_bytes!("data/auxtools/dup_exact.clstr"));
}

#[test]
fn dup_near_dup_e1_matches_reference() {
    let mut p = params("in");
    p.errors = 1;
    let out = cd_hit_dup(DUP_IN, &p).unwrap();
    assert_eq!(out.reps_r1, include_bytes!("data/auxtools/dup_e1.rep"));
    assert_eq!(out.clstr, include_bytes!("data/auxtools/dup_e1.clstr"));
}

#[test]
fn dup_chimera_matches_reference() {
    let mut p = params("in");
    p.nochimeric = true; // -f true
    let out = cd_hit_dup(CHIM_IN, &p).unwrap();
    assert_eq!(out.reps_r1, include_bytes!("data/auxtools/dup_chim.rep"));
    assert_eq!(out.clstr, include_bytes!("data/auxtools/dup_chim.clstr"));
    assert_eq!(
        out.clstr2,
        include_bytes!("data/auxtools/dup_chim2.clstr"),
        "chimeric-cluster file must be byte-identical to reference"
    );
}

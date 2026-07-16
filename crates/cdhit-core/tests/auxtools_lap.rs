//! Golden test for `cd-hit-lap`: assert the Rust port's representative FASTA
//! and `.clstr` output are byte-identical to the reference C++
//! (`cd-hit-auxtools/cd-hit-lap -m 20`). The fixture mixes overlapping
//! fragments of a template (prefix, suffix and reverse-complement overlaps),
//! an exact duplicate, and unrelated reads.

use cdhit_core::{cd_hit_lap, LapParams};

const INPUT: &[u8] = include_bytes!("data/auxtools/lap_in.fa");
const EXPECTED_REP: &[u8] = include_bytes!("data/auxtools/lap_m20.rep");
const EXPECTED_CLSTR: &[u8] = include_bytes!("data/auxtools/lap_m20.clstr");

#[test]
fn cd_hit_lap_matches_reference() {
    let out = cd_hit_lap(
        INPUT,
        &LapParams {
            minlen: 20,
            minper: 0.0,
            deslen: 0,
            seed: 0,
        },
    );
    assert_eq!(
        out.rep, EXPECTED_REP,
        "lap representative FASTA must be byte-identical to reference C++"
    );
    assert_eq!(
        out.clstr, EXPECTED_CLSTR,
        "lap .clstr must be byte-identical to reference C++"
    );
}

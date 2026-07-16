//! Golden test for `read-linker`: assert the Rust port's contig output is
//! byte-identical to the reference C++ (`cd-hit-auxtools/read-linker`) captured
//! output. Fixtures are 6 paired-end reads with a clean 20 bp overlap.

use cdhit_core::{read_linker, LinkerParams};

const R1: &[u8] = include_bytes!("data/auxtools/linker_R1.fq");
const R2: &[u8] = include_bytes!("data/auxtools/linker_R2.fq");
const EXPECTED_OUT: &[u8] = include_bytes!("data/auxtools/linker_l10_e1.out");

#[test]
fn read_linker_matches_reference() {
    let out = read_linker(
        R1,
        R2,
        &LinkerParams {
            min: 10,
            error: 1,
            maxlen: 0,
        },
    );
    assert_eq!(
        out.output,
        EXPECTED_OUT,
        "linker contig output must be byte-identical to reference C++"
    );
    // Summary log the C++ prints to stdout.
    let expected_log = "handled:         6\n\
        Total input pairs of read: 6\n\
        Total pairs of read used: 6\n\
        Total contigs: 6\n\
        \x200 mismatch:         6 pairs of reads =>         6 contigs\n\
        \x201 mismatch:         0 pairs of reads =>         0 contigs\n\
        Overlap range 20-20\n";
    assert_eq!(out.log, expected_log, "linker stdout log must match reference");
}

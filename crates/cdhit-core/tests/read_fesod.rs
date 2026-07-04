//! Phase 2 fidelity check: parsing the real FeSOD fixture must yield the same
//! (size, identifier) pairs the C++ reports in its `.clstr` output.

use cdhit_core::io::read_database;
use cdhit_core::Options;

#[test]
fn fesod_sizes_and_identifiers_match_cpp() {
    let fa = include_bytes!("data/fesod.fasta");
    let expected = include_str!("data/fesod_pairs.txt");

    let mut opt = Options::default();
    opt.des_len = 20; // cd-hit default
    let db = read_database(fa, &opt);
    assert_eq!(db.sequences.len(), 20);

    // Build "size\tidentifier(without '>')" pairs, sorted, to compare.
    let mut got: Vec<String> = db
        .sequences
        .iter()
        .map(|s| {
            let ident = String::from_utf8_lossy(&s.identifier[1..]); // drop '>'
            format!("{}\t{}", s.size, ident)
        })
        .collect();
    got.sort();

    let mut want: Vec<String> = expected.lines().map(|l| l.to_string()).collect();
    want.sort();

    assert_eq!(got, want);
}

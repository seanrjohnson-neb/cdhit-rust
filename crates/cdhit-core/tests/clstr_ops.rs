//! Golden tests for the `clstr_*` post-processing subcommands: assert each
//! matches the reference Perl script's stdout byte-for-byte. The fixture
//! `sample.clstr` is real `cd-hit-est` output with varied cluster sizes and
//! identities.

use cdhit_core::clstr::ops;

const SAMPLE: &[u8] = include_bytes!("data/clstr/sample.clstr");

#[test]
fn clstr_sort_by_no_matches_perl() {
    assert_eq!(ops::sort_by(SAMPLE, "no"), include_bytes!("data/clstr/sort_no.out"));
}

#[test]
fn clstr_sort_by_len_matches_perl() {
    assert_eq!(
        ops::sort_by(SAMPLE, "len"),
        include_bytes!("data/clstr/sort_len.out")
    );
}

#[test]
fn clstr_size_stat_matches_perl() {
    assert_eq!(ops::size_stat(SAMPLE), include_bytes!("data/clstr/size_stat.out"));
}

#[test]
fn clstr_size_histogram_matches_perl() {
    assert_eq!(
        ops::size_histogram(SAMPLE, 2),
        include_bytes!("data/clstr/size_hist_bin2.out")
    );
}

#[test]
fn clstr2txt_matches_perl() {
    assert_eq!(ops::to_txt(SAMPLE), include_bytes!("data/clstr/txt.out"));
}

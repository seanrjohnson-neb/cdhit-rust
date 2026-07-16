//! Golden tests for the `clstr_*` post-processing subcommands: assert each
//! matches the reference Perl script's stdout byte-for-byte. The fixture
//! `sample.clstr` is real `cd-hit-est` output with varied cluster sizes and
//! identities.

use cdhit_core::clstr::ops;

const SAMPLE: &[u8] = include_bytes!("data/clstr/sample.clstr");
const PROT: &[u8] = include_bytes!("data/clstr/prot.clstr");

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

#[test]
fn clstr_renumber_matches_perl() {
    assert_eq!(ops::renumber(SAMPLE), include_bytes!("data/clstr/renumber.out"));
}

#[test]
fn clstr_select_matches_perl() {
    assert_eq!(ops::select(SAMPLE, 2, 4), include_bytes!("data/clstr/select_2_4.out"));
}

#[test]
fn clstr_cut_matches_perl() {
    assert_eq!(ops::cut(SAMPLE, 3), include_bytes!("data/clstr/cut_3.out"));
}

#[test]
fn clstr_rep_matches_perl() {
    assert_eq!(ops::rep(PROT).unwrap(), include_bytes!("data/clstr/rep.out"));
}

#[test]
fn clstr2tree_matches_perl() {
    assert_eq!(ops::to_tree(SAMPLE, "0.9"), include_bytes!("data/clstr/tree_0.9.out"));
}

#[test]
fn clstr2blm8_nt_matches_perl() {
    assert_eq!(ops::to_blm8(SAMPLE), include_bytes!("data/clstr/blm8_nt.out"));
}

#[test]
fn clstr2blm8_aa_matches_perl() {
    assert_eq!(ops::to_blm8(PROT), include_bytes!("data/clstr/blm8_aa.out"));
}

#[test]
fn clstr_reduce_matches_perl() {
    assert_eq!(
        ops::reduce(SAMPLE, "1-2,3-10", 2),
        include_bytes!("data/clstr/reduce.out")
    );
}

#[test]
fn clstr_rev_matches_perl() {
    let fine = include_bytes!("data/clstr/rev_fine.clstr");
    let coarse = include_bytes!("data/clstr/rev_coarse.clstr");
    assert_eq!(ops::rev(fine, coarse), include_bytes!("data/clstr/rev.out"));
}

#[test]
fn clstr_merge_matches_perl() {
    // master merged with a copy of itself as the single div file.
    assert_eq!(
        ops::merge(SAMPLE, &[SAMPLE]),
        include_bytes!("data/clstr/merge_self.out")
    );
}

#[test]
fn clstr_reps_faa_rev_matches_perl() {
    let fasta = include_bytes!("data/clstr/src.fna");
    assert_eq!(
        ops::reps_faa_rev(SAMPLE, fasta, 2),
        include_bytes!("data/clstr/reps_faa_rev_c2.out")
    );
}

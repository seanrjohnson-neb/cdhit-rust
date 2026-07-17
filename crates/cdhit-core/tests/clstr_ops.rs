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

#[test]
fn clstr_select_rep_matches_perl() {
    assert_eq!(
        ops::select_rep(SAMPLE, 2, 6).unwrap(),
        include_bytes!("data/clstr/select_rep_2_6.out")
    );
}

#[test]
fn clstr_sort_prot_by_len_matches_perl() {
    assert_eq!(
        ops::sort_prot_by(PROT, "len"),
        include_bytes!("data/clstr/sort_prot_len.out")
    );
}

#[test]
fn clstr_sort_prot_by_id_matches_perl() {
    assert_eq!(
        ops::sort_prot_by(PROT, "id"),
        include_bytes!("data/clstr/sort_prot_id.out")
    );
}

#[test]
fn clstr_merge_noorder_matches_perl() {
    // master merged with a copy of itself as the single (unordered) div file.
    assert_eq!(
        ops::merge_noorder(SAMPLE, &[SAMPLE]),
        include_bytes!("data/clstr/merge_noorder_self.out")
    );
}

#[test]
fn clstr_quality_eval_by_link_matches_perl() {
    let bench = include_bytes!("data/clstr/bench.clstr");
    assert_eq!(
        ops::quality_eval_by_link(bench).unwrap(),
        include_bytes!("data/clstr/qual_eval_link.out")
    );
}

#[test]
fn plot_len1_matches_perl() {
    assert_eq!(
        ops::plot_len1(SAMPLE, "1,2-3,4-up", "1-149,150-250"),
        include_bytes!("data/clstr/plot_len1.out")
    );
}

#[test]
fn clstr_sql_tbl_sort_matches_perl() {
    let tbl = include_bytes!("data/clstr/sql_tbl.txt");
    assert_eq!(
        ops::sql_tbl_sort(tbl, 1).unwrap(),
        include_bytes!("data/clstr/sql_tbl_sort_1.out")
    );
}

#[test]
fn clstr_sql_tbl_sort_bad_level_errors() {
    // A 4-column table cannot satisfy level 2 (needs >= 6 columns).
    let tbl = include_bytes!("data/clstr/sql_tbl.txt");
    assert!(ops::sql_tbl_sort(tbl, 2).is_err());
}

#[test]
fn dup_pe_out_fastq_matches_perl() {
    let clstr = include_bytes!("data/clstr/pe.clstr");
    let r1 = include_bytes!("data/clstr/pe_R1.fq");
    let r2 = include_bytes!("data/clstr/pe_R2.fq");
    let out = ops::dup_pe_out(clstr, r1, r2);
    assert_eq!(out.out1, include_bytes!("data/clstr/pe_out_R1.fq"));
    assert_eq!(out.out2, include_bytes!("data/clstr/pe_out_R2.fq"));
}

#[test]
fn dup_pe_out_fasta_matches_perl() {
    let clstr = include_bytes!("data/clstr/pe.clstr");
    let r1 = include_bytes!("data/clstr/pe_R1.fa");
    let r2 = include_bytes!("data/clstr/pe_R2.fa");
    let out = ops::dup_pe_out(clstr, r1, r2);
    assert_eq!(out.out1, include_bytes!("data/clstr/pe_out_R1.fa"));
    assert_eq!(out.out2, include_bytes!("data/clstr/pe_out_R2.fa"));
}

#[test]
fn clstr_sql_tbl_create_matches_perl() {
    assert_eq!(
        ops::sql_tbl(PROT, None).unwrap(),
        include_bytes!("data/clstr/sqltbl_create.out")
    );
}

#[test]
fn clstr_sql_tbl_append_matches_perl() {
    let level1 = include_bytes!("data/clstr/sqltbl_create.out");
    let coarse = include_bytes!("data/clstr/sqltbl_coarse.clstr");
    assert_eq!(
        ops::sql_tbl(coarse, Some(level1)).unwrap(),
        include_bytes!("data/clstr/sqltbl_append.out")
    );
}

#[test]
fn make_multi_seq_matches_perl() {
    let fasta = include_bytes!("data/clstr/multiseq.faa");
    let files = ops::make_multi_seq(fasta, PROT, 3).unwrap();
    // Clusters 0,2,3,4,5,6 qualify (size >= 3); 1 and 7 are dropped.
    let cids: Vec<&[u8]> = files.iter().map(|f| f.cid.as_slice()).collect();
    assert_eq!(cids, vec![&b"0"[..], b"2", b"3", b"4", b"5", b"6"]);
    let by_cid = |cid: &[u8]| -> &[u8] {
        &files.iter().find(|f| f.cid == cid).unwrap().content
    };
    assert_eq!(by_cid(b"0"), include_bytes!("data/clstr/multiseq_c3_0.out"));
    assert_eq!(by_cid(b"2"), include_bytes!("data/clstr/multiseq_c3_2.out"));
    assert_eq!(by_cid(b"3"), include_bytes!("data/clstr/multiseq_c3_3.out"));
    assert_eq!(by_cid(b"4"), include_bytes!("data/clstr/multiseq_c3_4.out"));
    assert_eq!(by_cid(b"5"), include_bytes!("data/clstr/multiseq_c3_5.out"));
    assert_eq!(by_cid(b"6"), include_bytes!("data/clstr/multiseq_c3_6.out"));
}

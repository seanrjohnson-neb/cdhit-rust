//! Output writers: the representative FASTA (`-o`) and the `.clstr` membership
//! file.
//!
//! Port of `Sequence::PrintInfo` (cdhit-common.c++:1665-1682),
//! `WriteExtra1D` (:2592-2651), and `WriteClusters` (:2385-2421). The `.clstr`
//! text and the representative FASTA bytes must match the C++ exactly.

use crate::options::Options;
use crate::sequence::{Sequence, IS_MINUS_STRAND};

/// Format one `.clstr` member line (port of `PrintInfo`). `id` is the member's
/// position within its cluster.
fn print_info(out: &mut String, seq: &Sequence, id: i32, options: &Options) {
    let tag = if options.is_est { "nt" } else { "aa" };
    let strand = options.is_est;
    // identifier retains its leading '>'; skip it (C uses identifier+1).
    let ident = &seq.identifier[1.min(seq.identifier.len())..];
    let ident = String::from_utf8_lossy(ident);
    out.push_str(&format!("{}\t{}{}, >{}...", id, seq.size, tag, ident));
    if seq.identity != 0.0 {
        out.push_str(" at ");
        if options.print != 0 {
            let c = &seq.coverage;
            out.push_str(&format!("{}:{}:{}:{}/", c[0], c[1], c[2], c[3]));
        }
        if strand {
            let s = if seq.state & IS_MINUS_STRAND != 0 { '-' } else { '+' };
            out.push_str(&format!("{}/", s));
        }
        // identity*100 done in f32 then widened to f64 for the %.2f print.
        let pc = (seq.identity * 100.0f32) as f64;
        out.push_str(&format!("{:.2}%", pc));
        if options.use_distance {
            let d = (seq.distance * 100.0f32) as f64;
            out.push_str(&format!("/{:.2}%", d));
        }
        out.push('\n');
    } else {
        out.push_str(" *\n");
    }
}

/// Port of `quick_sort_idxr` (cdhit-common.c++:4052-4075): sort `a` descending,
/// carrying `idx` along. Ported verbatim (its exact partition/pivot decides the
/// order of equal-size clusters), so `-sc` output matches the C++.
pub fn quick_sort_idxr(a: &mut [i32], idx: &mut [i32], lo0: i32, hi0: i32) {
    if hi0 > lo0 {
        let mut lo = lo0;
        let mut hi = hi0;
        let mid = a[((lo0 + hi0) / 2) as usize];
        while lo <= hi {
            while lo < hi0 && a[lo as usize] > mid {
                lo += 1;
            }
            while hi > lo0 && a[hi as usize] < mid {
                hi -= 1;
            }
            if lo <= hi {
                a.swap(lo as usize, hi as usize);
                idx.swap(lo as usize, hi as usize);
                lo += 1;
                hi -= 1;
            }
        }
        if lo0 < hi {
            quick_sort_idxr(a, idx, lo0, hi);
        }
        if lo < hi0 {
            quick_sort_idxr(a, idx, lo, hi0);
        }
    }
}

/// Build the `.clstr` file contents (port of `WriteExtra1D`). Members within
/// each cluster are ordered by original input index; with `sort_output` (`-sc`)
/// clusters are emitted in descending size order and renumbered.
pub fn write_extra_1d(sequences: &[Sequence], rep_seqs: &[i32], options: &Options) -> String {
    let n = sequences.len();
    // sorting = (index << 32) | i, sorted -> sequences in original-index order.
    let mut sorting: Vec<u64> = (0..n)
        .map(|i| ((sequences[i].index as u64) << 32) | i as u64)
        .collect();
    sorting.sort_unstable();

    let m = rep_seqs.len();
    let mut clusters: Vec<Vec<usize>> = vec![Vec::new(); m];
    for &s in &sorting {
        let k = (s & 0xffff_ffff) as usize;
        let id = sequences[k].cluster_id;
        clusters[id as usize].push(k);
    }

    let mut out = String::new();
    if options.sort_output != 0 {
        let mut clstr_size: Vec<i32> = (0..m).map(|i| clusters[i].len() as i32).collect();
        let mut clstr_idx1: Vec<i32> = (0..m as i32).collect();
        quick_sort_idxr(&mut clstr_size, &mut clstr_idx1, 0, m as i32 - 1);
        for i in 0..m {
            let i0 = clstr_idx1[i] as usize;
            out.push_str(&format!(">Cluster {}\n", i));
            for (k, &member) in clusters[i0].iter().enumerate() {
                print_info(&mut out, &sequences[member], k as i32, options);
            }
        }
    } else {
        for i in 0..m {
            out.push_str(&format!(">Cluster {}\n", i));
            for (k, &member) in clusters[i].iter().enumerate() {
                print_info(&mut out, &sequences[member], k as i32, options);
            }
        }
    }
    out
}

/// Build the 2D `.clstr` file (port of `WriteExtra2D`, cdhit-common.c++:
/// 2652-2693). Each db1 sequence is its own cluster header/representative, with
/// matching (redundant) db2 sequences listed beneath it, ordered by db2 index.
pub fn write_extra_2d(
    db1: &[Sequence],
    db2: &[Sequence],
    options: &Options,
) -> String {
    let n = db1.len();
    let n2 = db2.len();
    // db2 grouped by cluster_id (= db1 index), in db2 original-index order.
    let mut sorting: Vec<u64> = (0..n2)
        .map(|i| ((db2[i].index as u64) << 32) | i as u64)
        .collect();
    sorting.sort_unstable();

    let mut clusters: Vec<Vec<usize>> = vec![Vec::new(); n];
    for &s in &sorting {
        let i = (s & 0xffff_ffff) as usize;
        if db2[i].state & crate::sequence::IS_REDUNDANT != 0 {
            let id = db2[i].cluster_id as usize;
            clusters[id].push(i);
        }
    }

    let mut out = String::new();
    for i in 0..n {
        out.push_str(&format!(">Cluster {}\n", i));
        print_info(&mut out, &db1[i], 0, options);
        for (k, &member) in clusters[i].iter().enumerate() {
            print_info(&mut out, &db2[member], k as i32 + 1, options);
        }
    }
    out
}

/// Build the representative FASTA (`-o`) by copying the original record bytes
/// for each representative, ordered by original input index (port of
/// `WriteClusters`, default `sort_outputf == 0`).
pub fn write_clusters(sequences: &[Sequence], rep_seqs: &[i32], raw: &[u8]) -> Vec<u8> {
    let n = rep_seqs.len();
    let mut sorting: Vec<u64> = (0..n)
        .map(|i| {
            let ri = rep_seqs[i] as usize;
            ((sequences[ri].index as u64) << 32) | ri as u64
        })
        .collect();
    sorting.sort_unstable();

    let mut out = Vec::new();
    for &s in &sorting {
        let seq = &sequences[(s & 0xffff_ffff) as usize];
        let begin = seq.des_begin as usize;
        let end = begin + seq.tot_length as usize;
        out.extend_from_slice(&raw[begin..end.min(raw.len())]);
    }
    out
}

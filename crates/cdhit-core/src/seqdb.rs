//! `SequenceDB` — the master container and clustering driver.
//!
//! This module grows across phases. For now it holds the sequence set, the
//! length statistics, and `sort_divide` (the stable counting sort). Later
//! phases add the word table, clustering, and output.

use crate::options::{Options, Scoring};
use crate::sequence::Sequence;

/// Matches `MAX_TABLE_SEQ` (cdhit-common.h). Used for `max_entries` and the
/// per-table representative cap.
pub const MAX_TABLE_SEQ: usize = 4_000_000;

pub struct SequenceDb {
    pub sequences: Vec<Sequence>,
    pub rep_seqs: Vec<i32>,
    /// Number of word-table rows = `NAAN_array[NAA]` (set by the caller).
    pub naan: i32,
    pub total_letter: i64,
    pub total_desc: i64,
    pub max_len: i32,
    pub min_len: i32,
    pub len_n50: i32,
}

impl SequenceDb {
    pub fn new(sequences: Vec<Sequence>) -> Self {
        SequenceDb {
            sequences,
            rep_seqs: Vec::new(),
            naan: 0,
            total_letter: 0,
            total_desc: 0,
            max_len: 0,
            min_len: 0,
            len_n50: 0,
        }
    }

    /// Port of `SequenceDB::SortDivide` (cdhit-common.c++:2228-2300).
    ///
    /// Computes length statistics, encodes each sequence to the integer
    /// alphabet, and (when `sort`) performs a **stable** counting sort from
    /// longest to shortest — sequences of equal length retain their input
    /// order, which is load-bearing for which sequence becomes a cluster
    /// representative. Also computes `len_n50` and `max_entries`.
    pub fn sort_divide(&mut self, options: &mut Options, scoring: &Scoring, sort: bool) {
        let n = self.sequences.len();
        self.total_letter = 0;
        self.total_desc = 0;
        self.max_len = 0;
        self.min_len = i32::MAX;

        for seq in &mut self.sequences {
            let len = seq.size;
            self.total_letter += len as i64;
            if len > self.max_len {
                self.max_len = len;
            }
            if len < self.min_len {
                self.min_len = len;
            }
            seq.convert_bases(&scoring.aa2idx);
            self.total_desc += seq.identifier.len() as i64;
        }
        if n == 0 {
            self.min_len = 0;
        }

        options.max_entries = (self.max_len as u64) * (MAX_TABLE_SEQ as u64);

        // len_n50 default before sorting (cdhit-common.c++:2257).
        self.len_n50 = (self.max_len + self.min_len) / 2;

        if sort {
            let m = (self.max_len - self.min_len + 1) as usize;
            let mut count = vec![0i32; m];
            let mut accum = vec![0i32; m];
            let offset_len = m;
            let mut offset = vec![0i32; offset_len];
            let mut sorting: Vec<usize> = vec![0; n];

            for seq in &self.sequences {
                count[(self.max_len - seq.size) as usize] += 1;
            }
            for i in 1..m {
                accum[i] = accum[i - 1] + count[i - 1];
            }
            let mut sum: i64 = 0;
            for i in 0..m {
                sum += (self.max_len - i as i32) as i64 * count[i] as i64;
                if sum >= self.total_letter >> 1 {
                    self.len_n50 = self.max_len - i as i32;
                    break;
                }
            }
            for (i, seq) in self.sequences.iter().enumerate() {
                let len = (self.max_len - seq.size) as usize;
                let id = accum[len] + offset[len];
                sorting[id as usize] = i;
                offset[len] += 1;
            }
            // Reorder sequences into sorted order.
            let old = std::mem::take(&mut self.sequences);
            let mut old_opt: Vec<Option<Sequence>> = old.into_iter().map(Some).collect();
            let mut sorted: Vec<Sequence> = Vec::with_capacity(n);
            for &src in &sorting {
                sorted.push(old_opt[src].take().unwrap());
            }
            self.sequences = sorted;

            options.max_entries = 0;
            for (i, seq) in self.sequences.iter().enumerate() {
                if i < MAX_TABLE_SEQ {
                    options.max_entries += seq.size as u64;
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::read_database;

    fn seq(id: &str, len: usize) -> Sequence {
        let mut s = Sequence::new();
        s.identifier = format!(">{id}").into_bytes();
        s.data = vec![b'A'; len];
        s.size = len as i32;
        s
    }

    #[test]
    fn counting_sort_is_stable_longest_first() {
        // Two length-5 (b, c in input order), one length-8 (a), one length-3 (d).
        let seqs = vec![seq("b", 5), seq("a", 8), seq("c", 5), seq("d", 3)];
        let mut db = SequenceDb::new(seqs);
        let mut opt = Options::default();
        let sc = Scoring::default();
        db.sort_divide(&mut opt, &sc, true);
        let ids: Vec<String> = db
            .sequences
            .iter()
            .map(|s| String::from_utf8_lossy(&s.identifier[1..]).into_owned())
            .collect();
        // longest first; equal-length b before c (stable).
        assert_eq!(ids, vec!["a", "b", "c", "d"]);
        assert_eq!(db.max_len, 8);
        assert_eq!(db.min_len, 3);
    }

    #[test]
    fn fesod_longest_first_matches_cpp_cluster_order() {
        // After sort, the first (longest) sequence is the rep of Cluster 0 at
        // -c 0.9: FeSOD_A0A060HP82 (215aa), per the C++ reference.
        let fa = include_bytes!("../tests/data/fesod.fasta");
        let mut opt = Options::default();
        let sc = Scoring::default();
        let db_in = read_database(fa, &opt);
        let mut db = SequenceDb::new(db_in.sequences);
        db.sort_divide(&mut opt, &sc, true);
        assert_eq!(db.max_len, 215);
        assert_eq!(
            String::from_utf8_lossy(&db.sequences[0].identifier[1..]),
            "FeSOD_A0A060HP82"
        );
        // Verify non-increasing length order.
        for w in db.sequences.windows(2) {
            assert!(w[0].size >= w[1].size);
        }
    }
}

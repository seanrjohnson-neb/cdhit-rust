//! `WordTable` — the short-word (k-mer) inverted index.
//!
//! Port of `WordTable` (cdhit-common.h:225-253) and its methods
//! (cdhit-common.c++:1254-1478). Each of the `NAAN` rows holds the list of
//! `(table-local sequence index, word count)` pairs for representatives that
//! contain that encoded word.

/// `{index, count}` pair (cdhit-common.h:215-221).
#[derive(Clone, Copy, Default, PartialEq, Debug)]
pub struct IndexCount {
    pub index: i32,
    pub count: i32,
}

impl IndexCount {
    pub fn new(index: i32, count: i32) -> Self {
        IndexCount { index, count }
    }
}

pub struct WordTable {
    pub index_counts: Vec<Vec<IndexCount>>,
    /// Global sequence indices of the representatives added to this table, in
    /// table order (the C++ stores `Sequence*`; we store the global index).
    pub sequences: Vec<usize>,
    pub naa: i32,
    pub naan: i32,
    pub is_aa: bool,
    pub size: usize,
    pub frag_count: i32,
}

impl WordTable {
    /// Port of `WordTable::WordTable` + `Init` (cdhit-common.c++:1254-1274).
    pub fn new(naa: i32, naan: i32) -> WordTable {
        WordTable {
            index_counts: vec![Vec::new(); naan.max(0) as usize],
            sequences: Vec::new(),
            naa,
            naan,
            is_aa: true,
            size: 0,
            frag_count: 0,
        }
    }

    pub fn set_dna(&mut self) {
        self.is_aa = false;
    }

    /// Port of `WordTable::Clear` (cdhit-common.c++:1276-1297).
    pub fn clear(&mut self) {
        self.size = 0;
        self.frag_count = 0;
        self.sequences.clear();
        for row in &mut self.index_counts {
            row.clear();
        }
    }

    /// Add a representative's word counts from a pre-built `IndexCount` list
    /// (port of the `NVector<IndexCount>` overload, cdhit-common.c++:1299-1316).
    /// `seq_global_idx` is recorded in `sequences`; the table-local index used
    /// inside rows is the current `sequences.len()`.
    pub fn add_word_counts_seq(
        &mut self,
        counts: &[IndexCount],
        seq_global_idx: usize,
        skip_n: bool,
    ) {
        let idx = self.sequences.len() as i32;
        for ic in counts {
            let k = ic.count;
            if k != 0 {
                let j = ic.index;
                if skip_n && j < 0 {
                    continue;
                }
                self.index_counts[j as usize].push(IndexCount::new(idx, k));
                self.size += 1;
            }
        }
        self.sequences.push(seq_global_idx);
    }

    /// Add a representative's word counts from encoded word arrays (port of the
    /// `word_encodes` overload, cdhit-common.c++:1321-1337). The caller supplies
    /// the table-local index `idx`.
    pub fn add_word_counts_encoded(
        &mut self,
        aan_no: usize,
        word_encodes: &[i32],
        word_encodes_no: &[u32],
        idx: i32,
        skip_n: bool,
    ) {
        for i in 0..aan_no {
            let k = word_encodes_no[i];
            if k != 0 {
                let j = word_encodes[i];
                if skip_n && j < 0 {
                    continue;
                }
                self.index_counts[j as usize].push(IndexCount::new(idx, k as i32));
                self.size += 1;
            }
        }
    }

    /// Port of `WordTable::CountWords` (cdhit-common.c++:1431-1478).
    ///
    /// Accumulates, for each candidate representative sharing words with the
    /// query, the summed `min(rep_word_count, query_word_count)` into
    /// `look_counts`. `index_mapping` (scratch, sized reps+2) tracks the
    /// position of each rep in `look_counts` (value = position+1, 0 = absent).
    /// A final sentinel entry with `count = 0` is written at `look_counts[len]`.
    /// Returns the number of populated entries.
    ///
    /// `prev_len` is the number of entries populated by the previous call on the
    /// same scratch buffers, used to reset `index_mapping`. This mirrors the C++
    /// where `lookCounts.size` persists on the per-thread buffer between calls.
    #[allow(clippy::too_many_arguments)]
    pub fn count_words(
        &self,
        aan_no: usize,
        word_encodes: &[i32],
        word_encodes_no: &[u32],
        look_counts: &mut [IndexCount],
        index_mapping: &mut [u32],
        est: bool,
        min: i32,
        prev_len: usize,
    ) -> usize {
        // Reset mapping for previously-populated entries.
        for j in 0..prev_len {
            let ix = look_counts[j].index;
            index_mapping[ix as usize] = 0;
        }
        let mut look_len = 0usize;

        // Skip leading N-words (est) that were marked with encode < 0.
        let mut j0 = 0usize;
        if est {
            while j0 < aan_no && word_encodes[j0] < 0 {
                j0 += 1;
            }
        }

        while j0 < aan_no {
            let j = word_encodes[j0];
            let j1 = word_encodes_no[j0];
            if j1 == 0 {
                j0 += 1;
                continue;
            }
            let row = &self.index_counts[j as usize];
            let rest = aan_no as i32 - j0 as i32 + 1;
            for ic in row {
                let c = if ic.count < j1 as i32 {
                    ic.count
                } else {
                    j1 as i32
                };
                let idm = &mut index_mapping[ic.index as usize];
                if *idm == 0 {
                    if rest < min {
                        continue;
                    }
                    look_counts[look_len].index = ic.index;
                    look_counts[look_len].count = c;
                    look_len += 1;
                    *idm = look_len as u32;
                } else {
                    look_counts[(*idm - 1) as usize].count += c;
                }
            }
            j0 += 1;
        }
        // Sentinel (cdhit-common.c++:1475).
        look_counts[look_len].count = 0;
        look_len
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet::Naa;
    use crate::buffer::WorkingBuffer;
    use crate::options::Options;

    // Encode a short protein sequence's words and add to the table, then count
    // shared words against a query — checked by hand.
    #[test]
    fn count_words_shared_kmers() {
        let naa_len = 2i32;
        let naa = Naa::init(21);
        let naan = naa.array[naa_len as usize];
        let opt = Options::default();
        let mut buf = WorkingBuffer::new(8, 64, &opt);

        // Two identical sequences of codes; all 2-mers shared.
        // codes: A=0, R=1 -> "ARAR" = [0,1,0,1]
        let rep: Vec<u8> = vec![0, 1, 0, 1];
        let query: Vec<u8> = vec![0, 1, 0, 1];

        let mut table = WordTable::new(naa_len, naan);

        // Encode rep and add.
        let aan_no = rep.len() - naa_len as usize + 1; // 3 words
        buf.encode_words(&rep, rep.len() as i32, naa_len, &naa, false);
        table.add_word_counts_encoded(
            aan_no,
            &buf.word_encodes,
            &buf.word_encodes_no,
            0,
            false,
        );
        table.sequences.push(0);

        // Encode query and count.
        buf.encode_words(&query, query.len() as i32, naa_len, &naa, false);
        let mut look = vec![IndexCount::default(); 8];
        let mut map = vec![0u32; 8];
        let n = table.count_words(
            aan_no,
            &buf.word_encodes,
            &buf.word_encodes_no,
            &mut look,
            &mut map,
            false,
            0,
            0,
        );
        // One candidate (rep 0). Shared word count: query has words AR,RA,AR ->
        // encodes: AR=0*21+1=1 (x2), RA=1*21+0=21 (x1). Rep same. min counts:
        // AR min(2,2)=2, RA min(1,1)=1 => total 3.
        assert_eq!(n, 1);
        assert_eq!(look[0].index, 0);
        assert_eq!(look[0].count, 3);
    }

    #[test]
    fn encode_words_run_length_counts() {
        let naa = Naa::init(21);
        let opt = Options::default();
        let mut buf = WorkingBuffer::new(4, 64, &opt);
        // "ARAR" 2-mers: AR(1), RA(21), AR(1) -> sorted [1,1,21], counts [2,0,1].
        let data: Vec<u8> = vec![0, 1, 0, 1];
        buf.encode_words(&data, 4, 2, &naa, false);
        assert_eq!(&buf.word_encodes[0..3], &[1, 1, 21]);
        assert_eq!(&buf.word_encodes_no[0..3], &[2, 0, 1]);
    }
}

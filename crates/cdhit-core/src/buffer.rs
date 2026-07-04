//! `WorkingBuffer` — per-thread scratch memory and the word/AAP encoders.
//!
//! Port of `WorkingBuffer` (cdhit-common.h:468-532) and its methods
//! `EncodeWords`, `ComputeAAP`, `ComputeAAP2` (cdhit-common.c++:2836-2909).
//! The DP matrices (`score_mat`/`back_mat`) are added in the alignment phase.

use crate::alphabet::{Naa, MAX_SEQ, MAX_UAA};
use crate::options::Options;
use crate::wordtable::IndexCount;

/// `MAX_DIAG = MAX_SEQ << 1` (cdhit-common.h:55).
pub const MAX_DIAG: usize = (MAX_SEQ as usize) << 1;

pub struct WorkingBuffer {
    pub taap: Vec<i32>,
    pub word_encodes: Vec<i32>,
    pub word_encodes_backup: Vec<i32>,
    pub word_encodes_no: Vec<u32>,
    pub aap_list: Vec<u32>,
    pub aap_begin: Vec<u32>,
    pub look_counts: Vec<IndexCount>,
    /// `indexMapping[i]` = (position in `look_counts`)+1, or 0 if absent.
    pub index_mapping: Vec<u32>,
    pub diag_score: Vec<i32>,
    pub diag_score2: Vec<i32>,
    pub aan_list_comp: Vec<i32>,
    pub seqi_comp: Vec<u8>,
}

impl WorkingBuffer {
    /// Port of `WorkingBuffer::Set` (cdhit-common.c++:491-527). `frag` is the
    /// max representatives per table batch; `maxlen` the longest sequence.
    pub fn new(frag: usize, maxlen: usize, options: &Options) -> WorkingBuffer {
        let est = options.is_est;
        let mut m = (MAX_UAA as usize) * (MAX_UAA as usize);
        if est {
            m *= m;
        }
        let max_len = maxlen;
        let mut frag = frag;
        if frag > crate::seqdb::MAX_TABLE_SEQ {
            frag = crate::seqdb::MAX_TABLE_SEQ;
        }
        WorkingBuffer {
            taap: vec![0; m],
            word_encodes: vec![0; max_len],
            word_encodes_backup: vec![0; max_len],
            word_encodes_no: vec![0; max_len],
            aap_list: vec![0; max_len],
            aap_begin: vec![0; m],
            look_counts: vec![IndexCount::default(); frag + 2],
            index_mapping: vec![0; frag + 2],
            diag_score: vec![0; MAX_DIAG],
            diag_score2: vec![0; MAX_DIAG],
            aan_list_comp: vec![0; max_len],
            seqi_comp: vec![0; MAX_SEQ as usize],
        }
    }

    /// Port of `WorkingBuffer::EncodeWords` (cdhit-common.c++:2836-2872).
    ///
    /// Encodes each length-`naa` word of `data` as a base-`NAA` integer, marks
    /// words containing `N` (code >= 4) as -1 in EST mode, then sorts the word
    /// list and run-length-counts duplicates into `word_encodes_no`. Returns the
    /// count of skipped (N-containing) words.
    pub fn encode_words(&mut self, data: &[u8], size: i32, naa: i32, naan: &Naa, est: bool) -> i32 {
        let len = size as usize;
        let aan_no = len as i32 - naa + 1;
        let aan_no = aan_no.max(0) as usize;
        let mut skip = 0i32;

        for j in 0..aan_no {
            let mut encode = 0i32;
            // encode += word[k] * NAAN_array[NAA-1-k]
            let mut k1 = (naa - 1) as usize;
            for k in 0..naa as usize {
                encode += (data[j + k] as i32) * naan.array[k1];
                if k1 > 0 {
                    k1 -= 1;
                }
            }
            self.word_encodes[j] = encode;
            self.word_encodes_backup[j] = encode;
        }

        if est {
            for j in 0..len {
                if data[j] >= 4 {
                    let i0 = if j as i32 - naa + 1 > 0 {
                        (j as i32 - naa + 1) as usize
                    } else {
                        0
                    };
                    let i1 = if (j as i32) < aan_no as i32 {
                        j
                    } else {
                        aan_no - 1
                    };
                    for i in i0..=i1 {
                        self.word_encodes[i] = -1;
                    }
                }
            }
            for j in 0..aan_no {
                skip += (self.word_encodes[j] == -1) as i32;
            }
        }

        self.word_encodes[0..aan_no].sort_unstable();
        for j in 0..aan_no {
            self.word_encodes_no[j] = 1;
        }
        // for(j=aan_no-1; j; j--): merge runs of equal encodings, right to left.
        for j in (1..aan_no).rev() {
            if self.word_encodes[j] == self.word_encodes[j - 1] {
                self.word_encodes_no[j - 1] += self.word_encodes_no[j];
                self.word_encodes_no[j] = 0;
            }
        }
        skip
    }

    /// Port of `WorkingBuffer::ComputeAAP` (protein amino-acid-pair index,
    /// cdhit-common.c++:2874-2890).
    pub fn compute_aap(&mut self, seqi: &[u8], size: i32, naa1: i32, naa2: i32) {
        let len1 = (size - 1) as usize;
        for sk in 0..naa2 as usize {
            self.taap[sk] = 0;
        }
        for j1 in 0..len1 {
            let c22 = (seqi[j1] as i32) * naa1 + seqi[j1 + 1] as i32;
            self.taap[c22 as usize] += 1;
        }
        let mut mm = 0i32;
        for sk in 0..naa2 as usize {
            self.aap_begin[sk] = mm as u32;
            mm += self.taap[sk];
            self.taap[sk] = 0;
        }
        for j1 in 0..len1 {
            let c22 = ((seqi[j1] as i32) * naa1 + seqi[j1 + 1] as i32) as usize;
            let pos = self.aap_begin[c22] + self.taap[c22] as u32;
            self.aap_list[pos as usize] = j1 as u32;
            self.taap[c22] += 1;
        }
    }

    /// Port of `WorkingBuffer::ComputeAAP2` (EST 4-mer index,
    /// cdhit-common.c++:2891-2909).
    pub fn compute_aap2(&mut self, seqi: &[u8], size: i32, naa: &Naa) {
        let (naa1, naa2, naa3, naa4) = (naa.array[1], naa.array[2], naa.array[3], naa.array[4]);
        let len1 = (size - 3) as usize;
        for sk in 0..naa4 as usize {
            self.taap[sk] = 0;
        }
        let is_n = |b: u8| b >= 4;
        for j1 in 0..len1 {
            if is_n(seqi[j1]) || is_n(seqi[j1 + 1]) || is_n(seqi[j1 + 2]) || is_n(seqi[j1 + 3]) {
                continue;
            }
            let c22 = (seqi[j1] as i32) * naa3
                + (seqi[j1 + 1] as i32) * naa2
                + (seqi[j1 + 2] as i32) * naa1
                + seqi[j1 + 3] as i32;
            self.taap[c22 as usize] += 1;
        }
        let mut mm = 0i32;
        for sk in 0..naa4 as usize {
            self.aap_begin[sk] = mm as u32;
            mm += self.taap[sk];
            self.taap[sk] = 0;
        }
        for j1 in 0..len1 {
            if is_n(seqi[j1]) || is_n(seqi[j1 + 1]) || is_n(seqi[j1 + 2]) || is_n(seqi[j1 + 3]) {
                continue;
            }
            let c22 = ((seqi[j1] as i32) * naa3
                + (seqi[j1 + 1] as i32) * naa2
                + (seqi[j1 + 2] as i32) * naa1
                + seqi[j1 + 3] as i32) as usize;
            let pos = self.aap_begin[c22] + self.taap[c22] as u32;
            self.aap_list[pos as usize] = j1 as u32;
            self.taap[c22] += 1;
        }
    }
}

//! `read-linker` — join paired-end reads by 3'/5' overlap.
//!
//! Port of `cd-hit-auxtools/read-linker.cxx`. For each read pair, the second
//! end is reverse-complemented and overlapped against the first; on a
//! sufficient overlap the two are merged into one (or, when the overlap has
//! real mismatches, 2^k) contig(s). Output is bit-for-bit compatible with the
//! reference, including the quirk that the emitted `mismatch_no=` uses the
//! Check error count (which also counts `N`/`N` pairs) rather than the number
//! of real substitutions used to enumerate contigs.
//!
//! Input is consumed in batches of 10000 pairs (as the C++ `SequenceCache`
//! does); this reproduces both the per-batch `handled:` log lines and the
//! behaviour of discarding the tail of a batch when the two files are of
//! unequal length.

use super::bioseq::{parse_fastaq, Sequence};

const BATCH: usize = 10000;

/// Parameters mirroring the `read-linker` command line.
pub struct LinkerParams {
    /// `-l`: minimum overlap length. Negative selects the fixed-selection
    /// ("Link2") mode.
    pub min: i32,
    /// `-e`: maximum number of mismatches tolerated in the overlap.
    pub error: i32,
    /// `-m`: cutoff length used only in the negative-`-l` mode.
    pub maxlen: i32,
}

/// Result of a `read-linker` run.
pub struct LinkerOutput {
    /// The linked contigs, as FASTA or FASTQ bytes.
    pub output: Vec<u8>,
    /// The text the C++ prints to stdout (progress and summary).
    pub log: String,
}

struct Linker {
    min_overlap: i32,
    tolerance: i32,
    maxlen: i32,
    overlap: i32,
    error: i32,
    count: usize,
    sequences: Vec<Sequence>,
}

impl Linker {
    fn new(min_overlap: i32, tolerance: i32) -> Self {
        Linker {
            min_overlap,
            tolerance,
            maxlen: 0,
            overlap: 0,
            error: 0,
            count: 0,
            sequences: Vec::new(),
        }
    }

    fn reset(&mut self) {
        self.count = 0;
        self.sequences.clear();
    }

    /// Port of `SequenceLinker::Check`. Finds the longest overlap (>= min) whose
    /// mismatch count is within tolerance. Note that an `N`/`N` pair is counted
    /// as an error here, matching the C++.
    fn check(&mut self, first: &[u8], second: &[u8]) -> i32 {
        let n1 = first.len() as i32;
        let n2 = second.len() as i32;
        let m = n1.min(n2);
        let mut k = m;
        while k >= self.min_overlap {
            self.error = 0;
            let mut i = n1 - k;
            let mut j = 0;
            while j < k {
                let a = first[i as usize];
                let b = second[j as usize];
                if a != b || (a == b'N' && b == b'N') {
                    self.error += 1;
                }
                if self.error > self.tolerance {
                    break;
                }
                i += 1;
                j += 1;
            }
            self.overlap = k;
            if self.error <= self.tolerance {
                return k;
            }
            k -= 1;
        }
        self.overlap = 0;
        0
    }

    /// Port of `SequenceLinker::Link` (the `-l > 0` path).
    fn link(&mut self, first: &Sequence, second: &Sequence) -> usize {
        let has_qs = !first.qs.is_empty();
        let s1 = &first.seq;
        let s2 = &second.seq;
        let qs1 = &first.qs;
        let qs2 = &second.qs;
        let n1 = s1.len();
        let n2 = s2.len();
        let o = self.check(s1, s2);
        if o == 0 {
            return 0;
        }
        let o = o as usize;

        // Enumerate real substitutions (N-tolerant) and record their positions.
        let mut error2 = 0u32;
        let mut locs: Vec<u8> = Vec::new();
        {
            let mut i = n1 - o;
            let mut j = 0usize;
            while j < o {
                let a = s1[i];
                let b = s2[j];
                if a != b && a != b'N' && b != b'N' {
                    error2 += 1;
                }
                if a != b {
                    if !locs.is_empty() {
                        locs.push(b',');
                    }
                    locs.extend_from_slice((i + 1).to_string().as_bytes());
                }
                i += 1;
                j += 1;
            }
        }
        let m = 1usize << error2;

        // Build the m contigs (identical except at the real-mismatch columns).
        let tag_error = self.error; // NB: Check error, incl. N/N pairs
        let base = self.count;
        let need = base + m;
        while self.sequences.len() < need {
            self.sequences.push(Sequence::new());
        }
        for i in 0..m {
            let seq = &mut self.sequences[base + i];
            seq.des.clear();
            seq.seq.clear();
            seq.qs.clear();
            let tag = format!(
                ".contig.{} length={} overlap={} mismatch_no={}",
                i + 1,
                n1 + n2 - o,
                o,
                tag_error
            );
            seq.des = first.get_description(0);
            seq.des.extend_from_slice(tag.as_bytes());
            if !locs.is_empty() {
                seq.des.extend_from_slice(b" mismatch_pos=");
                seq.des.extend_from_slice(&locs);
            }
            seq.seq.extend_from_slice(s1);
            seq.seq.extend_from_slice(&s2[o..]);
            if has_qs {
                seq.qs.extend_from_slice(qs1);
                seq.qs.extend_from_slice(&qs2[o..]);
            }
        }

        // Fill the overlap columns: consensus where one side is N, otherwise
        // split the 2^k contigs across the two alleles.
        let mut step = m;
        {
            let mut i = n1 - o;
            let mut j = 0usize;
            while j < o {
                let a = s1[i];
                let b = s2[j];
                if a == b {
                    i += 1;
                    j += 1;
                    continue;
                }
                if a == b'N' {
                    for k in 0..m {
                        let seq = &mut self.sequences[base + k];
                        seq.seq[i] = b;
                        if has_qs {
                            seq.qs[i] = qs2[j];
                        }
                    }
                } else if b == b'N' {
                    for k in 0..m {
                        let seq = &mut self.sequences[base + k];
                        seq.seq[i] = a;
                        if has_qs {
                            seq.qs[i] = qs1[i];
                        }
                    }
                } else {
                    step >>= 1;
                    for k in 0..m {
                        let seq = &mut self.sequences[base + k];
                        let odd = (k / step) % 2;
                        seq.seq[i] = if odd == 1 { b } else { a };
                        if has_qs {
                            seq.qs[i] = if odd == 1 { qs2[j] } else { qs1[i] };
                        }
                    }
                }
                i += 1;
                j += 1;
            }
        }
        self.count += m;
        m
    }

    /// Port of `SequenceLinker::Link2` (the `-l < 0` fixed-selection path).
    fn link2(&mut self, first: &Sequence, second: &Sequence) -> usize {
        let has_qs = !first.qs.is_empty();
        let s1 = &first.seq;
        let s2 = &second.seq;
        let qs2 = &second.qs;
        let n1 = s1.len() as i32;
        let n2 = s2.len() as i32;
        let min = self.min_overlap;
        let maxlen = self.maxlen;

        if n1 < min || n2 < min {
            return 0;
        }
        if n2 < maxlen {
            return 0;
        }
        let min = min as usize;
        let maxlen = maxlen as usize;

        let mut joint = Sequence::new();
        joint.des = first.get_description(0);
        joint.seq.extend_from_slice(&s1[..min]);
        joint.seq.extend_from_slice(&s2[(maxlen - min)..maxlen]);
        if has_qs {
            joint.qs.extend_from_slice(&first.qs[..min]);
            joint.qs.extend_from_slice(&qs2[(maxlen - min)..maxlen]);
        }
        self.sequences.push(joint);
        self.count += 1;
        1
    }

    fn write(&self, out: &mut Vec<u8>) {
        for i in 0..self.count {
            self.sequences[i].print(out);
        }
    }
}

/// Run `read-linker` over two in-memory paired-end files.
pub fn read_linker(first: &[u8], second: &[u8], params: &LinkerParams) -> LinkerOutput {
    let min = params.min;
    let error = params.error;
    let mut maxlen = params.maxlen;

    let seqs1 = parse_fastaq(first, false);
    let mut seqs2 = parse_fastaq(second, false);

    let mut linker = Linker::new(min.abs(), error);
    linker.maxlen = maxlen;

    let mut output = Vec::new();
    let mut log = String::new();

    let mut omin: i32 = 0x7fff_ffff;
    let mut omax: i32 = 0;
    let mut count = 0usize;
    let mut count_input_pairs = 0i64;
    let mut count_used_pairs = 0i64;
    let mut count_all_contig = 0i64;
    let ecap = (error.max(0) + 1) as usize;
    let mut count_pairs = vec![0i64; ecap];
    let mut count_contigs = vec![0i64; ecap];

    let mut off1 = 0usize;
    let mut off2 = 0usize;
    loop {
        let n1 = (seqs1.len() - off1).min(BATCH);
        let n2 = (seqs2.len() - off2).min(BATCH);
        let n = n1.min(n2);
        if n1 != n2 {
            log.push_str("Warning: the pair end files contain different number of reads!\n");
        }
        if n == 0 {
            break;
        }
        if min < 0 && maxlen == 0 {
            maxlen = seqs2[off2].len() as i32;
            linker.maxlen = maxlen;
        }
        for j in 0..n {
            seqs2[off2 + j].to_reverse_complement();
            let c = if min > 0 {
                let c = linker.link(&seqs1[off1 + j], &seqs2[off2 + j]);
                let o = linker.overlap;
                if o > omax {
                    omax = o;
                }
                if o != 0 && o < omin {
                    omin = o;
                }
                c
            } else {
                linker.link2(&seqs1[off1 + j], &seqs2[off2 + j])
            };
            if c != 0 {
                count_used_pairs += 1;
                let e = linker.error as usize;
                count_pairs[e] += 1;
                count_contigs[e] += c as i64;
            }
        }
        count_input_pairs += n as i64;
        count_all_contig += linker.count as i64;
        linker.write(&mut output);
        linker.reset();
        count += n;
        log.push_str(&format!("handled: {:9}\n", count));
        off1 += n1;
        off2 += n2;
    }

    log.push_str(&format!("Total input pairs of read: {}\n", count_input_pairs));
    log.push_str(&format!("Total pairs of read used: {}\n", count_used_pairs));
    log.push_str(&format!("Total contigs: {}\n", count_all_contig));
    if min > 0 {
        for i in 0..=(error as usize) {
            log.push_str(&format!(
                "{:2} mismatch: {:9} pairs of reads => {:9} contigs\n",
                i, count_pairs[i], count_contigs[i]
            ));
        }
        log.push_str(&format!("Overlap range {}-{}\n", omin, omax));
    }

    LinkerOutput { output, log }
}

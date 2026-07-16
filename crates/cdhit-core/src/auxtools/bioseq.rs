//! Sequence I/O shared by the `cd-hit-auxtools` programs (`cd-hit-dup`,
//! `cd-hit-lap`, `read-linker`).
//!
//! Port of `cd-hit-auxtools/bioSequence.{cxx,hxx}`. Unlike the clustering
//! reader in [`crate::io::reader`], these tools retain the FASTQ quality string
//! and re-emit FASTA/FASTQ records verbatim, so a dedicated `Sequence` type is
//! used here. The parser reproduces the quirks of `SequenceList::ReadFastAQ` /
//! `SequenceCache::Update` required for bit-for-bit output parity — in
//! particular the "upper-case every record except the file's last" behaviour
//! (the C++ omits `ToUpper()` on the final flush).

/// DNA reverse-complement table, indexed by `byte - b'A'`. Verbatim copy of
/// `dna_reverse_complimentary` in `bioSequence.cxx` (`T B G D E F C H I J K L M
/// N O P Q R S A A V W X Y Z`).
const DNA_RC: &[u8; 26] = b"TBGDEFCHIJKLMNOPQRSAAVWXYZ";

/// A biological sequence with an optional FASTQ quality string.
#[derive(Clone, Default)]
pub struct Sequence {
    /// Numeric id, used by `Print` only when the description is empty. The C++
    /// never assigns it during parsing, so it stays 0 for parsed records.
    pub id: u32,
    /// Description (header line with the leading `>`/`@` marker and trailing
    /// newline removed).
    pub des: Vec<u8>,
    /// Sequence residues.
    pub seq: Vec<u8>,
    /// Quality scores (empty for FASTA).
    pub qs: Vec<u8>,
}

impl Sequence {
    pub fn new() -> Self {
        Sequence::default()
    }

    pub fn len(&self) -> usize {
        self.seq.len()
    }

    pub fn is_empty(&self) -> bool {
        self.seq.is_empty()
    }

    /// Port of `Sequence::GetDescription`. Returns the identifier token: the
    /// description trimmed to the first whitespace, after skipping an optional
    /// leading `>`/`@`/`+` marker and a following space/tab, bounded by
    /// `deslen` (0 = full length).
    pub fn get_description(&self, deslen: usize) -> Vec<u8> {
        let des = &self.des;
        let n = des.len();
        let mut i = 0usize;
        if i < n && (des[i] == b'>' || des[i] == b'@' || des[i] == b'+') {
            i += 1;
        }
        if i < n && (des[i] == b' ' || des[i] == b'\t') {
            i += 1;
        }
        let mut deslen = deslen;
        if deslen == 0 || deslen > n {
            deslen = n;
        }
        while i < deslen && !is_space(des[i]) {
            i += 1;
        }
        des[..i].to_vec()
    }

    /// Port of `Sequence::ToReverseComplimentary`. Reverse-complements `seq`
    /// (using [`DNA_RC`]) and reverses `qs` when its length matches `seq`.
    pub fn to_reverse_complement(&mut self) {
        let n = self.seq.len();
        let m = n / 2 + n % 2;
        for i in 0..m {
            let j = n - 1 - i;
            let l = self.seq[i];
            let r = self.seq[j];
            self.seq[i] = rc(r);
            self.seq[j] = rc(l);
        }
        if self.qs.len() != n {
            return;
        }
        for i in 0..m {
            let j = n - 1 - i;
            self.qs.swap(i, j);
        }
    }

    /// Port of `Sequence::Print`. Appends a FASTQ record when a quality string
    /// is present, otherwise a FASTA record. The `width` argument of the C++ is
    /// unused there and omitted here.
    pub fn print(&self, out: &mut Vec<u8>) {
        if !self.qs.is_empty() {
            debug_assert_eq!(self.seq.len(), self.qs.len());
            out.push(b'@');
            self.push_id_or_des(out);
            out.push(b'\n');
            out.extend_from_slice(&self.seq);
            out.push(b'\n');
            out.push(b'+');
            self.push_id_or_des(out);
            out.push(b'\n');
            out.extend_from_slice(&self.qs);
            out.push(b'\n');
        } else {
            out.push(b'>');
            self.push_id_or_des(out);
            out.push(b'\n');
            out.extend_from_slice(&self.seq);
            out.push(b'\n');
        }
    }

    fn push_id_or_des(&self, out: &mut Vec<u8>) {
        if !self.des.is_empty() {
            out.extend_from_slice(&self.des);
        } else {
            out.extend_from_slice(self.id.to_string().as_bytes());
        }
    }
}

fn rc(c: u8) -> u8 {
    if (b'A'..=b'Z').contains(&c) {
        DNA_RC[(c - b'A') as usize]
    } else {
        c // non A-Z passes through (C++ would read out of table bounds)
    }
}

/// `isspace` for the C locale bytes we care about.
fn is_space(c: u8) -> bool {
    matches!(c, b' ' | b'\t' | b'\n' | b'\r' | 0x0b | 0x0c)
}

fn to_upper(seq: &mut [u8]) {
    for b in seq.iter_mut() {
        b.make_ascii_uppercase();
    }
}

/// Iterate logical lines, each slice including its trailing `\n` when present.
struct Lines<'a> {
    data: &'a [u8],
    pos: usize,
}

impl<'a> Iterator for Lines<'a> {
    type Item = &'a [u8];
    fn next(&mut self) -> Option<Self::Item> {
        if self.pos >= self.data.len() {
            return None;
        }
        let start = self.pos;
        match self.data[start..].iter().position(|&b| b == b'\n') {
            Some(rel) => {
                let end = start + rel + 1;
                self.pos = end;
                Some(&self.data[start..end])
            }
            None => {
                self.pos = self.data.len();
                Some(&self.data[start..])
            }
        }
    }
}

#[derive(Clone, Copy, PartialEq)]
enum Target {
    Seq,
    Qual,
}

/// Parse an in-memory FASTA/FASTQ buffer into a list of sequences, reproducing
/// `SequenceList::ReadFastAQ` / `SequenceCache::Update`.
///
/// `idonly` truncates each description at its first whitespace during parsing
/// (the `-idonly` behaviour). Every record's residues are upper-cased *except
/// the last one in the buffer*, matching the C++ which omits `ToUpper()` on the
/// final flush.
pub fn parse_fastaq(data: &[u8], idonly: bool) -> Vec<Sequence> {
    let mut sequences: Vec<Sequence> = Vec::new();
    let mut getdes = true;
    let mut target = Target::Seq;
    let mut cur = Sequence::new();

    for line in (Lines { data, pos: 0 }) {
        let first = line.first().copied();
        if getdes && matches!(first, Some(b'>') | Some(b'@') | Some(b'+')) {
            getdes = false;
            if first == Some(b'+') {
                target = Target::Qual;
            } else {
                if !cur.seq.is_empty() {
                    to_upper(&mut cur.seq); // upper-cased on flush (not the last record)
                    sequences.push(std::mem::take(&mut cur));
                }
                target = Target::Seq;
                // description = line without the leading marker, minus the
                // trailing '\n' (a single Chop in the C++).
                let mut ds = line[1..].to_vec();
                if ds.last() == Some(&b'\n') {
                    ds.pop();
                }
                if idonly {
                    let mut d = 0usize;
                    while d < ds.len() && !is_space(ds[d]) {
                        d += 1;
                    }
                    ds.truncate(d);
                }
                cur.des = ds;
            }
        } else {
            let buf = match target {
                Target::Seq => &mut cur.seq,
                Target::Qual => &mut cur.qs,
            };
            buf.extend_from_slice(line);
            while buf.last().map_or(false, |&c| is_space(c)) {
                buf.pop();
            }
            getdes =
                target == Target::Seq || (target == Target::Qual && cur.qs.len() == cur.seq.len());
        }
    }
    // Final flush — note: NOT upper-cased, matching the C++.
    if !cur.seq.is_empty() {
        sequences.push(cur);
    }
    sequences
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_fasta_uppercases_all_but_last() {
        let fa = b">a\nacgt\n>b\nttgca\n";
        let seqs = parse_fastaq(fa, false);
        assert_eq!(seqs.len(), 2);
        assert_eq!(seqs[0].des, b"a");
        assert_eq!(seqs[0].seq, b"ACGT"); // upper-cased on flush
        assert_eq!(seqs[1].des, b"b");
        assert_eq!(seqs[1].seq, b"ttgca"); // last record NOT upper-cased (C++ quirk)
    }

    #[test]
    fn parses_fastq_with_quality() {
        let fq = b"@r1 desc\nACGT\n+r1 desc\n!!!!\n@r2\nTTGC\n+\n####\n";
        let seqs = parse_fastaq(fq, false);
        assert_eq!(seqs.len(), 2);
        assert_eq!(seqs[0].des, b"r1 desc");
        assert_eq!(seqs[0].seq, b"ACGT");
        assert_eq!(seqs[0].qs, b"!!!!");
        assert_eq!(seqs[1].qs, b"####");
    }

    #[test]
    fn fastq_quality_line_starting_with_at() {
        // A quality line may legitimately start with '@'; getdes must be false
        // until qs.len() == seq.len().
        let fq = b"@r1\nACGT\n+\n@@@@\n";
        let seqs = parse_fastaq(fq, false);
        assert_eq!(seqs.len(), 1);
        assert_eq!(seqs[0].seq, b"ACGT");
        assert_eq!(seqs[0].qs, b"@@@@");
    }

    #[test]
    fn get_description_trims_at_whitespace() {
        let s = Sequence {
            des: b"read1 some comment".to_vec(),
            ..Default::default()
        };
        assert_eq!(s.get_description(0), b"read1");
    }

    #[test]
    fn reverse_complement_dna() {
        let mut s = Sequence {
            seq: b"ACGTN".to_vec(),
            ..Default::default()
        };
        s.to_reverse_complement();
        // revcomp of ACGTN -> N A C G T
        assert_eq!(s.seq, b"NACGT");
    }

    #[test]
    fn reverse_complement_reverses_quality() {
        let mut s = Sequence {
            seq: b"ACGT".to_vec(),
            qs: b"1234".to_vec(),
            ..Default::default()
        };
        s.to_reverse_complement();
        assert_eq!(s.seq, b"ACGT".iter().rev().map(|&c| rc(c)).collect::<Vec<_>>());
        assert_eq!(s.qs, b"4321");
    }
}

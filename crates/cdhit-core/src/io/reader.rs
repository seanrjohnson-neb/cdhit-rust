//! FASTA / FASTQ reader.
//!
//! Port of `SequenceDB::Read` (cdhit-common.c++:1784-1878). Reads the whole
//! input into memory and parses it line by line, reproducing the C++ record
//! flushing, byte-offset (`des_begin`) and byte-length (`tot_length`)
//! accounting, and identifier-token extraction exactly. Keeping the raw input
//! is required because representative FASTA output re-copies the original record
//! bytes verbatim (see `writer.rs`).
//!
//! The C++ streams via `fgets` with a fixed buffer, splitting over-long lines
//! into chunks; because every byte still flows through the same accumulation
//! and byte counters, processing whole logical lines is behaviourally identical.

use crate::options::Options;
use crate::sequence::{c_isspace, Sequence};

/// Result of reading a database: the parsed sequences plus the original raw
/// bytes (needed to reconstruct representative records for `-o` output).
pub struct Database {
    pub sequences: Vec<Sequence>,
    pub raw: Vec<u8>,
}

/// Iterate logical lines, each slice including its trailing `\n` when present.
/// Yields `(start_offset, line_bytes)`.
struct Lines<'a> {
    data: &'a [u8],
    pos: usize,
}

impl<'a> Iterator for Lines<'a> {
    type Item = (usize, &'a [u8]);
    fn next(&mut self) -> Option<Self::Item> {
        if self.pos >= self.data.len() {
            return None;
        }
        let start = self.pos;
        match self.data[start..].iter().position(|&b| b == b'\n') {
            Some(rel) => {
                let end = start + rel + 1; // include the '\n'
                self.pos = end;
                Some((start, &self.data[start..end]))
            }
            None => {
                self.pos = self.data.len();
                Some((start, &self.data[start..]))
            }
        }
    }
}

/// Extract the identifier token from a header line (port of the logic at
/// cdhit-common.c++:1857-1863). Returns the token bytes including the leading
/// `>`/`@`/`+` marker.
fn parse_identifier(line: &[u8], des_len: i32) -> Vec<u8> {
    let full = line.len();
    let mut i = 0usize;
    if i < full && (line[i] == b'>' || line[i] == b'@' || line[i] == b'+') {
        i += 1;
    }
    if i < full && (line[i] == b' ' || line[i] == b'\t') {
        i += 1;
    }
    let mut des_size = full;
    if des_len > 0 && (des_len as usize) < des_size {
        des_size = des_len as usize;
    }
    while i < des_size && !c_isspace(line[i]) {
        i += 1;
    }
    line[..i].to_vec()
}

/// Parse an in-memory FASTA/FASTQ database.
pub fn read_database(input: &[u8], options: &Options) -> Database {
    let option_l = options.min_length;
    let mut sequences: Vec<Sequence> = Vec::new();

    let mut one = Sequence::new();
    let mut have_ident = false; // whether `one` has a header (identifier) yet
    let mut has_data = false; // whether any data has been accumulated

    let mut lines = Lines { data: input, pos: 0 }.peekable();

    // Flush the pending record `one` (port of cdhit-common.c++:1826-1839).
    macro_rules! flush_pending {
        () => {
            if one.size > 0 {
                let invalid = !have_ident || one.format() != 0;
                if invalid {
                    // "Discarding invalid sequence ..." (warning omitted here).
                    one.size = 0;
                }
                one.index = sequences.len() as i32;
                if one.size > option_l {
                    if options.trim_len > 0 {
                        one.trim(options.trim_len);
                    }
                    sequences.push(one.clone());
                }
            }
        };
    }

    while let Some((start, line)) = lines.next() {
        let first = line.first().copied().unwrap_or(b'>');
        if first == b'+' {
            // FASTQ separator: this line plus the following quality line belong
            // to the current record (cdhit-common.c++:1805-1824).
            one.tot_length += line.len() as i32;
            if let Some((_, qual)) = lines.next() {
                one.tot_length += qual.len() as i32;
            }
        } else if first == b'>' || first == b'@' {
            // Flush previous record, then begin a new one.
            flush_pending!();
            one = Sequence::new();
            has_data = false;

            one.des_begin = start as u64;
            one.tot_length += line.len() as i32;
            one.identifier = parse_identifier(line, options.des_len);
            have_ident = true;
        } else {
            // Sequence data line; keep raw bytes (newlines stripped in Format).
            one.tot_length += line.len() as i32;
            one.data.extend_from_slice(line);
            one.size = one.data.len() as i32;
            has_data = true;
        }
    }
    // Flush the final record.
    let _ = has_data;
    flush_pending!();

    Database {
        sequences,
        raw: input.to_vec(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_simple_fasta() {
        let fa = b">seq1 description here\nACDEFGHIKLMNPQRSTVWY\n>seq2\nACDEFGHIKL\n";
        let mut opt = Options::default();
        opt.min_length = 5;
        let db = read_database(fa, &opt);
        assert_eq!(db.sequences.len(), 2);
        assert_eq!(db.sequences[0].identifier, b">seq1".to_vec());
        assert_eq!(db.sequences[0].size, 20);
        assert_eq!(db.sequences[0].index, 0);
        assert_eq!(db.sequences[1].identifier, b">seq2".to_vec());
        assert_eq!(db.sequences[1].size, 10);
        assert_eq!(db.sequences[1].index, 1);
    }

    #[test]
    fn min_length_filter_is_strict_greater() {
        // size must be > option_l; a length-10 seq with option_l=10 is dropped.
        let fa = b">a\nACDEFGHIKL\n>b\nACDEFGHIKLM\n";
        let mut opt = Options::default();
        opt.min_length = 10;
        let db = read_database(fa, &opt);
        assert_eq!(db.sequences.len(), 1);
        assert_eq!(db.sequences[0].identifier, b">b".to_vec());
        // index reflects the position it lands at (0), since `a` was skipped.
        assert_eq!(db.sequences[0].index, 0);
    }

    #[test]
    fn des_len_truncates_identifier() {
        // Long ID with no whitespace, des_len default 20 -> truncated to 20 bytes.
        let id = "A".repeat(40);
        let fa = format!(">{}\nACDEFGHIKLMNPQRST\n", id);
        let opt = Options::default();
        let db = read_database(fa.as_bytes(), &opt);
        assert_eq!(db.sequences.len(), 1);
        assert_eq!(db.sequences[0].identifier.len(), 20);
        assert_eq!(&db.sequences[0].identifier[..1], b">");
    }

    #[test]
    fn des_begin_and_tot_length_cover_record() {
        let fa = b">a\nACDEFGHIKLM\n>b\nACDEFGHIKLMN\n";
        let opt = Options::default();
        let db = read_database(fa, &opt);
        // Record 0 is ">a\n" (3) + "ACDEFGHIKLM\n" (12) = 15 bytes.
        assert_eq!(db.sequences[0].des_begin, 0);
        assert_eq!(db.sequences[0].tot_length, 15);
        // Record 1 starts right after.
        assert_eq!(db.sequences[1].des_begin, 15);
    }

    #[test]
    fn invalid_chars_discard_sequence() {
        // A digit in the sequence makes Format reject it.
        let fa = b">bad\nACDE123FGHIKLM\n>good\nACDEFGHIKLMN\n";
        let opt = Options::default();
        let db = read_database(fa, &opt);
        assert_eq!(db.sequences.len(), 1);
        assert_eq!(db.sequences[0].identifier, b">good".to_vec());
    }
}

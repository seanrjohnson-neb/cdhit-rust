//! `Sequence` — one database entry.
//!
//! Port of the `Sequence` struct (cdhit-common.h:374-429) and its methods
//! (cdhit-common.c++:1480-1682). The raw `char *data` becomes a `Vec<u8>` that
//! is encoded in place to the small integer alphabet by [`Sequence::convert_bases`]
//! (called from `SortDivide`, matching the C++). The identifier retains its
//! leading `>`/`@`, exactly as the C++ stores it (`PrintInfo` does `identifier+1`).

/// C-locale `isspace`: space, tab, newline, vertical tab, form feed, carriage
/// return. Rust's `is_ascii_whitespace` omits vertical tab (0x0B), so we match
/// C explicitly for parsing fidelity.
#[inline]
pub fn c_isspace(b: u8) -> bool {
    matches!(b, b' ' | b'\t' | b'\n' | 0x0B | 0x0C | b'\r')
}

// Sequence state bitflags (cdhit-common.h:67-70).
pub const IS_REP: i16 = 1;
pub const IS_REDUNDANT: i16 = 2;
pub const IS_PROCESSED: i16 = 16;
pub const IS_MINUS_STRAND: i16 = 32;

#[derive(Clone)]
pub struct Sequence {
    /// Sequence letters. After `Read`+`Format`, uppercase ASCII; after
    /// `convert_bases`, small integer codes.
    pub data: Vec<u8>,
    /// `size` in the C++ (kept explicit because it can differ from `data.len()`
    /// only transiently; we keep them in sync).
    pub size: i32,
    /// For back-to-back merged PE reads: length of the R2 portion.
    pub size_r2: i32,

    /// Byte offset of the record start in the source file (for re-reading the
    /// original record when writing representative FASTA).
    pub des_begin: u64,
    pub des_begin2: u64,
    /// Total byte length of the source record.
    pub tot_length: i32,
    pub tot_length2: i32,

    /// Identifier token, retaining the leading `>`/`@` marker. Stored as raw
    /// bytes (not `String`) so output is byte-identical for non-UTF-8 headers.
    pub identifier: Vec<u8>,

    /// Index of the sequence in the original (pre-sort) database.
    pub index: i32,
    pub state: i16,
    pub cluster_id: i32,
    pub identity: f32,
    pub distance: f32,
    pub coverage: [i32; 4],
}

impl Default for Sequence {
    fn default() -> Self {
        // Mirrors `Sequence()` (memset 0, distance = 2.0).
        Sequence {
            data: Vec::new(),
            size: 0,
            size_r2: 0,
            des_begin: 0,
            des_begin2: 0,
            tot_length: 0,
            tot_length2: 0,
            identifier: Vec::new(),
            index: 0,
            state: 0,
            // C++ `Sequence()` memsets to 0, so cluster_id starts at 0 (the EST
            // tie-breaks compare against it before assignment).
            cluster_id: 0,
            identity: 0.0,
            distance: 2.0,
            coverage: [0; 4],
        }
    }
}

impl Sequence {
    pub fn new() -> Self {
        Sequence::default()
    }

    /// Port of `Sequence::Format` (cdhit-common.c++:1628-1646).
    ///
    /// Trims trailing whitespace and a single trailing `*`, then verifies every
    /// remaining char is alphabetic or whitespace. Returns the count of invalid
    /// characters (0 = valid). On success, keeps only uppercased alphabetics.
    pub fn format(&mut self) -> i32 {
        let mut size = self.size as usize;
        while size > 0 && c_isspace(self.data[size - 1]) {
            size -= 1;
        }
        if size > 0 && self.data[size - 1] == b'*' {
            size -= 1;
        }
        self.data.truncate(size);
        self.size = size as i32;

        let mut m = 0i32;
        for i in 0..size {
            let ch = self.data[i];
            // C: m += !(isalpha(ch) | isspace(ch))
            let ok = ch.is_ascii_alphabetic() || c_isspace(ch);
            if !ok {
                m += 1;
            }
        }
        if m != 0 {
            return m;
        }
        let mut j = 0usize;
        for i in 0..size {
            let ch = self.data[i];
            if ch.is_ascii_alphabetic() {
                self.data[j] = ch.to_ascii_uppercase();
                j += 1;
            }
        }
        self.data.truncate(j);
        self.size = j as i32;
        0
    }

    /// Port of `Sequence::ConvertBases` (cdhit-common.c++:1614-1618): encode each
    /// letter in place through the active `aa2idx` table.
    pub fn convert_bases(&mut self, aa2idx: &[i32; 26]) {
        for i in 0..self.size as usize {
            let idx = (self.data[i] - b'A') as usize;
            self.data[i] = aa2idx[idx] as u8;
        }
    }

    /// Port of `Sequence::trim` (cdhit-common.c++:1609-1613).
    pub fn trim(&mut self, trim_len: i32) {
        if trim_len >= self.size {
            return;
        }
        self.size = trim_len;
        self.data.truncate(trim_len as usize);
    }
}

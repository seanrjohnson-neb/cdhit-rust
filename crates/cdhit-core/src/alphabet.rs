//! Alphabet encoding tables, score matrices, and word-count power tables.
//!
//! Everything here is a verbatim port of the corresponding constants and helper
//! routines in the original `cdhit-common.c++`. The values must match exactly,
//! because they feed the integer-alphabet encoding, the banded DP scoring, and
//! the k-mer word encoding — all of which are load-bearing for bit-for-bit
//! output fidelity.

// Compile-time size constants (cdhit-common.h:48-62).
pub const MAX_SEQ: i32 = 655360;
pub const MAX_AA: usize = 23;
pub const MAX_NA: usize = 6;
pub const MAX_UAA: i32 = 21;

/// Ordered residue alphabet: index -> character (cdhit-common.c++:50).
pub const AA: &[u8; MAX_AA] = b"ARNDCQEGHILKMFPSTWYVBZX";

/// `aa2idx[c - 'A']` -> residue index, for protein sequences
/// (cdhit-common.c++:52-53).
pub const AA2IDX: [i32; 26] = [
    0, 2, 4, 3, 6, 13, 7, 8, 9, 20, 11, 10, 12, 2, 20, 14, 5, 1, 15, 16, 20, 19, 17, 20, 18, 6,
];

/// `na2idx[c - 'A']` -> base index, for nucleotide sequences
/// (cdhit-common.c++:88-89).
pub const NA2IDX: [i32; 26] = [
    0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4,
];

/// BLOSUM62, packed lower-triangular in `AA` order (cdhit-common.c++:59-85).
#[rustfmt::skip]
pub const BLOSUM62: [i32; 276] = [
     4,
    -1, 5,
    -2, 0, 6,
    -2,-2, 1, 6,
     0,-3,-3,-3, 9,
    -1, 1, 0, 0,-3, 5,
    -1, 0, 0, 2,-4, 2, 5,
     0,-2, 0,-1,-3,-2,-2, 6,
    -2, 0, 1,-1,-3, 0, 0,-2, 8,
    -1,-3,-3,-3,-1,-3,-3,-4,-3, 4,
    -1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,
    -1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,
    -1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5,
    -2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,
    -1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,
     1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4,
     0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,
    -3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11,
    -2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,
     0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4,
    -2,-1, 3, 4,-3, 0, 1,-1, 0,-3,-4, 0,-3,-3,-2, 0,-1,-4,-3,-3, 4,
    -1, 0, 0, 1,-3, 3, 4,-2, 0,-3,-3, 1,-1,-3,-1, 0,-1,-3,-2,-2, 1, 4,
     0,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-2, 0, 0,-2,-1,-1,-1,-1,-1,
];

/// Nucleotide scoring matrix, packed lower-triangular (cdhit-common.c++:94-104).
#[rustfmt::skip]
pub const BLOSUM62_NA: [i32; 28] = [
     2,
    -2, 2,
    -2,-2, 2,
    -2,-2,-2, 2,
    -2,-2,-2, 1, 2,
    -2,-2,-2,-2,-2, 1,
     0, 0, 0, 0, 0, 0, 1,
];

/// Score matrix used by the banded DP. Mirrors `ScoreMatrix` in
/// `cdhit-common.c++`: every stored value is the raw score multiplied by
/// `MAX_SEQ`, and gaps are likewise pre-scaled.
#[derive(Clone)]
pub struct ScoreMatrix {
    pub matrix: [[i32; MAX_AA]; MAX_AA],
    pub gap: i32,
    pub ext_gap: i32,
}

impl Default for ScoreMatrix {
    fn default() -> Self {
        let mut m = ScoreMatrix {
            matrix: [[0; MAX_AA]; MAX_AA],
            gap: 0,
            ext_gap: 0,
        };
        m.init();
        m
    }
}

impl ScoreMatrix {
    /// Protein defaults: gap -11/-1, BLOSUM62 (cdhit-common.c++:1210-1214).
    pub fn init(&mut self) {
        self.set_gap(-11, -1);
        self.set_matrix(&BLOSUM62);
    }

    /// `gap = MAX_SEQ * gap1`, `ext_gap = MAX_SEQ * ext_gap1`
    /// (cdhit-common.c++:1216-1221).
    pub fn set_gap(&mut self, gap1: i32, ext_gap1: i32) {
        self.gap = MAX_SEQ * gap1;
        self.ext_gap = MAX_SEQ * ext_gap1;
    }

    /// Unpack a lower-triangular matrix into the full square, scaled by MAX_SEQ
    /// (cdhit-common.c++:1223-1230).
    ///
    /// The C++ loop always runs `MAX_AA` rows regardless of the input length,
    /// reading past the 28-element NA matrix into adjacent globals for the
    /// unused rows 7..22. Those cells are never touched during NA alignment (NA
    /// codes are 0..5), so we simply stop once the input is exhausted, leaving
    /// the higher rows at their prior (protein) values.
    pub fn set_matrix(&mut self, mat1: &[i32]) {
        let mut k = 0usize;
        'outer: for i in 0..MAX_AA {
            for j in 0..=i {
                if k >= mat1.len() {
                    break 'outer;
                }
                let v = MAX_SEQ * mat1[k];
                self.matrix[j][i] = v;
                self.matrix[i][j] = v;
                k += 1;
            }
        }
    }

    /// Nucleotide defaults: gap -6/-1, BLOSUM62_na (cdhit-common.c++:1232-1236).
    pub fn set_to_na(&mut self) {
        self.set_gap(-6, -1);
        self.set_matrix(&BLOSUM62_NA);
    }

    /// EST match score, diagonal entries 0..5 (cdhit-common.c++:1238-1243).
    pub fn set_match(&mut self, score: i32) {
        for i in 0..5 {
            self.matrix[i][i] = MAX_SEQ * score;
        }
    }

    /// EST mismatch score, off-diagonal entries; keep T/U (3/4) as a match
    /// (cdhit-common.c++:1245-1252).
    pub fn set_mismatch(&mut self, score: i32) {
        for i in 0..MAX_AA {
            for j in 0..i {
                let v = MAX_SEQ * score;
                self.matrix[j][i] = v;
                self.matrix[i][j] = v;
            }
        }
        self.matrix[3][4] = MAX_SEQ;
        self.matrix[4][3] = MAX_SEQ;
    }
}

/// Word-count power table: `naan_array[k] = NAA^k`, plus the individual `NAAk`
/// scalars the C++ keeps as globals (cdhit-common.c++:171-199).
///
/// `NAAN_array[0]` stays 1 (base case) exactly as the C++ initializer sets it.
#[derive(Clone, Copy)]
pub struct Naa {
    pub array: [i32; 13],
}

impl Naa {
    /// Port of `InitNAA(max)`. Uses wrapping multiplication because the C++
    /// computes all 12 powers as `int` even when the high ones overflow (e.g.
    /// base 21 for protein); those entries wrap but are never indexed
    /// (`NAA <= NAA_top_limit`), so only the wrap-free low powers are observed.
    pub fn init(max: i32) -> Naa {
        let mut a = [0i32; 13];
        a[0] = 1;
        a[1] = max; // NAA1
        a[2] = a[1].wrapping_mul(a[1]); // NAA2
        a[3] = a[1].wrapping_mul(a[2]); // NAA3
        a[4] = a[2].wrapping_mul(a[2]); // NAA4
        a[5] = a[2].wrapping_mul(a[3]); // NAA5
        a[6] = a[3].wrapping_mul(a[3]); // NAA6
        a[7] = a[3].wrapping_mul(a[4]); // NAA7
        a[8] = a[4].wrapping_mul(a[4]); // NAA8
        a[9] = a[4].wrapping_mul(a[5]); // NAA9
        a[10] = a[5].wrapping_mul(a[5]); // NAA10
        a[11] = a[5].wrapping_mul(a[6]); // NAA11
        a[12] = a[6].wrapping_mul(a[6]); // NAA12
        Naa { array: a }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn naan_array_matches_cpp() {
        // InitNAA(MAX_UAA) with MAX_UAA = 21.
        let naa = Naa::init(MAX_UAA);
        assert_eq!(naa.array[0], 1);
        assert_eq!(naa.array[1], 21);
        assert_eq!(naa.array[2], 441);
        assert_eq!(naa.array[3], 9261);
        assert_eq!(naa.array[5], 4084101);
    }

    #[test]
    fn blosum62_is_lower_triangular_276() {
        // 23 rows, triangular count = 23*24/2 = 276.
        assert_eq!(BLOSUM62.len(), 23 * 24 / 2);
    }

    #[test]
    fn score_matrix_symmetric_and_scaled() {
        let m = ScoreMatrix::default();
        // A-A (index 0) BLOSUM62 = 4, scaled by MAX_SEQ.
        assert_eq!(m.matrix[0][0], 4 * MAX_SEQ);
        // W-W (index 17) = 11.
        assert_eq!(m.matrix[17][17], 11 * MAX_SEQ);
        // symmetry
        assert_eq!(m.matrix[0][5], m.matrix[5][0]);
        assert_eq!(m.gap, -11 * MAX_SEQ);
        assert_eq!(m.ext_gap, -1 * MAX_SEQ);
    }

    #[test]
    fn na_matrix() {
        let mut m = ScoreMatrix::default();
        m.set_to_na();
        assert_eq!(m.matrix[0][0], 2 * MAX_SEQ);
        assert_eq!(m.gap, -6 * MAX_SEQ);
    }
}

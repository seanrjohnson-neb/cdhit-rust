//! Nucleotide (EST) helpers: reverse-complement sequence and word-index maps.
//!
//! Port of `make_comp_iseq` (cdhit-common.c++:205-209) and
//! `make_comp_short_word_index` (:3996-4017).

use crate::alphabet::Naa;

/// Reverse-complement an encoded nucleotide sequence into `out`
/// (`make_comp_iseq`). Complement map: A(0)->T(3), C(1)->G(2), G(2)->C(1),
/// T(3)->A(0), N(4)->4, mask(5)->5; then reversed.
pub fn make_comp_iseq(len: usize, out: &mut [u8], iseq: &[u8]) {
    let c = [3u8, 2, 1, 0, 4, 5];
    for i in 0..len {
        out[i] = c[iseq[len - i - 1] as usize];
    }
}

/// Build the reverse-complement word-index table `Comp_AAN_idx` of length
/// `NAAN` (`make_comp_short_word_index`): for each encoded word, the encoding of
/// its reverse complement.
pub fn make_comp_short_word_index(naa: i32, naan: &Naa) -> Vec<i32> {
    let c = [3usize, 2, 1, 0];
    let naa1 = naan.array[1];
    let naan_len = naan.array[naa as usize];
    let mut comp = vec![0i32; naan_len.max(0) as usize];
    let mut short_word = [0u8; 32];
    for i in 0..naan_len {
        // Decompose i back into its NAA base-NAA1 digits.
        let mut k = i;
        for j in 0..naa as usize {
            short_word[j] = (k % naa1) as u8;
            k /= naa1;
        }
        // Encode the complement (digits reversed and complemented).
        let mut icomp = 0i32;
        let mut k1 = (naa - 1) as usize;
        for kk in 0..naa as usize {
            icomp += (c[short_word[k1] as usize] as i32) * naan.array[kk];
            if k1 > 0 {
                k1 -= 1;
            }
        }
        comp[i as usize] = icomp;
    }
    comp
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn revcomp_of_acgt() {
        // ACGT = [0,1,2,3] -> revcomp = ACGT again ([0,1,2,3])? No:
        // reverse = TGCA = [3,2,1,0]; complement each: T->A(0),G->C(1),C->G(2),A->T(3)
        // = [0,1,2,3]. ACGT is its own reverse complement.
        let seq = [0u8, 1, 2, 3];
        let mut out = [0u8; 4];
        make_comp_iseq(4, &mut out, &seq);
        assert_eq!(out, [0, 1, 2, 3]);
    }

    #[test]
    fn revcomp_simple() {
        // AAAC = [0,0,0,1] -> revcomp: reverse=CAAA=[1,0,0,0], complement:
        // C->G(2),A->T(3),A->T(3),A->T(3) = [2,3,3,3] = GTTT
        let seq = [0u8, 0, 0, 1];
        let mut out = [0u8; 4];
        make_comp_iseq(4, &mut out, &seq);
        assert_eq!(out, [2, 3, 3, 3]);
    }
}

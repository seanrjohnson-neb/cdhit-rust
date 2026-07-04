//! Diagonal (short-word) pre-filter.
//!
//! Port of `diag_test_aapn` (protein, cdhit-common.c++:517-611) and
//! `diag_test_aapn_est` (nucleotide, :614-749). Both accumulate amino-acid-pair
//! / 4-mer hits per alignment diagonal, then slide a `band_width`-wide window to
//! find the best band and its bounds. The rep's AAP index must already be built
//! in `buffer` via `ComputeAAP`/`ComputeAAP2`.

use crate::buffer::{WorkingBuffer, MAX_DIAG};
use crate::options::Options;

/// Output of a diagonal test.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct DiagResult {
    pub best_sum: i32,
    pub band_left: i32,
    pub band_center: i32,
    pub band_right: i32,
}

/// Shared band-selection logic (the identical tail of both diag_test variants).
/// `band_b`/`band_e` are the pre-clamped diagonal search bounds.
fn select_band(
    buffer: &WorkingBuffer,
    nall: i32,
    len1: i32,
    band_width: i32,
    band_b: i32,
    band_e: i32,
    cluster_thd: f64,
) -> DiagResult {
    let diag_score = &buffer.diag_score;
    let diag_score2 = &buffer.diag_score2;

    let band_m = if band_b + band_width - 1 < band_e {
        band_b + band_width - 1
    } else {
        band_e
    };
    let mut best_score = 0i32;
    let mut best_score2 = 0i32;
    let mut max_diag2 = 0i32;
    let mut imax_diag = 0i32;
    for i in band_b..=band_m {
        best_score += diag_score[i as usize];
        best_score2 += diag_score2[i as usize];
        if diag_score2[i as usize] > max_diag2 {
            max_diag2 = diag_score2[i as usize];
            imax_diag = i;
        }
    }
    let mut from = band_b;
    let mut end = band_m;
    let mut score = best_score;
    let mut score2 = best_score2;
    let mut k = from;
    let mut j = band_m + 1;
    while j < band_e {
        score -= diag_score[k as usize];
        score += diag_score[j as usize];
        score2 -= diag_score2[k as usize];
        score2 += diag_score2[j as usize];
        if score2 > best_score2 {
            from = k + 1;
            end = j;
            best_score = score;
            best_score2 = score2;
            if diag_score2[j as usize] > max_diag2 {
                max_diag2 = diag_score2[j as usize];
                imax_diag = j;
            }
        }
        j += 1;
        k += 1;
    }
    let mut mlen = imax_diag;
    if imax_diag > len1 {
        mlen = nall - imax_diag;
    }
    let emax = ((1.0 - cluster_thd) * mlen as f64) as i32 + 1;
    let mut j = from;
    while j < imax_diag {
        if (imax_diag - j) > emax || diag_score[j as usize] < 1 {
            best_score -= diag_score[j as usize];
            from += 1;
        } else {
            break;
        }
        j += 1;
    }
    let mut j = end;
    while j > imax_diag {
        if (j - imax_diag) > emax || diag_score[j as usize] < 1 {
            best_score -= diag_score[j as usize];
            end -= 1;
        } else {
            break;
        }
        j -= 1;
    }

    DiagResult {
        best_sum: best_score,
        band_left: from - len1 + 1,
        band_right: end - len1 + 1,
        band_center: imax_diag - len1 + 1,
    }
}

/// Protein diagonal test (`diag_test_aapn`).
pub fn diag_test_aapn(
    naa1: i32,
    iseq2: &[u8],
    len1: i32,
    len2: i32,
    buffer: &mut WorkingBuffer,
    band_width: i32,
    required_aa1: i32,
    options: &Options,
) -> DiagResult {
    let nall = len1 + len2 - 1;
    assert!(nall as usize <= MAX_DIAG, "in diag_test_aapn, MAX_DIAG reached");
    for i in 0..nall as usize {
        buffer.diag_score[i] = 0;
        buffer.diag_score2[i] = 0;
    }

    let len22 = (len2 - 1) as usize;
    let mut i1 = (len1 - 1) as i32;
    for i in 0..len22 {
        let c22 = (iseq2[i] as i32) * naa1 + iseq2[i + 1] as i32;
        let cpx = 1 + (iseq2[i] != iseq2[i + 1]) as i32;
        let j = buffer.taap[c22 as usize];
        if j == 0 {
            i1 += 1;
            continue;
        }
        let m = buffer.aap_begin[c22 as usize] as usize;
        for k in 0..j as usize {
            let d = (i1 - buffer.aap_list[m + k] as i32) as usize;
            buffer.diag_score[d] += 1;
            buffer.diag_score2[d] += cpx;
        }
        i1 += 1;
    }

    let band_b = if required_aa1 - 1 >= 0 { required_aa1 - 1 } else { 0 };
    let band_e = nall - band_b;
    select_band(buffer, nall, len1, band_width, band_b, band_e, options.cluster_thd)
}

/// Nucleotide diagonal test (`diag_test_aapn_est`).
pub fn diag_test_aapn_est(
    naa1: i32,
    iseq2: &[u8],
    len1: i32,
    len2: i32,
    buffer: &mut WorkingBuffer,
    band_width: i32,
    required_aa1: i32,
    options: &Options,
) -> DiagResult {
    let nall = len1 + len2 - 1;
    let naa2 = naa1 * naa1;
    let naa3 = naa2 * naa1;
    assert!(
        nall as usize <= MAX_DIAG,
        "in diag_test_aapn_est, MAX_DIAG reached"
    );
    for i in 0..nall as usize {
        buffer.diag_score[i] = 0;
        buffer.diag_score2[i] = 0;
    }

    let len22 = len2 - 3;
    let mut i1 = (len1 - 1) as i32;
    let mut base = 0usize; // iseq2 is advanced by 1 each iteration in the C++
    for _i in 0..len22 {
        let c0 = iseq2[base];
        let c1 = iseq2[base + 1];
        let c2 = iseq2[base + 2];
        let c3 = iseq2[base + 3];
        if c0 >= 4 || c1 >= 4 || c2 >= 4 || c3 >= 4 {
            i1 += 1;
            base += 1;
            continue;
        }
        let c22 =
            (c0 as i32) * naa3 + (c1 as i32) * naa2 + (c2 as i32) * naa1 + c3 as i32;
        let j = buffer.taap[c22 as usize];
        if j == 0 {
            i1 += 1;
            base += 1;
            continue;
        }
        let cpx = 1 + (c0 != c1) as i32 + (c1 != c2) as i32 + (c2 != c3) as i32;
        let m = buffer.aap_begin[c22 as usize] as usize;
        for k in 0..j as usize {
            let d = (i1 - buffer.aap_list[m + k] as i32) as usize;
            buffer.diag_score[d] += 1;
            buffer.diag_score2[d] += cpx;
        }
        i1 += 1;
        base += 1;
    }

    let mut band_b = if required_aa1 - 1 >= 0 { required_aa1 - 1 } else { 0 };
    let mut band_e = nall - band_b;
    if options.is454 {
        band_b = len1 - band_width;
        band_e = len1 + band_width;
        if band_b < 0 {
            band_b = 0;
        }
        if band_e > nall {
            band_e = nall;
        }
    }

    let mut result = select_band(buffer, nall, len1, band_width, band_b, band_e, options.cluster_thd);
    if options.is454 {
        if result.band_left > 0 {
            result.best_sum = 0;
        }
        if result.band_right < 0 {
            result.best_sum = 0;
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet::Naa;

    #[test]
    fn identical_protein_seq_best_band_on_main_diagonal() {
        let naa = Naa::init(21);
        let opt = Options::default();
        let mut buf = WorkingBuffer::new(8, 64, &opt);
        // rep == query, codes [0,1,2,3,4,5]
        let seq: Vec<u8> = vec![0, 1, 2, 3, 4, 5];
        let len = seq.len() as i32;
        buf.compute_aap(&seq, len, naa.array[1], naa.array[2]);
        let r = diag_test_aapn(naa.array[1], &seq, len, len, &mut buf, 20, 1, &opt);
        // Best band centered on the main diagonal (band_center == 0).
        assert_eq!(r.band_center, 0);
        // Every adjacent pair matches: 5 pairs on the main diagonal.
        assert!(r.best_sum >= 5);
    }
}

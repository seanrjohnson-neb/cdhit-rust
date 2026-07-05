//! Banded local/global DP alignment.
//!
//! Port of `local_band_align` (cdhit-common.c++:789-1193), excluding the
//! `#ifdef PRINT` / `#ifdef MAKEALIGN` debug blocks (not compiled by default).
//! This is the most fidelity-critical routine: the traceback tie-breaks
//! (`>` vs `>=`, the `DP_BACK_*` preference order) and the exact integer/float
//! arithmetic determine identity, alignment length, and distance.
//!
//! `iseq1` is the query, `iseq2` the representative. Scores live in `i64`
//! (matrix values are pre-scaled by `MAX_SEQ`); `best_score` is truncated to
//! `i32` on assignment exactly as the C++ `int` assignment does.

use crate::alphabet::ScoreMatrix;
use crate::buffer::WorkingBuffer;
use crate::options::Options;

// Back-pointer directions (cdhit-common.h:466).
const DP_BACK_NONE: i32 = 0;
const DP_BACK_LEFT_TOP: i32 = 1;
const DP_BACK_LEFT: i32 = 2;
const DP_BACK_TOP: i32 = 3;

/// Alignment result (the out-parameters of `local_band_align`).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct BandAlign {
    pub best_score: i32,
    pub iden_no: i32,
    pub alnln: i32,
    pub dist: f32,
    pub alninfo: [i32; 5],
}

/// Run the banded alignment. Returns `None` for `FAILED_FUNC` (invalid band,
/// 454 indel overflow, or unmatched-region limits exceeded). `want_alninfo`
/// mirrors passing a non-null `alninfo` pointer in the C++ (it gates the
/// global-identity begin/end capture).
#[allow(clippy::too_many_arguments)]
pub fn local_band_align(
    iseq1: &[u8],
    iseq2: &[u8],
    len1: i32,
    len2: i32,
    mat: &ScoreMatrix,
    band_left: i32,
    band_center: i32,
    band_right: i32,
    options: &Options,
    want_alninfo: bool,
    buffer: &mut WorkingBuffer,
) -> Option<BandAlign> {
    if band_right >= len2 || band_left <= -len1 || band_left > band_right {
        return None;
    }

    let band_width = band_right - band_left + 1;
    let band_width1 = band_width + 1;

    // score_mat / back_mat [i][j1]: reused per-thread buffer (flat rows*w),
    // grown on demand. The C++ likewise reuses a per-thread buffer; the
    // recurrence only ever reads cells it wrote this call, so no zeroing between
    // calls is required (verified: this matches the previous fresh-zeroed
    // implementation, and the C++, bit-for-bit).
    let rows = (len1 + 1) as usize;
    let w = band_width1 as usize;
    let needed = rows * w;
    if buffer.score_mat.len() < needed {
        buffer.score_mat.resize(needed, 0);
        buffer.back_mat.resize(needed, 0);
    }
    let score_mat = &mut buffer.score_mat;
    let back_mat = &mut buffer.back_mat;
    let idx = |i: i32, j1: i32| (i as usize) * w + (j1 as usize);

    // Left border (band_left < 0): leading query hanging residues.
    if band_left < 0 {
        let tband = if band_right < 0 { band_right } else { 0 };
        let mut k = band_left;
        while k <= tband {
            let i = -k;
            let j1 = k - band_left;
            score_mat[idx(i, j1)] = (mat.ext_gap as i64) * i as i64;
            back_mat[idx(i, j1)] = DP_BACK_TOP;
            k += 1;
        }
        back_mat[idx(-tband, tband - band_left)] = DP_BACK_NONE;
    }

    // Top border (band_right >= 0).
    if band_right >= 0 {
        let tband = if band_left > 0 { band_left } else { 0 };
        for j in tband..=band_right {
            let j1 = j - band_left;
            score_mat[idx(0, j1)] = (mat.ext_gap as i64) * j as i64;
            back_mat[idx(0, j1)] = DP_BACK_LEFT;
        }
        back_mat[idx(0, tband - band_left)] = DP_BACK_NONE;
    }

    let gap_open = [mat.gap, mat.ext_gap];
    let max_diag = band_center - band_left;
    let extra_score = [4i32, 3, 2, 1];

    for i in 1..=len1 {
        let mut j0 = 1 - band_left - i;
        let mut j1_hi = len2 - band_left - i;
        if j0 < 0 {
            j0 = 0;
        }
        if j1_hi >= band_width {
            j1_hi = band_width;
        }
        // Hoist per-row values out of the inner loop and use unchecked indexing
        // in the hot DP recurrence (all indices are provably in bounds: codes
        // are < MAX_AA, and flat offsets are < rows*w). Arithmetic is identical
        // to the checked version above/below.
        let row_i = (i as usize) * w;
        let row_im1 = ((i - 1) as usize) * w;
        let ci = iseq1[(i - 1) as usize] as usize;
        let ci_row = unsafe { mat.matrix.get_unchecked(ci) };
        let at_last_row = i == len1;
        for j1 in j0..=j1_hi {
            let j = j1 + i + band_left;
            let j1u = j1 as usize;

            let cj = unsafe { *iseq2.get_unchecked((j - 1) as usize) } as usize;
            let mut sij = unsafe { *ci_row.get_unchecked(cj) };

            // Extra score by distance to the best diagonal (max distance 3).
            let extra = unsafe { *extra_score.get_unchecked(((j1 - max_diag).abs() & 3) as usize) };
            sij += extra * (sij > 0) as i32;

            let mut back = DP_BACK_LEFT_TOP;
            let mut best_score1 =
                unsafe { *score_mat.get_unchecked(row_im1 + j1u) } + sij as i64;
            let gap0 = unsafe {
                *gap_open.get_unchecked((at_last_row as usize) | ((j == len2) as usize))
            };

            if j1 > 0 {
                let mut gap = gap0;
                if unsafe { *back_mat.get_unchecked(row_i + j1u - 1) } == DP_BACK_LEFT {
                    gap = mat.ext_gap;
                }
                let score = unsafe { *score_mat.get_unchecked(row_i + j1u - 1) } + gap as i64;
                if score > best_score1 {
                    back = DP_BACK_LEFT;
                    best_score1 = score;
                }
            }
            if j1 + 1 < band_width {
                let mut gap = gap0;
                if unsafe { *back_mat.get_unchecked(row_im1 + j1u + 1) } == DP_BACK_TOP {
                    gap = mat.ext_gap;
                }
                let score = unsafe { *score_mat.get_unchecked(row_im1 + j1u + 1) } + gap as i64;
                if score > best_score1 {
                    back = DP_BACK_TOP;
                    best_score1 = score;
                }
            }
            unsafe {
                *score_mat.get_unchecked_mut(row_i + j1u) = best_score1;
                *back_mat.get_unchecked_mut(row_i + j1u) = back;
            }
        }
    }

    // Choose the endpoint (cdhit-common.c++:921-933).
    let (mut i, mut j);
    if len2 - band_left < len1 {
        i = len2 - band_left;
        j = len2;
    } else if len1 + band_right < len2 {
        i = len1;
        j = len1 + band_right;
    } else {
        i = len1;
        j = len2;
    }
    let mut j1 = j - i - band_left;
    let best_score_i64 = score_mat[idx(i, j1)];

    // Traceback.
    let mut back = back_mat[idx(i, j1)];
    let mut last = back;
    let mut count = 0i32;
    let mut count2 = 0i32;
    let mut count3 = 0i32;
    let (mut begin1, mut begin2, mut end1, mut end2) = (0i32, 0i32, 0i32, 0i32);
    let (mut gbegin1, mut gbegin2, mut gend1, mut gend2) = (0i32, 0i32, 0i32, 0i32);
    let mut smin = best_score_i64;
    let mut smax = best_score_i64 - 1;
    let mut posmin = 0i32;
    let mut posmax = 0i32;
    let mut pos = 0i32;
    let mut dlen = 0i32;
    let mut dcount = 0i32;
    let mut masked = 0i32;
    let mut indels = 0i32;
    let mut max_indels = 0i32;

    while back != DP_BACK_NONE {
        match back {
            DP_BACK_TOP => {
                let bl = ((last != back) & (j != 1) & (j != len2)) as i32;
                dlen += bl;
                dcount += bl;
                let score = score_mat[idx(i, j1)];
                if score < smin {
                    count2 = 0;
                    smin = score;
                    posmin = pos - 1;
                    begin1 = i;
                    begin2 = j;
                }
                i -= 1;
                j1 += 1;
            }
            DP_BACK_LEFT => {
                let bl = ((last != back) & (i != 1) & (i != len1)) as i32;
                dlen += bl;
                dcount += bl;
                let score = score_mat[idx(i, j1)];
                if score < smin {
                    count2 = 0;
                    smin = score;
                    posmin = pos - 1;
                    begin1 = i;
                    begin2 = j;
                }
                j1 -= 1;
                j -= 1;
            }
            DP_BACK_LEFT_TOP => {
                if want_alninfo && options.global_identity {
                    if i == 1 || j == 1 {
                        gbegin1 = i - 1;
                        gbegin2 = j - 1;
                    } else if i == len1 || j == len2 {
                        gend1 = i - 1;
                        gend2 = j - 1;
                    }
                }
                let score = score_mat[idx(i, j1)];
                i -= 1;
                j -= 1;
                let is_match = iseq1[i as usize] == iseq2[j as usize];
                if score > smax {
                    count = 0;
                    smax = score;
                    posmax = pos;
                    end1 = i;
                    end2 = j;
                }
                if options.is_est
                    && (iseq1[i as usize] > 4 || iseq2[j as usize] > 4)
                {
                    masked += 1;
                } else {
                    dlen += 1;
                    dcount += (!is_match) as i32;
                    count += is_match as i32;
                    count2 += is_match as i32;
                    count3 += is_match as i32;
                }
                if score < smin {
                    let mm = (is_match == false) as i32;
                    count2 = 0;
                    smin = score;
                    posmin = pos - mm;
                    begin1 = i + mm;
                    begin2 = j + mm;
                }
            }
            _ => {}
        }
        if options.is454 {
            if back == DP_BACK_LEFT_TOP {
                if indels > max_indels {
                    max_indels = indels;
                }
                indels = 0;
            } else if last == DP_BACK_LEFT_TOP {
                indels = 1;
            } else if indels != 0 {
                indels += 1;
            }
        }
        pos += 1;
        last = back;
        back = back_mat[idx(i, j1)];
    }

    if options.is454 && max_indels > options.max_indel {
        return None;
    }

    let iden_no = if options.global_identity {
        count3
    } else {
        count - count2
    };
    let alnln = posmin - posmax + 1 - masked;
    let dist = dcount as f32 / dlen as f32;

    let umtail1 = len1 - 1 - end1;
    let umtail2 = len2 - 1 - end2;
    let umhead = if begin1 < begin2 { begin1 } else { begin2 };
    let umtail = if umtail1 < umtail2 { umtail1 } else { umtail2 };
    let umlen = umhead + umtail;
    if umlen > options.unmatch_len {
        return None;
    }
    if (umlen as f64) > len1 as f64 * options.short_unmatch_per {
        return None;
    }
    if (umlen as f64) > len2 as f64 * options.long_unmatch_per {
        return None;
    }

    let mut alninfo = [0i32; 5];
    if want_alninfo {
        alninfo[0] = begin1;
        alninfo[1] = end1;
        alninfo[2] = begin2;
        alninfo[3] = end2;
        alninfo[4] = masked;
        if options.global_identity {
            alninfo[0] = gbegin1;
            alninfo[1] = gend1;
            alninfo[2] = gbegin2;
            alninfo[3] = gend2;
        }
    }

    Some(BandAlign {
        best_score: best_score_i64 as i32,
        iden_no,
        alnln,
        dist,
        alninfo,
    })
}

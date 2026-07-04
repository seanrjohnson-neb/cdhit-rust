//! Greedy incremental clustering: `CheckOneAA`, `ClusterOne`, `DoClustering`.
//!
//! Port of the serial clustering path (cdhit-common.c++:3228-3405, 2911-2954,
//! 3661-3784). The C++ organises work into memory-bounded blocks, but that is a
//! pure memory optimisation: each candidate is still compared against every
//! earlier representative, so a single growing word table produces identical
//! clusters. We use that single-table greedy loop.
//!
//! This module currently implements the protein (`CheckOneAA`) path for the
//! non-fragment, non-2D configuration used by default `cd-hit`. EST/2D/frag
//! paths are added in a later phase.

use crate::align::{diag_test_aapn, diag_test_aapn_est, local_band_align};
use crate::alphabet::{Naa, ScoreMatrix};
use crate::buffer::WorkingBuffer;
use crate::cutoff::{
    cal_aax_cutoff, update_aax_cutoff, upper_bound_length_rep, WorkingParam,
};
use crate::est::{make_comp_iseq, make_comp_short_word_index};
use crate::options::Options;
use crate::seqdb::SequenceDb;
use crate::sequence::{Sequence, IS_MINUS_STRAND, IS_REDUNDANT, IS_REP};
use crate::wordtable::WordTable;

const OK_FUNC: i32 = 0;

/// Outcome of comparing a candidate against the current word table.
struct CheckContext<'a> {
    reps: &'a [Sequence],
    table: &'a WordTable,
    naa: &'a Naa,
    mat: &'a ScoreMatrix,
    /// Reverse-complement word-index table (EST, `option_r`); empty otherwise.
    comp_aan_idx: &'a [i32],
}

/// Port of `CheckOneAA` (cdhit-common.c++:3235-3405) for the frag_size == 0,
/// has2D == false configuration. Returns 1 if `seq` is redundant against an
/// existing representative (fields on `seq` are updated), else 0.
fn check_one_aa(
    seq: &mut Sequence,
    ctx: &CheckContext,
    param: &mut WorkingParam,
    buf: &mut WorkingBuffer,
    options: &Options,
) -> i32 {
    let mut aa1_cutoff = param.aa1_cutoff;
    let mut aa2_cutoff = param.aas_cutoff;
    let mut aan_cutoff = param.aan_cutoff;

    let table = ctx.table;
    let naa1 = ctx.naa.array[1];
    let naa2 = ctx.naa.array[2];
    let naa = options.naa;
    let len = seq.size;
    let mut flag = 0;

    let s = table.sequences.len();
    let mut len_eff = len;

    if s != 0 {
        let min_size = ctx.reps[table.sequences[s - 1]].size;
        if min_size < len {
            let mut min = min_size;
            if (len as f64 * options.diff_cutoff2) > min as f64 {
                min = (len as f64 * options.diff_cutoff2) as i32;
            }
            if (len - options.diff_cutoff_aa2) > min {
                min = len - options.diff_cutoff_aa2;
            }
            len_eff = min;
        }
    }

    // -aL early skip (cdhit-common.c++:3274-3278).
    if s != 0 {
        let min_size = ctx.reps[table.sequences[s - 1]].size;
        let min_red = min_size as f64 * options.long_coverage - options.band_width as f64;
        if (len as f64) < min_red {
            return 0;
        }
    }

    param.control_short_coverage(len_eff, options);
    param.compute_required_bases(naa, 2, options);

    buf.encode_words(&seq.data, seq.size, naa, ctx.naa, false);

    if options.min_control > len {
        return 0;
    }

    let aan_no = (len - naa + 1) as usize;
    let required_aan = param.required_aan;
    let look_len = table.count_words(
        aan_no,
        &buf.word_encodes,
        &buf.word_encodes_no,
        &mut buf.look_counts,
        &mut buf.index_mapping,
        false,
        required_aan,
        0,
    );

    let len_upper_bound = param.len_upper_bound;
    let len_lower_bound = param.len_lower_bound;
    let mut has_aa2 = false;

    for t in 0..look_len {
        let ic = buf.look_counts[t];
        // Non-frag path: clear mapping and skip weak candidates.
        buf.index_mapping[ic.index as usize] = 0;
        if ic.count < required_aan {
            continue;
        }

        let rep = &ctx.reps[table.sequences[ic.index as usize]];
        let len2 = rep.size;
        if len2 > len_upper_bound {
            continue;
        }
        if options.has2d && len2 < len_lower_bound {
            continue;
        }

        param.control_long_coverage(len2, options);

        if !has_aa2 {
            buf.compute_aap(&seq.data, seq.size, naa1, naa2);
            has_aa2 = true;
        }

        let band_width1 = if options.band_width < len + len2 - 2 {
            options.band_width
        } else {
            len + len2 - 2
        };
        let diag = diag_test_aapn(
            naa1,
            &rep.data,
            len,
            len2,
            buf,
            band_width1,
            param.required_aa1,
            options,
        );
        if diag.best_sum < param.required_aas {
            continue;
        }

        let align = local_band_align(
            &seq.data,
            &rep.data,
            len,
            len2,
            ctx.mat,
            diag.band_left,
            diag.band_center,
            diag.band_right,
            options,
            true,
        );
        let align = match align {
            Some(a) => a,
            None => continue,
        };
        if align.iden_no < param.required_aa1 {
            continue;
        }
        let lens = if options.has2d && len > len2 { len2 } else { len };
        let len_eff1 = if !options.global_identity {
            align.alnln
        } else {
            lens - align.alninfo[4]
        };
        let tiden_pc = align.iden_no as f32 / len_eff1 as f32;
        if options.use_distance {
            if align.dist > options.distance_thd as f32 {
                continue;
            }
            if align.dist >= seq.distance {
                continue;
            }
        } else {
            if tiden_pc < options.cluster_thd as f32 {
                continue;
            }
            if tiden_pc <= seq.identity {
                continue;
            }
        }
        if param.aln_cover_flag != 0 {
            if align.alninfo[3] - align.alninfo[2] + 1 < param.min_aln_len_l {
                continue;
            }
            if align.alninfo[1] - align.alninfo[0] + 1 < param.min_aln_len_s {
                continue;
            }
        }
        if options.has2d {
            seq.state |= IS_REDUNDANT;
        }
        flag = 1;
        seq.identity = tiden_pc;
        seq.cluster_id = rep.cluster_id;
        seq.distance = align.dist;
        seq.coverage[0] = align.alninfo[0] + 1;
        seq.coverage[1] = align.alninfo[1] + 1;
        seq.coverage[2] = align.alninfo[2] + 1;
        seq.coverage[3] = align.alninfo[3] + 1;
        if !options.cluster_best {
            break;
        }
        // NOTE: the C++ updates only these *local* cutoff copies; it does not
        // write them back to `param`, so the subsequent ComputeRequiredBases
        // still uses the original param cutoffs. Replicate that exactly — do
        // NOT propagate to `param` (doing so diverges the -g 1 results).
        update_aax_cutoff(
            &mut aa1_cutoff,
            &mut aa2_cutoff,
            &mut aan_cutoff,
            options.tolerance,
            naa,
            tiden_pc as f64,
        );
        param.compute_required_bases(naa, 2, options);
    }

    // Clean up index_mapping for all touched entries so the buffer is clean for
    // the next call (matches the C++ post-loop clear).
    for t in 0..look_len {
        buf.index_mapping[buf.look_counts[t].index as usize] = 0;
    }

    if flag == 1 && !options.cluster_best {
        seq.data.clear();
        seq.state |= IS_REDUNDANT;
    }
    flag
}

/// Port of `CheckOneEST` (cdhit-common.c++:3406-3593) for frag_size == 0.
/// Runs a forward pass and (if `option_r`) a reverse-complement pass; returns
/// 1 for a forward match, -1 for a reverse-complement match, 0 for no match.
fn check_one_est(
    seq: &mut Sequence,
    ctx: &CheckContext,
    param: &mut WorkingParam,
    buf: &mut WorkingBuffer,
    options: &Options,
) -> i32 {
    let table = ctx.table;
    let naa1 = ctx.naa.array[1];
    let naa = options.naa;
    let len = seq.size;
    let mut flag = 0;

    let s = table.sequences.len();
    let mut len_eff = len;
    if s != 0 {
        let min_size = ctx.reps[table.sequences[s - 1]].size;
        if min_size < len {
            let mut min = min_size;
            if (len as f64 * options.diff_cutoff2) > min as f64 {
                min = (len as f64 * options.diff_cutoff2) as i32;
            }
            if (len - options.diff_cutoff_aa2) > min {
                min = len - options.diff_cutoff_aa2;
            }
            len_eff = min;
        }
    }
    if s != 0 {
        let min_size = ctx.reps[table.sequences[s - 1]].size;
        let min_red = min_size as f64 * options.long_coverage - options.band_width as f64;
        if (len as f64) < min_red {
            return 0;
        }
    }

    param.control_short_coverage(len_eff, options);
    param.compute_required_bases(naa, 4, options);
    let skip = buf.encode_words(&seq.data, seq.size, naa, ctx.naa, true);
    param.required_aan -= skip;
    param.required_aas -= skip;
    param.required_aa1 -= skip;
    if param.required_aan <= 0 {
        param.required_aan = 1;
    }
    if param.required_aas <= 0 {
        param.required_aas = 1;
    }
    if param.required_aa1 <= 0 {
        param.required_aa1 = 1;
    }

    if options.min_control > len {
        return 0;
    }

    let aan_no = (len - naa + 1) as usize;
    let len_upper_bound = param.len_upper_bound;
    let len_lower_bound = param.len_lower_bound;

    // Forward seqi = seq.data; reverse pass uses a local revcomp buffer.
    let comp_seq: Vec<u8> = Vec::new();
    let mut comp_seq = comp_seq;

    for comp in 0..2 {
        // Build the word list for this pass.
        let word_list: Vec<i32> = if comp != 0 {
            let mut list = vec![0i32; aan_no];
            for j0 in 0..aan_no {
                let j = buf.word_encodes[j0];
                list[j0] = if j < 0 { j } else { ctx.comp_aan_idx[j as usize] };
            }
            comp_seq = vec![0u8; len as usize];
            make_comp_iseq(len as usize, &mut comp_seq, &seq.data);
            list
        } else {
            Vec::new()
        };

        let required_aan = param.required_aan;
        let look_len = if comp != 0 {
            table.count_words(
                aan_no,
                &word_list,
                &buf.word_encodes_no,
                &mut buf.look_counts,
                &mut buf.index_mapping,
                true,
                required_aan,
                0,
            )
        } else {
            table.count_words(
                aan_no,
                &buf.word_encodes,
                &buf.word_encodes_no,
                &mut buf.look_counts,
                &mut buf.index_mapping,
                true,
                required_aan,
                0,
            )
        };

        // seqi for this pass.
        let seqi: &[u8] = if comp != 0 { &comp_seq } else { &seq.data };
        let mut has_aas = false;

        for t in 0..look_len {
            let ic = buf.look_counts[t];
            buf.index_mapping[ic.index as usize] = 0;
            if ic.count < required_aan {
                continue;
            }
            let rep = &ctx.reps[table.sequences[ic.index as usize]];
            let len2 = rep.size;
            if len2 > len_upper_bound {
                continue;
            }
            if options.has2d && len2 < len_lower_bound {
                continue;
            }

            param.control_long_coverage(len2, options);

            if !has_aas {
                buf.compute_aap2(seqi, seq.size, ctx.naa);
                has_aas = true;
            }

            let band_width1 = if options.band_width < len + len2 - 2 {
                options.band_width
            } else {
                len + len2 - 2
            };
            let diag = diag_test_aapn_est(
                naa1,
                &rep.data,
                len,
                len2,
                buf,
                band_width1,
                param.required_aa1,
                options,
            );
            if diag.best_sum < param.required_aas {
                continue;
            }

            let align = local_band_align(
                seqi,
                &rep.data,
                len,
                len2,
                ctx.mat,
                diag.band_left,
                diag.band_center,
                diag.band_right,
                options,
                true,
            );
            let mut align = match align {
                Some(a) => a,
                None => continue,
            };
            if (options.print != 0 || param.aln_cover_flag != 0) && comp != 0 {
                align.alninfo[0] = len - align.alninfo[0] - 1;
                align.alninfo[1] = len - align.alninfo[1] - 1;
            }
            if align.iden_no < param.required_aa1 {
                continue;
            }
            if options.is454 {
                if align.alninfo[2] != align.alninfo[0] {
                    continue;
                }
                if align.alninfo[0] > 1 {
                    continue;
                }
                if (len - align.alninfo[1]) > 2 {
                    continue;
                }
            }

            let lens = if options.has2d && len > len2 { len2 } else { len };
            let len_eff1 = if !options.global_identity {
                align.alnln
            } else {
                lens - align.alninfo[4]
            };
            let tiden_pc = align.iden_no as f32 / len_eff1 as f32;
            if options.use_distance {
                if align.dist > options.distance_thd as f32 {
                    continue;
                }
                if options.cluster_best && align.dist >= seq.distance {
                    continue;
                }
            } else {
                if tiden_pc < options.cluster_thd as f32 {
                    continue;
                }
                if options.cluster_best && tiden_pc < seq.identity {
                    continue;
                }
            }
            if param.aln_cover_flag != 0 {
                if align.alninfo[3] - align.alninfo[2] + 1 < param.min_aln_len_l {
                    continue;
                }
                if comp != 0 {
                    if align.alninfo[0] - align.alninfo[1] + 1 < param.min_aln_len_s {
                        continue;
                    }
                } else if align.alninfo[1] - align.alninfo[0] + 1 < param.min_aln_len_s {
                    continue;
                }
            }
            if options.cluster_best
                && (tiden_pc - seq.identity).abs() < 1e-9
                && rep.cluster_id >= seq.cluster_id
            {
                continue;
            }
            if !options.cluster_best && flag != 0 && rep.cluster_id >= seq.cluster_id {
                continue;
            }
            flag = if comp != 0 { -1 } else { 1 };
            seq.identity = tiden_pc;
            seq.distance = align.dist;
            seq.cluster_id = rep.cluster_id;
            seq.coverage[0] = align.alninfo[0] + 1;
            seq.coverage[1] = align.alninfo[1] + 1;
            seq.coverage[2] = align.alninfo[2] + 1;
            seq.coverage[3] = align.alninfo[3] + 1;
            if !options.cluster_best {
                break;
            }
        }
        for t in 0..look_len {
            buf.index_mapping[buf.look_counts[t].index as usize] = 0;
        }
        if !(options.option_r != 0) {
            break;
        }
    }

    if flag == 1 || flag == -1 {
        if !options.cluster_best {
            seq.data.clear();
            seq.state |= IS_REDUNDANT;
        }
        if flag == -1 {
            seq.state |= IS_MINUS_STRAND;
        } else {
            seq.state &= !IS_MINUS_STRAND;
        }
    }
    flag
}

impl crate::seqdb::SequenceDb {
    /// Port of the serial `DoClustering` (cdhit-common.c++:3661-3784) as a
    /// single-table greedy loop. Populates `rep_seqs` and per-sequence cluster
    /// assignments.
    pub fn do_clustering(&mut self, options: &Options, scoring: &Scoring, naa: &Naa) {
        let seq_no = self.sequences.len();
        let mut aa1_cutoff = options.cluster_thd;
        let mut aas_cutoff = 1.0 - (1.0 - options.cluster_thd) * 4.0;
        let mut aan_cutoff = 1.0 - (1.0 - options.cluster_thd) * options.naa as f64;
        if !options.is_est {
            let (a1, a2, an) =
                cal_aax_cutoff(options.cluster_thd, options.tolerance, options.naa);
            aa1_cutoff = a1;
            aas_cutoff = a2;
            aan_cutoff = an;
        }

        let mut param = WorkingParam::new(aa1_cutoff, aas_cutoff, aan_cutoff);
        let mut buf = WorkingBuffer::new(seq_no, self.max_len as usize, options);
        let mut table = WordTable::new(options.naa, self.naan);
        if options.is_est {
            table.set_dna();
        }

        // Reverse-complement word-index table for EST + option_r
        // (cdhit-est.c++:73-76).
        let comp_aan_idx: Vec<i32> = if options.is_est && options.option_r != 0 {
            make_comp_short_word_index(options.naa, naa)
        } else {
            Vec::new()
        };

        for ks in 0..seq_no {
            if self.sequences[ks].state & IS_REDUNDANT != 0 {
                continue;
            }
            cluster_one(
                &mut self.sequences,
                ks,
                &mut table,
                &mut param,
                &mut buf,
                options,
                scoring,
                naa,
                &comp_aan_idx,
                &mut self.rep_seqs,
            );
        }
    }
}

use crate::options::Scoring;

impl crate::seqdb::SequenceDb {
    /// Port of `ClusterTo` (cdhit-common.c++:3786-3962): cluster this database
    /// (db2) against a fixed reference `other` (db1). db2 sequences that match a
    /// db1 representative are marked redundant; `rep_seqs` collects the novel
    /// db2 sequences. Both databases must already be sorted/encoded.
    ///
    /// Uses a single word table over all of db1 (the C++ block-batching is a
    /// memory optimisation that does not change the result: each db2 sequence is
    /// still compared against every db1 representative within its length bounds).
    pub fn cluster_to(
        &mut self,
        other: &SequenceDb,
        options: &Options,
        scoring: &Scoring,
        naa: &Naa,
        comp_aan_idx: &[i32],
    ) {
        let mut aa1_cutoff = options.cluster_thd;
        let mut aas_cutoff = 1.0 - (1.0 - options.cluster_thd) * 4.0;
        let mut aan_cutoff = 1.0 - (1.0 - options.cluster_thd) * options.naa as f64;
        if !options.is_est {
            let (a1, a2, an) =
                cal_aax_cutoff(options.cluster_thd, options.tolerance, options.naa);
            aa1_cutoff = a1;
            aas_cutoff = a2;
            aan_cutoff = an;
        }
        let mut param = WorkingParam::new(aa1_cutoff, aas_cutoff, aan_cutoff);
        let n = other.sequences.len();
        let mut buf = WorkingBuffer::new(n, self.max_len.max(other.max_len) as usize, options);
        let mut table = WordTable::new(options.naa, self.naan);
        if options.is_est {
            table.set_dna();
        }

        // Build the word table from all of db1 (other).
        for ks in 0..n {
            let seq = &other.sequences[ks];
            let aan_no = (seq.size - options.naa + 1) as usize;
            buf.encode_words(&seq.data, seq.size, options.naa, naa, options.is_est);
            let idx = table.sequences.len() as i32;
            table.add_word_counts_encoded(
                aan_no,
                &buf.word_encodes,
                &buf.word_encodes_no,
                idx,
                options.is_est,
            );
            table.sequences.push(ks);
        }

        let ctx = CheckContext {
            reps: &other.sequences,
            table: &table,
            naa,
            mat: &scoring.mat,
            comp_aan_idx,
        };

        let m = self.sequences.len();
        for j in 0..m {
            if self.sequences[j].state & IS_REDUNDANT != 0 {
                continue;
            }
            let len = self.sequences[j].size;
            let len_upper_bound = upper_bound_length_rep(len, options);
            let mut len_lower_bound = len - options.diff_cutoff_aa2;
            let len_tmp = (len as f64 * options.diff_cutoff2) as i32;
            if len_tmp < len_lower_bound {
                len_lower_bound = len_tmp;
            }
            param.len_upper_bound = len_upper_bound;
            param.len_lower_bound = len_lower_bound;

            let seq = &mut self.sequences[j];
            let flag = if options.is_est {
                check_one_est(seq, &ctx, &mut param, &mut buf, options)
            } else {
                check_one_aa(seq, &ctx, &mut param, &mut buf, options)
            };
            if flag == 1 || flag == -1 {
                if !options.cluster_best {
                    seq.data.clear();
                    seq.state |= IS_REDUNDANT;
                }
                if flag == -1 {
                    seq.state |= IS_MINUS_STRAND;
                }
            }
        }

        // cluster_best: finalize redundancy by identity.
        if options.cluster_best {
            for seq in &mut self.sequences {
                if seq.identity > 0.0 {
                    seq.state |= IS_REDUNDANT;
                }
            }
        }
        for i in 0..self.sequences.len() {
            if self.sequences[i].identity < 0.0 {
                self.sequences[i].identity *= -1.0;
            }
            if self.sequences[i].state & IS_REDUNDANT == 0 {
                self.rep_seqs.push(i as i32);
            }
        }
    }
}

/// Port of `ClusterOne` (cdhit-common.c++:2911-2954), frag_size == 0.
#[allow(clippy::too_many_arguments)]
fn cluster_one(
    sequences: &mut [Sequence],
    id: usize,
    table: &mut WordTable,
    param: &mut WorkingParam,
    buf: &mut WorkingBuffer,
    options: &Options,
    scoring: &Scoring,
    naa: &Naa,
    comp_aan_idx: &[i32],
    rep_seqs: &mut Vec<i32>,
) {
    if sequences[id].state & IS_REDUNDANT != 0 {
        return;
    }
    let len = sequences[id].size;
    param.len_upper_bound = upper_bound_length_rep(len, options);

    // All representatives have global index < id, so split there and treat the
    // left part as read-only reps, the candidate as the single mutable element.
    let (left, right) = sequences.split_at_mut(id);
    let seq = &mut right[0];
    param.len_upper_bound = upper_bound_length_rep(len, options);
    let ctx = CheckContext {
        reps: left,
        table,
        naa,
        mat: &scoring.mat,
        comp_aan_idx,
    };
    let flag = if options.is_est {
        check_one_est(seq, &ctx, param, buf, options)
    } else {
        check_one_aa(seq, &ctx, param, buf, options)
    };

    if flag == 0 {
        if seq.identity > 0.0 && options.cluster_best {
            seq.state |= IS_REDUNDANT;
            seq.data.clear();
        } else {
            let aan_no = (len - options.naa + 1) as usize;
            let size = rep_seqs.len();
            rep_seqs.push(id as i32);
            seq.cluster_id = size as i32;
            seq.identity = 0.0;
            seq.state |= IS_REP;
            let idx = table.sequences.len() as i32;
            table.add_word_counts_encoded(
                aan_no,
                &buf.word_encodes,
                &buf.word_encodes_no,
                idx,
                options.is_est,
            );
            table.sequences.push(id);
        }
    }
    let _ = OK_FUNC;
}

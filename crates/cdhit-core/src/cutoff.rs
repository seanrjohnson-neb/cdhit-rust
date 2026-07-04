//! `WorkingParam` and the cutoff computations.
//!
//! Port of `WorkingParam` methods, `cal_aax_cutoff`, `update_aax_cutoff`,
//! `upper_bound_length_rep`, and `ComputeRequiredBases`
//! (cdhit-common.c++:2694-2835).

use crate::naa_stat::{NAA_STAT, NAA_STAT_START_PERCENT};
use crate::options::Options;

#[derive(Clone, Debug)]
pub struct WorkingParam {
    pub aa1_cutoff: f64,
    pub aas_cutoff: f64,
    pub aan_cutoff: f64,
    pub len_upper_bound: i32,
    pub len_lower_bound: i32,

    pub len_eff: i32,
    pub aln_cover_flag: i32,
    pub min_aln_len_s: i32,
    pub min_aln_len_l: i32,
    pub required_aa1: i32,
    pub required_aas: i32,
    pub required_aan: i32,
}

impl WorkingParam {
    pub fn new(a1: f64, a2: f64, an: f64) -> Self {
        WorkingParam {
            aa1_cutoff: a1,
            aas_cutoff: a2,
            aan_cutoff: an,
            len_upper_bound: 0,
            len_lower_bound: 0,
            len_eff: 0,
            aln_cover_flag: 0,
            min_aln_len_s: 0,
            min_aln_len_l: 0,
            required_aa1: 0,
            required_aas: 0,
            required_aan: 0,
        }
    }

    /// Port of `ControlShortCoverage` (cdhit-common.c++:2694-2705).
    pub fn control_short_coverage(&mut self, len: i32, options: &Options) {
        self.len_eff = len;
        self.aln_cover_flag = 0;
        if options.short_coverage > 0.0 || options.min_control > 0 {
            self.aln_cover_flag = 1;
            self.min_aln_len_s = (len as f64 * options.short_coverage) as i32;
            if len - options.short_control > self.min_aln_len_s {
                self.min_aln_len_s = len - options.short_control;
            }
            if options.min_control > self.min_aln_len_s {
                self.min_aln_len_s = options.min_control;
            }
        }
        if !options.global_identity {
            self.len_eff = self.min_aln_len_s;
        }
    }

    /// Port of `ControlLongCoverage` (cdhit-common.c++:2706-2713).
    pub fn control_long_coverage(&mut self, len2: i32, options: &Options) {
        if self.aln_cover_flag != 0 {
            self.min_aln_len_l = (len2 as f64 * options.long_coverage) as i32;
            if len2 - options.long_control > self.min_aln_len_l {
                self.min_aln_len_l = len2 - options.long_control;
            }
            if options.min_control > self.min_aln_len_l {
                self.min_aln_len_l = options.min_control;
            }
        }
    }

    /// Port of `ComputeRequiredBases` (cdhit-common.c++:2772-2835).
    pub fn compute_required_bases(&mut self, naa: i32, ss: i32, options: &Options) {
        if options.use_distance {
            let band = options.band_width + 1;
            let invd = (1.0 / (options.distance_thd + 1e-9)) as i32;
            let k = if self.len_eff < invd { self.len_eff } else { invd };
            let kn = self.len_eff - naa + 1;
            let ks2 = invd - ss;
            let kn2 = invd - naa;
            let ks = self.len_eff - ss + 1;
            let _ks3 = ((self.len_eff - band + 1) as f64 * (1.0 - options.distance_thd * ss as f64)) as i32;
            let _kn3 = ((self.len_eff - band + 1) as f64 * (1.0 - options.distance_thd * naa as f64)) as i32;
            self.required_aa1 = if ks2 < ks { ks2 } else { ks };
            self.required_aas = self.required_aa1;
            self.required_aan = if kn2 < kn { kn2 } else { kn };
            if self.required_aa1 <= 0 {
                self.required_aa1 = 1;
                self.required_aas = 1;
            }
            if self.required_aan <= 0 {
                self.required_aan = 1;
            }
            let _ = k;
            return;
        }
        let len_eff = self.len_eff;
        self.required_aa1 =
            (len_eff - ss) - (ss as f64 * ((1.0 - self.aa1_cutoff) * len_eff as f64).ceil()) as i32;
        if self.required_aa1 < 0 {
            self.required_aa1 = 0;
        }
        self.required_aas = self.required_aa1;
        self.required_aan =
            (len_eff - naa) - (naa as f64 * ((1.0 - self.aa1_cutoff) * len_eff as f64).ceil()) as i32;
        if self.required_aan < 0 {
            self.required_aan = 0;
        }

        let aa1_old = (self.aa1_cutoff * len_eff as f64) as i32 - ss + 1;
        let aas_old = (self.aas_cutoff * len_eff as f64) as i32;
        let aan_old = (self.aan_cutoff * len_eff as f64) as i32;

        let thd = options.cluster_thd;
        let rest = (len_eff - naa) as f64 / (len_eff * naa) as f64;
        let thd0 = 1.0 - rest;
        let mut fnew = 0.0;
        let mut fold = 1.0;
        if thd > thd0 {
            fnew = (thd - thd0) / rest;
            fold = 1.0 - fnew;
        }
        self.required_aa1 = (fnew * self.required_aa1 as f64 + fold * aa1_old as f64) as i32;
        self.required_aas = (fnew * self.required_aas as f64 + fold * aas_old as f64) as i32;
        self.required_aan = (fnew * self.required_aan as f64 + fold * aan_old as f64) as i32;
    }
}

/// Port of `upper_bound_length_rep` (cdhit-common.c++:2718-2735).
pub fn upper_bound_length_rep(len: i32, options: &Options) -> i32 {
    let opt_s = options.diff_cutoff;
    let opt_s_aa = options.diff_cutoff_aa;
    let opt_al = options.long_coverage;
    let opt_al_ctrl = options.long_control;
    let mut len_upper_bound = 99_999_999i32;
    let r1 = if opt_s > opt_al { opt_s } else { opt_al };
    let a2 = if opt_s_aa < opt_al_ctrl { opt_s_aa } else { opt_al_ctrl };
    if r1 > 0.0 {
        len_upper_bound = (len as f32 / r1 as f32) as i32;
    }
    if (len + a2) < len_upper_bound {
        len_upper_bound = len + a2;
    }
    len_upper_bound
}

/// Port of `cal_aax_cutoff` (cdhit-common.c++:2738-2754).
pub fn cal_aax_cutoff(
    cluster_thd: f64,
    tolerance: i32,
    naa: i32,
) -> (f64, f64, f64) {
    let aa1_cutoff = cluster_thd;
    let mut aa2_cutoff = 1.0 - (1.0 - cluster_thd) * 2.0;
    let mut aan_cutoff = 1.0 - (1.0 - cluster_thd) * naa as f64;
    if tolerance == 0 {
        return (aa1_cutoff, aa2_cutoff, aan_cutoff);
    }
    let mut clstr_idx = (cluster_thd * 100.0) as i32 - NAA_STAT_START_PERCENT;
    if clstr_idx < 0 {
        clstr_idx = 0;
    }
    let d2 = NAA_STAT[(tolerance - 1) as usize][clstr_idx as usize][3] as f64 / 100.0;
    let dn = NAA_STAT[(tolerance - 1) as usize][clstr_idx as usize][(5 - naa) as usize] as f64 / 100.0;
    if d2 > aa2_cutoff {
        aa2_cutoff = d2;
    }
    if dn > aan_cutoff {
        aan_cutoff = dn;
    }
    (aa1_cutoff, aa2_cutoff, aan_cutoff)
}

/// Port of `update_aax_cutoff` (cdhit-common.c++:2757-2770). Raises the cutoffs
/// toward those implied by a newly-observed identity (used only with `-g 1`).
pub fn update_aax_cutoff(
    aa1_cutoff: &mut f64,
    aa2_cutoff: &mut f64,
    aan_cutoff: &mut f64,
    tolerance: i32,
    naa: i32,
    mut cluster_thd: f64,
) {
    if cluster_thd > 1.0 {
        cluster_thd = 1.0;
    }
    let (aa1_t, aa2_t, aan_t) = cal_aax_cutoff(cluster_thd, tolerance, naa);
    if aa1_t > *aa1_cutoff {
        *aa1_cutoff = aa1_t;
    }
    if aa2_t > *aa2_cutoff {
        *aa2_cutoff = aa2_t;
    }
    if aan_t > *aan_cutoff {
        *aan_cutoff = aan_t;
    }
}

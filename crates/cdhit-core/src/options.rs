//! Command-line options, defaults, parsing, and validation.
//!
//! Port of `Options` (cdhit-common.h:254-370) and the `SetOption*` / `Validate`
//! routines (cdhit-common.c++:211-429). Field names and default values mirror
//! the C++ exactly so parsed configurations are byte-identical.

use crate::alphabet::{ScoreMatrix, AA2IDX, MAX_SEQ, NA2IDX};
use crate::error::{CdError, Result};
use crate::naa_stat::{NAA_STAT, NAA_STAT_START_PERCENT};

/// Runtime alphabet/scoring state mutated by EST options (`-gap`, `-match`,
/// `-mismatch`, `-mask`). The C++ keeps these as globals (`mat`, `aa2idx`,
/// `na2idx`); we thread them explicitly so the engine is reentrant.
#[derive(Clone)]
pub struct Scoring {
    pub mat: ScoreMatrix,
    /// Active encoding table used by `Sequence::convert_bases` (the C++ global
    /// `aa2idx`). For EST it is overwritten with `na2idx` via `setaa_to_na`.
    pub aa2idx: [i32; 26],
    /// Nucleotide code table, a mutable copy of `NA2IDX` (only `-mask` edits it).
    pub na2idx: [i32; 26],
}

impl Default for Scoring {
    fn default() -> Self {
        Scoring {
            mat: ScoreMatrix::default(),
            aa2idx: AA2IDX,
            na2idx: NA2IDX,
        }
    }
}

/// Port of `setaa_to_na` (cdhit-common.c++:1197-1201): copy `na2idx` into the
/// active encoding table. The C++ EST driver separately calls
/// `mat.set_to_na()`; that is done during EST setup, not here.
pub fn setaa_to_na(scoring: &mut Scoring) {
    scoring.aa2idx = scoring.na2idx;
}

#[derive(Clone, Debug, PartialEq)]
pub struct Options {
    pub naa: i32,
    pub naan: i32,
    pub naa_top_limit: i32,

    pub max_memory: u64,
    pub min_length: i32,
    pub cluster_best: bool,
    pub global_identity: bool,
    pub store_disk: bool,
    pub band_width: i32,
    pub cluster_thd: f64,
    pub distance_thd: f64,
    pub diff_cutoff: f64,
    pub diff_cutoff2: f64,
    pub diff_cutoff_aa: i32,
    pub diff_cutoff_aa2: i32,
    pub tolerance: i32,
    pub long_coverage: f64,
    pub long_control: i32,
    pub short_coverage: f64,
    pub short_control: i32,
    pub min_control: i32,
    pub long_unmatch_per: f64,
    pub short_unmatch_per: f64,
    pub unmatch_len: i32,
    pub max_indel: i32,
    pub print: i32,
    pub des_len: i32,
    pub frag_size: i32,
    pub option_r: i32,
    pub threads: i32,
    pub pe_mode: i32,
    pub trim_len: i32,
    pub trim_len_r2: i32,
    pub align_pos: i32,

    pub max_entries: u64,
    pub max_sequences: u64,
    pub mem_limit: u64,

    pub has2d: bool,
    pub is_est: bool,
    pub is454: bool,
    pub use_identity: bool,
    pub use_distance: bool,
    pub backup_file: bool,

    pub input: String,
    pub input_pe: String,
    pub input2: String,
    pub input2_pe: String,
    pub output: String,
    pub output_pe: String,

    pub sort_output: i32,
    pub sort_outputf: i32,
}

impl Default for Options {
    /// Port of the `Options()` constructor (cdhit-common.h:313-358).
    fn default() -> Self {
        Options {
            naa: 5,
            naan: 0,
            naa_top_limit: 5,
            max_memory: 800_000_000,
            min_length: 10,
            cluster_best: false,
            global_identity: true,
            store_disk: false,
            band_width: 20,
            cluster_thd: 0.9,
            distance_thd: 0.0,
            diff_cutoff: 0.0,
            diff_cutoff2: 1.0,
            diff_cutoff_aa: 99_999_999,
            diff_cutoff_aa2: 0,
            tolerance: 2,
            long_coverage: 0.0,
            long_control: 99_999_999,
            short_coverage: 0.0,
            short_control: 99_999_999,
            min_control: 0,
            long_unmatch_per: 1.0,
            short_unmatch_per: 1.0,
            unmatch_len: 99_999_999,
            max_indel: 1,
            print: 0,
            des_len: 20,
            frag_size: 0,
            option_r: 1,
            threads: 1,
            pe_mode: 0,
            trim_len: 0,
            trim_len_r2: 0,
            align_pos: 0,
            max_entries: 0,
            max_sequences: 1 << 20,
            mem_limit: 100_000_000,
            has2d: false,
            is_est: false,
            is454: false,
            use_identity: false,
            use_distance: false,
            backup_file: false,
            input: String::new(),
            input_pe: String::new(),
            input2: String::new(),
            input2_pe: String::new(),
            output: String::new(),
            output_pe: String::new(),
            sort_output: 0,
            sort_outputf: 0,
        }
    }
}

// Faithful ports of C's atoi/atof/atoll: parse a leading numeric prefix,
// returning 0 on no match (never erroring). This matters because CD-HIT relies
// on that lax behaviour for some flags.
fn atoi(s: &str) -> i32 {
    parse_int_prefix(s) as i32
}
fn atoll(s: &str) -> i64 {
    parse_int_prefix(s)
}
fn parse_int_prefix(s: &str) -> i64 {
    let s = s.trim_start();
    let bytes = s.as_bytes();
    let mut i = 0;
    let mut sign = 1i64;
    if i < bytes.len() && (bytes[i] == b'+' || bytes[i] == b'-') {
        if bytes[i] == b'-' {
            sign = -1;
        }
        i += 1;
    }
    let start = i;
    let mut val: i64 = 0;
    while i < bytes.len() && bytes[i].is_ascii_digit() {
        val = val.saturating_mul(10).saturating_add((bytes[i] - b'0') as i64);
        i += 1;
    }
    if i == start {
        return 0;
    }
    sign * val
}
fn atof(s: &str) -> f64 {
    let s = s.trim_start();
    // Find the longest valid f64 prefix.
    let mut end = 0;
    let bytes = s.as_bytes();
    let mut seen_dot = false;
    let mut seen_e = false;
    let mut i = 0;
    if i < bytes.len() && (bytes[i] == b'+' || bytes[i] == b'-') {
        i += 1;
    }
    while i < bytes.len() {
        let c = bytes[i];
        if c.is_ascii_digit() {
            end = i + 1;
        } else if c == b'.' && !seen_dot && !seen_e {
            seen_dot = true;
        } else if (c == b'e' || c == b'E') && !seen_e && end > 0 {
            seen_e = true;
            if i + 1 < bytes.len() && (bytes[i + 1] == b'+' || bytes[i + 1] == b'-') {
                i += 1;
            }
        } else {
            break;
        }
        i += 1;
    }
    if end == 0 {
        return 0.0;
    }
    s[..end].parse::<f64>().unwrap_or(0.0)
}

impl Options {
    /// Port of `SetOptionCommon` (cdhit-common.c++:211-268). Returns true if the
    /// flag was recognised. `-T` is accepted but not clamped to CPU count here
    /// (that is a native-only concern handled by the CLI layer).
    fn set_option_common(&mut self, flag: &str, value: &str) -> bool {
        let intval = atoi(value);
        match flag {
            "-i" => self.input = value.to_string(),
            "-j" => self.input_pe = value.to_string(),
            "-o" => self.output = value.to_string(),
            "-op" => self.output_pe = value.to_string(),
            "-M" => self.max_memory = (atoll(value) * 1_000_000) as u64,
            "-l" => self.min_length = intval,
            "-c" => {
                self.cluster_thd = atof(value);
                self.use_identity = true;
            }
            "-D" => {
                self.distance_thd = atof(value);
                self.use_distance = true;
            }
            "-b" => self.band_width = intval,
            "-n" => self.naa = intval,
            "-d" => self.des_len = intval,
            "-s" => self.diff_cutoff = atof(value),
            "-S" => self.diff_cutoff_aa = intval,
            "-B" => self.store_disk = intval != 0,
            "-P" => self.pe_mode = intval,
            "-cx" => self.trim_len = intval,
            "-cy" => self.trim_len_r2 = intval,
            "-ap" => self.align_pos = intval,
            "-sc" => self.sort_output = intval,
            "-sf" => self.sort_outputf = intval,
            "-p" => self.print = intval,
            "-g" => self.cluster_best = intval != 0,
            "-G" => self.global_identity = intval != 0,
            "-aL" => self.long_coverage = atof(value),
            "-AL" => self.long_control = intval,
            "-aS" => self.short_coverage = atof(value),
            "-AS" => self.short_control = intval,
            "-A" => self.min_control = intval,
            "-uL" => self.long_unmatch_per = atof(value),
            "-uS" => self.short_unmatch_per = atof(value),
            "-U" => self.unmatch_len = intval,
            "-tmp" => { /* temp dir: disk-swap dropped in v1, accepted and ignored */ }
            "-bak" => self.backup_file = intval != 0,
            "-T" => self.threads = intval,
            _ => return false,
        }
        true
    }

    /// Port of `SetOption` (cdhit-common.c++:269-286).
    fn set_option(&mut self, flag: &str, value: &str, scoring: &mut Scoring) -> bool {
        if self.is454 {
            match flag {
                "-s" | "-S" | "-G" | "-A" | "-r" => return false,
                "-D" => {
                    self.max_indel = atoi(value);
                    return true;
                }
                _ => {}
            }
        }
        if self.set_option_common(flag, value) {
            return true;
        }
        match flag {
            "-t" => self.tolerance = atoi(value),
            "-F" => self.frag_size = atoi(value),
            _ => {
                if self.has2d && self.set_option_2d(flag, value) {
                    return true;
                }
                if self.is_est && self.set_option_est(flag, value, scoring) {
                    return true;
                }
                return false;
            }
        }
        true
    }

    /// Port of `SetOption2D` (cdhit-common.c++:287-296).
    fn set_option_2d(&mut self, flag: &str, value: &str) -> bool {
        if self.set_option_common(flag, value) {
            return true;
        }
        match flag {
            "-i2" => self.input2 = value.to_string(),
            "-j2" => self.input2_pe = value.to_string(),
            "-s2" => self.diff_cutoff2 = atof(value),
            "-S2" => self.diff_cutoff_aa2 = atoi(value),
            _ => return false,
        }
        true
    }

    /// Port of `SetOptionEST` (cdhit-common.c++:297-317).
    fn set_option_est(&mut self, flag: &str, value: &str, scoring: &mut Scoring) -> bool {
        self.naa_top_limit = 12;
        if self.set_option_common(flag, value) {
            return true;
        }
        match flag {
            "-r" => self.option_r = atoi(value),
            "-gap" => scoring.mat.gap = MAX_SEQ * atoi(value),
            "-gap-ext" => scoring.mat.ext_gap = MAX_SEQ * atoi(value),
            "-match" => scoring.mat.set_match(atoi(value)),
            "-mismatch" => scoring.mat.set_mismatch(atoi(value)),
            "-mask" => {
                for ch in value.chars() {
                    let ch = ch.to_ascii_uppercase();
                    if !ch.is_ascii_uppercase() {
                        continue;
                    }
                    scoring.na2idx[(ch as u8 - b'A') as usize] = 5;
                }
                setaa_to_na(scoring);
            }
            _ => return false,
        }
        true
    }

    /// Parse a full argument vector (excluding argv[0]). Mirrors the pairwise
    /// `SetOptions` loop (cdhit-common.c++:346-347): flags and values alternate,
    /// and a trailing lone flag is an error.
    pub fn parse(
        args: &[String],
        twod: bool,
        est: bool,
        scoring: &mut Scoring,
    ) -> Result<Options> {
        let mut opt = Options::default();
        opt.has2d = twod;
        opt.is_est = est;
        if est {
            // cd-hit-est program defaults, applied before parsing so user
            // arguments override them (cdhit-est.c++:49-53).
            opt.cluster_thd = 0.95;
            opt.naa = 10;
            opt.naa_top_limit = 12;
        }
        opt.parse_args(args, scoring)?;
        Ok(opt)
    }

    /// Parse `args` (flag/value pairs) onto an already-seeded `Options`. Use
    /// this when a program needs custom defaults (e.g. `cd-hit-454`) applied
    /// before parsing so user arguments still override them.
    pub fn parse_args(&mut self, args: &[String], scoring: &mut Scoring) -> Result<()> {
        let mut i = 0;
        while i + 1 < args.len() {
            if !self.set_option(&args[i], &args[i + 1], scoring) {
                return Err(CdError::BadOption(args[i].clone()));
            }
            i += 2;
        }
        if i < args.len() {
            return Err(CdError::BadOption(args[i].clone()));
        }
        Ok(())
    }

    /// Port of `Options::Validate` (cdhit-common.c++:352-429). Returns an error
    /// instead of calling `bomb_error`/`exit`. Word-length "faster if" hints are
    /// emitted to `hints` for the caller to print (they go to stdout in C++).
    pub fn validate(&mut self, hints: &mut Vec<String>) -> Result<()> {
        if self.use_identity && self.use_distance {
            return Err(CdError::Validation(
                "can not use both identity cutoff and distance cutoff".into(),
            ));
        }
        if self.use_distance {
            if self.distance_thd > 1.0 || self.distance_thd < 0.0 {
                return Err(CdError::Validation("invalid distance threshold".into()));
            }
        } else if self.is_est {
            if self.cluster_thd > 1.0 || self.cluster_thd < 0.8 {
                return Err(CdError::Validation(
                    "invalid clstr threshold, should >=0.8".into(),
                ));
            }
        } else if self.cluster_thd > 1.0 || self.cluster_thd < 0.4 {
            return Err(CdError::Validation("invalid clstr".into()));
        }

        if self.input.is_empty() {
            return Err(CdError::Validation("no input file".into()));
        }
        if self.output.is_empty() {
            return Err(CdError::Validation("no output file".into()));
        }
        if self.pe_mode != 0 {
            if self.input_pe.is_empty() {
                return Err(CdError::Validation(
                    "no input file for R2 sequences in PE mode".into(),
                ));
            }
            if self.output_pe.is_empty() {
                return Err(CdError::Validation(
                    "no output file for R2 sequences in PE mode".into(),
                ));
            }
        }
        if self.is_est && self.align_pos == 1 {
            self.option_r = 0;
        }

        if self.band_width < 1 {
            return Err(CdError::Validation("invalid band width".into()));
        }
        if self.naa < 2 || self.naa > self.naa_top_limit {
            return Err(CdError::Validation("invalid word length".into()));
        }
        if self.des_len < 0 {
            return Err(CdError::Validation(
                "too short description, not enough to identify sequences".into(),
            ));
        }
        if !self.is_est && (self.tolerance < 0 || self.tolerance > 5) {
            return Err(CdError::Validation("invalid tolerance".into()));
        }
        if self.diff_cutoff < 0.0 || self.diff_cutoff > 1.0 {
            return Err(CdError::Validation("invalid value for -s".into()));
        }
        if self.diff_cutoff_aa < 0 {
            return Err(CdError::Validation("invalid value for -S".into()));
        }
        if self.has2d {
            if self.diff_cutoff2 < 0.0 || self.diff_cutoff2 > 1.0 {
                return Err(CdError::Validation("invalid value for -s2".into()));
            }
            if self.diff_cutoff_aa2 < 0 {
                return Err(CdError::Validation("invalid value for -S2".into()));
            }
            if self.pe_mode != 0 && self.input2_pe.is_empty() {
                return Err(CdError::Validation(
                    "no input file for R2 sequences for 2nd db in PE mode".into(),
                ));
            }
        }
        if !self.global_identity {
            self.print = 1;
        }
        if self.short_coverage < self.long_coverage {
            self.short_coverage = self.long_coverage;
        }
        if self.short_control > self.long_control {
            self.short_control = self.long_control;
        }
        if !self.global_identity && self.short_coverage == 0.0 && self.min_control == 0 {
            return Err(CdError::Validation(
                "You are using local identity, but no -aS -aL -A option".into(),
            ));
        }
        if self.frag_size < 0 {
            return Err(CdError::Validation("invalid fragment size".into()));
        }

        let message = |naa: i32, i: i32| format!("Your word length is {naa}, using {i} may be faster!");
        if !self.is_est && self.tolerance != 0 {
            let clstr_idx = (self.cluster_thd * 100.0) as i32 - NAA_STAT_START_PERCENT;
            let tcutoff = NAA_STAT[(self.tolerance - 1) as usize][clstr_idx as usize]
                [(5 - self.naa) as usize];
            if tcutoff < 5 {
                return Err(CdError::Validation(
                    "Too low cluster threshold for the word length.\nIncrease the threshold or the tolerance, or decrease the word length.".into(),
                ));
            }
            let mut i = 5;
            while i > self.naa {
                if NAA_STAT[(self.tolerance - 1) as usize][clstr_idx as usize][(5 - i) as usize] > 10
                {
                    hints.push(message(self.naa, i));
                    break;
                }
                i -= 1;
            }
        } else if self.is_est {
            if self.cluster_thd > 0.9 && self.naa < 8 {
                hints.push(message(self.naa, 8));
            } else if self.cluster_thd > 0.87 && self.naa < 5 {
                hints.push(message(self.naa, 5));
            } else if self.cluster_thd > 0.80 && self.naa < 4 {
                hints.push(message(self.naa, 4));
            } else if self.cluster_thd > 0.75 && self.naa < 3 {
                hints.push(message(self.naa, 3));
            }
        } else if self.cluster_thd > 0.85 && self.naa < 5 {
            hints.push(message(self.naa, 5));
        } else if self.cluster_thd > 0.80 && self.naa < 4 {
            hints.push(message(self.naa, 4));
        } else if self.cluster_thd > 0.75 && self.naa < 3 {
            hints.push(message(self.naa, 3));
        }

        if (self.min_length + 1) < self.naa {
            return Err(CdError::Validation("Too short -l, redefine it".into()));
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn args(v: &[&str]) -> Vec<String> {
        v.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn defaults_match_cpp() {
        let o = Options::default();
        assert_eq!(o.naa, 5);
        assert_eq!(o.cluster_thd, 0.9);
        assert_eq!(o.band_width, 20);
        assert_eq!(o.max_memory, 800_000_000);
        assert!(o.global_identity);
        assert_eq!(o.diff_cutoff_aa, 99_999_999);
    }

    #[test]
    fn parse_basic_protein() {
        let mut sc = Scoring::default();
        let o = Options::parse(
            &args(&["-i", "in.fa", "-o", "out", "-c", "0.9", "-n", "5"]),
            false,
            false,
            &mut sc,
        )
        .unwrap();
        assert_eq!(o.input, "in.fa");
        assert_eq!(o.output, "out");
        assert_eq!(o.cluster_thd, 0.9);
        assert!(o.use_identity);
        assert_eq!(o.naa, 5);
    }

    #[test]
    fn trailing_flag_is_error() {
        let mut sc = Scoring::default();
        let r = Options::parse(&args(&["-i", "in.fa", "-o"]), false, false, &mut sc);
        assert!(r.is_err());
    }

    #[test]
    fn validate_rejects_low_threshold_protein() {
        let mut sc = Scoring::default();
        let mut o = Options::parse(
            &args(&["-i", "in.fa", "-o", "out", "-c", "0.3"]),
            false,
            false,
            &mut sc,
        )
        .unwrap();
        let mut hints = vec![];
        assert!(o.validate(&mut hints).is_err());
    }

    #[test]
    fn validate_accepts_default_protein() {
        let mut sc = Scoring::default();
        let mut o = Options::parse(
            &args(&["-i", "in.fa", "-o", "out", "-c", "0.9", "-n", "5"]),
            false,
            false,
            &mut sc,
        )
        .unwrap();
        let mut hints = vec![];
        o.validate(&mut hints).unwrap();
    }

    #[test]
    fn atoi_atof_c_semantics() {
        assert_eq!(atoi("12abc"), 12);
        assert_eq!(atoi("abc"), 0);
        assert_eq!(atoi("-7"), -7);
        assert_eq!(atof("0.9x"), 0.9);
        assert_eq!(atof("abc"), 0.0);
    }
}

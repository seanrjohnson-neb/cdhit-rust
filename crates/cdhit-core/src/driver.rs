//! High-level driver tying the pipeline together for the 1D clustering
//! programs (`cd-hit`, `cd-hit-est`). Mirrors the `main()` flow of `cdhit.c++`:
//! parse options → InitNAA → read → sort_divide → do_clustering → write.

use crate::alphabet::{Naa, MAX_UAA};
use crate::error::Result;
use crate::est::make_comp_short_word_index;
use crate::io::{read_database, write_clusters, write_extra_1d, write_extra_2d};
use crate::options::{setaa_to_na, Options, Scoring};
use crate::seqdb::SequenceDb;
use crate::sequence::{IS_REP};

/// The 1D clustering programs.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Program {
    /// `cd-hit` — protein.
    CdHit,
    /// `cd-hit-est` — nucleotide.
    CdHitEst,
    /// `cd-hit-454` — 454 pyrosequencing duplicate detection.
    CdHit454,
}

/// 454 nucleotide scoring matrix `myBLOSUM62_na2` (cdhit-454.c++:40-49).
#[rustfmt::skip]
const MY_BLOSUM62_NA2: [i32; 21] = [
     2,
    -1, 2,
    -1,-1, 2,
    -1,-1,-1, 2,
    -1,-1,-1, 2, 2,
    -1,-1,-1,-1,-1, 2,
];

/// Output of a clustering run.
pub struct ClusterOutput {
    /// Representative FASTA (the `-o` file contents).
    pub rep_fasta: Vec<u8>,
    /// `.clstr` membership file contents.
    pub clstr: String,
    /// Number of clusters.
    pub num_clusters: usize,
    /// Informational hints from option validation (printed to stdout in C++).
    pub hints: Vec<String>,
}

/// Run a 1D clustering job on in-memory input.
///
/// `args` are the CLI arguments after argv[0] (flag/value pairs). Backwards
/// compatible helper: `est == true` selects `cd-hit-est`, else `cd-hit`.
pub fn cluster_1d(input: &[u8], args: &[String], est: bool) -> Result<ClusterOutput> {
    cluster_1d_program(
        input,
        args,
        if est { Program::CdHitEst } else { Program::CdHit },
    )
}

/// Run a 1D clustering job for a specific program.
pub fn cluster_1d_program(
    input: &[u8],
    args: &[String],
    program: Program,
) -> Result<ClusterOutput> {
    let est = program != Program::CdHit;
    let mut scoring = Scoring::default();
    if est {
        // Nucleotide alphabet + matrix (cd-hit-est.c++:54-55).
        setaa_to_na(&mut scoring);
        scoring.mat.set_to_na();
    }

    // Program-specific defaults applied before parsing (so args override).
    let mut options = if program == Program::CdHit454 {
        // cd-hit-454.c++:59-74 program defaults.
        let mut o = Options::default();
        o.is_est = true;
        o.is454 = true;
        o.naa = 10;
        o.naa_top_limit = 12;
        o.cluster_thd = 0.98;
        o.band_width = 10;
        o.print = 1;
        o.des_len = 0;
        o.option_r = 0;
        scoring.mat.set_gap(-3, -1);
        scoring.mat.set_matrix(&MY_BLOSUM62_NA2);
        o.parse_args(args, &mut scoring)?;
        o
    } else {
        Options::parse(args, false, est, &mut scoring)?
    };

    let mut hints = Vec::new();
    options.validate(&mut hints)?;

    // InitNAA: protein uses base MAX_UAA (21); EST/454 use base 4.
    let naa = if est { Naa::init(4) } else { Naa::init(MAX_UAA) };
    options.naan = naa.array[options.naa as usize];

    let database = read_database(input, &options);
    let raw = database.raw;
    let mut db = SequenceDb::new(database.sequences);
    db.naan = options.naan;

    db.sort_divide(&mut options, &scoring, true);
    db.do_clustering(&options, &scoring, &naa);

    let rep_fasta = write_clusters(&db.sequences, &db.rep_seqs, &raw);
    let clstr = write_extra_1d(&db.sequences, &db.rep_seqs, &options);
    let num_clusters = db.rep_seqs.len();

    Ok(ClusterOutput {
        rep_fasta,
        clstr,
        num_clusters,
        hints,
    })
}

/// `MAX_BIN_SWAP` (cdhit-common.h:61).
const MAX_BIN_SWAP: u64 = 2_000_000_000;

/// Run `cd-hit-div`: sort the input longest-first and split it into `n`
/// roughly-equal segments by letter count (port of `DivideSave`,
/// cdhit-common.c++:2302+). Returns the segment file contents in order
/// (the C++ writes them as `<output>-0`, `<output>-1`, ...).
pub fn divide(input: &[u8], n: i32) -> Vec<Vec<u8>> {
    // cd-hit-div uses default Options (no SetOptions), so default min_length.
    let mut options = Options::default();
    let scoring = Scoring::default();
    let naa = Naa::init(MAX_UAA);
    options.naan = naa.array[options.naa as usize];

    let d = read_database(input, &options);
    let raw = d.raw;
    let mut db = SequenceDb::new(d.sequences);
    db.naan = options.naan;
    db.sort_divide(&mut options, &scoring, true);

    if n == 0 || db.sequences.is_empty() {
        return Vec::new();
    }
    let mut max_seg = db.total_letter as u64 / n as u64 + db.sequences[0].size as u64;
    if max_seg >= MAX_BIN_SWAP {
        max_seg = MAX_BIN_SWAP;
    }

    let mut segments: Vec<Vec<u8>> = vec![Vec::new()];
    let mut seg_size: u64 = 0;
    for seq in &db.sequences {
        seg_size += seq.size as u64;
        if seg_size >= max_seg {
            segments.push(Vec::new());
            seg_size = seq.size as u64;
        }
        let begin = seq.des_begin as usize;
        let end = (begin + seq.tot_length as usize).min(raw.len());
        segments.last_mut().unwrap().extend_from_slice(&raw[begin..end]);
    }
    segments
}

/// Run a 2D clustering job (`cd-hit-2d` / `cd-hit-est-2d`): cluster db2 against
/// a fixed reference db1. `input1` is db1 (`-i`), `input2` is db2 (`-i2`). The
/// `-o` representatives are the novel db2 sequences.
pub fn cluster_2d(
    input1: &[u8],
    input2: &[u8],
    args: &[String],
    est: bool,
) -> Result<ClusterOutput> {
    let mut scoring = Scoring::default();
    if est {
        setaa_to_na(&mut scoring);
        scoring.mat.set_to_na();
    }
    let mut options = Options::parse(args, true, est, &mut scoring)?;
    let mut hints = Vec::new();
    options.validate(&mut hints)?;

    let naa = if est { Naa::init(4) } else { Naa::init(MAX_UAA) };
    options.naan = naa.array[options.naa as usize];

    let comp_aan_idx: Vec<i32> = if est && options.option_r != 0 {
        make_comp_short_word_index(options.naa, &naa)
    } else {
        Vec::new()
    };

    // db1 (reference) and db2 (query).
    let d1 = read_database(input1, &options);
    let d2 = read_database(input2, &options);
    let raw2 = d2.raw;
    let mut db1 = SequenceDb::new(d1.sequences);
    let mut db2 = SequenceDb::new(d2.sequences);
    db1.naan = options.naan;
    db2.naan = options.naan;

    db1.sort_divide(&mut options, &scoring, true);
    db2.sort_divide(&mut options, &scoring, true);

    // db1 sequences are the fixed representatives: cluster_id = sorted index.
    for (k, seq) in db1.sequences.iter_mut().enumerate() {
        seq.cluster_id = k as i32;
        seq.state |= IS_REP;
    }

    db2.cluster_to(&db1, &options, &scoring, &naa, &comp_aan_idx);

    let rep_fasta = write_clusters(&db2.sequences, &db2.rep_seqs, &raw2);
    let clstr = write_extra_2d(&db1.sequences, &db2.sequences, &options);
    let num_clusters = db1.sequences.len();

    Ok(ClusterOutput {
        rep_fasta,
        clstr,
        num_clusters,
        hints,
    })
}

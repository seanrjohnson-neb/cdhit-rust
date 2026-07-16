//! `cdhit-core` — a Rust port of the CD-HIT sequence-clustering engine.
//!
//! The port targets bit-for-bit output parity with the original C++
//! (`original_code/cdhit`). Modules are added in dependency order; see the
//! project plan for the phased milestones.

pub mod align;
pub mod alphabet;
pub mod auxtools;
pub mod buffer;
pub mod clstr;
pub mod cluster;
pub mod cutoff;
pub mod driver;
pub mod est;
pub mod error;
pub mod io;
pub mod naa_stat;
pub mod options;
pub mod seqdb;
pub mod sequence;
pub mod wordtable;

pub use auxtools::{
    cd_hit_dup, cd_hit_lap, read_linker, DupOutput, DupParams, LapOutput, LapParams, LinkerOutput,
    LinkerParams,
};
pub use driver::{cluster_1d, cluster_1d_program, cluster_2d, divide, ClusterOutput, Program};
pub use error::{CdError, Result};
pub use options::{Options, Scoring};
pub use sequence::Sequence;

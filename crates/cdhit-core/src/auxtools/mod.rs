//! Ports of the `cd-hit-auxtools` package (`cd-hit-dup`, `cd-hit-lap`,
//! `read-linker`).
//!
//! These tools are self-contained in the C++ (they use their own `mintlib`
//! containers and `bioSequence` I/O rather than `cdhit-common`), so the Rust
//! port lives in its own module tree and shares only the FASTA/FASTQ handling
//! in [`bioseq`].

pub mod bioseq;
pub mod dup;
pub mod dup_chimera;
pub mod lap;
pub mod linker;

pub use bioseq::Sequence;
pub use dup::{cd_hit_dup, DupOutput, DupParams};
pub use lap::{cd_hit_lap, LapOutput, LapParams};
pub use linker::{read_linker, LinkerOutput, LinkerParams};

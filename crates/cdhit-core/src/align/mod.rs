//! Alignment: the diagonal pre-filter and the banded DP.

pub mod band;
pub mod diag;

pub use band::{local_band_align, BandAlign};
pub use diag::{diag_test_aapn, diag_test_aapn_est, DiagResult};

//! Error type for the engine. The C++ calls `bomb_error`/`exit(1)`; we return
//! `Result` so the library is embeddable (CLI maps to an exit code, WASM to a
//! JS exception).

use std::fmt;

#[derive(Debug, Clone, PartialEq)]
pub enum CdError {
    /// Unrecognised or malformed command-line flag.
    BadOption(String),
    /// Option validation failure (message mirrors the C++ `bomb_error` text).
    Validation(String),
    /// Input parsing / format error.
    Parse(String),
    /// I/O error (message from the underlying `std::io::Error`).
    Io(String),
}

impl fmt::Display for CdError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            CdError::BadOption(s) => write!(f, "unknown or malformed option: {s}"),
            CdError::Validation(s) => write!(f, "{s}"),
            CdError::Parse(s) => write!(f, "parse error: {s}"),
            CdError::Io(s) => write!(f, "io error: {s}"),
        }
    }
}

impl std::error::Error for CdError {}

impl From<std::io::Error> for CdError {
    fn from(e: std::io::Error) -> Self {
        CdError::Io(e.to_string())
    }
}

pub type Result<T> = std::result::Result<T, CdError>;

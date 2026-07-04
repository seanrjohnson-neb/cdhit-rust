//! `cdhit` CLI — a thin dispatcher over the `cdhit-core` engine.
//!
//! CD-HIT uses non-GNU flags (`-aL`, `-s2`, `-cx`, …), so rather than force
//! them through a standard argument parser we pass the argument list verbatim
//! to the core option parser, preserving exact upstream semantics. The first
//! argument selects the program:
//!
//! ```text
//! cdhit cd-hit     -i in.fa -o out -c 0.9 ...
//! cdhit cd-hit-est -i in.fa -o out -c 0.95 ...
//! cdhit cd-hit-454 -i in.fa -o out ...
//! ```
//!
//! For drop-in compatibility, invoking the binary under the historical names
//! (`cd-hit`, `cd-hit-est`, `cd-hit-454`, e.g. via a symlink) also works.

use std::io::Read;
use std::process::exit;

use cdhit_core::{cluster_1d_program, cluster_2d, divide, Program};

/// Which program the CLI is running.
#[derive(Clone, Copy, PartialEq)]
enum Prog {
    OneD(Program),
    TwoD { est: bool },
    Div,
}

fn usage() -> ! {
    eprintln!(
        "Usage: cdhit <program> [options]\n\
         \n\
         Programs:\n\
         \x20 cd-hit         protein clustering\n\
         \x20 cd-hit-est     nucleotide clustering\n\
         \x20 cd-hit-2d      protein: compare db2 (-i2) against db1 (-i)\n\
         \x20 cd-hit-est-2d  nucleotide 2D comparison\n\
         \x20 cd-hit-454     454 read duplicate detection\n\
         \x20 cd-hit-div     split a database into -div segments\n\
         \n\
         Run a clustering program with -i <input> -o <output> -c <threshold> ...\n\
         Output is written to <output> (representatives) and <output>.clstr."
    );
    exit(1);
}

/// Map a program/subcommand token to a `Prog`. Note `cdhit` (the dispatcher
/// binary name) is intentionally NOT mapped here — it requires an explicit
/// subcommand argument.
fn program_from_name(name: &str) -> Option<Prog> {
    match name {
        "cd-hit" => Some(Prog::OneD(Program::CdHit)),
        "cd-hit-est" => Some(Prog::OneD(Program::CdHitEst)),
        "cd-hit-454" => Some(Prog::OneD(Program::CdHit454)),
        "cd-hit-2d" => Some(Prog::TwoD { est: false }),
        "cd-hit-est-2d" => Some(Prog::TwoD { est: true }),
        "cd-hit-div" => Some(Prog::Div),
        _ => None,
    }
}

/// Read a possibly-gzipped input file into memory.
fn read_input(path: &str) -> std::io::Result<Vec<u8>> {
    let bytes = std::fs::read(path)?;
    if path.ends_with(".gz") {
        let mut d = flate2::read::GzDecoder::new(&bytes[..]);
        let mut out = Vec::new();
        d.read_to_end(&mut out)?;
        Ok(out)
    } else {
        Ok(bytes)
    }
}

fn main() {
    let argv: Vec<String> = std::env::args().collect();

    // Determine the program from argv[0]'s basename (drop-in names) or argv[1].
    let exe = argv
        .first()
        .and_then(|p| p.rsplit(['/', '\\']).next())
        .unwrap_or("");
    let (program, rest) = if let Some(p) = program_from_name(exe) {
        (p, &argv[1..])
    } else if argv.len() >= 2 {
        match program_from_name(&argv[1]) {
            Some(p) => (p, &argv[2..]),
            None => usage(),
        }
    } else {
        usage();
    };

    // Extract -i / -i2 / -o / -div for file handling.
    let mut input = String::new();
    let mut input2 = String::new();
    let mut output = String::new();
    let mut div_n = 0i32;
    let mut i = 0;
    while i + 1 < rest.len() {
        match rest[i].as_str() {
            "-i" => input = rest[i + 1].clone(),
            "-i2" => input2 = rest[i + 1].clone(),
            "-o" => output = rest[i + 1].clone(),
            "-div" => div_n = rest[i + 1].parse().unwrap_or(0),
            _ => {}
        }
        i += 2;
    }
    if input.is_empty() || output.is_empty() {
        usage();
    }
    let read = |p: &str| -> Vec<u8> {
        match read_input(p) {
            Ok(d) => d,
            Err(e) => {
                eprintln!("Failed to read input '{p}': {e}");
                exit(1);
            }
        }
    };

    // cd-hit-div: split into <output>-0, <output>-1, ...
    if program == Prog::Div {
        if div_n <= 1 {
            eprintln!("Warning: -div must be greater than 1.");
        }
        let data = read(&input);
        for (k, seg) in divide(&data, div_n).into_iter().enumerate() {
            let path = format!("{output}-{k}");
            if let Err(e) = std::fs::write(&path, seg) {
                eprintln!("Failed to write '{path}': {e}");
                exit(1);
            }
        }
        return;
    }

    let data = read(&input);
    let result = match program {
        Prog::OneD(p) => cluster_1d_program(&data, rest, p),
        Prog::TwoD { est } => {
            if input2.is_empty() {
                eprintln!("cd-hit-2d requires -i2 <second database>");
                exit(1);
            }
            let data2 = read(&input2);
            cluster_2d(&data, &data2, rest, est)
        }
        Prog::Div => unreachable!(),
    };

    match result {
        Ok(out) => {
            for h in &out.hints {
                eprintln!("{h}");
            }
            if let Err(e) = std::fs::write(&output, &out.rep_fasta) {
                eprintln!("Failed to write '{output}': {e}");
                exit(1);
            }
            if let Err(e) = std::fs::write(format!("{output}.clstr"), out.clstr) {
                eprintln!("Failed to write '{output}.clstr': {e}");
                exit(1);
            }
            eprintln!("{} clusters", out.num_clusters);
        }
        Err(e) => {
            eprintln!("\nFatal Error:\n{e}\nProgram halted !!\n");
            exit(1);
        }
    }
}

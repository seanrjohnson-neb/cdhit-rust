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

use std::io::{Read, Write};
use std::process::exit;

use cdhit_core::clstr::ops as clstr_ops;
use cdhit_core::{
    cd_hit_dup, cd_hit_lap, cluster_1d_program, cluster_2d, divide, read_linker, DupParams,
    LapParams, LinkerParams, Program,
};

/// Which program the CLI is running.
#[derive(Clone, Copy, PartialEq)]
enum Prog {
    OneD(Program),
    TwoD { est: bool },
    Div,
    Linker,
    Lap,
    Dup,
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
         \x20 cd-hit-dup     duplicate / near-duplicate read detection (auxtools)\n\
         \x20 cd-hit-lap     cluster overlapping reads (auxtools)\n\
         \x20 read-linker    join paired-end reads by overlap (auxtools)\n\
         \n\
         Run a clustering program with -i <input> -o <output> -c <threshold> ...\n\
         Output is written to <output> (representatives) and <output>.clstr.\n\
         The auxtools programs use their own flags; run one with no args for help."
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
        "read-linker" => Some(Prog::Linker),
        "cd-hit-lap" => Some(Prog::Lap),
        "cd-hit-dup" => Some(Prog::Dup),
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

const READ_LINKER_HELP: &str = "Options:\n\
    \x20   -1 file       Input file, first end;\n\
    \x20   -2 file       Input file, second end;\n\
    \x20   -o file       Output file;\n\
    \x20   -l number     Minimum overlapping length (default 10);\n\
    \x20   -e number     Maximum number of errors (mismatches, default 1);\n";

/// `read-linker`: join paired-end reads by overlap. Argument parsing mirrors
/// the C++ (`argv` scanned in pairs; at least three options required).
fn run_read_linker(rest: &[String]) {
    if rest.len() < 6 {
        print!("{READ_LINKER_HELP}\n");
        exit(1);
    }
    let mut first = String::new();
    let mut second = String::new();
    let mut output = String::new();
    let mut min = 10i32;
    let mut error = 1i32;
    let mut maxlen = 0i32;
    let mut i = 0;
    while i < rest.len() {
        if i + 1 == rest.len() {
            print!("Incomplete argument {}\n\n{READ_LINKER_HELP}\n", rest[i]);
            exit(1);
        }
        let val = &rest[i + 1];
        match rest[i].as_str() {
            "-1" => first = val.clone(),
            "-2" => second = val.clone(),
            "-o" => output = val.clone(),
            "-l" => min = val.parse().unwrap_or(0),
            "-e" => error = val.parse().unwrap_or(0),
            "-m" => maxlen = val.parse().unwrap_or(0),
            other => {
                print!("Unknown argument {other}\n\n{READ_LINKER_HELP}\n");
                exit(1);
            }
        }
        i += 2;
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
    let d1 = read(&first);
    let d2 = read(&second);
    let out = read_linker(&d1, &d2, &LinkerParams { min, error, maxlen });
    if let Err(e) = std::fs::write(&output, &out.output) {
        eprintln!("Failed to write '{output}': {e}");
        exit(1);
    }
    print!("{}", out.log);
}

const CD_HIT_LAP_HELP: &str = "Options:\n\
    \x20   -i        Input file;\n\
    \x20   -o        Output file;\n\
    \x20   -m        Minimum length of overlapping part (default 20);\n\
    \x20   -p        Minimum percentage of overlapping part (default 0, any percentage);\n\
    \x20   -d        Description length (default 0, truncate at the first whitespace character)\n\
    \x20   -s        Random number seed for shuffling (default 0, no shuffling);\n\
    \x20   -stdout   Standard output type (default \"log\", other options \"rep\", \"clstr\");\n";

/// `cd-hit-lap`: cluster overlapping reads. Argument parsing mirrors the C++.
fn run_cd_hit_lap(rest: &[String]) {
    if rest.len() < 4 {
        print!("{CD_HIT_LAP_HELP}\n");
        exit(1);
    }
    let mut input = String::new();
    let mut output = String::new();
    let mut stdout_type = String::from("log");
    let mut minlen = 20i32;
    let mut minper = 0f32;
    let mut deslen = 0i32;
    let mut seed = 0u32;
    let mut i = 0;
    while i < rest.len() {
        if i + 1 == rest.len() {
            print!("Incomplete argument {}\n\n{CD_HIT_LAP_HELP}\n", rest[i]);
            exit(1);
        }
        let val = &rest[i + 1];
        match rest[i].as_str() {
            "-i" => input = val.clone(),
            "-o" => output = val.clone(),
            "-stdout" => stdout_type = val.clone(),
            "-m" => minlen = val.parse().unwrap_or(0),
            "-p" => minper = val.parse().unwrap_or(0.0),
            "-d" => deslen = val.parse().unwrap_or(0),
            "-s" => seed = val.parse().unwrap_or(0),
            other => {
                print!("Unknown argument {other}\n\n{CD_HIT_LAP_HELP}\n");
                exit(1);
            }
        }
        i += 2;
    }
    if !matches!(stdout_type.as_str(), "log" | "rep" | "clstr") {
        print!("Unknown standard output type {stdout_type}\n");
        exit(1);
    }

    let data = match read_input(&input) {
        Ok(d) => d,
        Err(_) => {
            print!("File openning failed: {input}\n");
            exit(1);
        }
    };
    let out = cd_hit_lap(
        &data,
        &LapParams {
            minlen,
            minper,
            deslen,
            seed,
        },
    );

    // Route rep / clstr / log to files or stdout per -stdout, matching the C++.
    let write_file = |path: String, bytes: &[u8]| {
        if let Err(e) = std::fs::write(&path, bytes) {
            eprintln!("Failed to write '{path}': {e}");
            exit(1);
        }
    };
    match stdout_type.as_str() {
        "log" => {
            write_file(output.clone(), &out.rep);
            write_file(format!("{output}.clstr"), &out.clstr);
            print!("{}", out.log);
        }
        "rep" => {
            write_file(format!("{output}.log"), out.log.as_bytes());
            write_file(format!("{output}.clstr"), &out.clstr);
            std::io::stdout().write_all(&out.rep).ok();
        }
        "clstr" => {
            write_file(format!("{output}.log"), out.log.as_bytes());
            write_file(output.clone(), &out.rep);
            std::io::stdout().write_all(&out.clstr).ok();
        }
        _ => unreachable!(),
    }
}

const CD_HIT_DUP_HELP: &str = "CD-HIT-DUP\n\n\
Usage:\n\
cd-hit-dup -i input.fa -o output [other options] (for single reads FASTQ)\n\
cd-hit-dup -i input.fq -o output [other options] (for single reads FASTA)\n\
cd-hit-dup -i R1.fq -i2 R2.fq -o output -o2 output-R2 [other options] (for PE reads FASTQ)\n\
cd-hit-dup -i R1.fa -i2 R2.fa -o output -o2 output-R2 [other options] (for PE reads FASTA)\n\n\
Options:\n\
    \x20   -i        Input file (FASTQ or FASTA);\n\
    \x20   -i2       Second input file (FASTQ or FASTA);\n\
    \x20   -o        Output file;\n\
    \x20   -o2       Output file for R2;\n\
    \x20   -d        Description length (default 0, truncate at the first whitespace character)\n\
    \x20   -u        Length of prefix to be used in the analysis (default 0, for full/maximum length);\n\
    \x20   -m        Match length (true/false, default true);\n\
    \x20   -e        Maximum number of mismatches allowd;\n\
    \x20   -f        Filter out chimeric clusters (true/false, default false);\n\
    \x20   -s        Minimum length of common sequence shared between a chimeric read\n\
    \x20             and each of its parents (default 30, minimum 20);\n\
    \x20   -a        Abundance cutoff (default 1 without chimeric filtering, 2 with chimeric filtering);\n\
    \x20   -b        Abundance ratio between a parent read and chimeric read (default 1);\n\
    \x20   -p        Dissimilarity control for chimeric filtering (default 1);\n";

/// `cd-hit-dup`: duplicate / near-duplicate read detection. Argument parsing
/// mirrors the C++.
fn run_cd_hit_dup(rest: &[String]) {
    if rest.len() < 4 {
        print!("{CD_HIT_DUP_HELP}\n");
        exit(1);
    }
    let mut input = String::new();
    let mut input2 = String::new();
    let mut output = String::new();
    let mut output2 = String::new();
    let mut match_length = true;
    let mut nochimeric = false;
    let mut abundance = -1i32;
    let mut deslen = 0i32;
    let mut shared = 30i32;
    let mut uselen = 0i32;
    let mut errors = 0i32;
    let mut errors2 = 0f32;
    let mut abratio = 1f32;
    let mut percent = 1f32;
    let mut i = 0;
    while i < rest.len() {
        if i + 1 == rest.len() {
            print!("Incomplete argument {}\n\n{CD_HIT_DUP_HELP}\n", rest[i]);
            exit(1);
        }
        let v = &rest[i + 1];
        match rest[i].as_str() {
            "-i" => input = v.clone(),
            "-i2" => input2 = v.clone(),
            "-o" => output = v.clone(),
            "-o2" => output2 = v.clone(),
            "-a" => abundance = v.parse().unwrap_or(0),
            "-d" => deslen = v.parse().unwrap_or(0),
            "-u" => uselen = v.parse().unwrap_or(0),
            "-b" => abratio = v.parse().unwrap_or(0.0),
            "-p" => percent = v.parse().unwrap_or(0.0),
            "-m" => match_length = v == "true",
            "-f" => nochimeric = v == "true",
            "-e" => {
                errors = v.parse().unwrap_or(0);
                errors2 = v.parse().unwrap_or(0.0);
            }
            "-s" => {
                shared = v.parse().unwrap_or(0);
                nochimeric = true;
                print!("Chimeric cluster filtering is automatically enabled by \"-s\" parameter!\n");
            }
            other => {
                print!("Unknown argument {other}\n\n{CD_HIT_DUP_HELP}\n");
                exit(1);
            }
        }
        i += 2;
    }
    if !input2.is_empty() && output2.is_empty() {
        output2 = format!("{output}.2");
    }

    let read = |p: &str| -> Vec<u8> {
        match read_input(p) {
            Ok(d) => d,
            Err(_) => {
                print!("File openning failed: {p}\n");
                exit(1);
            }
        }
    };
    let d1 = read(&input);
    let d2 = if input2.is_empty() { None } else { Some(read(&input2)) };

    let params = DupParams {
        input_name: &input,
        input2: d2.as_deref().map(|b| (b, input2.as_str())),
        match_length,
        abundance,
        deslen,
        uselen,
        errors,
        errors2,
        nochimeric,
        shared,
        abratio,
        percent,
    };
    match cd_hit_dup(&d1, &params) {
        Ok(out) => {
            let write_file = |path: String, bytes: &[u8]| {
                if let Err(e) = std::fs::write(&path, bytes) {
                    eprintln!("Failed to write '{path}': {e}");
                    exit(1);
                }
            };
            write_file(output.clone(), &out.reps_r1);
            write_file(format!("{output}.clstr"), &out.clstr);
            write_file(format!("{output}2.clstr"), &out.clstr2);
            if let Some(r2) = &out.reps_r2 {
                write_file(output2.clone(), r2);
            }
            print!("{}", out.log);
        }
        Err(msg) => {
            print!("{msg}\n");
            exit(1);
        }
    }
}

/// Read `.clstr`-style text from a file path, or from stdin when `path` is None.
fn read_clstr_input(path: Option<&str>) -> Vec<u8> {
    match path {
        Some(p) => match std::fs::read(p) {
            Ok(d) => d,
            Err(_) => {
                eprintln!("Can not open file {p}");
                exit(1);
            }
        },
        None => {
            let mut buf = Vec::new();
            std::io::stdin().read_to_end(&mut buf).ok();
            buf
        }
    }
}

/// Handle the `clstr_*` post-processing subcommands (ports of the Perl
/// scripts). Returns true if `name` matched one of them. These read `.clstr`
/// from a file argument (or stdin) and write the result to stdout.
fn try_run_clstr(name: &str, args: &[String]) -> bool {
    let mut stdout = std::io::stdout();
    match name {
        "clstr_sort_by" | "clstr-sort-by" => {
            // clstr_sort_by.pl <no|len> < file   (key is the first positional)
            let key = args.first().map(|s| s.as_str()).unwrap_or("no");
            let file = args.get(1).map(|s| s.as_str());
            let input = read_clstr_input(file);
            stdout.write_all(&clstr_ops::sort_by(&input, key)).ok();
        }
        "clstr_size_stat" | "clstr-size-stat" => {
            if args.is_empty() {
                print!("Usage:\n\tclstr_size_stat.pl clstr_file\n");
                exit(1);
            }
            let input = read_clstr_input(Some(&args[0]));
            stdout.write_all(&clstr_ops::size_stat(&input)).ok();
        }
        "clstr_size_histogram" | "clstr-size-histogram" => {
            if args.is_empty() {
                print!("Usage:\n\tclstr_size_histogram.pl [-bin N] clstr_file\n");
                exit(1);
            }
            let (step, file) = if args[0] == "-bin" {
                (args.get(1).and_then(|s| s.parse().ok()).unwrap_or(100), args.get(2))
            } else {
                (100i64, args.first())
            };
            let input = read_clstr_input(file.map(|s| s.as_str()));
            stdout.write_all(&clstr_ops::size_histogram(&input, step)).ok();
        }
        "clstr2txt" | "clstr2txt.pl" => {
            let input = read_clstr_input(args.first().map(|s| s.as_str()));
            stdout.write_all(&clstr_ops::to_txt(&input)).ok();
        }
        _ => return false,
    }
    true
}

fn main() {
    let argv: Vec<String> = std::env::args().collect();

    // clstr_* post-processing subcommands (dispatched by argv[1]).
    if argv.len() >= 2 && try_run_clstr(&argv[1], &argv[2..]) {
        return;
    }

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

    // Aux tools have their own flag conventions; dispatch them before the
    // clustering-oriented -i/-o extraction below.
    if program == Prog::Linker {
        run_read_linker(rest);
        return;
    }
    if program == Prog::Lap {
        run_cd_hit_lap(rest);
        return;
    }
    if program == Prog::Dup {
        run_cd_hit_dup(rest);
        return;
    }

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
        Prog::Div | Prog::Linker | Prog::Lap | Prog::Dup => unreachable!(),
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

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
         The auxtools programs use their own flags; run one with no args for help.\n\
         \n\
         .clstr post-processing (read a file arg or stdin, write stdout):\n\
         \x20 clstr_sort_by, clstr_sort_prot_by, clstr_size_stat, clstr_size_histogram,\n\
         \x20 clstr2txt, clstr_renumber, clstr_select, clstr_select_rep, clstr_cut,\n\
         \x20 clstr_rep, clstr2tree, cd-hit-clstr_2_blm8, clstr_reduce, clstr_rev,\n\
         \x20 clstr_merge, clstr_merge_noorder, clstr_reps_faa_rev, plot_len1,\n\
         \x20 clstr_quality_eval_by_link, clstr_quality_eval, clstr2xml,\n\
         \x20 clstr_sql_tbl_sort.\n\
         \n\
         .clstr post-processing that writes files (native only):\n\
         \x20 clstr_sql_tbl <clstr> <tbl>, make_multi_seq <fasta> <clstr> <dir> [size],\n\
         \x20 cd-hit-dup-PE-out -i R1 -j R2 -c clstr -o out1 -p out2."
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

const DUP_PE_OUT_USAGE: &str = "This script exports the representative PE reads into two seperate files after running\n\
cd-hit-dup\n \n\
    \x20    -i fasta or fastq file of PE read 1\n\
    \x20    -j fasta or fastq file of PE read 2\n\
    \x20    -c .clstr file produced by cd-hit-dup\n\
    \x20    -o output file of representative reads, PE read 1\n\
    \x20    -p output file of representative reads, PE read 2\n";

/// `cd-hit-dup-PE-out`: export representative paired-end reads to two files,
/// given the `.clstr` from `cd-hit-dup`. Flags mirror the Perl `getopts`.
fn run_dup_pe_out(args: &[String]) {
    let mut fastq1 = String::new();
    let mut fastq2 = String::new();
    let mut clstr_file = String::new();
    let mut out1 = String::new();
    let mut out2 = String::new();
    let mut i = 0;
    while i + 1 < args.len() {
        let v = &args[i + 1];
        match args[i].as_str() {
            "-i" => fastq1 = v.clone(),
            "-j" => fastq2 = v.clone(),
            "-c" => clstr_file = v.clone(),
            "-o" => out1 = v.clone(),
            "-p" => out2 = v.clone(),
            _ => {}
        }
        i += 2;
    }
    if fastq1.is_empty()
        || fastq2.is_empty()
        || out1.is_empty()
        || out2.is_empty()
        || clstr_file.is_empty()
    {
        print!("{DUP_PE_OUT_USAGE}");
        exit(255);
    }

    let read = |p: &str| -> Vec<u8> {
        match read_input(p) {
            Ok(d) => d,
            Err(_) => {
                eprintln!("can not open {p}");
                exit(1);
            }
        }
    };
    let clstr = read(&clstr_file);
    let in1 = read(&fastq1);
    let in2 = read(&fastq2);
    let out = clstr_ops::dup_pe_out(&clstr, &in1, &in2);
    if let Err(e) = std::fs::write(&out1, &out.out1) {
        eprintln!("can not write to {out1}: {e}");
        exit(1);
    }
    if let Err(e) = std::fs::write(&out2, &out.out2) {
        eprintln!("can not write to {out2}: {e}");
        exit(1);
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
        "clstr_renumber" | "clstr-renumber" => {
            let input = read_clstr_input(args.first().map(|s| s.as_str()));
            stdout.write_all(&clstr_ops::renumber(&input)).ok();
        }
        "clstr_select" | "clstr-select" => {
            let min: usize = args.first().and_then(|s| s.parse().ok()).unwrap_or(0);
            let max: usize = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(0);
            let input = read_clstr_input(args.get(2).map(|s| s.as_str()));
            stdout.write_all(&clstr_ops::select(&input, min, max)).ok();
        }
        "clstr_cut" | "clstr-cut" => {
            let n: i64 = match args.first().and_then(|s| s.parse().ok()) {
                Some(n) if n != 0 => n,
                _ => {
                    print!("no number\n");
                    exit(1);
                }
            };
            let input = read_clstr_input(args.get(1).map(|s| s.as_str()));
            stdout.write_all(&clstr_ops::cut(&input, n)).ok();
        }
        "clstr_rep" | "clstr-rep" => {
            let input = read_clstr_input(args.first().map(|s| s.as_str()));
            match clstr_ops::rep(&input) {
                Ok(o) => {
                    stdout.write_all(&o).ok();
                }
                Err(e) => {
                    eprintln!("{e}");
                    exit(1);
                }
            }
        }
        "clstr2tree" | "clstr2tree.pl" => {
            // clstr2tree.pl <clstr> <fr>
            let input = read_clstr_input(args.first().map(|s| s.as_str()));
            let fr = args.get(1).map(|s| s.as_str()).unwrap_or("");
            stdout.write_all(&clstr_ops::to_tree(&input, fr)).ok();
        }
        "cd-hit-clstr_2_blm8" | "cd-hit-clstr_2_blm8.pl" => {
            let input = read_clstr_input(args.first().map(|s| s.as_str()));
            stdout.write_all(&clstr_ops::to_blm8(&input)).ok();
        }
        "clstr_reduce" | "clstr-reduce" => {
            // clstr_reduce.pl <clstr> <segs> <rate>
            let input = read_clstr_input(args.first().map(|s| s.as_str()));
            let segs = args.get(1).map(|s| s.as_str()).unwrap_or("");
            let rate: i64 = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(1);
            stdout.write_all(&clstr_ops::reduce(&input, segs, rate.max(1))).ok();
        }
        "clstr_rev" | "clstr-rev" => {
            // clstr_rev.pl <file90> <file80>
            let f90 = read_clstr_input(args.first().map(|s| s.as_str()));
            let f80 = read_clstr_input(args.get(1).map(|s| s.as_str()));
            stdout.write_all(&clstr_ops::rev(&f90, &f80)).ok();
        }
        "clstr_merge" | "clstr-merge" => {
            // clstr_merge.pl <master> <div1> [div2 ...]
            if args.len() < 2 {
                eprintln!("Usage: clstr_merge <master.clstr> <div1.clstr> [div2.clstr ...]");
                exit(1);
            }
            let master = read_clstr_input(Some(&args[0]));
            let divs: Vec<Vec<u8>> = args[1..]
                .iter()
                .map(|p| read_clstr_input(Some(p)))
                .collect();
            let div_refs: Vec<&[u8]> = divs.iter().map(|d| d.as_slice()).collect();
            stdout.write_all(&clstr_ops::merge(&master, &div_refs)).ok();
        }
        "clstr_reps_faa_rev" | "clstr_reps_faa_rev.pl" => {
            // clstr_reps_faa_rev.pl <clstr> <fasta> <cutoff>
            if args.len() < 3 {
                exit(1);
            }
            let clstr = read_clstr_input(Some(&args[0]));
            let fasta = read_clstr_input(Some(&args[1]));
            let cutoff: usize = args[2].parse().unwrap_or(0);
            stdout
                .write_all(&clstr_ops::reps_faa_rev(&clstr, &fasta, cutoff))
                .ok();
        }
        "clstr_select_rep" | "clstr-select-rep" => {
            // clstr_select_rep.pl <min> <max> [file]
            let min: usize = args.first().and_then(|s| s.parse().ok()).unwrap_or(0);
            let max: usize = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(0);
            let input = read_clstr_input(args.get(2).map(|s| s.as_str()));
            match clstr_ops::select_rep(&input, min, max) {
                Ok(o) => {
                    stdout.write_all(&o).ok();
                }
                Err(e) => {
                    eprintln!("{e}");
                    exit(1);
                }
            }
        }
        "clstr_sort_prot_by" | "clstr-sort-prot-by" => {
            // clstr_sort_prot_by.pl [len|id] < file   (key is the first positional)
            let key = args.first().map(|s| s.as_str()).unwrap_or("len");
            let file = args.get(1).map(|s| s.as_str());
            let input = read_clstr_input(file);
            stdout.write_all(&clstr_ops::sort_prot_by(&input, key)).ok();
        }
        "clstr_merge_noorder" | "clstr-merge-noorder" => {
            // clstr_merge_noorder.pl <master> <div1> [div2 ...]
            if args.len() < 2 {
                eprintln!("Usage: clstr_merge_noorder <master.clstr> <div1.clstr> [div2.clstr ...]");
                exit(1);
            }
            let master = read_clstr_input(Some(&args[0]));
            let divs: Vec<Vec<u8>> = args[1..].iter().map(|p| read_clstr_input(Some(p))).collect();
            let div_refs: Vec<&[u8]> = divs.iter().map(|d| d.as_slice()).collect();
            stdout
                .write_all(&clstr_ops::merge_noorder(&master, &div_refs))
                .ok();
        }
        "clstr_quality_eval_by_link" | "clstr-quality-eval-by-link" => {
            let input = read_clstr_input(args.first().map(|s| s.as_str()));
            match clstr_ops::quality_eval_by_link(&input) {
                Ok(o) => {
                    stdout.write_all(&o).ok();
                }
                Err(e) => {
                    eprintln!("{e}");
                    exit(1);
                }
            }
        }
        "plot_len1" | "plot_len1.pl" => {
            // plot_len1.pl <clstr> <segs> <len_segs>
            if args.len() < 3 {
                eprintln!("Usage: plot_len1 <clstr> <size_segs> <length_segs>");
                eprintln!("  e.g. plot_len1 in.clstr 1,2-5,6-up 1-100,101-200,201-up");
                exit(1);
            }
            let input = read_clstr_input(Some(&args[0]));
            stdout
                .write_all(&clstr_ops::plot_len1(&input, &args[1], &args[2]))
                .ok();
        }
        "clstr_sql_tbl_sort" | "clstr-sql-tbl-sort" => {
            // clstr_sql_tbl_sort.pl <table_file> <level>
            if args.is_empty() {
                print!("Usage:\n\tclstr_sql_tbl_sort.pl table_file level\n");
                exit(1);
            }
            let input = read_clstr_input(Some(&args[0]));
            let level: i64 = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(0);
            eprintln!("done reading {}", args[0]);
            match clstr_ops::sql_tbl_sort(&input, level) {
                Ok(o) => {
                    stdout.write_all(&o).ok();
                }
                Err(e) => {
                    print!("{e}\n");
                    exit(1);
                }
            }
        }
        "cd-hit-dup-PE-out" | "cd-hit-dup-PE-out.pl" => {
            run_dup_pe_out(args);
        }
        "clstr_sql_tbl" | "clstr-sql-tbl" => {
            // clstr_sql_tbl.pl <clstr_file> <tbl_file>: create the table if it
            // does not exist, otherwise append two columns for this level.
            if args.len() < 2 {
                print!("Usage:\n\tclstr_sql_tbl.pl clstr_file tbl_file\n");
                exit(1);
            }
            let clstr = read_clstr_input(Some(&args[0]));
            let tbl_path = &args[1];
            let result = if std::path::Path::new(tbl_path).exists() {
                let existing = read_clstr_input(Some(tbl_path));
                eprintln!("done reading {}", args[0]);
                clstr_ops::sql_tbl(&clstr, Some(&existing))
            } else {
                clstr_ops::sql_tbl(&clstr, None)
            };
            match result {
                Ok(o) => {
                    if let Err(e) = std::fs::write(tbl_path, &o) {
                        eprintln!("Failed to write '{tbl_path}': {e}");
                        exit(1);
                    }
                }
                Err(e) => {
                    eprintln!("{e}");
                    exit(1);
                }
            }
        }
        "clstr_quality_eval" | "clstr-quality-eval" => {
            let input = read_clstr_input(args.first().map(|s| s.as_str()));
            match clstr_ops::quality_eval(&input) {
                Ok(o) => {
                    stdout.write_all(&o).ok();
                }
                Err(e) => {
                    eprintln!("{e}");
                    exit(1);
                }
            }
        }
        "clstr2xml" | "clstr2xml.pl" => {
            // clstr2xml.pl [-len|-size] input1.clstr [input2.clstr ...]
            let mut option = "-len";
            let mut rest: &[String] = args;
            if let Some(first) = args.first() {
                if first.starts_with('-') {
                    option = first;
                    rest = &args[1..];
                }
            }
            if rest.is_empty() {
                print!("Usage:\n\tclstr2xml.pl [-len|-size] input1.clstr [input2.clstr input3.clstr ...]\n");
                exit(0);
            }
            let files: Vec<Vec<u8>> = rest.iter().map(|p| read_clstr_input(Some(p))).collect();
            let refs: Vec<&[u8]> = files.iter().map(|f| f.as_slice()).collect();
            stdout.write_all(&clstr_ops::to_xml(option, &refs)).ok();
        }
        "make_multi_seq" | "make_multi_seq.pl" => {
            // make_multi_seq.pl <fasta> <clstr> <out_dir> <size_cutoff>
            if args.len() < 3 {
                eprintln!("Usage: make_multi_seq <fasta> <clstr> <out_dir> [size_cutoff]");
                exit(1);
            }
            let fasta = read_clstr_input(Some(&args[0]));
            let clstr = read_clstr_input(Some(&args[1]));
            let out_dir = &args[2];
            // Perl: size_cutoff defaults to 1; a 0/non-numeric value also -> 1.
            let cutoff = args
                .get(3)
                .and_then(|s| s.parse::<usize>().ok())
                .filter(|&n| n != 0)
                .unwrap_or(1);
            if let Err(e) = std::fs::create_dir_all(out_dir) {
                eprintln!("can not create {out_dir}: {e}");
                exit(1);
            }
            match clstr_ops::make_multi_seq(&fasta, &clstr, cutoff) {
                Ok(files) => {
                    for f in files {
                        let name = String::from_utf8_lossy(&f.cid);
                        let path = format!("{out_dir}/{name}");
                        if let Err(e) = std::fs::write(&path, &f.content) {
                            eprintln!("can not open file to write {path}: {e}");
                            exit(1);
                        }
                    }
                }
                Err(e) => {
                    eprintln!("{e}");
                    exit(1);
                }
            }
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

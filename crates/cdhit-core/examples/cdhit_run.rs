//! Minimal runner used by the regression harness: reads a FASTA via `-i`,
//! writes `<out>` (rep FASTA) and `<out>.clstr` via `-o`, mirroring cd-hit's
//! file interface for a subset of protein options.

use std::process::exit;

fn main() {
    let argv: Vec<String> = std::env::args().skip(1).collect();
    // Extract -i and -o for file handling; pass the whole arg list to the core.
    let mut input = String::new();
    let mut input2 = String::new();
    let mut output = String::new();
    let mut i = 0;
    while i + 1 < argv.len() {
        match argv[i].as_str() {
            "-i" => input = argv[i + 1].clone(),
            "-i2" => input2 = argv[i + 1].clone(),
            "-o" => output = argv[i + 1].clone(),
            _ => {}
        }
        i += 2;
    }
    let program = match std::env::var("CDHIT_PROG").as_deref() {
        Ok("est") => cdhit_core::Program::CdHitEst,
        Ok("454") => cdhit_core::Program::CdHit454,
        _ if std::env::var("CDHIT_EST").is_ok() => cdhit_core::Program::CdHitEst,
        _ => cdhit_core::Program::CdHit,
    };
    let data = std::fs::read(&input).expect("read input");
    let result = if !input2.is_empty() {
        let data2 = std::fs::read(&input2).expect("read input2");
        let est = program != cdhit_core::Program::CdHit;
        cdhit_core::cluster_2d(&data, &data2, &argv, est)
    } else {
        cdhit_core::cluster_1d_program(&data, &argv, program)
    };
    match result {
        Ok(out) => {
            std::fs::write(&output, &out.rep_fasta).expect("write rep fasta");
            std::fs::write(format!("{output}.clstr"), out.clstr).expect("write clstr");
        }
        Err(e) => {
            eprintln!("error: {e}");
            exit(1);
        }
    }
}

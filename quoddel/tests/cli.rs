use assert_cmd::Command;
use predicates::prelude::*;
use rand::{distributions::Alphanumeric, Rng};
use std::fs;

type TestResult = Result<(), Box<dyn std::error::Error>>;

const PRG: &str = "quoddel";
const FA1: &str = "tests/inputs/fasta1.fa";
//const FA2: &str = "tests/inputs/fasta2.fa";

// --------------------------------------------------
fn gen_bad_file() -> String {
    loop {
        let filename = rand::thread_rng()
            .sample_iter(&Alphanumeric)
            .take(7)
            .map(char::from)
            .collect();

        if fs::metadata(&filename).is_err() {
            return filename;
        }
    }
}

// --------------------------------------------------
fn run(args: &[&str], expected_file: &str) -> TestResult {
    let expected = fs::read_to_string(expected_file)?;
    Command::cargo_bin(PRG)?
        .args(args)
        .assert()
        .success()
        .stdout(expected);
    Ok(())
}

#[test]
fn skips_bad_file() -> TestResult {
    let bad = gen_bad_file();
    let expected = format!(".* [(]os error 2[)]");
    Command::cargo_bin(PRG)?
        .arg(bad)
        .assert()
        .success()
        .stderr(predicate::str::is_match(expected)?);
    Ok(())
}

#[test]
fn fasta1_gc() -> TestResult {
    //todo
    Command::cargo_bin(PRG)?
        .arg(FA1)
        .arg("-m")
        .arg("10")
        .assert()
        .success()
        .stdout(
            predicate::str::contains("min_contig_length: 10")
                .and(predicate::str::contains("num_contigs_gr0: 4"))
                .and(predicate::str::contains("num_contigs_ge_cutoff: 2"))
                .and(predicate::str::contains("total_length_ge_cutoff: 84"))
                .and(predicate::str::contains("total_length_gr0: 90")) //todo greater equal..
                .and(predicate::str::contains("gc_percent_ge_cutoff: 0.2857143")),
        );
    Ok(())
}
#[test]
fn fasta1_gc_m0() -> TestResult {
    //todo
    Command::cargo_bin(PRG)?
        .arg(FA1)
        .arg("-m")
        .arg("0")
        .assert()
        .success()
        .stdout(
            predicate::str::contains("min_contig_length: 0")
                .and(predicate::str::contains("num_contigs_gr0: 4"))
                .and(predicate::str::contains("num_contigs_ge_cutoff: 4"))
                .and(predicate::str::contains("total_length_ge_cutoff: 90"))
                .and(predicate::str::contains("total_length_gr0: 90")) //todo greater equal..
                .and(predicate::str::contains("gc_percent_ge_cutoff: 0.26666668")),
        );
    Ok(())
}

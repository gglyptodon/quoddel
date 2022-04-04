use assert_cmd::Command;
use predicates::prelude::*;
use rand::{distributions::Alphanumeric, Rng};
use std::fs;

type TestResult = Result<(), Box<dyn std::error::Error>>;

const PRG: &str = "quoddel";
const FA1: &str = "tests/inputs/fasta1.fa";
const FA2: &str = "tests/inputs/fasta2.fa";

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

#[test]
fn dies_bad_file() -> TestResult {
    let bad = gen_bad_file();
    let expected = format!(".* [(]os error 2[)]");
    Command::cargo_bin(PRG)?
        .arg(bad)
        .assert()
        .failure()
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
        .arg("--debug")
        .assert()
        .success()
        .stdout(
            predicate::str::contains("min_contig_length: 10")
                .and(predicate::str::contains("num_contigs_ge0: 4"))
                .and(predicate::str::contains("num_contigs_ge_cutoff: 2"))
                .and(predicate::str::contains("total_length_ge_cutoff: 84"))
                .and(predicate::str::contains("total_length_ge0: 90")) //todo greater equal..
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
        .arg("--debug")
        .assert()
        .success()
        .stdout(
            predicate::str::contains("min_contig_length: 0")
                .and(predicate::str::contains("num_contigs_ge0: 4"))
                .and(predicate::str::contains("num_contigs_ge_cutoff: 4"))
                .and(predicate::str::contains("total_length_ge_cutoff: 90"))
                .and(predicate::str::contains("total_length_ge0: 90")) //todo greater equal..
                .and(predicate::str::contains("gc_percent_ge_cutoff: 0.26666668")),
        );
    Ok(())
}

#[test]
fn fasta1_gc_m1() -> TestResult {
    //todo
    Command::cargo_bin(PRG)?
        .arg(FA1)
        .arg("-m")
        .arg("1")
        .arg("--debug")
        .assert()
        .success()
        .stdout(
            predicate::str::contains("min_contig_length: 1")
                .and(predicate::str::contains("num_contigs_ge0: 4"))
                .and(predicate::str::contains("num_contigs_ge_cutoff: 3"))
                .and(predicate::str::contains("total_length_ge_cutoff: 90"))
                .and(predicate::str::contains("total_length_ge0: 90")) //todo greater equal..
                .and(predicate::str::contains("gc_percent_ge_cutoff: 0.26666668")),
        );
    Ok(())
}

#[test]
fn fasta2_gc_m1_n90_l90_n50_l50() -> TestResult {
    //todo
    Command::cargo_bin(PRG)?
        .arg(FA2)
        .arg("-m")
        .arg("1")
        .arg("--debug")
        .assert()
        .success()
        .stdout(
            predicate::str::contains("min_contig_length: 1")
                .and(predicate::str::contains("num_contigs_ge0: 7"))
                .and(predicate::str::contains("num_contigs_ge_cutoff: 7"))
                .and(predicate::str::contains("total_length_ge_cutoff: 525"))
                .and(predicate::str::contains("total_length_ge0: 525")) //todo greater equal..
                .and(predicate::str::contains("gc_percent_ge_cutoff: 0.422"))
                .and(predicate::str::contains("n50_ge_cutoff: 70"))
                .and(predicate::str::contains("l50: 2"))
                .and(predicate::str::contains("n90_ge_cutoff: 35"))
                .and(predicate::str::contains("l90: 6")),
        );
    Ok(())
}

#[test]
pub fn fasta2_gc_m1_n90_l90_n50_l50_stdin() -> TestResult {
    let input = fs::read_to_string(FA1)?;

    Command::cargo_bin(PRG)?
        .write_stdin(input)
        .arg(FA2)
        .arg("-m")
        .arg("1")
        .arg("--debug")
        .assert()
        .success()
        .stdout(
            predicate::str::contains("min_contig_length: 1")
                .and(predicate::str::contains("num_contigs_ge0: 7"))
                .and(predicate::str::contains("num_contigs_ge_cutoff: 7"))
                .and(predicate::str::contains("total_length_ge_cutoff: 525"))
                .and(predicate::str::contains("total_length_ge0: 525"))
                .and(predicate::str::contains("gc_percent_ge_cutoff: 0.422"))
                .and(predicate::str::contains("n50_ge_cutoff: 70"))
                .and(predicate::str::contains("l50: 2"))
                .and(predicate::str::contains("n90_ge_cutoff: 35"))
                .and(predicate::str::contains("l90: 6")),
        );
    Ok(())
}
#[test]
fn fasta1_gc_m0_stdin() -> TestResult {
    let input = fs::read_to_string(FA1)?;
    Command::cargo_bin(PRG)?
        .write_stdin(input)
        .arg("-m")
        .arg("0")
        .arg("--debug")
        .assert()
        .success()
        .stdout(
            predicate::str::contains("min_contig_length: 0")
                .and(predicate::str::contains("num_contigs_ge0: 4"))
                .and(predicate::str::contains("num_contigs_ge_cutoff: 4"))
                .and(predicate::str::contains("total_length_ge_cutoff: 90"))
                .and(predicate::str::contains("total_length_ge0: 90"))
                .and(predicate::str::contains("gc_percent_ge_cutoff: 0.26666668")),
        );
    Ok(())
}

#[test]
fn fasta1_gc_m2_stdin() -> TestResult {
    let input = fs::read_to_string(FA1)?;
    Command::cargo_bin(PRG)?
        .write_stdin(input)
        .arg("-m")
        .arg("1")
        .arg("--debug")
        .assert()
        .success()
        .stdout(
            predicate::str::contains("min_contig_length: 1")
                .and(predicate::str::contains("num_contigs_ge0: 4"))
                .and(predicate::str::contains("num_contigs_ge_cutoff: 3"))
                .and(predicate::str::contains("total_length_ge_cutoff: 90"))
                .and(predicate::str::contains("total_length_ge0: 90"))
                .and(predicate::str::contains("gc_percent_ge_cutoff: 0.26666668")),
        );
    Ok(())
}
//todo test n per 100k bp

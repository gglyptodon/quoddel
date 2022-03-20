pub mod calc;
pub mod output;

use std::error::Error;
use std::io;
use std::io::{BufRead, Read};

use seq_io::fasta::{Reader, Record};

use clap::{Arg, Command};

use crate::calc::*;

type QuoddelResult<T> = Result<T, Box<dyn Error>>;

#[derive(Debug)]
pub struct Config {
    files: Vec<String>, //...
    min_contig_length: usize,
}

#[derive(Debug)]
pub struct FastaInfo {
    num_contigs_gr0: usize,
    num_contigs_gr1000: usize,
    num_contigs_gr5000: usize,
    num_contigs_gr10000: usize,
    num_contigs_gr25000: usize,
    num_contigs_gr50000: usize,

    total_length_gr0: usize,
    total_length_gr1000: usize,
    total_length_gr5000: usize,
    total_length_gr10000: usize,
    total_length_gr25000: usize,
    total_length_gr50000: usize,

    num_contigs_ge_cutoff: usize,

    largest_contig_ge_cutoff: usize,
    gc_percent_ge_cutoff: f32,
    n50_ge_cutoff: usize,
    n90_ge_cutoff: usize,
    //auN: f64, ?
    l50: usize,
    l90: usize,
    //num_n_per_100_kbp: f64,
    total_length_ge_cutoff: usize,
}

pub fn get_args() -> QuoddelResult<Config> {
    let matches = Command::new("quoddel")
        .about("Shows some stats for fasta files. Default minimum contig length is 500bp.")
        .arg(
            Arg::new("files")
                .allow_invalid_utf8(true)
                .default_value("-") //todo
                .min_values(1)
                .help("fasta files (contigs)"),
        ).arg(
        Arg::new("min_contig_length").short('m')
            .long("min-contig")
            .help("minimum contig length to be considered for some stats (to be compatible w quast output)")
            .default_value("500")
    )
        .get_matches();
    let files = matches.values_of_lossy("files").unwrap();
    let min_contig_length = matches.value_of("min_contig_length").unwrap().parse()?;
    Ok(Config {
        files,
        min_contig_length,
    })
}

pub fn run(config: Config) -> QuoddelResult<()> {
    println!("{:#?}", config);
    for file in config.files {
        println!(
            "{:#?}",
            if file == "-" {
                //stdin
                read_fasta_sequences_stdin(config.min_contig_length)?
            } else {
                read_fasta_sequences_file(&file, config.min_contig_length)?
            }
        );
    }
    Ok(())
}
pub fn read_fasta_sequences_file(file: &str, min_contig_length: usize) -> QuoddelResult<FastaInfo> {
    let mut reader = Reader::from_path(file)?;
    let mut info = FastaInfo {
        num_contigs_gr0: 0,
        num_contigs_gr1000: 0,
        num_contigs_gr5000: 0,
        num_contigs_gr10000: 0,
        num_contigs_gr25000: 0,
        num_contigs_gr50000: 0,
        total_length_gr0: 0,
        total_length_gr1000: 0,
        total_length_gr5000: 0,
        total_length_gr10000: 0,
        total_length_gr25000: 0,
        total_length_gr50000: 0,
        num_contigs_ge_cutoff: 0,
        largest_contig_ge_cutoff: 0,
        gc_percent_ge_cutoff: 0.0,
        n50_ge_cutoff: 0,
        n90_ge_cutoff: 0,
        l50: 0,
        l90: 0,
        total_length_ge_cutoff: 0,
    };
    let mut largest_contig_ge_cutoff: usize = 0;
    let mut seq_lengths: Vec<usize> = Vec::new();
    let mut seq_gc: Vec<usize> = Vec::new();
    let mut seq_at: Vec<usize> = Vec::new();

    while let Some(result) = reader.next() {
        let record = result?;
        let seqlen = record.seq_lines().fold(0, |l, seq| l + seq.len());
        info.num_contigs_gr0 += 1;
        info.total_length_gr0 += seqlen;
        if seqlen >= min_contig_length {
            let seqgc = get_gc_num(&record.owned_seq());
            seq_gc.push(seqgc);
            let seqat = get_at_num(&record.owned_seq());
            seq_at.push(seqat);

            seq_lengths.push(seqlen);
            info.num_contigs_ge_cutoff += 1;
            info.total_length_ge_cutoff += seqlen;
            //todo count for all categories
            if seqlen >= 1000 {
                info.num_contigs_gr1000 += 1;
                info.total_length_gr1000 += seqlen;
            }
            if seqlen >= 5000 {
                info.num_contigs_gr5000 += 1;
                info.total_length_gr5000 += seqlen;
            }
            if seqlen >= 10000 {
                info.num_contigs_gr10000 += 1;
                info.total_length_gr10000 = seqlen;
            }
            if seqlen >= 25000 {
                info.num_contigs_gr25000 += 1;
                info.total_length_gr25000 += seqlen
            }
            if seqlen >= 50000 {
                info.num_contigs_gr50000 += 1;
                info.total_length_gr50000 += seqlen;
            }
            // largest
            if seqlen > largest_contig_ge_cutoff {
                largest_contig_ge_cutoff = seqlen;
            }
        }
    }

    let nl_stats = calc_stats(&seq_lengths);
    info.n50_ge_cutoff = nl_stats.n50;
    info.l50 = nl_stats.l50;
    info.n90_ge_cutoff = nl_stats.n90;
    info.l90 = nl_stats.l90;
    info.largest_contig_ge_cutoff = largest_contig_ge_cutoff;

    let gcsum: usize = seq_gc.iter().sum(); //todo clean up
    info.gc_percent_ge_cutoff = gcsum as f32 / (seq_at.iter().sum::<usize>() + gcsum) as f32;

    Ok(info)
}

pub fn read_fasta_sequences_stdin(min_contig_len: usize) -> QuoddelResult<FastaInfo> {
    let mut info = FastaInfo {
        num_contigs_gr0: 0,
        num_contigs_gr1000: 0,
        num_contigs_gr5000: 0,
        num_contigs_gr10000: 0,
        num_contigs_gr25000: 0,
        num_contigs_gr50000: 0,
        total_length_gr0: 0,
        total_length_gr1000: 0,
        total_length_gr5000: 0,
        total_length_gr10000: 0,
        total_length_gr25000: 0,
        total_length_gr50000: 0,
        num_contigs_ge_cutoff: 0,
        largest_contig_ge_cutoff: 0,
        gc_percent_ge_cutoff: 0.0,
        n50_ge_cutoff: 0,
        n90_ge_cutoff: 0,
        l50: 0,
        l90: 0,
        total_length_ge_cutoff: 0,
    };
    let mut reader = Reader::new(io::stdin());
    let mut largest_contig_ge_cutoff: usize = 0;
    let mut seq_lengths: Vec<usize> = Vec::new();
    let mut seq_gc: Vec<usize> = Vec::new();
    let mut seq_at: Vec<usize> = Vec::new();

    while let Some(result) = reader.next() {
        let record = result?;
        let seqlen = record.seq_lines().fold(0, |l, seq| l + seq.len());
        info.num_contigs_gr0 += 1;
        info.total_length_gr0 += seqlen;
        if seqlen >= min_contig_len {
            let seqgc = get_gc_num(&record.owned_seq());
            seq_gc.push(seqgc);
            let seqat = get_at_num(&record.owned_seq());
            seq_at.push(seqat);

            seq_lengths.push(seqlen);
            info.num_contigs_ge_cutoff += 1;
            info.total_length_ge_cutoff += seqlen;
            //todo count for all categories
            if seqlen >= 1000 {
                info.num_contigs_gr1000 += 1;
                info.total_length_gr1000 += seqlen;
            }
            if seqlen >= 5000 {
                info.num_contigs_gr5000 += 1;
                info.total_length_gr5000 += seqlen;
            }
            if seqlen >= 10000 {
                info.num_contigs_gr10000 += 1;
                info.total_length_gr10000 += seqlen;
            }
            if seqlen >= 25000 {
                info.num_contigs_gr25000 += 1;
                info.total_length_gr25000 += seqlen
            }
            if seqlen >= 50000 {
                info.num_contigs_gr50000 += 1;
                info.total_length_gr50000 += seqlen;
            }
            // largest
            if seqlen > largest_contig_ge_cutoff {
                largest_contig_ge_cutoff = seqlen;
            }
        }
    }

    let nl_stats = calc_stats(&seq_lengths);
    info.n50_ge_cutoff = nl_stats.n50;
    info.l50 = nl_stats.l50;
    info.n90_ge_cutoff = nl_stats.n90;
    info.l90 = nl_stats.l90;
    info.largest_contig_ge_cutoff = largest_contig_ge_cutoff;

    let gcsum: usize = seq_gc.iter().sum(); //todo clean up
    info.gc_percent_ge_cutoff = gcsum as f32 / (seq_at.iter().sum::<usize>() + gcsum) as f32;

    Ok(info)
}

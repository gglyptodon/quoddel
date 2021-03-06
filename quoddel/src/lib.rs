pub mod calc;
pub mod output;

use clap::{Arg, Command};
use seq_io::fasta::Reader;
use std::error::Error;

use crate::calc::*;
use crate::output::FastaInfo;

type QuoddelResult<T> = Result<T, Box<dyn Error>>;

#[derive(Debug)]
pub struct Config {
    files: Vec<String>, //...
    min_contig_length: usize,
    debug: bool,
}

pub fn get_args() -> QuoddelResult<Config> {
    let matches = Command::new("quoddel")
        .about("Shows some stats for nucleotide fasta files, e.g. genome assemblies.")
        .arg(
            Arg::new("files")
                .allow_invalid_utf8(true)
                .default_value("-") //todo
                .min_values(1)
                .help("fasta files (contigs)"),
        ).arg(
        Arg::new("min_contig_length").short('m')
            .long("min-contig")
            .help("minimum contig length to be considered for some stats (to be compatible with QUAST output)")
            .default_value("500")
    ).arg(Arg::new("debug").long("--debug").takes_value(false).help("print debug output to stdout"))
        .get_matches();
    let files = matches.values_of_lossy("files").unwrap();
    let min_contig_length = matches.value_of("min_contig_length").unwrap().parse()?;
    let debug = matches.is_present("debug");
    Ok(Config {
        files,
        min_contig_length,
        debug,
    })
}

pub fn run(config: Config) -> QuoddelResult<()> {
    if config.debug {
        println!("{:#?}", config);
    }
    for file in config.files {
        let result = if file == "-" {
            //stdin
            let reader = Reader::new(std::io::stdin());
            read_fasta_sequences(String::from("STDIN"), config.min_contig_length, reader)
        } else {
            let reader = Reader::from_path(&file);
            match reader {
                Err(e) => {
                    eprintln!("{}: {}", file, e);
                    std::process::exit(1)
                }
                Ok(r) => read_fasta_sequences(file.to_string(), config.min_contig_length, r),
            }
        };
        if config.debug {
            println!("{:#?}", result);
        } else {
            match result {
                Ok(r) => print!("{}", r),
                Err(e) => {
                    eprintln!("{}", e);
                    std::process::exit(1)
                }
            }
        }
    }

    Ok(())
}

pub fn read_fasta_sequences<T: std::io::Read>(
    name: String,
    min_contig_length: usize,
    mut reader: Reader<T>,
) -> QuoddelResult<FastaInfo> {
    let mut info = FastaInfo {
        name,
        min_contig_length_cutoff_used: min_contig_length,
        ..Default::default()
    };
    let mut largest_contig_ge_cutoff: usize = 0;
    let mut seq_lengths: Vec<usize> = Vec::new();
    let mut atgcn_vec: Vec<NucCount> = Vec::new();

    while let Some(result) = reader.next() {
        let record = result?;
        let seqlen = record.seq_lines().fold(0, |l, seq| l + seq.len());

        info.num_contigs_ge0 += 1;
        info.total_length_ge0 += seqlen;
        if seqlen >= min_contig_length {
            let atgcn = get_atgcn_num(&record.owned_seq());
            atgcn_vec.push(atgcn);
            seq_lengths.push(seqlen);
            info.num_contigs_ge_cutoff += 1;
            info.total_length_ge_cutoff += seqlen;

            if seqlen >= 1000 {
                info.num_contigs_ge1000 += 1;
                info.total_length_ge1000 += seqlen;
            }
            if seqlen >= 5000 {
                info.num_contigs_ge5000 += 1;
                info.total_length_ge5000 += seqlen;
            }
            if seqlen >= 10000 {
                info.num_contigs_ge10000 += 1;
                info.total_length_ge10000 += seqlen;
            }
            if seqlen >= 25000 {
                info.num_contigs_ge25000 += 1;
                info.total_length_ge25000 += seqlen
            }
            if seqlen >= 50000 {
                info.num_contigs_ge50000 += 1;
                info.total_length_ge50000 += seqlen;
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

    let nucsum: NucCount = atgcn_vec.iter().fold(
        NucCount {
            num_a: 0,
            num_t: 0,
            num_c: 0,
            num_g: 0,
            num_n: 0,
        },
        |count, &new| count + new,
    );

    let gcsum = nucsum.num_g + nucsum.num_c;
    let atsum = nucsum.num_a + nucsum.num_t;

    info.gc_percent_ge_cutoff = gcsum as f32 / (gcsum + atsum) as f32;
    info.num_n_per_100_kbp =
        nucsum.num_n as f32 * 100.0 / (info.total_length_ge_cutoff as f32 / 1000.0);
    Ok(info)
}

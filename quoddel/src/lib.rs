pub mod calc;
pub mod output;

use std::error::Error;
use std::io::BufRead;

use seq_io::fasta::{Reader};

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
    //l50: usize,
    //l90: usize,
    //num_n_per_100_kbp: f64,
    total_length_ge_cutoff: usize,
}

pub fn get_args() -> QuoddelResult<Config> {
    let matches = Command::new("quoddel")
        .about("Shows some stats for fasta files")
        .arg(
            Arg::new("files")
                .allow_invalid_utf8(true)
                //.default_value("-") //todo
                .required(true)
                .min_values(1)
                .help("fasta files (contigs)"),
        ).arg(
        Arg::new("min_contig_length").short('m')
            .long("min-contig")
            .help("minimum contig length to be considered for some stats (to be compatible w quast output)")
            .default_value("500")
    )
        .get_matches();
    // add min contig length to be comparable to quast default 500 -m, --min-contig
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
            read_fasta_sequences(&file, config.min_contig_length)?
        );
    }
    Ok(())
}

pub fn read_fasta_sequences(file: &str, min_contig_length: usize) -> QuoddelResult<FastaInfo> {
    //todo refactor
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
        total_length_ge_cutoff: 0,
    };
    // let mut num_seqs:usize = 0;
    let mut reader = Reader::from_path(file)?;
    let mut reader_largest_contig = Reader::from_path(file)?;
    let mut reader_n50 = Reader::from_path(file)?;

    let mut gr0: usize = 0;
    let mut gr1000: usize = 0;
    let mut gr5000: usize = 0;
    let mut gr10000: usize = 0;
    let mut gr25000: usize = 0;
    let mut gr50000: usize = 0;

    let mut tlen0: usize = 0;
    let mut tlen1000: usize = 0;
    let mut tlen5000: usize = 0;
    let mut tlen10000: usize = 0;
    let mut tlen25000: usize = 0;
    let mut tlen50000: usize = 0;

    let mut count_ge_cutoff = 0;
    for x in reader.records().map(|x| x.unwrap().seq.len()) {
        if x >= min_contig_length {
            count_ge_cutoff += 1;
        }
        match x {
            //0 => {}
            0..=999 => {
                gr0 += 1;
                tlen0 += x
            }
            1000..=4999 => {
                gr1000 += 1;
                tlen1000 += x
            }
            5000..=9999 => {
                gr5000 += 1;
                tlen5000 += x
            }
            10000..=24999 => {
                gr10000 += 1;
                tlen10000 += x
            }
            25000..=49999 => {
                gr25000 += 1;
                tlen25000 += x
            }
            _ => {
                gr50000 += 1;
                tlen50000 += x
            }
        }
        //println!("{}",x);
    }
    // num_seqs= reader.records().count();
    //while let Some(_record) = reader.next() {
    //let record = record?;
    //.expect("Error reading record");
    //     num_seqs+=1;
    //println!("{}", record.id().unwrap());
    //
    // }
    info.num_contigs_gr0 = [gr50000, gr25000, gr10000, gr5000, gr1000, gr0]
        .iter()
        .sum();
    info.num_contigs_gr1000 = [gr50000, gr25000, gr10000, gr5000, gr1000].iter().sum();
    info.num_contigs_gr5000 = [gr50000, gr25000, gr10000, gr5000].iter().sum();
    info.num_contigs_gr10000 = [gr50000, gr25000, gr10000].iter().sum();
    info.num_contigs_gr25000 = [gr50000, gr25000].iter().sum();
    info.num_contigs_gr50000 = gr50000;

    info.total_length_ge_cutoff = [tlen50000, tlen25000, tlen10000, tlen5000, tlen1000, tlen0]
        .iter()
        .sum();

    info.total_length_gr0 = info.total_length_ge_cutoff;
    info.total_length_gr1000 = info.total_length_ge_cutoff - tlen0;
    info.total_length_gr5000 = info.total_length_gr1000 - tlen1000;
    info.total_length_gr10000 = info.total_length_gr5000 - tlen5000;
    info.total_length_gr25000 = info.total_length_gr10000 - tlen10000;
    info.total_length_gr50000 = info.total_length_gr25000 - tlen25000;

    info.num_contigs_ge_cutoff = count_ge_cutoff; //todo

    if let Some(largest) = reader_largest_contig
        .records()
        .map(|x| x.unwrap().seq.len())
        .max()
    {
        info.largest_contig_ge_cutoff = largest;
    }
    let lengths: Vec<usize> = reader_n50
        .records()
        .map(|x| x.unwrap().seq.len())
        .filter(|&y| y >= min_contig_length)
        .collect();
    info.total_length_ge_cutoff = lengths.iter().sum(); //todo refactor
    info.n50_ge_cutoff = n50(&lengths);
    info.n90_ge_cutoff = n90(&lengths);

    // let tmp = calc::gc_total(file)?;//TODO!
    //  println!("gc total num: {}, at: {}", tmp.0, tmp.1);
    info.gc_percent_ge_cutoff = calc_gc(file, min_contig_length)?; //tmp.0 as f32 /(tmp.1+tmp.0) as f32; //todo in quast this if for contigs : > cutoff only!
    Ok(info)
}

pub fn get_info(mut buff: impl BufRead) -> QuoddelResult<FastaInfo> {
    let info = FastaInfo {
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
        total_length_ge_cutoff: 0,
    };

    let mut res_buff = String::new();
    loop {
        let num_bytes_read = buff.read_line(&mut res_buff)?;
        if num_bytes_read == 0 {
            break;
        }
        // num_lines += 1;
        //num_bytes += num_bytes_read;
        //num_chars += res_buff.chars().count();
        //num_words += res_buff.split_ascii_whitespace().count();
        res_buff.clear();
    }
    Ok(info)
}
/*
pub fn read_gz(mut buff: impl BufRead) -> QuoddelResult<()> {
    let mut d = GzDecoder::new(buff);
    let mut s = [08; 10];
    let read = d.read(&mut s)?;
    println!("{:?}", s);
    Ok(())
}
*/

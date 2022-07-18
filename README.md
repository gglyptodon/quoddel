# quoddel
[![Rust](https://github.com/gglyptodon/quoddel/actions/workflows/rust.yml/badge.svg)](https://github.com/gglyptodon/quoddel/actions/workflows/rust.yml)
[![Lint Code Base](https://github.com/gglyptodon/quoddel/actions/workflows/linter.yml/badge.svg)](https://github.com/gglyptodon/quoddel/actions/workflows/linter.yml)

Display key quality metrics about genome assemblies.

Input: 
- Fasta file with contigs. (Can also be read from stdin, e.g. for on-the-fly decompression)

Output: 
- Table of metrics (same output format as QUAST's output for this subset of metrics; tab-separated; printed to stdout).

## Build / Install
Pre-built binaries for Linux and MacOS can be found here: https://github.com/gglyptodon/quoddel/releases

To build yourself:
```
# install rust as described here: https://www.rust-lang.org/tools/install
git clone git@github.com:gglyptodon/quoddel.git
cd quoddel
cargo build --release
# then run the binary from the quoddel/target/release directory or copy to a directory in your PATH
```

## Usage

```text
./quoddel -h
quoddel 
Shows some stats for nucleotide fasta files, e.g. genome assemblies.

USAGE:
    quoddel [OPTIONS] [files]...

ARGS:
    <files>...    fasta files (contigs) [default: -]

OPTIONS:
        --debug
            print debug output to stdout

    -h, --help
            Print help information

    -m, --min-contig <min_contig_length>
            minimum contig length to be considered for some stats (to be compatible with QUAST
            output) [default: 500]

```

## Example

Use gunzip to pipe gzipped fasta to quoddel:

```text
gunzip -kc LargeAssembly.fna.gz | quoddel > LargeAssembly.tsv
```

LargeAssembly.tsv:

```text
Assembly        STDIN
num contigs (>= 0 bp)   12668
num contigs (>= 1000 bp)        12668
num contigs (>= 5000 bp)        11893
num contigs (>= 10000 bp)       4709
num contigs (>= 25000 bp)       312
num contigs (>= 50000 bp)       137
total length (>= 0 bp)  40054341269
total length (>= 1000 bp)       40054341269
total length (>= 5000 bp)       40050941494
total length (>= 10000 bp)      39999622610
total length (>= 25000 bp)      39943414367
total length (>= 50000 bp)      39937391001
minimum contig length cutoff    500
num contigs     12668
largest contig  2000000000
total length    40054341269
GC (%)  40.212
N50     1974114035
N90     863342500
L50     11
L90     22
num N's per 100 kbp     2479.399
```

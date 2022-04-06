[![Rust](https://github.com/gglyptodon/quoddel/actions/workflows/rust.yml/badge.svg)](https://github.com/gglyptodon/quoddel/actions/workflows/rust.yml)

# quoddel

Display stats about fasta files.
Output is printed to stdout.


## Usage

```
./quoddel -h
quoddel 
Shows some stats for fasta files. Default minimum contig length is 500bp.

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
            minimum contig length to be considered for some stats (to be compatible w quast output)
            [default: 500]
```

## Example

Use gunzip to pipe gzipped fasta to quoddel:

```
gunzip -kc LargeAssembly.fna.gz | ./quoddel-v0.1.0 > LargeAssembly.tsv
```

LargeAssembly.tsv:

```
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
```

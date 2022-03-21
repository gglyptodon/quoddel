use std::fmt;

#[derive(Default, Debug)]
pub struct FastaInfo {
    pub(crate) name: String,
    pub(crate) min_contig_length_cutoff_used: usize,
    pub(crate) num_contigs_ge0: usize,
    pub(crate) num_contigs_ge1000: usize,
    pub(crate) num_contigs_ge5000: usize,
    pub(crate) num_contigs_ge10000: usize,
    pub(crate) num_contigs_ge25000: usize,
    pub(crate) num_contigs_ge50000: usize,

    pub(crate) total_length_ge0: usize,
    pub(crate) total_length_ge1000: usize,
    pub(crate) total_length_ge5000: usize,
    pub(crate) total_length_ge10000: usize,
    pub(crate) total_length_ge25000: usize,
    pub(crate) total_length_ge50000: usize,

    pub(crate) num_contigs_ge_cutoff: usize,

    pub(crate) largest_contig_ge_cutoff: usize,
    pub(crate) gc_percent_ge_cutoff: f32,
    pub(crate) n50_ge_cutoff: usize,
    pub(crate) n90_ge_cutoff: usize,
    //auN: f32, ?
    pub(crate) l50: usize,
    pub(crate) l90: usize,
    pub(crate) num_n_per_100_kbp: f32,
    pub(crate) total_length_ge_cutoff: usize,
}

impl fmt::Display for FastaInfo {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "Assembly\t{}
num contigs (>= 0 bp)\t{}
num contigs (>= 1000 bp)\t{}
num contigs (>= 5000 bp)\t{}
num contigs (>= 10000 bp)\t{}
num contigs (>= 25000 bp)\t{}
num contigs (>= 50000 bp)\t{}
total length (>= 0 bp)\t{}
total length (>= 1000 bp)\t{}
total length (>= 5000 bp)\t{}
total length (>= 10000 bp)\t{}
total length (>= 25000 bp)\t{}
total length (>= 50000 bp)\t{}
minimum contig length cutoff\t{}
num contigs\t{}
largest contig\t{}
total length\t{}
GC (%)\t{:.3}
N50\t{}
N90\t{}
L50\t{}
L90\t{}
num N's per 100 kbp\t{:.3}
        ",
            self.name,
            self.num_contigs_ge0,
            self.num_contigs_ge1000,
            self.num_contigs_ge5000,
            self.num_contigs_ge10000,
            self.num_contigs_ge25000,
            self.num_contigs_ge50000,
            self.total_length_ge0,
            self.total_length_ge1000,
            self.total_length_ge5000,
            self.total_length_ge10000,
            self.total_length_ge25000,
            self.total_length_ge50000,
            self.min_contig_length_cutoff_used,
            self.num_contigs_ge_cutoff,
            self.largest_contig_ge_cutoff,
            self.total_length_ge_cutoff,
            self.gc_percent_ge_cutoff * 100.0,
            self.n50_ge_cutoff,
            self.n90_ge_cutoff,
            self.l50,
            self.l90,
            self.num_n_per_100_kbp,
        )
    }
}

use crate::QuoddelResult;
use seq_io::fasta::Reader;
use std::ops::Add;

#[derive(Debug)]
pub struct NLStats {
    pub n50: usize,
    pub n90: usize,
    pub l50: usize,
    pub l90: usize,
}
pub fn calc_stats(lengths: &[usize]) -> NLStats {
    let fifty: f32 = (lengths.iter().map(|&x| x as f32).sum::<f32>()) * 0.5;
    let ninety: f32 = (lengths.iter().map(|&x| x as f32).sum::<f32>()) * 0.9;
    let mut tmp = lengths.to_owned();
    tmp.sort_by(|a, b| b.cmp(a));
    let mut n90sum: usize = 0;
    let mut n50sum: usize = 0;
    let mut current_50: usize = 0;
    let mut current_90: usize = 0;
    let mut count_50: usize = 0;
    let mut count_90: usize = 0;
    for decreasing in tmp {
        if n90sum >= ninety as usize {
            break;
        }
        if n50sum < fifty as usize {
            //todo test on edge case
            n50sum += decreasing;
            current_50 = decreasing;
            count_50 += 1;
        }
        n90sum += decreasing;
        current_90 = decreasing;
        count_90 += 1;
    }
    NLStats {
        n50: current_50,
        n90: current_90,
        l50: count_50,
        l90: count_90,
    }
}

pub fn n90(lengths: &Vec<usize>) -> usize {
    //todo de-uglify
    let ninety: f32 = (lengths.iter().map(|&x| x as f32).sum::<f32>()) * 0.9;
    let mut tmp = lengths.to_owned(); //clone();
    tmp.sort_by(|a, b| b.cmp(a));
    let mut n90sum: usize = 0;
    let mut current: usize = 0;
    for decreasing in tmp {
        //println!("decreasing: {}", decreasing);
        if n90sum >= ninety as usize {
            break;
        }
        n90sum += decreasing;
        current = decreasing;
    }
    current
}

pub fn n50(lengths: &Vec<usize>) -> usize {
    //todo oof ugly ey
    //longest to shortest
    let half: f32 = (lengths.iter().map(|&x| x as f32).sum::<f32>()) * 0.5;
    let mut tmp = lengths.to_owned(); //.clone();
    tmp.sort_by(|a, b| b.cmp(a));
    let mut n50sum: usize = 0;
    let mut current: usize = 0;
    for decreasing in tmp {
        if n50sum >= half as usize {
            //todo...
            break;
        }
        n50sum += decreasing;
        current = decreasing;
    }
    current
}

pub fn get_gc_num(seq: &Vec<u8>) -> usize {
    seq.iter()
        .filter(|&&c| c == b'C' || c == b'G' || c == b'g' || c == b'c')
        .count()
}
pub fn get_at_num(seq: &Vec<u8>) -> usize {
    seq.iter()
        .filter(|&&c| c == b'A' || c == b'T' || c == b'a' || c == b't')
        .count()
}

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct NucCount {
    pub(crate) num_a: usize,
    pub(crate) num_t: usize,
    pub(crate) num_c: usize,
    pub(crate) num_g: usize,
    pub(crate) num_n: usize,
}
impl Add for NucCount {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self {
            num_a: self.num_a + other.num_a,
            num_t: self.num_t + other.num_t,
            num_g: self.num_g + other.num_g,
            num_c: self.num_c + other.num_c,
            num_n: self.num_n + other.num_n,
        }
    }
}

pub fn count_nuc(nuc: &u8) -> NucCount {
    let mut nuc_count = NucCount {
        num_a: 0,
        num_t: 0,
        num_c: 0,
        num_g: 0,
        num_n: 0,
    };
    match nuc {
        b'A' | b'a' => nuc_count.num_a += 1,
        b'T' | b't' => nuc_count.num_t += 1,
        b'G' | b'g' => nuc_count.num_g += 1,
        b'C' | b'c' => nuc_count.num_c += 1,
        b'N' | b'n' => nuc_count.num_n += 1,
        _ => {}
    }
    nuc_count
}
pub fn get_atgcn_num(seq: &Vec<u8>) -> NucCount {
    seq.iter()
        //.filter(|&&c| c in &[b'a',b'A',b't',b'T', b'g',b'G',b'c',b'C', b'N', b'n'])
        .map(count_nuc)
        .fold(
            NucCount {
                num_a: 0,
                num_t: 0,
                num_c: 0,
                num_g: 0,
                num_n: 0,
            },
            |count, new| count + new,
        )
}

pub fn calc_gc(file: &str, length_cutoff: usize) -> QuoddelResult<f32> {
    let mut reader_gc = Reader::from_path(file)?;
    //let mut reader_at = Reader::from_path(file)?;

    let filtered: Vec<Vec<u8>> = reader_gc
        .records()
        .map(|x| x.unwrap().seq)
        .filter(|x| x.len() >= length_cutoff)
        .collect();
    let gcnum = filtered.iter().map(get_gc_num).sum::<usize>();
    let atnum = filtered.iter().map(get_at_num).sum::<usize>();

    Ok(gcnum as f32 / (atnum + gcnum) as f32)
}

#[cfg(test)]
mod tests {
    use crate::{
        calc::{n50, n90},
        get_at_num, get_atgcn_num, get_gc_num, FastaInfo, NucCount,
    };

    #[test]
    fn test_n50_odd() {
        let len_vec: Vec<usize> = vec![15, 11, 12, 16, 14, 10, 1];
        let expected: usize = 14;
        let result = n50(&len_vec);
        assert_eq!(result, expected);
    }
    #[test]
    fn test_n50_even() {
        let len_vec: Vec<usize> = vec![1, 3, 15, 16, 8, 4];
        // 47.0*0.5 -> 23.5
        // 16 < 23.5
        // 16+15 >= 23.5
        let expected: usize = 15;
        let result = n50(&len_vec);
        assert_eq!(result, expected);
    }

    #[test]
    fn test_n50_odd_first() {
        let len_vec: Vec<usize> = vec![150, 1, 2];
        // 76.5
        // 150 > 76.5
        let expected: usize = 150;
        let result = n50(&len_vec);
        assert_eq!(result, expected);
    }

    #[test]
    fn test_n50_even_first() {
        let len_vec: Vec<usize> = vec![1, 3, 150, 2];
        let expected: usize = 150;
        let result = n50(&len_vec);
        assert_eq!(result, expected);
    }

    #[test]
    fn test_n90_odd() {
        let len_vec: Vec<usize> = vec![15, 11, 12, 16, 14, 10, 1];
        let expected: usize = 10;
        let result = n90(&len_vec);
        assert_eq!(result, expected);
    }
    #[test]
    fn test_n90_even() {
        let len_vec: Vec<usize> = vec![1, 3, 15, 16, 8, 4];
        // 47.0*0.5 -> 23.5
        // 16 < 23.5
        // 16+15 >= 23.5
        let expected: usize = 4;
        let result = n90(&len_vec);
        assert_eq!(result, expected);
    }

    #[test]
    fn test_n90_odd_first() {
        let len_vec: Vec<usize> = vec![150, 1, 2];
        // 76.5
        // 150 > 76.5
        let expected: usize = 150;
        let result = n90(&len_vec);
        assert_eq!(result, expected);
    }

    #[test]
    fn test_n90_even_first() {
        let len_vec: Vec<usize> = vec![1, 3, 140, 2];
        let expected: usize = 140;
        let result = n90(&len_vec);
        assert_eq!(result, expected);
    }
    #[test]
    fn test_get_gc_num() {
        let seq: Vec<u8> = vec!['A', 'T', 'G', 'C'].iter().map(|&x| x as u8).collect();
        let expected: usize = 2;
        let result = get_gc_num(&seq);
        assert_eq!(result, expected);
    }
    #[test]
    fn test_get_at_num() {
        let seq: Vec<u8> = vec!['A', 'T', 'G', 'C'].iter().map(|&x| x as u8).collect();
        let expected: usize = 2;
        let result = get_gc_num(&seq);
        assert_eq!(result, expected);
    }
    #[test]
    fn test_get_gc_num_n() {
        let seq: Vec<u8> = vec!['A', 'N', 'T', 'G', 'C', 'a', 'N']
            .iter()
            .map(|&x| x as u8)
            .collect();
        let expected: usize = 2;
        let result = get_gc_num(&seq);
        assert_eq!(result, expected);
    }
    #[test]
    fn test_get_gc_num_lc() {
        let seq: Vec<u8> = vec!['A', 'N', 'T', 'g', 'C', 'N', 'G']
            .iter()
            .map(|&x| x as u8)
            .collect();
        let expected: usize = 3;
        let result = get_gc_num(&seq);
        assert_eq!(result, expected);
    }
    #[test]
    fn test_get_at_num_lc() {
        let seq: Vec<u8> = vec!['a', 'N', 't', 'g', 'C', 'N', 'G']
            .iter()
            .map(|&x| x as u8)
            .collect();
        let expected: usize = 2;
        let result = get_at_num(&seq);
        assert_eq!(result, expected);
    }
    #[test]
    fn test_gc_no_minlength_lc() {
        let seqs = ["atggggggggggtt", "atg", "atggggggggggttatg", "tt"]
            .iter()
            .map(|x| x.chars().map(|y| y as u8).collect::<Vec<u8>>())
            .collect::<Vec<Vec<u8>>>();
        let expected: [usize; 4] = [10, 1, 11, 0];
        for (seq, expected) in seqs.iter().zip(expected) {
            assert_eq!(get_gc_num(&seq), expected);
        }
    }
    #[test]
    fn test_gc_no_minlength_mc() {
        let seqs = ["atggGgggggggtt", "atg", "atGGGggggcCgggttatg", "tt"]
            .iter()
            .map(|x| x.chars().map(|y| y as u8).collect::<Vec<u8>>())
            .collect::<Vec<Vec<u8>>>();
        let expected: [usize; 4] = [10, 1, 13, 0];
        for (seq, expected) in seqs.iter().zip(expected) {
            assert_eq!(get_gc_num(&seq), expected);
        }
    }
    #[test]
    fn test_gc_minlength_3_mc() {
        let seqs = ["atggGgggggggtt", "atg", "atGGGggggcCgggttatg", "tt"]
            .iter()
            .map(|x| x.chars().map(|y| y as u8).collect::<Vec<u8>>())
            .collect::<Vec<Vec<u8>>>();
        let expected: [usize; 4] = [10, 1, 13, 0];
        for (seq, expected) in seqs.iter().zip(expected) {
            assert_eq!(get_gc_num(&seq), expected);
        }
    }

    #[test]
    fn test_nuc_count() {
        let seqs = ["atggGggngggggtt", "atg", "antGGGggggcCgggttatgN", "tt"]
            .iter()
            .map(|x| x.chars().map(|y| y as u8).collect::<Vec<u8>>())
            .collect::<Vec<Vec<u8>>>();
        let expected = [
            NucCount {
                num_a: 1,
                num_t: 3,
                num_c: 0,
                num_g: 10,
                num_n: 1,
            },
            NucCount {
                num_a: 1,
                num_t: 1,
                num_c: 0,
                num_g: 1,
                num_n: 0,
            },
            NucCount {
                num_a: 2,
                num_t: 4,
                num_c: 2,
                num_g: 11,
                num_n: 2,
            },
            NucCount {
                num_a: 0,
                num_t: 2,
                num_c: 0,
                num_g: 0,
                num_n: 0,
            },
        ];
        let mut res: Vec<NucCount> = Vec::new();
        let expected_total: NucCount = NucCount {
            num_a: 4,
            num_t: 10,
            num_g: 22,
            num_c: 2,
            num_n: 3,
        };
        for (s, ex) in seqs.iter().zip(expected) {
            let count = get_atgcn_num(&s);
            res.push(count);
            assert_eq!(count, ex)
        }
        assert_eq!(
            res.iter().fold(
                NucCount {
                    num_a: 0,
                    num_t: 0,
                    num_c: 0,
                    num_g: 0,
                    num_n: 0,
                },
                |count, &new| count + new,
            ),
            expected_total
        )
    }
}

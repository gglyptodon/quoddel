use crate::QuoddelResult;
use seq_io::fasta::{Reader};

pub fn n90(lengths: &Vec<usize>) -> usize {
    //todo de-uglify
    let ninety: f32 = (lengths.iter().map(|&x| x as f32).sum::<f32>()) * 0.9;
    let mut tmp = lengths.to_owned();//clone();
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
    let mut tmp = lengths.to_owned();//.clone();
    tmp.sort_by(|a, b| b.cmp(a));
    let mut n50sum: usize = 0;
    let mut current: usize = 0;
    for decreasing in tmp {
        //println!("decreasing: {}", decreasing);
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
        .filter(|&&c| c == b'A'|| c == b'T' || c == b'a'|| c == b't')
        .count()
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
    //.map(|x| get_gc_num(x.unwrap().seq))
    //.sum();
    //println!("gcnum: {}", gcnum);
    let atnum = filtered.iter().map(get_at_num).sum::<usize>();
    //.map(|x| get_gc_num(x.unwrap().seq))
    //.sum();
    //println!("atnum: {}", atnum);

    // let atnum: usize = reader_at
    //     .records()
    //.filter(|x| x.unwrap().seq.len() >= length_cutoff)
    //     .map(|x| get_at_num(x.unwrap().seq))
    //     .sum();
    // println!("atnum: {}", atnum);

    Ok(gcnum as f32 / (atnum + gcnum) as f32)
}

/*pub fn gc_total(file: &str) ->QuoddelResult<(u64,u64)>{ //todo min length cutoff needs to be regarded
    let mut reader = Reader::from_path(file)?;
    let mut reader_at = Reader::from_path(file)?;
    let gc_num =reader.records().map(|x| x.unwrap()
        .seq.iter()
        .filter(|&&c| c =='C' as u8 ||c=='G' as u8 ||c=='g' as u8 ||c=='c' as u8 ).map(|&x|x as u64).sum::<u64>()//todo
    ).sum::<u64>();
    let at_num = reader_at.records().map(|x| x.unwrap()
        .seq.iter()
        .filter(|&&c| c =='A' as u8 ||c=='T' as u8 ||c=='a' as u8 ||c=='t' as u8 ).map(|&x|x as u64).sum::<u64>()//todo
    ).sum::<u64>();
    // {
        //println!("{}",x);
    //}
    Ok((gc_num,at_num))
}
*/

#[cfg(test)]
mod tests {
    use crate::{
        calc::{n50, n90},
        get_at_num, get_gc_num, FastaInfo,
    };
    use std::io::Cursor;
    //use std::io::Cursor;

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
}

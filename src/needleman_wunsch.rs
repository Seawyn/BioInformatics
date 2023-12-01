use std::clone::Clone;

#[derive(Clone)]
#[derive(Debug)]
struct Mcell {
    score: i32,
    // Binary represents encoding of trace combinations
    // 000 (origin), 100 (left), 010 (diagonal), 001 (up)
    trace: i8
}

fn setup_matrix(strand_1: &str, strand_2: &str) -> Vec<Vec<Mcell>> {
    let size_1 = strand_1.len();
    let size_2 = strand_2.len();

    let mut matrix = vec![vec![Mcell{score: 0, trace: 0}; size_2 + 1]; size_1 + 1];

    for i in 1..size_1 + 1 {
        matrix[0][i] = Mcell{score: i as i32 * -1, trace: 0 | 0b100};
    }

    for i in 1..size_2 + 1 {
        matrix[i][0] = Mcell{score: i as i32 * -1, trace: 0 | 0b001};
    }

    return matrix;
}

fn match_score(a: char, b: char, s_match: i32, s_mismatch: i32) -> i32 {
    let r = if a == b {s_match} else {s_mismatch};
    return r;
}

fn fill_table(
    strand_1: &str,
    strand_2: &str,
    mut m: Vec<Vec<Mcell>>, 
    s_match: i32,
    s_mismatch: i32,
    s_indel: i32,
) -> Vec<Vec<Mcell>> {
    for i in 1..strand_2.len() + 1 {
        for j in 1..strand_1.len() + 1 {
            let curr_c1 = strand_1.chars().nth(j - 1).unwrap();
            let curr_c2 = strand_2.chars().nth(i - 1).unwrap();

            let score_up = m[i - 1][j].score + s_indel;
            let score_left = m[i][j - 1].score + s_indel;
            let score_diag = m[i - 1][j - 1].score + match_score(curr_c1, curr_c2, s_match, s_mismatch);

            let curr_max = *[score_up, score_left, score_diag].iter().max().unwrap();
            m[i][j] = Mcell{score: curr_max, trace: 0};
            
            if score_up == curr_max { m[i][j].trace |= 0b001 }
            if score_left == curr_max { m[i][j].trace |= 0b100 }
            if score_diag == curr_max { m[i][j].trace |= 0b010 }
        }
    }

    m
}

fn align(
    strand_1: &str, strand_2: &str,
    alignment: Vec<i8>
) -> [String; 2] {
    let mut cursor_1 = (strand_1.len() - 1) as i8;
    let mut cursor_2 = (strand_2.len() - 1) as i8;

    let mut aligned_1 = String::new();
    let mut aligned_2 = String::new();

    for option in alignment.iter() {
        match option {
            // Diagonal, alignment
            0b010 => {
                aligned_1.insert(0, strand_1.chars().nth(cursor_1 as usize).unwrap());
                aligned_2.insert(0, strand_2.chars().nth(cursor_2 as usize).unwrap());
                cursor_1 -= 1;
                cursor_2 -= 1;
            }
            // Left
            0b100 => {
                aligned_1.insert(0, strand_1.chars().nth(cursor_1 as usize).unwrap());
                aligned_2.insert(0, '-');
                cursor_1 -= 1;
            }
            // Up
            0b001 => {
                aligned_1.insert(0, '-');
                aligned_2.insert(0, strand_2.chars().nth(cursor_2 as usize).unwrap());
                cursor_2 -= 1;
            }
            0b000 => break,

            _ => panic!("Unknown value"),
        }
    }

    [aligned_1, aligned_2]

}

pub fn needleman_wunsch(
    strand_1: &str, strand_2: &str,
    s_match: i32, s_mismatch: i32, s_indel: i32,
    debug: bool
) -> i32 {
    let mut m = setup_matrix(strand_1, strand_2);
    m = fill_table(strand_1, strand_2, m, s_match, s_mismatch, s_indel);

    let last_pos = m.last().unwrap().last().unwrap();

    if debug {
        for i in 0..m.len() {
            println!("{:?}", m[i]);
        }
    }
    
    last_pos.score
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn match_score_differ() {
        let a: char = 'A';
        let b: char = 'C';

        let s_match = 1;
        let s_mismatch = -1;

        let res = match_score(a, b, s_match, s_mismatch);

        assert_eq!(res, s_mismatch);
    }

    #[test]
    fn match_score_equal() {
        let a: char = 'G';
        let b: char = 'G';

        let s_match = 1;
        let s_mismatch = -1;

        let res = match_score(a, b, s_match, s_mismatch);

        assert_eq!(res, s_match);
    }

    #[test]
    fn setup_trace_line() {
        let strand_1 = "AAAAAAAA";
        let strand_2 = "AAAAAAAA";

        let match_score = 1;
        let mismatch_score = -1;
        let indel_score = -1;

        let mut m = setup_matrix(strand_1, strand_2);
        m = fill_table(strand_1, strand_2, m, 
            match_score, mismatch_score, indel_score);
        
        assert_eq!(m[0][0].trace, 0);

        for i in 1..strand_1.len() + 1 {
            assert_ne!(m[0][i].trace & 0b100, 0);
        }
    }

    #[test]
    fn setup_trace_col() {
        let strand_1 = "AAAAAAAA";
        let strand_2 = "AAAAAAAA";

        let match_score = 1;
        let mismatch_score = -1;
        let indel_score = -1;

        let mut m = setup_matrix(strand_1, strand_2);
        m = fill_table(strand_1, strand_2, m, 
            match_score, mismatch_score, indel_score);
        
        assert_eq!(m[0][0].trace, 0);

        for i in 1..strand_2.len() + 1 {
            assert_ne!(m[i][0].trace & 0b001, 0);
        }
    }

    #[test]
    fn setup_score_line() {
        let strand_1 = "AAAAAAAA";
        let strand_2 = "AAAAAAAA";

        let match_score = 1;
        let mismatch_score = -1;
        let indel_score = -1;

        let mut m = setup_matrix(strand_1, strand_2);
        m = fill_table(strand_1, strand_2, m, 
            match_score, mismatch_score, indel_score);
        
        for i in 1..strand_1.len() + 1 {
            assert_eq!(m[0][i].score, i as i32 * -1);
        }
    }

    #[test]
    fn setup_score_col() {
        let strand_1 = "AAAAAAAA";
        let strand_2 = "AAAAAAAA";

        let match_score = 1;
        let mismatch_score = -1;
        let indel_score = -1;

        let mut m = setup_matrix(strand_1, strand_2);
        m = fill_table(strand_1, strand_2, m, 
            match_score, mismatch_score, indel_score);
        
        for i in 1..strand_2.len() + 1 {
            assert_eq!(m[i][0].score, i as i32 * -1);
        }
    }

    #[test]
    fn all_traces() {
        let strand_1 = "ACGTACGT";
        let strand_2 = "TGCATGCA";

        let match_score = 1;
        let mismatch_score = -1;
        let indel_score = -1;

        let mut m = setup_matrix(strand_1, strand_2);
        m = fill_table(strand_1, strand_2, m, 
            match_score, mismatch_score, indel_score);
        
        for i in 1..strand_1.len() + 1 {
            for j in 1..strand_2.len() + 1 {
                assert_ne!(m[i][j].trace | 0b111, 0);
            }
        }
    }

    #[test]
    fn alignment_1() {
        let strand_1 = "AAGTCGGAT";
        let strand_2 = "AGGCGTATT";

        let alignment = vec![0b010, 0b001, 0b010, 0b010, 0b010, 0b010, 0b010, 0b010, 0b010, 0b100];

        let res = align(strand_1, strand_2, alignment);

        assert_eq!(res[0], "AAGTCGGA-T");
        assert_eq!(res[1], "-AGGCGTATT");
    }
}
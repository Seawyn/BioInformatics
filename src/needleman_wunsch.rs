use std::cmp;

fn setup_matrix(strand_1: &str, strand_2: &str) -> Vec<Vec<i32>> {
    let size_1 = strand_1.len();
    let size_2 = strand_2.len();

    let mut matrix = vec![vec![0 as i32; size_2 + 1]; size_1 + 1];

    for i in 1..size_1 + 1 {
        matrix[0][i] = i as i32 * -1;
    }

    for i in 1..size_2 + 1 {
        matrix[i][0] = i as i32 * -1;
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
    mut m: Vec<Vec<i32>>, 
    s_match: i32,
    s_mismatch: i32,
    s_indel: i32,
) -> Vec<Vec<i32>> {
    for i in 1..strand_2.len() + 1 {
        for j in 1..strand_1.len() + 1 {
            let curr_c1 = strand_1.chars().nth(j - 1).unwrap();
            let curr_c2 = strand_2.chars().nth(i - 1).unwrap();

            let score_up = m[i - 1][j] + s_indel;
            let score_left = m[i][j - 1] + s_indel;
            let score_diag = m[i - 1][j - 1] + match_score(curr_c1, curr_c2, s_match, s_mismatch);

            // TODO: change behaviour to check multiple scores
            m[i][j] = cmp::max(score_diag, cmp::max(score_up, score_left));
        }
    }

    m
}

fn print_arr(m: Vec<Vec<i32>>) {
    for i in 0..m.len() {
        println!("{:?}", m[i]);
    }
}

pub fn needleman_wunsch(
    strand_1: &str, strand_2: &str,
    s_match: i32, s_mismatch: i32, s_indel: i32,
    debug: bool
) -> i32 {
    let mut m = setup_matrix(strand_1, strand_2);
    m = fill_table(strand_1, strand_2, m, s_match, s_mismatch, s_indel);

    let last_pos = *m.last().unwrap().last().unwrap();

    if debug {
        print_arr(m);
    }
    
    last_pos
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
}
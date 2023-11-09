use bioinformatics::needleman_wunsch::needleman_wunsch;

#[test]
fn simple() {
    let strand_1 = "GCATGCG";
    let strand_2 = "GATTACA";

    let match_score = 1;
    let mismatch_score = -1;
    let indel_score = -1;

    let debug = false;

    let score = needleman_wunsch(strand_1, strand_2, 
        match_score, mismatch_score, indel_score, debug);

    assert_eq!(score, 0);
}
#include "cass/cass.h"
#include "imm.h"
#include "nmm.h"

void test_hmm_frame_state_0eps(void);
void test_hmm_frame_state_len1(void);
void test_hmm_frame_state_len2(void);
void test_hmm_frame_state_len3(void);
void test_hmm_frame_state_len4(void);
void test_hmm_frame_state_len5(void);

int main(void)
{
    test_hmm_frame_state_0eps();
    test_hmm_frame_state_len1();
    test_hmm_frame_state_len2();
    test_hmm_frame_state_len3();
    test_hmm_frame_state_len4();
    test_hmm_frame_state_len5();
    return cass_status();
}

void test_hmm_frame_state_0eps(void)
{
    struct imm_abc *bases = imm_abc_create("ACGT");
    struct imm_hmm *hmm = imm_hmm_create(bases);
    const int A = imm_abc_symbol_idx(bases, 'A');
    const int C = imm_abc_symbol_idx(bases, 'C');
    const int G = imm_abc_symbol_idx(bases, 'G');
    const int T = imm_abc_symbol_idx(bases, 'T');

    struct nmm_codon *codon = nmm_codon_create();
    nmm_codon_set_ninfs(codon);

    nmm_codon_set_lprob(codon, A, T, G, log(0.8));
    nmm_codon_set_lprob(codon, A, T, T, log(0.1));
    nmm_codon_set_lprob(codon, C, C, C, log(0.1));

    double base_lprobs[] = {log(0.25), log(0.25), log(0.5), LOG0};
    struct nmm_frame_state *state =
        nmm_frame_state_create("M", bases, base_lprobs, codon, 0.0);

    imm_hmm_add_state(hmm, imm_state_cast_c(state), log(1.0));

    struct imm_path *path = imm_path_create();
    imm_path_add(path, imm_state_cast_c(state), 3);
    cass_close(imm_hmm_likelihood(hmm, "ATT", path), -2.3025850930);
    cass_close(imm_hmm_likelihood(hmm, "ATG", path), -0.2231435513142097);
    cass_condition(imm_isnan(imm_hmm_likelihood(hmm, "AT", path)));
    imm_path_destroy(path);

    path = imm_path_create();
    imm_path_add(path, imm_state_cast_c(state), 2);
    cass_condition(imm_isninf(imm_hmm_likelihood(hmm, "AT", path)));
    imm_path_destroy(path);

    imm_hmm_destroy(hmm);
    nmm_frame_state_destroy(state);
    nmm_codon_destroy(codon);
    imm_abc_destroy(bases);
}

void test_hmm_frame_state_len1(void)
{
    struct imm_abc *bases = imm_abc_create("ACGT");
    struct imm_hmm *hmm = imm_hmm_create(bases);
    const int A = imm_abc_symbol_idx(bases, 'A');
    const int C = imm_abc_symbol_idx(bases, 'C');
    const int G = imm_abc_symbol_idx(bases, 'G');
    const int T = imm_abc_symbol_idx(bases, 'T');

    struct nmm_codon *codon = nmm_codon_create();
    nmm_codon_set_ninfs(codon);

    nmm_codon_set_lprob(codon, A, T, G, log(0.8));
    nmm_codon_set_lprob(codon, A, T, T, log(0.1));
    nmm_codon_set_lprob(codon, C, C, C, log(0.1));

    double base_lprobs[] = {log(0.25), log(0.25), log(0.5), LOG0};
    struct nmm_frame_state *state =
        nmm_frame_state_create("M", bases, base_lprobs, codon, 0.1);

    imm_hmm_add_state(hmm, imm_state_cast_c(state), log(1.0));

    struct imm_path *path = imm_path_create();
    imm_path_add(path, imm_state_cast_c(state), 1);
    cass_close(imm_hmm_likelihood(hmm, "A", path), -6.0198640216);
    imm_path_destroy(path);

    path = imm_path_create();
    imm_path_add(path, imm_state_cast_c(state), 1);
    cass_close(imm_hmm_likelihood(hmm, "C", path), -7.118476310297789);
    imm_path_destroy(path);

    cass_close(imm_hmm_viterbi(hmm, "A", imm_state_cast_c(state)), -6.0198640216);
    cass_close(imm_hmm_viterbi(hmm, "C", imm_state_cast_c(state)), -7.1184763103);

    imm_hmm_destroy(hmm);
    nmm_frame_state_destroy(state);
    nmm_codon_destroy(codon);
    imm_abc_destroy(bases);
}

void test_hmm_frame_state_len2(void)
{
    struct imm_abc *bases = imm_abc_create("ACGT");
    struct imm_hmm *hmm = imm_hmm_create(bases);
    const int A = imm_abc_symbol_idx(bases, 'A');
    const int C = imm_abc_symbol_idx(bases, 'C');
    const int G = imm_abc_symbol_idx(bases, 'G');
    const int T = imm_abc_symbol_idx(bases, 'T');

    struct nmm_codon *codon = nmm_codon_create();
    nmm_codon_set_ninfs(codon);

    nmm_codon_set_lprob(codon, A, T, G, log(0.8));
    nmm_codon_set_lprob(codon, A, T, T, log(0.1));
    nmm_codon_set_lprob(codon, C, C, C, log(0.1));

    double base_lprobs[] = {log(0.25), log(0.25), log(0.5), LOG0};
    struct nmm_frame_state *state =
        nmm_frame_state_create("M", bases, base_lprobs, codon, 0.1);

    imm_hmm_add_state(hmm, imm_state_cast_c(state), log(1.0));

    struct imm_path *path = imm_path_create();
    imm_path_add(path, imm_state_cast_c(state), 2);
    cass_close(imm_hmm_likelihood(hmm, "AA", path), -8.910235779525845);
    imm_path_destroy(path);

    path = imm_path_create();
    imm_path_add(path, imm_state_cast_c(state), 2);
    cass_close(imm_hmm_likelihood(hmm, "TG", path), -3.2434246977896133);
    imm_path_destroy(path);

    path = imm_path_create();
    imm_path_add(path, imm_state_cast_c(state), 2);
    cass_close(imm_hmm_likelihood(hmm, "CC", path), -4.225022885864217);
    imm_path_destroy(path);

    path = imm_path_create();
    imm_path_add(path, imm_state_cast_c(state), 2);
    cass_close(imm_hmm_likelihood(hmm, "TT", path), -5.326716841069734);
    imm_path_destroy(path);

    cass_close(imm_hmm_viterbi(hmm, "AA", imm_state_cast_c(state)), -8.910235779525845);
    cass_close(imm_hmm_viterbi(hmm, "TG", imm_state_cast_c(state)), -3.2434246977896133);
    cass_close(imm_hmm_viterbi(hmm, "CC", imm_state_cast_c(state)), -4.225022885864217);
    cass_close(imm_hmm_viterbi(hmm, "TT", imm_state_cast_c(state)), -5.326716841069734);


    imm_hmm_destroy(hmm);
    nmm_frame_state_destroy(state);
    nmm_codon_destroy(codon);
    imm_abc_destroy(bases);
}

void test_hmm_frame_state_len3(void)
{
    struct imm_abc *bases = imm_abc_create("ACGT");
    struct imm_hmm *hmm = imm_hmm_create(bases);
    const int A = imm_abc_symbol_idx(bases, 'A');
    const int C = imm_abc_symbol_idx(bases, 'C');
    const int G = imm_abc_symbol_idx(bases, 'G');
    const int T = imm_abc_symbol_idx(bases, 'T');

    struct nmm_codon *codon = nmm_codon_create();
    nmm_codon_set_ninfs(codon);

    nmm_codon_set_lprob(codon, A, T, G, log(0.8));
    nmm_codon_set_lprob(codon, A, T, T, log(0.1));
    nmm_codon_set_lprob(codon, C, C, C, log(0.1));

    double base_lprobs[] = {log(0.25), log(0.25), log(0.5), LOG0};
    struct nmm_frame_state *state =
        nmm_frame_state_create("M", bases, base_lprobs, codon, 0.1);

    imm_hmm_add_state(hmm, imm_state_cast_c(state), log(1.0));

    struct imm_path *path = imm_path_create();
    imm_path_add(path, imm_state_cast_c(state), 3);
    cass_close(imm_hmm_likelihood(hmm, "ATC", path), -7.012344487235739);
    imm_path_destroy(path);

    path = imm_path_create();
    imm_path_add(path, imm_state_cast_c(state), 3);
    cass_close(imm_hmm_likelihood(hmm, "ATG", path), -0.639793371602465);
    imm_path_destroy(path);

    cass_close(imm_hmm_viterbi(hmm, "ATC", imm_state_cast_c(state)), -7.012344487235739);
    cass_close(imm_hmm_viterbi(hmm, "ATG", imm_state_cast_c(state)), -0.639793371602465);

    imm_hmm_destroy(hmm);
    nmm_frame_state_destroy(state);
    nmm_codon_destroy(codon);
    imm_abc_destroy(bases);
}

void test_hmm_frame_state_len4(void)
{
    struct imm_abc *bases = imm_abc_create("ACGT");
    struct imm_hmm *hmm = imm_hmm_create(bases);
    const int A = imm_abc_symbol_idx(bases, 'A');
    const int C = imm_abc_symbol_idx(bases, 'C');
    const int G = imm_abc_symbol_idx(bases, 'G');
    const int T = imm_abc_symbol_idx(bases, 'T');

    struct nmm_codon *codon = nmm_codon_create();
    nmm_codon_set_ninfs(codon);

    nmm_codon_set_lprob(codon, A, T, G, log(0.8));
    nmm_codon_set_lprob(codon, A, T, T, log(0.1));
    nmm_codon_set_lprob(codon, C, C, C, log(0.1));

    double base_lprobs[] = {log(0.25), log(0.25), log(0.5), LOG0};
    struct nmm_frame_state *state =
        nmm_frame_state_create("M", bases, base_lprobs, codon, 0.1);

    imm_hmm_add_state(hmm, imm_state_cast_c(state), log(1.0));

    struct imm_path *path = imm_path_create();
    imm_path_add(path, imm_state_cast_c(state), 4);
    cass_close(imm_hmm_likelihood(hmm, "ATCC", path), -11.982929094215963);
    imm_path_destroy(path);

    cass_close(imm_hmm_viterbi(hmm, "ATCC", imm_state_cast_c(state)), -11.982929094215963);

    imm_hmm_destroy(hmm);
    nmm_frame_state_destroy(state);
    nmm_codon_destroy(codon);
    imm_abc_destroy(bases);
}

void test_hmm_frame_state_len5(void)
{
    struct imm_abc *bases = imm_abc_create("ACGT");
    struct imm_hmm *hmm = imm_hmm_create(bases);
    const int A = imm_abc_symbol_idx(bases, 'A');
    const int C = imm_abc_symbol_idx(bases, 'C');
    const int G = imm_abc_symbol_idx(bases, 'G');
    const int T = imm_abc_symbol_idx(bases, 'T');

    struct nmm_codon *codon = nmm_codon_create();
    nmm_codon_set_ninfs(codon);

    nmm_codon_set_lprob(codon, A, T, G, log(0.8));
    nmm_codon_set_lprob(codon, A, T, T, log(0.1));
    nmm_codon_set_lprob(codon, C, C, C, log(0.1));

    double base_lprobs[] = {log(0.25), log(0.25), log(0.5), LOG0};
    struct nmm_frame_state *state =
        nmm_frame_state_create("M", bases, base_lprobs, codon, 0.1);

    imm_hmm_add_state(hmm, imm_state_cast_c(state), log(1.0));

    struct imm_path *path = imm_path_create();
    imm_path_add(path, imm_state_cast_c(state), 5);
    cass_condition(imm_isninf(imm_hmm_likelihood(hmm, "ACGTA", path)));
    imm_path_destroy(path);

    path = imm_path_create();
    imm_path_add(path, imm_state_cast_c(state), 5);
    cass_close(imm_hmm_likelihood(hmm, "ACTAG", path), -10.11420858385178);
    imm_path_destroy(path);
    
    cass_condition(imm_isninf(imm_hmm_viterbi(hmm, "ACGTA", imm_state_cast_c(state))));
    cass_close(imm_hmm_viterbi(hmm, "ACTAG", imm_state_cast_c(state)), -10.11420858385178);

    imm_hmm_destroy(hmm);
    nmm_frame_state_destroy(state);
    nmm_codon_destroy(codon);
    imm_abc_destroy(bases);
}

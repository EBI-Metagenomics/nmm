#include "cass/cass.h"
#include "imm/imm.h"
#include "nmm/nmm.h"

void test_hmm_frame_state_0eps(void);
void test_hmm_frame_state_len1(void);
void test_hmm_frame_state_len2(void);
void test_hmm_frame_state_len3(void);
void test_hmm_frame_state_len4(void);
void test_hmm_frame_state_len5(void);

imm_float single_viterbi(struct imm_hmm const* hmm, struct imm_seq const* seq,
                         struct imm_state const* end_state, struct imm_path* path);

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
    struct nmm_base_abc const*   base = nmm_base_abc_create("ACGT", 'X');
    struct imm_abc const*        abc = nmm_base_abc_super(base);
    double const                 zero = imm_lprob_zero();
    struct nmm_base_table const* baset =
        nmm_base_table_create(base, log(0.25), log(0.25), log(0.5), zero);

    struct nmm_codon_lprob* codonp = nmm_codon_lprob_create(base);
    struct nmm_codon*       codon = nmm_codon_create(base);
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'T', 'G')) == 0);
    nmm_codon_lprob_set(codonp, codon, log(0.8));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'T', 'T')) == 0);
    nmm_codon_lprob_set(codonp, codon, log(0.1));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('C', 'C', 'C')) == 0);
    nmm_codon_lprob_set(codonp, codon, log(0.1));
    struct nmm_codon_table const* codont = nmm_codon_table_create(codonp);
    nmm_codon_lprob_destroy(codonp);

    struct imm_hmm* hmm = imm_hmm_create(abc);

    struct nmm_frame_state const* state = nmm_frame_state_create("M", baset, codont, 0.0);

    imm_hmm_add_state(hmm, nmm_frame_state_super(state), log(1.0));

    struct imm_path* path = imm_path_create();
    imm_path_append(path, imm_step_create(nmm_frame_state_super(state), 3));
    struct imm_seq const* seq = imm_seq_create("ATT", abc);
    cass_close(imm_hmm_likelihood(hmm, seq, path), -2.3025850930);
    imm_seq_destroy(seq);
    seq = imm_seq_create("ATG", abc);
    cass_close(imm_hmm_likelihood(hmm, seq, path), -0.2231435513142097);
    imm_seq_destroy(seq);
    seq = imm_seq_create("AT", abc);
    cass_cond(!imm_lprob_is_valid(imm_hmm_likelihood(hmm, seq, path)));
    imm_seq_destroy(seq);
    imm_path_destroy(path);

    path = imm_path_create();
    imm_path_append(path, imm_step_create(nmm_frame_state_super(state), 2));
    seq = imm_seq_create("AT", abc);
    cass_cond(imm_lprob_is_zero(imm_hmm_likelihood(hmm, seq, path)));
    imm_seq_destroy(seq);
    imm_path_destroy(path);

    nmm_base_abc_destroy(base);
    nmm_codon_destroy(codon);
    imm_hmm_destroy(hmm);
    nmm_frame_state_destroy(state);
    nmm_codon_table_destroy(codont);
    nmm_base_table_destroy(baset);
}

void test_hmm_frame_state_len1(void)
{
    struct nmm_base_abc const*   base = nmm_base_abc_create("ACGT", 'X');
    struct imm_abc const*        abc = nmm_base_abc_super(base);
    double const                 zero = imm_lprob_zero();
    struct nmm_base_table const* baset =
        nmm_base_table_create(base, log(0.25), log(0.25), log(0.5), zero);

    struct nmm_codon_lprob* codonp = nmm_codon_lprob_create(base);
    struct nmm_codon*       codon = nmm_codon_create(base);
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'T', 'G')) == 0);
    nmm_codon_lprob_set(codonp, codon, log(0.8));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'T', 'T')) == 0);
    nmm_codon_lprob_set(codonp, codon, log(0.1));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('C', 'C', 'C')) == 0);
    nmm_codon_lprob_set(codonp, codon, log(0.1));
    struct nmm_codon_table const* codont = nmm_codon_table_create(codonp);
    nmm_codon_lprob_destroy(codonp);

    struct imm_hmm* hmm = imm_hmm_create(abc);

    struct nmm_frame_state const* state = nmm_frame_state_create("M", baset, codont, 0.1);

    imm_hmm_add_state(hmm, nmm_frame_state_super(state), log(1.0));

    struct imm_path* path = imm_path_create();
    imm_path_append(path, imm_step_create(nmm_frame_state_super(state), 1));
    struct imm_seq const* seq = imm_seq_create("A", abc);
    cass_close(imm_hmm_likelihood(hmm, seq, path), -6.0198640216);
    imm_seq_destroy(seq);
    imm_path_destroy(path);

    path = imm_path_create();
    imm_path_append(path, imm_step_create(nmm_frame_state_super(state), 1));
    seq = imm_seq_create("C", abc);
    cass_close(imm_hmm_likelihood(hmm, seq, path), -7.118476310297789);
    imm_seq_destroy(seq);
    imm_path_destroy(path);

    path = imm_path_create();
    seq = imm_seq_create("A", abc);
    cass_close(single_viterbi(hmm, seq, nmm_frame_state_super(state), path), -6.0198640216);
    cass_close(imm_hmm_likelihood(hmm, seq, path), -6.0198640216);
    imm_seq_destroy(seq);
    imm_path_destroy(path);

    path = imm_path_create();
    seq = imm_seq_create("C", abc);
    cass_close(single_viterbi(hmm, seq, nmm_frame_state_super(state), path), -7.1184763103);
    cass_close(imm_hmm_likelihood(hmm, seq, path), -7.1184763103);
    imm_seq_destroy(seq);
    imm_path_destroy(path);

    nmm_base_abc_destroy(base);
    nmm_codon_destroy(codon);
    imm_hmm_destroy(hmm);
    nmm_frame_state_destroy(state);
    nmm_codon_table_destroy(codont);
    nmm_base_table_destroy(baset);
}

void test_hmm_frame_state_len2(void)
{
    struct nmm_base_abc const*   base = nmm_base_abc_create("ACGT", 'X');
    struct imm_abc const*        abc = nmm_base_abc_super(base);
    double const                 zero = imm_lprob_zero();
    struct nmm_base_table const* baset =
        nmm_base_table_create(base, log(0.25), log(0.25), log(0.5), zero);

    struct nmm_codon_lprob* codonp = nmm_codon_lprob_create(base);
    struct nmm_codon*       codon = nmm_codon_create(base);
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'T', 'G')) == 0);
    nmm_codon_lprob_set(codonp, codon, log(0.8));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'T', 'T')) == 0);
    nmm_codon_lprob_set(codonp, codon, log(0.1));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('C', 'C', 'C')) == 0);
    nmm_codon_lprob_set(codonp, codon, log(0.1));
    struct nmm_codon_table const* codont = nmm_codon_table_create(codonp);
    nmm_codon_lprob_destroy(codonp);

    struct imm_hmm* hmm = imm_hmm_create(abc);

    struct nmm_frame_state const* state = nmm_frame_state_create("M", baset, codont, 0.1);

    imm_hmm_add_state(hmm, nmm_frame_state_super(state), log(1.0));

    struct imm_path* path = imm_path_create();
    imm_path_append(path, imm_step_create(nmm_frame_state_super(state), 2));
    struct imm_seq const* seq = imm_seq_create("AA", abc);
    cass_close(imm_hmm_likelihood(hmm, seq, path), -8.910235779525845);
    imm_seq_destroy(seq);
    imm_path_destroy(path);

    path = imm_path_create();
    imm_path_append(path, imm_step_create(nmm_frame_state_super(state), 2));
    seq = imm_seq_create("TG", abc);
    cass_close(imm_hmm_likelihood(hmm, seq, path), -3.2434246977896133);
    imm_seq_destroy(seq);
    imm_path_destroy(path);

    path = imm_path_create();
    imm_path_append(path, imm_step_create(nmm_frame_state_super(state), 2));
    seq = imm_seq_create("CC", abc);
    cass_close(imm_hmm_likelihood(hmm, seq, path), -4.225022885864217);
    imm_seq_destroy(seq);
    imm_path_destroy(path);

    path = imm_path_create();
    imm_path_append(path, imm_step_create(nmm_frame_state_super(state), 2));
    seq = imm_seq_create("TT", abc);
    cass_close(imm_hmm_likelihood(hmm, seq, path), -5.326716841069734);
    imm_seq_destroy(seq);
    imm_path_destroy(path);

    path = imm_path_create();
    seq = imm_seq_create("AA", abc);
    cass_close(single_viterbi(hmm, seq, nmm_frame_state_super(state), path), -8.910235779525845);
    cass_close(imm_hmm_likelihood(hmm, seq, path), -8.910235779525845);
    imm_seq_destroy(seq);
    imm_path_destroy(path);

    path = imm_path_create();
    seq = imm_seq_create("TG", abc);
    cass_close(single_viterbi(hmm, seq, nmm_frame_state_super(state), path), -3.2434246977896133);
    cass_close(imm_hmm_likelihood(hmm, seq, path), -3.2434246977896133);
    imm_seq_destroy(seq);
    imm_path_destroy(path);

    path = imm_path_create();
    seq = imm_seq_create("CC", abc);
    cass_close(single_viterbi(hmm, seq, nmm_frame_state_super(state), path), -4.225022885864217);
    cass_close(imm_hmm_likelihood(hmm, seq, path), -4.225022885864217);
    imm_seq_destroy(seq);
    imm_path_destroy(path);

    path = imm_path_create();
    seq = imm_seq_create("TT", abc);
    cass_close(single_viterbi(hmm, seq, nmm_frame_state_super(state), path), -5.326716841069734);
    cass_close(imm_hmm_likelihood(hmm, seq, path), -5.326716841069734);
    imm_seq_destroy(seq);
    imm_path_destroy(path);

    nmm_base_abc_destroy(base);
    nmm_codon_destroy(codon);
    imm_hmm_destroy(hmm);
    nmm_frame_state_destroy(state);
    nmm_codon_table_destroy(codont);
    nmm_base_table_destroy(baset);
}

void test_hmm_frame_state_len3(void)
{
    struct nmm_base_abc const*   base = nmm_base_abc_create("ACGT", 'X');
    struct imm_abc const*        abc = nmm_base_abc_super(base);
    double const                 zero = imm_lprob_zero();
    struct nmm_base_table const* baset =
        nmm_base_table_create(base, log(0.25), log(0.25), log(0.5), zero);

    struct nmm_codon_lprob* codonp = nmm_codon_lprob_create(base);
    struct nmm_codon*       codon = nmm_codon_create(base);
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'T', 'G')) == 0);
    nmm_codon_lprob_set(codonp, codon, log(0.8));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'T', 'T')) == 0);
    nmm_codon_lprob_set(codonp, codon, log(0.1));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('C', 'C', 'C')) == 0);
    nmm_codon_lprob_set(codonp, codon, log(0.1));
    struct nmm_codon_table const* codont = nmm_codon_table_create(codonp);
    nmm_codon_lprob_destroy(codonp);

    struct imm_hmm* hmm = imm_hmm_create(abc);

    struct nmm_frame_state const* state = nmm_frame_state_create("M", baset, codont, 0.1);

    imm_hmm_add_state(hmm, nmm_frame_state_super(state), log(1.0));

    struct imm_path* path = imm_path_create();
    imm_path_append(path, imm_step_create(nmm_frame_state_super(state), 3));
    struct imm_seq const* seq = imm_seq_create("ATC", abc);
    cass_close(imm_hmm_likelihood(hmm, seq, path), -7.012344487235739);
    imm_seq_destroy(seq);
    imm_path_destroy(path);

    path = imm_path_create();
    imm_path_append(path, imm_step_create(nmm_frame_state_super(state), 3));
    seq = imm_seq_create("ATG", abc);
    cass_close(imm_hmm_likelihood(hmm, seq, path), -0.639793371602465);
    imm_seq_destroy(seq);
    imm_path_destroy(path);

    path = imm_path_create();
    seq = imm_seq_create("ATC", abc);
    cass_close(single_viterbi(hmm, seq, nmm_frame_state_super(state), path), -7.012344487235739);
    cass_close(imm_hmm_likelihood(hmm, seq, path), -7.012344487235739);
    imm_seq_destroy(seq);
    imm_path_destroy(path);

    path = imm_path_create();
    seq = imm_seq_create("ATG", abc);
    cass_close(single_viterbi(hmm, seq, nmm_frame_state_super(state), path), -0.639793371602465);
    cass_close(imm_hmm_likelihood(hmm, seq, path), -0.639793371602465);
    imm_seq_destroy(seq);
    imm_path_destroy(path);

    nmm_base_abc_destroy(base);
    nmm_codon_destroy(codon);
    imm_hmm_destroy(hmm);
    nmm_frame_state_destroy(state);
    nmm_codon_table_destroy(codont);
    nmm_base_table_destroy(baset);
}

void test_hmm_frame_state_len4(void)
{
    struct nmm_base_abc const*   base = nmm_base_abc_create("ACGT", 'X');
    struct imm_abc const*        abc = nmm_base_abc_super(base);
    double const                 zero = imm_lprob_zero();
    struct nmm_base_table const* baset =
        nmm_base_table_create(base, log(0.25), log(0.25), log(0.5), zero);

    struct nmm_codon_lprob* codonp = nmm_codon_lprob_create(base);
    struct nmm_codon*       codon = nmm_codon_create(base);
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'T', 'G')) == 0);
    nmm_codon_lprob_set(codonp, codon, log(0.8));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'T', 'T')) == 0);
    nmm_codon_lprob_set(codonp, codon, log(0.1));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('C', 'C', 'C')) == 0);
    nmm_codon_lprob_set(codonp, codon, log(0.1));
    struct nmm_codon_table const* codont = nmm_codon_table_create(codonp);
    nmm_codon_lprob_destroy(codonp);

    struct imm_hmm* hmm = imm_hmm_create(abc);

    struct nmm_frame_state const* state = nmm_frame_state_create("M", baset, codont, 0.1);

    imm_hmm_add_state(hmm, nmm_frame_state_super(state), log(1.0));

    struct imm_path* path = imm_path_create();
    imm_path_append(path, imm_step_create(nmm_frame_state_super(state), 4));
    struct imm_seq const* seq = imm_seq_create("ATCC", abc);
    cass_close(imm_hmm_likelihood(hmm, seq, path), -11.982929094215963);
    imm_seq_destroy(seq);
    imm_path_destroy(path);

    path = imm_path_create();
    seq = imm_seq_create("ATCC", abc);
    cass_close(single_viterbi(hmm, seq, nmm_frame_state_super(state), path), -11.982929094215963);
    cass_close(imm_hmm_likelihood(hmm, seq, path), -11.982929094215963);
    imm_seq_destroy(seq);
    imm_path_destroy(path);

    nmm_base_abc_destroy(base);
    nmm_codon_destroy(codon);
    imm_hmm_destroy(hmm);
    nmm_frame_state_destroy(state);
    nmm_codon_table_destroy(codont);
    nmm_base_table_destroy(baset);
}

void test_hmm_frame_state_len5(void)
{
    struct nmm_base_abc const*   base = nmm_base_abc_create("ACGT", 'X');
    struct imm_abc const*        abc = nmm_base_abc_super(base);
    double const                 zero = imm_lprob_zero();
    struct nmm_base_table const* baset =
        nmm_base_table_create(base, log(0.25), log(0.25), log(0.5), zero);

    struct nmm_codon_lprob* codonp = nmm_codon_lprob_create(base);
    struct nmm_codon*       codon = nmm_codon_create(base);
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'T', 'G')) == 0);
    nmm_codon_lprob_set(codonp, codon, log(0.8));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'T', 'T')) == 0);
    nmm_codon_lprob_set(codonp, codon, log(0.1));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('C', 'C', 'C')) == 0);
    nmm_codon_lprob_set(codonp, codon, log(0.1));
    struct nmm_codon_table const* codont = nmm_codon_table_create(codonp);
    nmm_codon_lprob_destroy(codonp);

    struct imm_hmm* hmm = imm_hmm_create(abc);

    struct nmm_frame_state const* state = nmm_frame_state_create("M", baset, codont, 0.1);

    imm_hmm_add_state(hmm, nmm_frame_state_super(state), log(1.0));

    struct imm_path* path = imm_path_create();
    imm_path_append(path, imm_step_create(nmm_frame_state_super(state), 5));
    struct imm_seq const* seq = imm_seq_create("ACGTA", abc);
    cass_cond(imm_lprob_is_zero(imm_hmm_likelihood(hmm, seq, path)));
    imm_seq_destroy(seq);
    imm_path_destroy(path);

    path = imm_path_create();
    imm_path_append(path, imm_step_create(nmm_frame_state_super(state), 5));
    seq = imm_seq_create("ACTAG", abc);
    cass_close(imm_hmm_likelihood(hmm, seq, path), -10.11420858385178);
    imm_seq_destroy(seq);
    imm_path_destroy(path);

    path = imm_path_create();
    seq = imm_seq_create("ACGTA", abc);
    cass_cond(!imm_lprob_is_valid(single_viterbi(hmm, seq, nmm_frame_state_super(state), path)));
    cass_cond(!imm_lprob_is_valid(imm_hmm_likelihood(hmm, seq, path)));
    imm_seq_destroy(seq);
    imm_path_destroy(path);

    path = imm_path_create();
    seq = imm_seq_create("ACTAG", abc);
    cass_close(single_viterbi(hmm, seq, nmm_frame_state_super(state), path), -10.11420858385178);
    cass_close(imm_hmm_likelihood(hmm, seq, path), -10.11420858385178);
    imm_seq_destroy(seq);
    imm_path_destroy(path);

    nmm_base_abc_destroy(base);
    nmm_codon_destroy(codon);
    imm_hmm_destroy(hmm);
    nmm_frame_state_destroy(state);
    nmm_codon_table_destroy(codont);
    nmm_base_table_destroy(baset);
}

imm_float single_viterbi(struct imm_hmm const* hmm, struct imm_seq const* seq,
                         struct imm_state const* end_state, struct imm_path* path)
{
    struct imm_dp const*      dp = imm_hmm_create_dp(hmm, end_state);
    struct imm_results const* results = imm_dp_viterbi(dp, seq, 0);
    if (results == NULL)
        return imm_lprob_invalid();
    struct imm_result const* r = imm_results_get(results, 0);

    struct imm_path const* src = imm_result_path(r);
    for (struct imm_step const* step = imm_path_first(src); step; step = imm_path_next(src, step))
        imm_path_append(path, imm_step_create(imm_step_state(step), imm_step_seq_len(step)));

    struct imm_subseq subseq = imm_result_subseq(r);
    imm_float         loglik = imm_hmm_likelihood(hmm, imm_subseq_cast(&subseq), path);
    imm_results_destroy(results);
    imm_dp_destroy(dp);

    return loglik;
}

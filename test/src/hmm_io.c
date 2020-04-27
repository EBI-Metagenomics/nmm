#include "cass.h"
#include "imm/imm.h"
#include "nmm/nmm.h"

void test_hmm_io(void);

int main(void)
{
    test_hmm_io();
    return cass_status();
}

void test_hmm_io(void)
{
    struct nmm_base_abc const*   base = nmm_base_abc_create("ACGT", 'X');
    struct imm_abc const*        abc = nmm_base_abc_parent(base);
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

    imm_hmm_add_state(hmm, nmm_frame_state_parent(state), log(1.0));

    struct imm_path* path = imm_path_create();
    imm_path_append(path, imm_step_create(nmm_frame_state_parent(state), 1));
    struct imm_seq const* seq = imm_seq_create("A", abc);
    cass_close(imm_hmm_likelihood(hmm, seq, path), -6.0198640216);
    imm_seq_destroy(seq);
    imm_path_destroy(path);

    path = imm_path_create();
    imm_path_append(path, imm_step_create(nmm_frame_state_parent(state), 1));
    seq = imm_seq_create("C", abc);
    cass_close(imm_hmm_likelihood(hmm, seq, path), -7.118476310297789);
    imm_seq_destroy(seq);
    imm_path_destroy(path);

    struct imm_dp const* dp = imm_hmm_create_dp(hmm, nmm_frame_state_parent(state));
    cass_cond(dp != NULL);

    seq = imm_seq_create("A", abc);
    struct imm_results const* results = imm_dp_viterbi(dp, seq, 0);
    cass_cond(results != NULL);
    cass_cond(imm_results_size(results) == 1);
    struct imm_result const* r = imm_results_get(results, 0);
    cass_close(imm_result_loglik(r), -6.0198640216);
    cass_close(imm_hmm_likelihood(hmm, seq, imm_result_path(r)), -6.0198640216);
    imm_seq_destroy(seq);
    imm_results_destroy(results);

    nmm_base_abc_destroy(base);
    nmm_codon_destroy(codon);
    imm_hmm_destroy(hmm);
    nmm_frame_state_destroy(state);
    nmm_codon_table_destroy(codont);
    nmm_base_table_destroy(baset);
    imm_dp_destroy(dp);
}

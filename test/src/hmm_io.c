#include "cass.h"
#include "imm/imm.h"
#include "nmm/nmm.h"

#define TMP_FOLDER "test_hmm_io.tmp"

void test_hmm_io(void);

int main(void)
{
    test_hmm_io();
    return cass_status();
}

void test_hmm_io(void)
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

    struct imm_hmm* hmm = imm_hmm_create(abc);

    struct nmm_frame_state const* state1 = nmm_frame_state_create("S0", baset, codont, 0.1);
    struct nmm_codon_state const* state2 = nmm_codon_state_create("S1", codonp);

    imm_hmm_add_state(hmm, nmm_frame_state_super(state1), log(1.0));
    imm_hmm_add_state(hmm, nmm_codon_state_super(state2), log(0.0001));
    imm_hmm_set_trans(hmm, nmm_frame_state_super(state1), nmm_codon_state_super(state2), log(0.2));

    struct imm_path* path = imm_path_create();
    imm_path_append(path, imm_step_create(nmm_frame_state_super(state1), 1));
    struct imm_seq const* seq = imm_seq_create("A", abc);
    cass_close(imm_hmm_likelihood(hmm, seq, path), -6.0198640216);
    imm_seq_destroy(seq);
    imm_path_destroy(path);

    path = imm_path_create();
    imm_path_append(path, imm_step_create(nmm_frame_state_super(state1), 1));
    seq = imm_seq_create("C", abc);
    cass_close(imm_hmm_likelihood(hmm, seq, path), -7.118476310297789);
    imm_seq_destroy(seq);
    imm_path_destroy(path);

    struct imm_dp const* dp = imm_hmm_create_dp(hmm, nmm_frame_state_super(state1));
    cass_cond(dp != NULL);

    seq = imm_seq_create("A", abc);
    struct imm_results const* results = imm_dp_viterbi(dp, seq, 0);
    cass_cond(results != NULL);
    cass_cond(imm_results_size(results) == 1);
    struct imm_result const* r = imm_results_get(results, 0);
    struct imm_subseq        subseq = imm_result_subseq(r);
    struct imm_path const*   p = imm_result_path(r);
    imm_float                loglik = imm_hmm_likelihood(hmm, imm_subseq_cast(&subseq), p);
    cass_close(loglik, -6.0198640216);
    cass_close(imm_hmm_likelihood(hmm, seq, imm_result_path(r)), -6.0198640216);
    imm_seq_destroy(seq);
    imm_results_destroy(results);

    struct nmm_output* output = nmm_output_create(TMP_FOLDER "/two_states.nmm");
    cass_cond(output != NULL);
    struct nmm_model const* model = nmm_model_create(hmm, dp);
    cass_equal_int(nmm_output_write(output, model), 0);
    nmm_model_destroy(model);
    cass_equal_int(nmm_output_destroy(output), 0);

    nmm_base_abc_destroy(base);
    nmm_codon_destroy(codon);
    imm_hmm_destroy(hmm);
    nmm_frame_state_destroy(state1);
    nmm_codon_state_destroy(state2);
    nmm_codon_table_destroy(codont);
    nmm_base_table_destroy(baset);
    imm_dp_destroy(dp);
    nmm_codon_lprob_destroy(codonp);

    struct nmm_input* input = nmm_input_create(TMP_FOLDER "/two_states.nmm");
    cass_cond(input != NULL);
    cass_cond(!nmm_input_eof(input));
    model = nmm_input_read(input);
    cass_cond(!nmm_input_eof(input));
    cass_cond(model != NULL);
    nmm_input_destroy(input);

    cass_equal_uint64(nmm_model_nstates(model), 2);

    abc = nmm_model_abc(model);
    hmm = nmm_model_hmm(model);
    dp = nmm_model_dp(model);

    seq = imm_seq_create("A", abc);
    results = imm_dp_viterbi(dp, seq, 0);
    cass_cond(results != NULL);
    cass_cond(imm_results_size(results) == 1);
    r = imm_results_get(results, 0);
    subseq = imm_result_subseq(r);
    p = imm_result_path(r);
    loglik = imm_hmm_likelihood(hmm, imm_subseq_cast(&subseq), p);
    cass_close(loglik, -6.0198640216);
    cass_close(imm_hmm_likelihood(hmm, seq, imm_result_path(r)), -6.0198640216);
    imm_results_destroy(results);

    for (uint16_t i = 0; i < nmm_model_nstates(model); ++i)
        imm_state_destroy(nmm_model_state(model, i));

    for (uint32_t i = 0; i < nmm_model_nbase_tables(model); ++i)
        nmm_base_table_destroy(nmm_model_base_table(model, i));

    for (uint32_t i = 0; i < nmm_model_ncodon_tables(model); ++i)
        nmm_codon_table_destroy(nmm_model_codon_table(model, i));

    for (uint32_t i = 0; i < nmm_model_ncodon_lprobs(model); ++i)
        nmm_codon_lprob_destroy(nmm_model_codon_lprob(model, i));

    imm_seq_destroy(seq);
    imm_abc_destroy(abc);
    imm_hmm_destroy(hmm);
    imm_dp_destroy(dp);
    nmm_model_destroy(model);
}

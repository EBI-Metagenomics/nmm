#include "cass/cass.h"
#include "helper.h"
#include "nmm/nmm.h"

#ifndef TMPDIR
#define TMPDIR ""
#endif

void test_hmm_io_two_states(void);
void test_hmm_io_two_hmm_models(void);

int main(void)
{
    test_hmm_io_two_states();
    test_hmm_io_two_hmm_models();
    return cass_status();
}

void test_hmm_io_two_states(void)
{
    struct nmm_base_abc const*   base = nmm_base_abc_create("ACGT", 'X');
    struct imm_abc const*        abc = nmm_base_abc_super(base);
    struct nmm_base_lprob const* basep =
        nmm_base_lprob_create(base, imm_log(0.25), imm_log(0.25), imm_log(0.5), zero());

    struct nmm_codon_lprob* codonp = nmm_codon_lprob_create(base);
    struct nmm_codon*       codon = nmm_codon_create(base);
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'T', 'G')) == 0);
    nmm_codon_lprob_set(codonp, codon, imm_log(0.8));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'T', 'T')) == 0);
    nmm_codon_lprob_set(codonp, codon, imm_log(0.1));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('C', 'C', 'C')) == 0);
    nmm_codon_lprob_set(codonp, codon, imm_log(0.1));
    struct nmm_codon_marg const* codont = nmm_codon_marg_create(codonp);

    struct imm_hmm* hmm = imm_hmm_create(abc);

    struct nmm_frame_state const* state1 = nmm_frame_state_create("S0", basep, codont, (imm_float)0.1);
    struct nmm_codon_state const* state2 = nmm_codon_state_create("S1", codonp);

    imm_hmm_add_state(hmm, nmm_frame_state_super(state1), imm_log(1.0));
    imm_hmm_add_state(hmm, nmm_codon_state_super(state2), imm_log(0.0001));
    imm_hmm_set_trans(hmm, nmm_frame_state_super(state1), nmm_codon_state_super(state2), imm_log(0.2));

    struct imm_path* path = imm_path_create();
    imm_path_append(path, imm_step_create(nmm_frame_state_super(state1), 1));
    struct imm_seq const* seq = imm_seq_create("A", abc);
    cass_close(imm_hmm_loglikelihood(hmm, seq, path), -6.0198640216);
    imm_seq_destroy(seq);
    imm_path_destroy(path);

    path = imm_path_create();
    imm_path_append(path, imm_step_create(nmm_frame_state_super(state1), 1));
    seq = imm_seq_create("C", abc);
    cass_close(imm_hmm_loglikelihood(hmm, seq, path), -7.118476310297789);
    imm_seq_destroy(seq);
    imm_path_destroy(path);

    struct imm_dp const* dp = imm_hmm_create_dp(hmm, nmm_frame_state_super(state1));
    cass_cond(dp != NULL);

    seq = imm_seq_create("A", abc);
    struct imm_dp_task* task = imm_dp_task_create(dp);
    imm_dp_task_setup(task, seq, 0);
    struct imm_results const* results = imm_dp_viterbi(dp, task);
    imm_dp_task_destroy(task);
    cass_cond(results != NULL);
    cass_cond(imm_results_size(results) == 1);
    struct imm_result const* r = imm_results_get(results, 0);
    struct imm_subseq        subseq = imm_result_subseq(r);
    struct imm_path const*   p_c = imm_result_path(r);
    imm_float                loglik = imm_hmm_loglikelihood(hmm, imm_subseq_cast(&subseq), p_c);
    cass_close(loglik, -6.0198640216);
    cass_close(imm_hmm_loglikelihood(hmm, seq, imm_result_path(r)), -6.0198640216);
    imm_seq_destroy(seq);
    imm_results_destroy(results);

    struct nmm_output* output = nmm_output_create(TMPDIR "/two_states.nmm");
    cass_cond(output != NULL);
    struct nmm_profile* p = nmm_profile_create(abc);
    nmm_profile_append_model(p, imm_model_create(hmm, dp));
    cass_equal(nmm_output_write(output, p), 0);
    nmm_profile_destroy(p, false);
    cass_equal(nmm_output_destroy(output), 0);

    nmm_base_abc_destroy(base);
    nmm_codon_destroy(codon);
    imm_hmm_destroy(hmm);
    nmm_frame_state_destroy(state1);
    nmm_codon_state_destroy(state2);
    nmm_codon_marg_destroy(codont);
    nmm_base_lprob_destroy(basep);
    imm_dp_destroy(dp);
    nmm_codon_lprob_destroy(codonp);

    struct nmm_input* input = nmm_input_create(TMPDIR "/two_states.nmm");
    cass_cond(input != NULL);
    cass_cond(!nmm_input_eof(input));
    struct nmm_profile const* prof = nmm_input_read(input);
    cass_cond(!nmm_input_eof(input));
    cass_cond(prof != NULL);
    nmm_input_destroy(input);

    struct imm_model* model = nmm_profile_get_model(prof, 0);

    cass_equal(imm_model_nstates(model), 2);

    abc = nmm_profile_abc(prof);
    hmm = imm_model_hmm(model);
    dp = imm_model_dp(model);

    seq = imm_seq_create("A", abc);
    task = imm_dp_task_create(dp);
    imm_dp_task_setup(task, seq, 0);
    results = imm_dp_viterbi(dp, task);
    cass_cond(results != NULL);
    cass_cond(imm_results_size(results) == 1);
    r = imm_results_get(results, 0);
    subseq = imm_result_subseq(r);
    p_c = imm_result_path(r);
    loglik = imm_hmm_loglikelihood(hmm, imm_subseq_cast(&subseq), p_c);
    cass_close(loglik, -6.0198640216);
    cass_close(imm_hmm_loglikelihood(hmm, seq, imm_result_path(r)), -6.0198640216);
    imm_results_destroy(results);

    imm_seq_destroy(seq);
    imm_dp_task_destroy(task);
    nmm_profile_destroy(prof, true);
}

void test_hmm_io_two_hmm_models(void)
{
    struct nmm_base_abc const*   base = nmm_base_abc_create("ACGT", 'X');
    struct imm_abc const*        abc = nmm_base_abc_super(base);
    struct nmm_base_lprob const* basep =
        nmm_base_lprob_create(base, imm_log(0.25), imm_log(0.25), imm_log(0.5), zero());

    struct nmm_codon_lprob* codonp = nmm_codon_lprob_create(base);
    struct nmm_codon*       codon = nmm_codon_create(base);
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'T', 'G')) == 0);
    nmm_codon_lprob_set(codonp, codon, imm_log(0.8));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'T', 'T')) == 0);
    nmm_codon_lprob_set(codonp, codon, imm_log(0.1));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('C', 'C', 'C')) == 0);
    nmm_codon_lprob_set(codonp, codon, imm_log(0.1));
    struct nmm_codon_marg const* codont = nmm_codon_marg_create(codonp);

    struct imm_hmm* hmm0 = imm_hmm_create(abc);
    struct imm_hmm* hmm1 = imm_hmm_create(abc);

    struct nmm_frame_state const* state1 = nmm_frame_state_create("S0", basep, codont, (imm_float)0.1);
    struct nmm_codon_state const* state2 = nmm_codon_state_create("S1", codonp);

    imm_hmm_add_state(hmm0, nmm_frame_state_super(state1), imm_log(1.0));
    imm_hmm_add_state(hmm0, nmm_codon_state_super(state2), imm_log(0.0001));
    imm_hmm_set_trans(hmm0, nmm_frame_state_super(state1), nmm_codon_state_super(state2), imm_log(0.2));

    imm_hmm_add_state(hmm1, nmm_frame_state_super(state1), imm_log(0.5));
    imm_hmm_add_state(hmm1, nmm_codon_state_super(state2), imm_log(0.01));
    imm_hmm_set_trans(hmm1, nmm_frame_state_super(state1), nmm_codon_state_super(state2), imm_log(0.5));

    struct imm_path* path = imm_path_create();
    imm_path_append(path, imm_step_create(nmm_frame_state_super(state1), 1));
    struct imm_seq const* seq = imm_seq_create("A", abc);
    cass_close(imm_hmm_loglikelihood(hmm0, seq, path), -6.0198640216);
    imm_seq_destroy(seq);
    imm_path_destroy(path);

    path = imm_path_create();
    imm_path_append(path, imm_step_create(nmm_frame_state_super(state1), 1));
    seq = imm_seq_create("A", abc);
    cass_close(imm_hmm_loglikelihood(hmm1, seq, path), -6.7130112648);
    imm_seq_destroy(seq);
    imm_path_destroy(path);

    path = imm_path_create();
    imm_path_append(path, imm_step_create(nmm_frame_state_super(state1), 1));
    seq = imm_seq_create("C", abc);
    cass_close(imm_hmm_loglikelihood(hmm0, seq, path), -7.118476310297789);
    imm_seq_destroy(seq);
    imm_path_destroy(path);

    struct imm_dp const* dp0 = imm_hmm_create_dp(hmm0, nmm_frame_state_super(state1));
    cass_cond(dp0 != NULL);
    struct imm_dp const* dp1 = imm_hmm_create_dp(hmm1, nmm_frame_state_super(state1));
    cass_cond(dp1 != NULL);

    seq = imm_seq_create("A", abc);
    struct imm_dp_task* task0 = imm_dp_task_create(dp0);
    imm_dp_task_setup(task0, seq, 0);
    struct imm_results const* results = imm_dp_viterbi(dp0, task0);
    imm_dp_task_destroy(task0);
    cass_cond(results != NULL);
    cass_cond(imm_results_size(results) == 1);
    struct imm_result const* r = imm_results_get(results, 0);
    struct imm_subseq        subseq = imm_result_subseq(r);
    struct imm_path const*   p = imm_result_path(r);
    imm_float                loglik = imm_hmm_loglikelihood(hmm0, imm_subseq_cast(&subseq), p);
    cass_close(loglik, -6.0198640216);
    cass_close(imm_hmm_loglikelihood(hmm0, seq, imm_result_path(r)), -6.0198640216);
    imm_seq_destroy(seq);
    imm_results_destroy(results);

    seq = imm_seq_create("A", abc);
    struct imm_dp_task* task1 = imm_dp_task_create(dp1);
    imm_dp_task_setup(task1, seq, 0);
    results = imm_dp_viterbi(dp1, task1);
    imm_dp_task_destroy(task1);
    cass_cond(results != NULL);
    cass_cond(imm_results_size(results) == 1);
    r = imm_results_get(results, 0);
    subseq = imm_result_subseq(r);
    p = imm_result_path(r);
    loglik = imm_hmm_loglikelihood(hmm1, imm_subseq_cast(&subseq), p);
    cass_close(loglik, -6.7130112648);
    cass_close(imm_hmm_loglikelihood(hmm1, seq, imm_result_path(r)), -6.7130112648);
    imm_seq_destroy(seq);
    imm_results_destroy(results);

    struct nmm_output* output = nmm_output_create(TMPDIR "/two_hmms.nmm");
    cass_cond(output != NULL);
    struct nmm_profile* m = nmm_profile_create(abc);
    nmm_profile_append_model(m, imm_model_create(hmm0, dp0));
    nmm_profile_append_model(m, imm_model_create(hmm1, dp1));
    cass_equal(nmm_output_write(output, m), 0);
    nmm_profile_destroy(m, false);
    cass_equal(nmm_output_destroy(output), 0);

    nmm_base_abc_destroy(base);
    nmm_codon_destroy(codon);
    imm_hmm_destroy(hmm0);
    imm_hmm_destroy(hmm1);
    nmm_frame_state_destroy(state1);
    nmm_codon_state_destroy(state2);
    nmm_codon_marg_destroy(codont);
    nmm_base_lprob_destroy(basep);
    imm_dp_destroy(dp0);
    imm_dp_destroy(dp1);
    nmm_codon_lprob_destroy(codonp);

    struct nmm_input* input = nmm_input_create(TMPDIR "/two_hmms.nmm");
    cass_cond(input != NULL);
    cass_cond(!nmm_input_eof(input));
    struct nmm_profile const* prof = nmm_input_read(input);
    cass_cond(!nmm_input_eof(input));
    cass_cond(prof != NULL);
    nmm_input_destroy(input);

    cass_equal(nmm_profile_nmodels(prof), 2);
    struct imm_model* model0 = nmm_profile_get_model(prof, 0);
    struct imm_model* model1 = nmm_profile_get_model(prof, 1);

    cass_equal(imm_model_nstates(model0), 2);
    cass_equal(imm_model_nstates(model1), 2);

    abc = nmm_profile_abc(prof);
    hmm0 = imm_model_hmm(model0);
    dp0 = imm_model_dp(model0);
    hmm1 = imm_model_hmm(model1);
    dp1 = imm_model_dp(model1);

    seq = imm_seq_create("A", abc);
    task0 = imm_dp_task_create(dp0);
    imm_dp_task_setup(task0, seq, 0);
    results = imm_dp_viterbi(dp0, task0);
    cass_cond(results != NULL);
    cass_cond(imm_results_size(results) == 1);
    r = imm_results_get(results, 0);
    subseq = imm_result_subseq(r);
    p = imm_result_path(r);
    loglik = imm_hmm_loglikelihood(hmm0, imm_subseq_cast(&subseq), p);
    cass_close(loglik, -6.0198640216);
    cass_close(imm_hmm_loglikelihood(hmm0, seq, imm_result_path(r)), -6.0198640216);
    imm_results_destroy(results);
    imm_seq_destroy(seq);

    seq = imm_seq_create("A", abc);
    task1 = imm_dp_task_create(dp1);
    imm_dp_task_setup(task1, seq, 0);
    results = imm_dp_viterbi(dp1, task1);
    cass_cond(results != NULL);
    cass_cond(imm_results_size(results) == 1);
    r = imm_results_get(results, 0);
    subseq = imm_result_subseq(r);
    p = imm_result_path(r);
    loglik = imm_hmm_loglikelihood(hmm1, imm_subseq_cast(&subseq), p);
    cass_close(loglik, -6.7130112648);
    cass_close(imm_hmm_loglikelihood(hmm1, seq, imm_result_path(r)), -6.7130112648);
    imm_results_destroy(results);

    imm_seq_destroy(seq);
    imm_dp_task_destroy(task0);
    imm_dp_task_destroy(task1);
    nmm_profile_destroy(prof, true);
}

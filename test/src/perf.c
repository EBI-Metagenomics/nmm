#include "cass/cass.h"
#include "helper.h"
#include "nmm/nmm.h"

#ifndef TMPDIR
#define TMPDIR ""
#endif

void test_perf_viterbi(void);

int main(void)
{
    test_perf_viterbi();
    return cass_status();
}

static inline struct imm_state const* frame_super(struct nmm_frame_state const* state)
{
    return nmm_frame_state_super(state);
}
static inline struct imm_state const* mute_super(struct imm_mute_state const* state)
{
    return imm_mute_state_super(state);
}

static struct nmm_codon_lprob* create_codonp(struct nmm_base_abc const* base);

static struct nmm_codon_marg const* create_codont(struct nmm_base_abc const* base, struct nmm_triplet const triplet,
                                                  imm_float const lprob)
{
    struct nmm_codon_lprob* codonp = create_codonp(base);
    struct nmm_codon*       codon = nmm_codon_create(base);
    cass_cond(nmm_codon_set_triplet(codon, triplet) == 0);
    nmm_codon_lprob_set(codonp, codon, lprob);
    struct nmm_codon_marg const* codont = nmm_codon_marg_create(codonp);
    nmm_codon_destroy(codon);
    nmm_codon_lprob_destroy(codonp);
    return codont;
}

static char const __seq[] = "AAAACGCGTGTCACGACAACGCGTACGTTTCGACGAGTACGACGCCCGGG"
                            "AAAACGCGTGTCGACGACGAACGCGTACGTTTACGACGAGTACGACGCCC"
                            "AAAACGCGTGTCACGACAACGCGTACGTTTCGACGAGTACGACGCCCGGG"
                            "AAAACGCGTGTCGACGACGAACGCGTACGTTTACGACGAGTACGACGCCC"
                            "AAAACGCGTGTCACGACAACGCGTACGTTTCGACGAGTACGACGCCCGGG"
                            "AAAACGCGTGTCGACGACGAACGCGTACGTTTACGACGAGTACGACGCCC"
                            "AAAACGCGTGTCACGACAACGCGTACGTTTCGACGAGTACGACGCCCGGG"
                            "AAAACGCGTGTCGACGACGAACGCGTACGTTTACGACGAGTACGACGCCC"
                            "AAAACGCGTGTCACGACAACGCGTACGTTTCGACGAGTACGACGCCCGGG"
                            "AAAACGCGTGTCGACGACGAACGCGTACGTTTACGACGAGTACGACGCCC"
                            "AAAACGCGTGTCACGACAACGCGTACGTTTCGACGAGTACGACGCCCGGG"
                            "AAAACGCGTGTCGACGACGAACGCGTACGTTTACGACGAGTACGACGCCC"
                            "AAAACGCGTGTCACGACAACGCGTACGTTTCGACGAGTACGACGCCCGGG"
                            "AAAACGCGTGTCGACGACGAACGCGTACGTTTACGACGAGTACGACGCCC"
                            "AAAACGCGTGTCACGACAACGCGTACGTTTCGACGAGTACGACGCCCGGG"
                            "AAAACGCGTGTCGACGACGAACGCGTACGTTTACGACGAGTACGACGCCC"
                            "AAAACGCGTGTCACGACAACGCGTACGTTTCGACGAGTACGACGCCCGGG"
                            "AAAACGCGTGTCGACGACGAACGCGTACGTTTACGACGAGTACGACGCCC"
                            "AAAACGCGTGTCACGACAACGCGTACGTTTCGACGAGTACGACGCCCGGG"
                            "AAAACGCGTGTCGACGACGAACGCGTACGTTTACGACGAGTACGACGCCC"
                            "AAAACGCGTGTCACGACAACGCGTACGTTTCGACGAGTACGACGCCCGGG"
                            "AAAACGCGTGTCGACGACGAACGCGTACGTTTACGACGAGTACGACGCCC"
                            "AAAACGCGTGTCACGACAACGCGTACGTTTCGACGAGTACGACGCCCGGG"
                            "AAAACGCGTGTCGACGACGAACGCGTACGTTTACGACGAGTACGACGCCC"
                            "AAAACGCGTGTCACGACAACGCGTACGTTTCGACGAGTACGACGCCCGGG"
                            "AAAACGCGTGTCGACGACGAACGCGTACGTTTACGACGAGTACGACGCCC"
                            "AAAACGCGTGTCACGACAACGCGTACGTTTCGACGAGTACGACGCCCGGG"
                            "AAAACGCGTGTCGACGACGAACGCGTACGTTTACGACGAGTACGACGCCC"
                            "AAAACGCGTGTCACGACAACGCGTACGTTTCGACGAGTACGACGCCCGGG"
                            "AAAACGCGTGTCGACGACGAACGCGTACGTTTACGACGAGTACGACGCCC"
                            "AAAACGCGTGTCACGACAACGCGTACGTTTCGACGAGTACGACGCCCGGG"
                            "AAAACGCGTGTCGACGACGAACGCGTACGTTTACGACGAGTACGACGCCC"
                            "AAAACGCGTGTCACGACAACGCGTACGTTTCGACGAGTACGACGCCCGGG"
                            "AAAACGCGTGTCGACGACGAACGCGTACGTTTACGACGAGTACGACGCCC"
                            "AAAACGCGTGTCACGACAACGCGTACGTTTCGACGAGTACGACGCCCGGG"
                            "AAAACGCGTGTCGACGACGAACGCGTACGTTTACGACGAGTACGACGCCC"
                            "AAAACGCGTGTCACGACAACGCGTACGTTTCGACGAGTACGACGCCCGGG"
                            "AAAACGCGTGTCGACGACGAACGCGTACGTTTACGACGAGTACGACGCCC"
                            "AAAACGCGTGTCACGACAACGCGTACGTTTCGACGAGTACGACGCCCGGG"
                            "AAAACGCGTGTCGACGACGAACGCGTACGTTTACGACGAGTACGACGCCC";

void test_perf_viterbi(void)
{
    unsigned const               ncore_nodes = 1000;
    imm_float const              epsilon = (imm_float)0.01;
    struct nmm_base_abc const*   base = nmm_base_abc_create("ACGT", 'X');
    struct imm_abc const*        abc = nmm_base_abc_super(base);
    struct nmm_base_lprob const* basep =
        nmm_base_lprob_create(base, imm_log(0.25), imm_log(0.25), imm_log(0.45), imm_log(0.05));

    struct nmm_codon_marg const* M_codont = create_codont(base, NMM_TRIPLET('A', 'C', 'G'), imm_log(100.0));
    struct nmm_codon_marg const* I_codont = create_codont(base, NMM_TRIPLET('C', 'G', 'T'), imm_log(100.0));
    struct nmm_codon_marg const* B_codont = create_codont(base, NMM_TRIPLET('A', 'A', 'A'), imm_log(100.0));
    struct nmm_codon_marg const* E_codont = create_codont(base, NMM_TRIPLET('C', 'C', 'C'), imm_log(100.0));
    struct nmm_codon_marg const* J_codont = create_codont(base, NMM_TRIPLET('G', 'G', 'G'), imm_log(100.0));

    struct imm_hmm* hmm = imm_hmm_create(abc);

    struct imm_mute_state const* start = imm_mute_state_create("START", abc);
    imm_hmm_add_state(hmm, mute_super(start), imm_log(1.0));

    struct imm_mute_state const* end = imm_mute_state_create("END", abc);
    imm_hmm_add_state(hmm, mute_super(end), zero());

    struct nmm_frame_state const* B = nmm_frame_state_create("B", basep, B_codont, epsilon);
    imm_hmm_add_state(hmm, frame_super(B), zero());

    struct nmm_frame_state const* E = nmm_frame_state_create("E", basep, B_codont, epsilon);
    imm_hmm_add_state(hmm, frame_super(E), zero());

    struct nmm_frame_state const* J = nmm_frame_state_create("J", basep, B_codont, epsilon);
    imm_hmm_add_state(hmm, frame_super(J), zero());

    imm_hmm_set_trans(hmm, mute_super(start), frame_super(B), imm_log(0.2));
    imm_hmm_set_trans(hmm, frame_super(B), frame_super(B), imm_log(0.2));
    imm_hmm_set_trans(hmm, frame_super(E), frame_super(E), imm_log(0.2));
    imm_hmm_set_trans(hmm, frame_super(J), frame_super(J), imm_log(0.2));
    imm_hmm_set_trans(hmm, frame_super(E), frame_super(J), imm_log(0.2));
    imm_hmm_set_trans(hmm, frame_super(J), frame_super(B), imm_log(0.2));
    imm_hmm_set_trans(hmm, frame_super(E), mute_super(end), imm_log(0.2));

    struct nmm_frame_state const** M = malloc(sizeof(struct nmm_frame_state*) * ncore_nodes);
    struct nmm_frame_state const** I = malloc(sizeof(struct nmm_frame_state*) * ncore_nodes);
    struct imm_mute_state const**  D = malloc(sizeof(struct nmm_frame_state*) * ncore_nodes);

    char name[10] = "\0";
    for (unsigned i = 0; i < ncore_nodes; ++i) {
        M[i] = nmm_frame_state_create(fmt_name(name, "M", i), basep, M_codont, epsilon);
        I[i] = nmm_frame_state_create(fmt_name(name, "I", i), basep, I_codont, epsilon);
        D[i] = imm_mute_state_create(fmt_name(name, "D", i), abc);

        imm_hmm_add_state(hmm, frame_super(M[i]), zero());
        imm_hmm_add_state(hmm, frame_super(I[i]), zero());
        imm_hmm_add_state(hmm, mute_super(D[i]), zero());

        if (i == 0)
            imm_hmm_set_trans(hmm, frame_super(B), frame_super(M[0]), imm_log(0.2));

        imm_hmm_set_trans(hmm, frame_super(M[i]), frame_super(I[i]), imm_log(0.2));
        imm_hmm_set_trans(hmm, frame_super(I[i]), frame_super(I[i]), imm_log(0.2));

        if (i > 0) {
            imm_hmm_set_trans(hmm, frame_super(M[i - 1]), frame_super(M[i]), imm_log(0.2));
            imm_hmm_set_trans(hmm, mute_super(D[i - 1]), frame_super(M[i]), imm_log(0.2));
            imm_hmm_set_trans(hmm, frame_super(I[i - 1]), frame_super(M[i]), imm_log(0.2));

            imm_hmm_set_trans(hmm, frame_super(M[i - 1]), mute_super(D[i]), imm_log(0.2));
            imm_hmm_set_trans(hmm, mute_super(D[i - 1]), mute_super(D[i]), imm_log(0.2));
        }

        if (i == ncore_nodes - 1) {
            imm_hmm_set_trans(hmm, frame_super(M[i]), frame_super(E), imm_log(0.2));
            imm_hmm_set_trans(hmm, mute_super(D[i]), frame_super(E), imm_log(0.2));
            imm_hmm_set_trans(hmm, frame_super(I[i]), frame_super(E), imm_log(0.2));
        }
    }

    struct imm_seq const* seq = imm_seq_create(__seq, abc);

    cass_cond(imm_seq_length(seq) == 2000);

    struct imm_dp const* dp = imm_hmm_create_dp(hmm, mute_super(end));

    struct imm_dp_task* task = imm_dp_task_create(dp);
    imm_dp_task_setup(task, seq, 0);
    struct imm_results const* results = imm_dp_viterbi(dp, task);
    imm_dp_task_destroy(task);

    cass_cond(imm_results_size(results) == 1);
    struct imm_result const* r = imm_results_get(results, 0);
    struct imm_subseq        subseq = imm_result_subseq(r);
    struct imm_path const*   path = imm_result_path(r);
    imm_float                loglik = imm_hmm_loglikelihood(hmm, imm_subseq_cast(&subseq), path);
    cass_cond(imm_lprob_is_valid(loglik));
    cass_cond(!imm_lprob_is_zero(loglik));
    cass_close(loglik, -1641.970511421383435);
    imm_results_destroy(results);
    imm_seq_destroy(seq);

    struct nmm_output* output = nmm_output_create(TMPDIR "/perf.nmm");
    cass_cond(output != NULL);
    struct nmm_profile* m = nmm_profile_create(abc);
    nmm_profile_append_model(m, imm_model_create(hmm, dp));
    cass_equal(nmm_output_write(output, m), 0);
    nmm_profile_destroy(m, false);
    cass_equal(nmm_output_destroy(output), 0);

    imm_hmm_destroy(hmm);
    imm_mute_state_destroy(start);
    imm_mute_state_destroy(end);
    nmm_frame_state_destroy(B);
    nmm_frame_state_destroy(E);
    nmm_frame_state_destroy(J);
    for (unsigned i = 0; i < ncore_nodes; ++i) {
        nmm_frame_state_destroy(M[i]);
        nmm_frame_state_destroy(I[i]);
        imm_mute_state_destroy(D[i]);
    }
    nmm_codon_marg_destroy(M_codont);
    nmm_codon_marg_destroy(I_codont);
    nmm_codon_marg_destroy(B_codont);
    nmm_codon_marg_destroy(E_codont);
    nmm_codon_marg_destroy(J_codont);
    nmm_base_lprob_destroy(basep);
    nmm_base_abc_destroy(base);
    imm_dp_destroy(dp);
    free(M);
    free(I);
    free(D);

    struct nmm_input* input = nmm_input_create(TMPDIR "/perf.nmm");
    cass_cond(input != NULL);
    cass_cond(!nmm_input_eof(input));
    struct nmm_profile const* prof = nmm_input_read(input);
    cass_cond(!nmm_input_eof(input));
    cass_cond(prof != NULL);
    nmm_input_destroy(input);

    struct imm_model* model = nmm_profile_get_model(prof, 0);

    cass_equal(imm_model_nstates(model), 3 * ncore_nodes + 5);

    abc = nmm_profile_abc(prof);
    hmm = imm_model_hmm(model);
    dp = imm_model_dp(model);

    seq = imm_seq_create(__seq, abc);
    task = imm_dp_task_create(dp);
    imm_dp_task_setup(task, seq, 0);
    results = imm_dp_viterbi(dp, task);
    imm_dp_task_destroy(task);
    cass_cond(imm_results_size(results) == 1);
    r = imm_results_get(results, 0);
    subseq = imm_result_subseq(r);
    path = imm_result_path(r);
    loglik = imm_hmm_loglikelihood(hmm, imm_subseq_cast(&subseq), path);
    cass_cond(imm_lprob_is_valid(loglik));
    cass_cond(!imm_lprob_is_zero(loglik));
    cass_close(loglik, -1641.970511421383435);
    imm_results_destroy(results);
    imm_seq_destroy(seq);

    for (uint16_t i = 0; i < imm_model_nstates(model); ++i)
        imm_state_destroy(imm_model_state(model, i));

    for (uint16_t i = 0; i < nmm_profile_nbase_lprobs(prof); ++i)
        nmm_base_lprob_destroy(nmm_profile_base_lprob(prof, i));

    for (uint16_t i = 0; i < nmm_profile_ncodon_margs(prof); ++i)
        nmm_codon_marg_destroy(nmm_profile_codon_marg(prof, i));

    for (uint16_t i = 0; i < nmm_profile_ncodon_lprobs(prof); ++i)
        nmm_codon_lprob_destroy(nmm_profile_codon_lprob(prof, i));

    imm_abc_destroy(abc);
    imm_hmm_destroy(hmm);
    imm_dp_destroy(dp);
    nmm_profile_destroy(prof, false);
}

static struct nmm_codon_lprob* create_codonp(struct nmm_base_abc const* base)
{
    struct nmm_codon_lprob* codonp = nmm_codon_lprob_create(base);
    struct nmm_codon*       codon = nmm_codon_create(base);

    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'A', 'A')) == 0);
    nmm_codon_lprob_set(codonp, codon, imm_log(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'A', 'C')) == 0);
    nmm_codon_lprob_set(codonp, codon, imm_log(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'A', 'G')) == 0);
    nmm_codon_lprob_set(codonp, codon, imm_log(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'A', 'T')) == 0);
    nmm_codon_lprob_set(codonp, codon, imm_log(0.0023173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'C', 'A')) == 0);
    nmm_codon_lprob_set(codonp, codon, imm_log(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'C', 'G')) == 0);
    nmm_codon_lprob_set(codonp, codon, imm_log(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'C', 'T')) == 0);
    nmm_codon_lprob_set(codonp, codon, imm_log(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'G', 'C')) == 0);
    nmm_codon_lprob_set(codonp, codon, imm_log(0.0029114));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'G', 'G')) == 0);
    nmm_codon_lprob_set(codonp, codon, imm_log(0.0029003));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'G', 'T')) == 0);
    nmm_codon_lprob_set(codonp, codon, imm_log(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'T', 'A')) == 0);
    nmm_codon_lprob_set(codonp, codon, imm_log(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'T', 'G')) == 0);
    nmm_codon_lprob_set(codonp, codon, imm_log(0.0049178));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'T', 'T')) == 0);
    nmm_codon_lprob_set(codonp, codon, imm_log(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('C', 'A', 'A')) == 0);
    nmm_codon_lprob_set(codonp, codon, imm_log(0.0029478));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('C', 'A', 'C')) == 0);
    nmm_codon_lprob_set(codonp, codon, imm_log(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('C', 'A', 'G')) == 0);
    nmm_codon_lprob_set(codonp, codon, imm_log(0.0029123));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('C', 'A', 'T')) == 0);
    nmm_codon_lprob_set(codonp, codon, imm_log(0.0009133));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('C', 'C', 'A')) == 0);
    nmm_codon_lprob_set(codonp, codon, imm_log(0.0029179));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('C', 'C', 'G')) == 0);
    nmm_codon_lprob_set(codonp, codon, imm_log(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('C', 'C', 'T')) == 0);
    nmm_codon_lprob_set(codonp, codon, imm_log(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('C', 'G', 'C')) == 0);
    nmm_codon_lprob_set(codonp, codon, imm_log(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('C', 'G', 'G')) == 0);
    nmm_codon_lprob_set(codonp, codon, imm_log(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('C', 'G', 'T')) == 0);
    nmm_codon_lprob_set(codonp, codon, imm_log(0.0041183));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('C', 'T', 'A')) == 0);
    nmm_codon_lprob_set(codonp, codon, imm_log(0.0038173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('C', 'T', 'G')) == 0);
    nmm_codon_lprob_set(codonp, codon, imm_log(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('C', 'T', 'T')) == 0);
    nmm_codon_lprob_set(codonp, codon, imm_log(0.0029111));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('G', 'A', 'A')) == 0);
    nmm_codon_lprob_set(codonp, codon, imm_log(0.0019173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('G', 'A', 'C')) == 0);
    nmm_codon_lprob_set(codonp, codon, imm_log(0.0029103));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('G', 'A', 'G')) == 0);
    nmm_codon_lprob_set(codonp, codon, imm_log(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('G', 'A', 'T')) == 0);
    nmm_codon_lprob_set(codonp, codon, imm_log(0.0029103));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('G', 'C', 'A')) == 0);
    nmm_codon_lprob_set(codonp, codon, imm_log(0.0003138));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('G', 'C', 'G')) == 0);
    nmm_codon_lprob_set(codonp, codon, imm_log(0.0039173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('G', 'C', 'T')) == 0);
    nmm_codon_lprob_set(codonp, codon, imm_log(0.0029103));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('G', 'G', 'C')) == 0);
    nmm_codon_lprob_set(codonp, codon, imm_log(0.0019173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('G', 'G', 'G')) == 0);
    nmm_codon_lprob_set(codonp, codon, imm_log(0.0009173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('G', 'G', 'T')) == 0);
    nmm_codon_lprob_set(codonp, codon, imm_log(0.0029143));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('G', 'T', 'A')) == 0);
    nmm_codon_lprob_set(codonp, codon, imm_log(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('G', 'T', 'G')) == 0);
    nmm_codon_lprob_set(codonp, codon, imm_log(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('G', 'T', 'T')) == 0);
    nmm_codon_lprob_set(codonp, codon, imm_log(0.0029171));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('T', 'A', 'A')) == 0);
    nmm_codon_lprob_set(codonp, codon, imm_log(0.0099113));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('T', 'A', 'C')) == 0);
    nmm_codon_lprob_set(codonp, codon, imm_log(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('T', 'A', 'G')) == 0);
    nmm_codon_lprob_set(codonp, codon, imm_log(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('T', 'A', 'T')) == 0);
    nmm_codon_lprob_set(codonp, codon, imm_log(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('T', 'C', 'A')) == 0);
    nmm_codon_lprob_set(codonp, codon, imm_log(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('T', 'C', 'G')) == 0);
    nmm_codon_lprob_set(codonp, codon, imm_log(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('T', 'C', 'T')) == 0);
    nmm_codon_lprob_set(codonp, codon, imm_log(0.0020173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('T', 'G', 'C')) == 0);
    nmm_codon_lprob_set(codonp, codon, imm_log(0.0029193));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('T', 'G', 'G')) == 0);
    nmm_codon_lprob_set(codonp, codon, imm_log(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('T', 'G', 'T')) == 0);
    nmm_codon_lprob_set(codonp, codon, imm_log(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('T', 'T', 'A')) == 0);
    nmm_codon_lprob_set(codonp, codon, imm_log(0.0089193));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('T', 'T', 'G')) == 0);
    nmm_codon_lprob_set(codonp, codon, imm_log(0.0029138));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('T', 'T', 'T')) == 0);
    nmm_codon_lprob_set(codonp, codon, imm_log(0.0089183));

    nmm_codon_destroy(codon);
    return codonp;
}

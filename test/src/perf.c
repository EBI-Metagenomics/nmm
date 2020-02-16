#include "cass.h"
#include "elapsed.h"
#include "imm/imm.h"
#include "nmm/nmm.h"

void test_perf_viterbi(void);

int main(void)
{
    test_perf_viterbi();
    return cass_status();
}

static inline double                  zero(void) { return imm_lprob_zero(); }
static inline int                     is_valid(double a) { return imm_lprob_is_valid(a); }
static inline int                     is_zero(double a) { return imm_lprob_is_zero(a); }
static inline struct imm_state const* cast_c(void const* s) { return imm_state_cast_c(s); }
static inline char*                   fmt_name(char* restrict buffer, char const* name, int i)
{
    sprintf(buffer, "%s%d", name, i);
    return buffer;
}

struct nmm_codon_lprob* create_codonp(struct nmm_base const* base);

static struct nmm_codon_table const* create_codont(struct nmm_base const*   base,
                                                   struct nmm_triplet const triplet,
                                                   double const             lprob)
{
    struct nmm_codon_lprob* codonp = create_codonp(base);
    struct nmm_codon*       codon = nmm_codon_create(base);
    cass_cond(nmm_codon_set_triplet(codon, triplet) == 0);
    nmm_codon_lprob_set_lprob(codonp, codon, lprob);
    struct nmm_codon_table const* codont = nmm_codon_table_create(codonp);
    nmm_codon_destroy(codon);
    nmm_codon_lprob_destroy(codonp);
    return codont;
}

void test_perf_viterbi(void)
{
    unsigned const          ncore_nodes = 1000;
    double const            epsilon = 0.01;
    struct imm_abc const*   abc = imm_abc_create("ACGT", 'X');
    struct nmm_base const*  base = nmm_base_create(abc);
    struct nmm_baset const* baset =
        nmm_baset_create(base, log(0.25), log(0.25), log(0.45), log(0.05));

    struct nmm_codon_table const* M_codont =
        create_codont(base, NMM_TRIPLET('A', 'C', 'G'), log(100.0));
    struct nmm_codon_table const* I_codont =
        create_codont(base, NMM_TRIPLET('C', 'G', 'T'), log(100.0));
    struct nmm_codon_table const* B_codont =
        create_codont(base, NMM_TRIPLET('A', 'A', 'A'), log(100.0));
    struct nmm_codon_table const* E_codont =
        create_codont(base, NMM_TRIPLET('C', 'C', 'C'), log(100.0));
    struct nmm_codon_table const* J_codont =
        create_codont(base, NMM_TRIPLET('G', 'G', 'G'), log(100.0));

    struct imm_hmm* hmm = imm_hmm_create(abc);

    struct imm_mute_state const* start = imm_mute_state_create("START", abc);
    imm_hmm_add_state(hmm, cast_c(start), log(1.0));

    struct imm_mute_state const* end = imm_mute_state_create("END", abc);
    imm_hmm_add_state(hmm, cast_c(end), zero());

    struct nmm_frame_state const* B = nmm_frame_state_create("B", baset, B_codont, epsilon);
    imm_hmm_add_state(hmm, cast_c(B), zero());

    struct nmm_frame_state const* E = nmm_frame_state_create("E", baset, B_codont, epsilon);
    imm_hmm_add_state(hmm, cast_c(E), zero());

    struct nmm_frame_state const* J = nmm_frame_state_create("J", baset, B_codont, epsilon);
    imm_hmm_add_state(hmm, cast_c(J), zero());

    imm_hmm_set_trans(hmm, cast_c(start), cast_c(B), log(0.2));
    imm_hmm_set_trans(hmm, cast_c(B), cast_c(B), log(0.2));
    imm_hmm_set_trans(hmm, cast_c(E), cast_c(E), log(0.2));
    imm_hmm_set_trans(hmm, cast_c(J), cast_c(J), log(0.2));
    imm_hmm_set_trans(hmm, cast_c(E), cast_c(J), log(0.2));
    imm_hmm_set_trans(hmm, cast_c(J), cast_c(B), log(0.2));
    imm_hmm_set_trans(hmm, cast_c(E), cast_c(end), log(0.2));

    struct nmm_frame_state const* M[ncore_nodes];
    struct nmm_frame_state const* I[ncore_nodes];
    struct imm_mute_state const*  D[ncore_nodes];

    char name[10] = "\0";
    for (unsigned i = 0; i < ncore_nodes; ++i) {
        M[i] = nmm_frame_state_create(fmt_name(name, "M", i), baset, M_codont, epsilon);
        I[i] = nmm_frame_state_create(fmt_name(name, "I", i), baset, I_codont, epsilon);
        D[i] = imm_mute_state_create(fmt_name(name, "D", i), abc);

        imm_hmm_add_state(hmm, cast_c(M[i]), zero());
        imm_hmm_add_state(hmm, cast_c(I[i]), zero());
        imm_hmm_add_state(hmm, cast_c(D[i]), zero());

        if (i == 0)
            imm_hmm_set_trans(hmm, cast_c(B), cast_c(M[0]), log(0.2));

        imm_hmm_set_trans(hmm, cast_c(M[i]), cast_c(I[i]), log(0.2));
        imm_hmm_set_trans(hmm, cast_c(I[i]), cast_c(I[i]), log(0.2));

        if (i > 0) {
            imm_hmm_set_trans(hmm, cast_c(M[i - 1]), cast_c(M[i]), log(0.2));
            imm_hmm_set_trans(hmm, cast_c(D[i - 1]), cast_c(M[i]), log(0.2));
            imm_hmm_set_trans(hmm, cast_c(I[i - 1]), cast_c(M[i]), log(0.2));

            imm_hmm_set_trans(hmm, cast_c(M[i - 1]), cast_c(D[i]), log(0.2));
            imm_hmm_set_trans(hmm, cast_c(D[i - 1]), cast_c(D[i]), log(0.2));
        }

        if (i == ncore_nodes - 1) {
            imm_hmm_set_trans(hmm, cast_c(M[i]), cast_c(E), log(0.2));
            imm_hmm_set_trans(hmm, cast_c(D[i]), cast_c(E), log(0.2));
            imm_hmm_set_trans(hmm, cast_c(I[i]), cast_c(E), log(0.2));
        }
    }

    struct imm_seq const* seq =
        imm_seq_create("AAAACGCGTGTCACGACAACGCGTACGTTTCGACGAGTACGACGCCCGGG"
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
                       "AAAACGCGTGTCGACGACGAACGCGTACGTTTACGACGAGTACGACGCCC",
                       abc);

    cass_cond(imm_seq_length(seq) == 2000);

    struct elapsed* elapsed = elapsed_create();
    elapsed_start(elapsed);
    struct imm_results const* results = imm_hmm_viterbi(hmm, seq, cast_c(end), 0);
    elapsed_end(elapsed);
    cass_cond(imm_results_size(results) == 1);
    struct imm_result const* r = imm_results_get(results, 0);
    struct imm_path const*   path = imm_result_path(r);
    double                   score = imm_result_loglik(r);
    cass_cond(is_valid(score) && !is_zero(score));
    cass_close(score, -1641.970511421383435);
    imm_results_destroy(results);
    imm_seq_destroy(seq);

#ifdef NDEBUG
    cass_cond(elapsed_seconds(elapsed) < 20.0);
#endif

    elapsed_destroy(elapsed);
    imm_hmm_destroy(hmm);
    imm_mute_state_destroy(start);
    imm_mute_state_destroy(end);
    nmm_frame_state_destroy(B);
    nmm_frame_state_destroy(E);
    nmm_frame_state_destroy(J);
    for (int i = 0; i < ncore_nodes; ++i) {
        nmm_frame_state_destroy(M[i]);
        nmm_frame_state_destroy(I[i]);
        imm_mute_state_destroy(D[i]);
    }
    nmm_codon_table_destroy(M_codont);
    nmm_codon_table_destroy(I_codont);
    nmm_codon_table_destroy(B_codont);
    nmm_codon_table_destroy(E_codont);
    nmm_codon_table_destroy(J_codont);
    nmm_baset_destroy(baset);
    nmm_base_destroy(base);
    imm_abc_destroy(abc);
}

struct nmm_codon_lprob* create_codonp(struct nmm_base const* base)
{
    struct nmm_codon_lprob* codonp = nmm_codon_lprob_create(base);
    struct nmm_codon*       codon = nmm_codon_create(base);

    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'A', 'A')) == 0);
    nmm_codon_lprob_set_lprob(codonp, codon, log(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'A', 'C')) == 0);
    nmm_codon_lprob_set_lprob(codonp, codon, log(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'A', 'G')) == 0);
    nmm_codon_lprob_set_lprob(codonp, codon, log(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'A', 'T')) == 0);
    nmm_codon_lprob_set_lprob(codonp, codon, log(0.0023173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'C', 'A')) == 0);
    nmm_codon_lprob_set_lprob(codonp, codon, log(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'C', 'G')) == 0);
    nmm_codon_lprob_set_lprob(codonp, codon, log(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'C', 'T')) == 0);
    nmm_codon_lprob_set_lprob(codonp, codon, log(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'G', 'C')) == 0);
    nmm_codon_lprob_set_lprob(codonp, codon, log(0.0029114));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'G', 'G')) == 0);
    nmm_codon_lprob_set_lprob(codonp, codon, log(0.0029003));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'G', 'T')) == 0);
    nmm_codon_lprob_set_lprob(codonp, codon, log(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'T', 'A')) == 0);
    nmm_codon_lprob_set_lprob(codonp, codon, log(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'T', 'G')) == 0);
    nmm_codon_lprob_set_lprob(codonp, codon, log(0.0049178));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'T', 'T')) == 0);
    nmm_codon_lprob_set_lprob(codonp, codon, log(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('C', 'A', 'A')) == 0);
    nmm_codon_lprob_set_lprob(codonp, codon, log(0.0029478));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('C', 'A', 'C')) == 0);
    nmm_codon_lprob_set_lprob(codonp, codon, log(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('C', 'A', 'G')) == 0);
    nmm_codon_lprob_set_lprob(codonp, codon, log(0.0029123));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('C', 'A', 'T')) == 0);
    nmm_codon_lprob_set_lprob(codonp, codon, log(0.0009133));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('C', 'C', 'A')) == 0);
    nmm_codon_lprob_set_lprob(codonp, codon, log(0.0029179));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('C', 'C', 'G')) == 0);
    nmm_codon_lprob_set_lprob(codonp, codon, log(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('C', 'C', 'T')) == 0);
    nmm_codon_lprob_set_lprob(codonp, codon, log(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('C', 'G', 'C')) == 0);
    nmm_codon_lprob_set_lprob(codonp, codon, log(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('C', 'G', 'G')) == 0);
    nmm_codon_lprob_set_lprob(codonp, codon, log(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('C', 'G', 'T')) == 0);
    nmm_codon_lprob_set_lprob(codonp, codon, log(0.0041183));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('C', 'T', 'A')) == 0);
    nmm_codon_lprob_set_lprob(codonp, codon, log(0.0038173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('C', 'T', 'G')) == 0);
    nmm_codon_lprob_set_lprob(codonp, codon, log(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('C', 'T', 'T')) == 0);
    nmm_codon_lprob_set_lprob(codonp, codon, log(0.0029111));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('G', 'A', 'A')) == 0);
    nmm_codon_lprob_set_lprob(codonp, codon, log(0.0019173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('G', 'A', 'C')) == 0);
    nmm_codon_lprob_set_lprob(codonp, codon, log(0.0029103));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('G', 'A', 'G')) == 0);
    nmm_codon_lprob_set_lprob(codonp, codon, log(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('G', 'A', 'T')) == 0);
    nmm_codon_lprob_set_lprob(codonp, codon, log(0.0029103));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('G', 'C', 'A')) == 0);
    nmm_codon_lprob_set_lprob(codonp, codon, log(0.0003138));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('G', 'C', 'G')) == 0);
    nmm_codon_lprob_set_lprob(codonp, codon, log(0.0039173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('G', 'C', 'T')) == 0);
    nmm_codon_lprob_set_lprob(codonp, codon, log(0.0029103));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('G', 'G', 'C')) == 0);
    nmm_codon_lprob_set_lprob(codonp, codon, log(0.0019173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('G', 'G', 'G')) == 0);
    nmm_codon_lprob_set_lprob(codonp, codon, log(0.0009173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('G', 'G', 'T')) == 0);
    nmm_codon_lprob_set_lprob(codonp, codon, log(0.0029143));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('G', 'T', 'A')) == 0);
    nmm_codon_lprob_set_lprob(codonp, codon, log(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('G', 'T', 'G')) == 0);
    nmm_codon_lprob_set_lprob(codonp, codon, log(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('G', 'T', 'T')) == 0);
    nmm_codon_lprob_set_lprob(codonp, codon, log(0.0029171));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('T', 'A', 'A')) == 0);
    nmm_codon_lprob_set_lprob(codonp, codon, log(0.0099113));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('T', 'A', 'C')) == 0);
    nmm_codon_lprob_set_lprob(codonp, codon, log(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('T', 'A', 'G')) == 0);
    nmm_codon_lprob_set_lprob(codonp, codon, log(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('T', 'A', 'T')) == 0);
    nmm_codon_lprob_set_lprob(codonp, codon, log(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('T', 'C', 'A')) == 0);
    nmm_codon_lprob_set_lprob(codonp, codon, log(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('T', 'C', 'G')) == 0);
    nmm_codon_lprob_set_lprob(codonp, codon, log(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('T', 'C', 'T')) == 0);
    nmm_codon_lprob_set_lprob(codonp, codon, log(0.0020173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('T', 'G', 'C')) == 0);
    nmm_codon_lprob_set_lprob(codonp, codon, log(0.0029193));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('T', 'G', 'G')) == 0);
    nmm_codon_lprob_set_lprob(codonp, codon, log(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('T', 'G', 'T')) == 0);
    nmm_codon_lprob_set_lprob(codonp, codon, log(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('T', 'T', 'A')) == 0);
    nmm_codon_lprob_set_lprob(codonp, codon, log(0.0089193));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('T', 'T', 'G')) == 0);
    nmm_codon_lprob_set_lprob(codonp, codon, log(0.0029138));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('T', 'T', 'T')) == 0);
    nmm_codon_lprob_set_lprob(codonp, codon, log(0.0089183));

    nmm_codon_destroy(codon);
    return codonp;
}

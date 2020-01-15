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

void set_codont(struct nmm_codont* codont);

void test_perf_viterbi(void)
{
    int const       ncore_nodes = 1000;
    double const    epsilon = 0.01;
    struct imm_abc* abc = imm_abc_create("ACGT", 'X');

    struct nmm_baset* baset = nmm_baset_create(abc);
    nmm_baset_set_lprob(baset, 'A', log(0.25));
    nmm_baset_set_lprob(baset, 'C', log(0.25));
    nmm_baset_set_lprob(baset, 'G', log(0.45));
    nmm_baset_set_lprob(baset, 'T', log(0.05));

    struct nmm_codont* M_codont = nmm_codont_create(abc);
    set_codont(M_codont);
    nmm_codont_set_lprob(M_codont, &NMM_CODON('A', 'C', 'G'), log(100.0));

    struct nmm_codont* I_codont = nmm_codont_create(abc);
    set_codont(I_codont);
    nmm_codont_set_lprob(I_codont, &NMM_CODON('C', 'G', 'T'), log(100.0));

    struct nmm_codont* B_codont = nmm_codont_create(abc);
    set_codont(B_codont);
    nmm_codont_set_lprob(B_codont, &NMM_CODON('A', 'A', 'A'), log(100.0));

    struct nmm_codont* E_codont = nmm_codont_create(abc);
    set_codont(E_codont);
    nmm_codont_set_lprob(E_codont, &NMM_CODON('C', 'C', 'C'), log(100.0));

    struct nmm_codont* J_codont = nmm_codont_create(abc);
    set_codont(J_codont);
    nmm_codont_set_lprob(J_codont, &NMM_CODON('G', 'G', 'G'), log(100.0));

    struct imm_hmm* hmm = imm_hmm_create(abc);

    struct imm_mute_state* start = imm_mute_state_create("START", abc);
    imm_hmm_add_state(hmm, cast_c(start), log(1.0));

    struct imm_mute_state* end = imm_mute_state_create("END", abc);
    imm_hmm_add_state(hmm, cast_c(end), zero());

    struct nmm_frame_state* B = nmm_frame_state_create("B", baset, B_codont, epsilon);
    imm_hmm_add_state(hmm, cast_c(B), zero());

    struct nmm_frame_state* E = nmm_frame_state_create("E", baset, B_codont, epsilon);
    imm_hmm_add_state(hmm, cast_c(E), zero());

    struct nmm_frame_state* J = nmm_frame_state_create("J", baset, B_codont, epsilon);
    imm_hmm_add_state(hmm, cast_c(J), zero());

    imm_hmm_set_trans(hmm, cast_c(start), cast_c(B), log(0.2));
    imm_hmm_set_trans(hmm, cast_c(B), cast_c(B), log(0.2));
    imm_hmm_set_trans(hmm, cast_c(E), cast_c(E), log(0.2));
    imm_hmm_set_trans(hmm, cast_c(J), cast_c(J), log(0.2));
    imm_hmm_set_trans(hmm, cast_c(E), cast_c(J), log(0.2));
    imm_hmm_set_trans(hmm, cast_c(J), cast_c(B), log(0.2));
    imm_hmm_set_trans(hmm, cast_c(E), cast_c(end), log(0.2));

    struct nmm_frame_state* M[ncore_nodes];
    struct nmm_frame_state* I[ncore_nodes];
    struct imm_mute_state*  D[ncore_nodes];

    char name[10] = "\0";
    for (int i = 0; i < ncore_nodes; ++i) {
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

    struct elapsed* elapsed = elapsed_create();
    struct imm_path* path = imm_path_create();

    char const seq[] = "AAAACGCGTGTCACGACAACGCGTACGTTTCGACGAGTACGACGCCCGGG"
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

    cass_cond(strlen(seq) == 2000);

    elapsed_start(elapsed);
    double score = imm_hmm_viterbi(hmm, seq, cast_c(end), path);
    cass_cond(is_valid(score) && !is_zero(score));
    cass_close(score, -1641.970511421383435);
    elapsed_end(elapsed);
    imm_path_destroy(path);

    cass_cond(elapsed_seconds(elapsed) < 20.0);
    /* Elapsed: 8.93 seconds */
    printf("Elapsed: %.10f seconds\n", elapsed_seconds(elapsed));

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
    nmm_codont_destroy(M_codont);
    nmm_codont_destroy(I_codont);
    nmm_codont_destroy(B_codont);
    nmm_codont_destroy(E_codont);
    nmm_codont_destroy(J_codont);
    nmm_baset_destroy(baset);
    imm_abc_destroy(abc);
}

void set_codont(struct nmm_codont* codont)
{
    nmm_codont_set_lprob(codont, &NMM_CODON('A', 'A', 'A'), log(0.0029173));
    nmm_codont_set_lprob(codont, &NMM_CODON('A', 'A', 'C'), log(0.0029173));
    nmm_codont_set_lprob(codont, &NMM_CODON('A', 'A', 'G'), log(0.0029173));
    nmm_codont_set_lprob(codont, &NMM_CODON('A', 'A', 'T'), log(0.0023173));
    nmm_codont_set_lprob(codont, &NMM_CODON('A', 'C', 'A'), log(0.0029173));
    nmm_codont_set_lprob(codont, &NMM_CODON('A', 'C', 'G'), log(0.0029173));
    nmm_codont_set_lprob(codont, &NMM_CODON('A', 'C', 'T'), log(0.0029173));
    nmm_codont_set_lprob(codont, &NMM_CODON('A', 'G', 'C'), log(0.0029114));
    nmm_codont_set_lprob(codont, &NMM_CODON('A', 'G', 'G'), log(0.0029003));
    nmm_codont_set_lprob(codont, &NMM_CODON('A', 'G', 'T'), log(0.0029173));
    nmm_codont_set_lprob(codont, &NMM_CODON('A', 'T', 'A'), log(0.0029173));
    nmm_codont_set_lprob(codont, &NMM_CODON('A', 'T', 'G'), log(0.0049178));
    nmm_codont_set_lprob(codont, &NMM_CODON('A', 'T', 'T'), log(0.0029173));
    nmm_codont_set_lprob(codont, &NMM_CODON('C', 'A', 'A'), log(0.0029478));
    nmm_codont_set_lprob(codont, &NMM_CODON('C', 'A', 'C'), log(0.0029173));
    nmm_codont_set_lprob(codont, &NMM_CODON('C', 'A', 'G'), log(0.0029123));
    nmm_codont_set_lprob(codont, &NMM_CODON('C', 'A', 'T'), log(0.0009133));
    nmm_codont_set_lprob(codont, &NMM_CODON('C', 'C', 'A'), log(0.0029179));
    nmm_codont_set_lprob(codont, &NMM_CODON('C', 'C', 'G'), log(0.0029173));
    nmm_codont_set_lprob(codont, &NMM_CODON('C', 'C', 'T'), log(0.0029173));
    nmm_codont_set_lprob(codont, &NMM_CODON('C', 'G', 'C'), log(0.0029173));
    nmm_codont_set_lprob(codont, &NMM_CODON('C', 'G', 'G'), log(0.0029173));
    nmm_codont_set_lprob(codont, &NMM_CODON('C', 'G', 'T'), log(0.0041183));
    nmm_codont_set_lprob(codont, &NMM_CODON('C', 'T', 'A'), log(0.0038173));
    nmm_codont_set_lprob(codont, &NMM_CODON('C', 'T', 'G'), log(0.0029173));
    nmm_codont_set_lprob(codont, &NMM_CODON('C', 'T', 'T'), log(0.0029111));
    nmm_codont_set_lprob(codont, &NMM_CODON('G', 'A', 'A'), log(0.0019173));
    nmm_codont_set_lprob(codont, &NMM_CODON('G', 'A', 'C'), log(0.0029103));
    nmm_codont_set_lprob(codont, &NMM_CODON('G', 'A', 'G'), log(0.0029173));
    nmm_codont_set_lprob(codont, &NMM_CODON('G', 'A', 'T'), log(0.0029103));
    nmm_codont_set_lprob(codont, &NMM_CODON('G', 'C', 'A'), log(0.0003138));
    nmm_codont_set_lprob(codont, &NMM_CODON('G', 'C', 'G'), log(0.0039173));
    nmm_codont_set_lprob(codont, &NMM_CODON('G', 'C', 'T'), log(0.0029103));
    nmm_codont_set_lprob(codont, &NMM_CODON('G', 'G', 'C'), log(0.0019173));
    nmm_codont_set_lprob(codont, &NMM_CODON('G', 'G', 'G'), log(0.0009173));
    nmm_codont_set_lprob(codont, &NMM_CODON('G', 'G', 'T'), log(0.0029143));
    nmm_codont_set_lprob(codont, &NMM_CODON('G', 'T', 'A'), log(0.0029173));
    nmm_codont_set_lprob(codont, &NMM_CODON('G', 'T', 'G'), log(0.0029173));
    nmm_codont_set_lprob(codont, &NMM_CODON('G', 'T', 'T'), log(0.0029171));
    nmm_codont_set_lprob(codont, &NMM_CODON('T', 'A', 'A'), log(0.0099113));
    nmm_codont_set_lprob(codont, &NMM_CODON('T', 'A', 'C'), log(0.0029173));
    nmm_codont_set_lprob(codont, &NMM_CODON('T', 'A', 'G'), log(0.0029173));
    nmm_codont_set_lprob(codont, &NMM_CODON('T', 'A', 'T'), log(0.0029173));
    nmm_codont_set_lprob(codont, &NMM_CODON('T', 'C', 'A'), log(0.0029173));
    nmm_codont_set_lprob(codont, &NMM_CODON('T', 'C', 'G'), log(0.0029173));
    nmm_codont_set_lprob(codont, &NMM_CODON('T', 'C', 'T'), log(0.0020173));
    nmm_codont_set_lprob(codont, &NMM_CODON('T', 'G', 'C'), log(0.0029193));
    nmm_codont_set_lprob(codont, &NMM_CODON('T', 'G', 'G'), log(0.0029173));
    nmm_codont_set_lprob(codont, &NMM_CODON('T', 'G', 'T'), log(0.0029173));
    nmm_codont_set_lprob(codont, &NMM_CODON('T', 'T', 'A'), log(0.0089193));
    nmm_codont_set_lprob(codont, &NMM_CODON('T', 'T', 'G'), log(0.0029138));
    nmm_codont_set_lprob(codont, &NMM_CODON('T', 'T', 'T'), log(0.0089183));
}

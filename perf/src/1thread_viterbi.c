#include "cass/cass.h"
#include "elapsed/elapsed.h"
#include "imm/imm.h"
#include "nmm/nmm.h"
#include "stats.h"

#define CLOSE(a, b) cass_close2(a, b, 1e-6, 0.0)
#define NSAMPLES 100
#define LOG(x) ((imm_float)log(x))

imm_float perf_1thread_viterbi(imm_float* seconds, uint16_t ncore_nodes, uint16_t seq_100length);

static struct nmm_codon_table const* create_codont(struct nmm_base_abc const* base,
                                                   struct nmm_triplet const   triplet,
                                                   imm_float const            lprob);
static struct imm_state const*       frame_super(struct nmm_frame_state const* state);
static struct imm_state const*       mute_super(struct imm_mute_state const* state);

static inline imm_float zero(void) { return imm_lprob_zero(); }
static inline int       is_valid(imm_float a) { return imm_lprob_is_valid(a); }
static inline int       is_zero(imm_float a) { return imm_lprob_is_zero(a); }
static char*            fmt_name(char* restrict buffer, char const* name, unsigned i);

int main(void)
{
    imm_float const logliks[] = {
        (imm_float)-157.9239472248127925,  (imm_float)-469.4680058306055912,
        (imm_float)-930.7839644835630679,  (imm_float)-1395.1589152663593723,
        (imm_float)-1856.4749222363013814, (imm_float)-2320.8498972262464122,
        (imm_float)-2782.1659041962057017, (imm_float)-3246.5408791861568716,
        (imm_float)-3707.8568861561161611, (imm_float)-4172.2318611460850661,
        (imm_float)-4633.5478681161512213, (imm_float)-801.6991121984581241,
        (imm_float)-818.5397430959862959,  (imm_float)-966.6884890100161556,
        (imm_float)-1459.2779880933458116, (imm_float)-1931.6654243218631564,
        (imm_float)-2419.0655189571257324, (imm_float)-2896.6423659701340512,
        (imm_float)-3384.0424606054048127, (imm_float)-3861.6193165051490723,
        (imm_float)-4349.0194111404780415, (imm_float)-4826.5962670403350785,
        (imm_float)-1606.4180684154907794, (imm_float)-1623.2586993130216797,
        (imm_float)-1641.9705114213859360, (imm_float)-1660.6823235297551946,
        (imm_float)-1941.4869290813458065, (imm_float)-2428.8870067061666305,
        (imm_float)-2926.7316511467915916, (imm_float)-3439.8649362890646444,
        (imm_float)-3881.2623126855455666, (imm_float)-4368.6624044120408143,
        (imm_float)-4856.0626420328935637,
    };

    printf("ncore_nodes,seq_length,median,std_err_mean,score,perf_name\n");

    uint16_t ncore_nodes[] = {100, 500, 1000};
    uint16_t seq_100len[] = {1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
    for (unsigned i = 0; i < IMM_ARRAY_SIZE(ncore_nodes); ++i) {
        for (unsigned j = 0; j < IMM_ARRAY_SIZE(seq_100len); ++j) {
            uint16_t  len = seq_100len[j];
            imm_float seconds[NSAMPLES] = {0.};
            imm_float loglik = perf_1thread_viterbi(seconds, ncore_nodes[i], len);
            printf("%.16f, ", loglik);
            cass_close(loglik, logliks[i * IMM_ARRAY_SIZE(seq_100len) + j]);
            struct stats stats = compute_stats(seconds, NSAMPLES);
            char const   fmt[] = "%d,%d,%.6f,%.6f,%.6f,1thread_viterbi\n";
            printf(fmt, ncore_nodes[i], len * 100, stats.median, stats.sem, loglik);
        }
    }

    return 0;
}

imm_float perf_1thread_viterbi(imm_float* seconds, uint16_t ncore_nodes, uint16_t seq_100length)
{
    imm_float const              epsilon = (imm_float)0.01;
    struct nmm_base_abc const*   base = nmm_base_abc_create("ACGT", 'X');
    struct imm_abc const*        abc = nmm_base_abc_super(base);
    struct nmm_base_lprob const* basep =
        nmm_base_lprob_create(base, LOG(0.25), LOG(0.25), LOG(0.45), LOG(0.05));

    struct nmm_codon_table const* M_codont =
        create_codont(base, NMM_TRIPLET('A', 'C', 'G'), LOG(100.0));
    struct nmm_codon_table const* I_codont =
        create_codont(base, NMM_TRIPLET('C', 'G', 'T'), LOG(100.0));
    struct nmm_codon_table const* B_codont =
        create_codont(base, NMM_TRIPLET('A', 'A', 'A'), LOG(100.0));
    struct nmm_codon_table const* E_codont =
        create_codont(base, NMM_TRIPLET('C', 'C', 'C'), LOG(100.0));
    struct nmm_codon_table const* J_codont =
        create_codont(base, NMM_TRIPLET('G', 'G', 'G'), LOG(100.0));

    struct imm_hmm* hmm = imm_hmm_create(abc);

    struct imm_mute_state const* start = imm_mute_state_create("START", abc);
    imm_hmm_add_state(hmm, mute_super(start), LOG(1.0));

    struct imm_mute_state const* end = imm_mute_state_create("END", abc);
    imm_hmm_add_state(hmm, mute_super(end), zero());

    struct nmm_frame_state const* B = nmm_frame_state_create("B", basep, B_codont, epsilon);
    imm_hmm_add_state(hmm, frame_super(B), zero());

    struct nmm_frame_state const* E = nmm_frame_state_create("E", basep, B_codont, epsilon);
    imm_hmm_add_state(hmm, frame_super(E), zero());

    struct nmm_frame_state const* J = nmm_frame_state_create("J", basep, B_codont, epsilon);
    imm_hmm_add_state(hmm, frame_super(J), zero());

    imm_hmm_set_trans(hmm, mute_super(start), frame_super(B), LOG(0.2));
    imm_hmm_set_trans(hmm, frame_super(B), frame_super(B), LOG(0.2));
    imm_hmm_set_trans(hmm, frame_super(E), frame_super(E), LOG(0.2));
    imm_hmm_set_trans(hmm, frame_super(J), frame_super(J), LOG(0.2));
    imm_hmm_set_trans(hmm, frame_super(E), frame_super(J), LOG(0.2));
    imm_hmm_set_trans(hmm, frame_super(J), frame_super(B), LOG(0.2));
    imm_hmm_set_trans(hmm, frame_super(E), mute_super(end), LOG(0.2));

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
            imm_hmm_set_trans(hmm, frame_super(B), frame_super(M[0]), LOG(0.2));

        imm_hmm_set_trans(hmm, frame_super(M[i]), frame_super(I[i]), LOG(0.2));
        imm_hmm_set_trans(hmm, frame_super(I[i]), frame_super(I[i]), LOG(0.2));

        if (i > 0) {
            imm_hmm_set_trans(hmm, frame_super(M[i - 1]), frame_super(M[i]), LOG(0.2));
            imm_hmm_set_trans(hmm, mute_super(D[i - 1]), frame_super(M[i]), LOG(0.2));
            imm_hmm_set_trans(hmm, frame_super(I[i - 1]), frame_super(M[i]), LOG(0.2));

            imm_hmm_set_trans(hmm, frame_super(M[i - 1]), mute_super(D[i]), LOG(0.2));
            imm_hmm_set_trans(hmm, mute_super(D[i - 1]), mute_super(D[i]), LOG(0.2));
        }

        if (i == ncore_nodes - 1) {
            imm_hmm_set_trans(hmm, frame_super(M[i]), frame_super(E), LOG(0.2));
            imm_hmm_set_trans(hmm, mute_super(D[i]), frame_super(E), LOG(0.2));
            imm_hmm_set_trans(hmm, frame_super(I[i]), frame_super(E), LOG(0.2));
        }
    }

    char* str = malloc(sizeof(*str) * (unsigned)(100 * seq_100length + 1));
    str[100 * seq_100length] = '\0';
    char const str0[] = "AAAACGCGTGTCACGACAACGCGTACGTTTCGACGAGTACGACGCCCGGG";
    char const str1[] = "AAAACGCGTGTCGACGACGAACGCGTACGTTTACGACGAGTACGACGCCC";
    for (uint16_t i = 0; i < seq_100length; ++i) {
        for (uint16_t j = 0; j < 50; ++j) {
            str[i * 100 + j] = str0[j];
            str[i * 100 + 50 + j] = str1[j];
        }
    }
    cass_cond(strlen(str) == 100 * seq_100length);

    struct imm_seq const* seq = imm_seq_create(str, abc);
    struct imm_dp const*  dp = imm_hmm_create_dp(hmm, mute_super(end));

    imm_float           loglik = 0.0;
    struct imm_dp_task* task = imm_dp_task_create(dp);
    imm_dp_task_setup(task, seq, 0);
    for (unsigned i = 0; i < NSAMPLES; ++i) {
        struct imm_results const* results = imm_dp_viterbi(dp, task);
        cass_cond(imm_results_size(results) == 1);
        struct imm_result const* r = imm_results_get(results, 0);
        struct imm_subseq        subseq = imm_result_subseq(r);
        struct imm_seq const*    s = imm_subseq_cast(&subseq);
        loglik = imm_hmm_likelihood(hmm, s, imm_result_path(r));
        cass_cond(is_valid(loglik) && !is_zero(loglik));
        seconds[i] = imm_result_seconds(r);
        imm_results_destroy(results);
    }

    imm_dp_task_destroy(task);
    imm_mute_state_destroy(start);
    imm_mute_state_destroy(end);
    nmm_frame_state_destroy(B);
    nmm_frame_state_destroy(E);
    nmm_frame_state_destroy(J);
    for (uint16_t i = 0; i < ncore_nodes; ++i) {
        nmm_frame_state_destroy(M[i]);
        nmm_frame_state_destroy(I[i]);
        imm_mute_state_destroy(D[i]);
    }
    nmm_codon_table_destroy(M_codont);
    nmm_codon_table_destroy(I_codont);
    nmm_codon_table_destroy(B_codont);
    nmm_codon_table_destroy(E_codont);
    nmm_codon_table_destroy(J_codont);
    nmm_base_lprob_destroy(basep);
    nmm_base_abc_destroy(base);
    free(M);
    free(I);
    free(D);

    imm_hmm_destroy(hmm);
    imm_dp_destroy(dp);
    imm_seq_destroy(seq);
    free((void*)str);

    return loglik;
}

static char* fmt_name(char* restrict buffer, char const* name, unsigned i)
{
    sprintf(buffer, "%s%u", name, i);
    return buffer;
}

static struct imm_state const* frame_super(struct nmm_frame_state const* state)
{
    return nmm_frame_state_super(state);
}

static struct imm_state const* mute_super(struct imm_mute_state const* state)
{
    return imm_mute_state_super(state);
}

static struct nmm_codon_lprob* create_codonp(struct nmm_base_abc const* base);

static struct nmm_codon_table const* create_codont(struct nmm_base_abc const* base,
                                                   struct nmm_triplet const   triplet,
                                                   imm_float const            lprob)
{
    struct nmm_codon_lprob* codonp = create_codonp(base);
    struct nmm_codon*       codon = nmm_codon_create(base);
    cass_cond(nmm_codon_set_triplet(codon, triplet) == 0);
    nmm_codon_lprob_set(codonp, codon, lprob);
    struct nmm_codon_table const* codont = nmm_codon_table_create(codonp);
    nmm_codon_destroy(codon);
    nmm_codon_lprob_destroy(codonp);
    return codont;
}

static struct nmm_codon_lprob* create_codonp(struct nmm_base_abc const* base)
{
    struct nmm_codon_lprob* codonp = nmm_codon_lprob_create(base);
    struct nmm_codon*       codon = nmm_codon_create(base);

    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'A', 'A')) == 0);
    nmm_codon_lprob_set(codonp, codon, LOG(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'A', 'C')) == 0);
    nmm_codon_lprob_set(codonp, codon, LOG(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'A', 'G')) == 0);
    nmm_codon_lprob_set(codonp, codon, LOG(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'A', 'T')) == 0);
    nmm_codon_lprob_set(codonp, codon, LOG(0.0023173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'C', 'A')) == 0);
    nmm_codon_lprob_set(codonp, codon, LOG(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'C', 'G')) == 0);
    nmm_codon_lprob_set(codonp, codon, LOG(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'C', 'T')) == 0);
    nmm_codon_lprob_set(codonp, codon, LOG(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'G', 'C')) == 0);
    nmm_codon_lprob_set(codonp, codon, LOG(0.0029114));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'G', 'G')) == 0);
    nmm_codon_lprob_set(codonp, codon, LOG(0.0029003));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'G', 'T')) == 0);
    nmm_codon_lprob_set(codonp, codon, LOG(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'T', 'A')) == 0);
    nmm_codon_lprob_set(codonp, codon, LOG(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'T', 'G')) == 0);
    nmm_codon_lprob_set(codonp, codon, LOG(0.0049178));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'T', 'T')) == 0);
    nmm_codon_lprob_set(codonp, codon, LOG(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('C', 'A', 'A')) == 0);
    nmm_codon_lprob_set(codonp, codon, LOG(0.0029478));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('C', 'A', 'C')) == 0);
    nmm_codon_lprob_set(codonp, codon, LOG(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('C', 'A', 'G')) == 0);
    nmm_codon_lprob_set(codonp, codon, LOG(0.0029123));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('C', 'A', 'T')) == 0);
    nmm_codon_lprob_set(codonp, codon, LOG(0.0009133));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('C', 'C', 'A')) == 0);
    nmm_codon_lprob_set(codonp, codon, LOG(0.0029179));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('C', 'C', 'G')) == 0);
    nmm_codon_lprob_set(codonp, codon, LOG(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('C', 'C', 'T')) == 0);
    nmm_codon_lprob_set(codonp, codon, LOG(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('C', 'G', 'C')) == 0);
    nmm_codon_lprob_set(codonp, codon, LOG(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('C', 'G', 'G')) == 0);
    nmm_codon_lprob_set(codonp, codon, LOG(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('C', 'G', 'T')) == 0);
    nmm_codon_lprob_set(codonp, codon, LOG(0.0041183));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('C', 'T', 'A')) == 0);
    nmm_codon_lprob_set(codonp, codon, LOG(0.0038173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('C', 'T', 'G')) == 0);
    nmm_codon_lprob_set(codonp, codon, LOG(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('C', 'T', 'T')) == 0);
    nmm_codon_lprob_set(codonp, codon, LOG(0.0029111));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('G', 'A', 'A')) == 0);
    nmm_codon_lprob_set(codonp, codon, LOG(0.0019173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('G', 'A', 'C')) == 0);
    nmm_codon_lprob_set(codonp, codon, LOG(0.0029103));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('G', 'A', 'G')) == 0);
    nmm_codon_lprob_set(codonp, codon, LOG(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('G', 'A', 'T')) == 0);
    nmm_codon_lprob_set(codonp, codon, LOG(0.0029103));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('G', 'C', 'A')) == 0);
    nmm_codon_lprob_set(codonp, codon, LOG(0.0003138));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('G', 'C', 'G')) == 0);
    nmm_codon_lprob_set(codonp, codon, LOG(0.0039173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('G', 'C', 'T')) == 0);
    nmm_codon_lprob_set(codonp, codon, LOG(0.0029103));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('G', 'G', 'C')) == 0);
    nmm_codon_lprob_set(codonp, codon, LOG(0.0019173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('G', 'G', 'G')) == 0);
    nmm_codon_lprob_set(codonp, codon, LOG(0.0009173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('G', 'G', 'T')) == 0);
    nmm_codon_lprob_set(codonp, codon, LOG(0.0029143));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('G', 'T', 'A')) == 0);
    nmm_codon_lprob_set(codonp, codon, LOG(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('G', 'T', 'G')) == 0);
    nmm_codon_lprob_set(codonp, codon, LOG(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('G', 'T', 'T')) == 0);
    nmm_codon_lprob_set(codonp, codon, LOG(0.0029171));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('T', 'A', 'A')) == 0);
    nmm_codon_lprob_set(codonp, codon, LOG(0.0099113));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('T', 'A', 'C')) == 0);
    nmm_codon_lprob_set(codonp, codon, LOG(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('T', 'A', 'G')) == 0);
    nmm_codon_lprob_set(codonp, codon, LOG(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('T', 'A', 'T')) == 0);
    nmm_codon_lprob_set(codonp, codon, LOG(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('T', 'C', 'A')) == 0);
    nmm_codon_lprob_set(codonp, codon, LOG(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('T', 'C', 'G')) == 0);
    nmm_codon_lprob_set(codonp, codon, LOG(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('T', 'C', 'T')) == 0);
    nmm_codon_lprob_set(codonp, codon, LOG(0.0020173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('T', 'G', 'C')) == 0);
    nmm_codon_lprob_set(codonp, codon, LOG(0.0029193));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('T', 'G', 'G')) == 0);
    nmm_codon_lprob_set(codonp, codon, LOG(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('T', 'G', 'T')) == 0);
    nmm_codon_lprob_set(codonp, codon, LOG(0.0029173));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('T', 'T', 'A')) == 0);
    nmm_codon_lprob_set(codonp, codon, LOG(0.0089193));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('T', 'T', 'G')) == 0);
    nmm_codon_lprob_set(codonp, codon, LOG(0.0029138));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('T', 'T', 'T')) == 0);
    nmm_codon_lprob_set(codonp, codon, LOG(0.0089183));

    nmm_codon_destroy(codon);
    return codonp;
}

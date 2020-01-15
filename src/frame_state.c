#include "imm/imm.h"
#include "logaddexp.h"
#include "nmm/nmm.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define ASCII_LAST_STD 127
#define NSYMBOLS (NMM_CODON_NBASES + 1)

struct codon_lprob
{
    double lprob[NSYMBOLS * NSYMBOLS * NSYMBOLS];
    int    symbol_idx[ASCII_LAST_STD + 1];
};

static void codon_lprob_init(struct codon_lprob*      codon_lprob,
                             struct nmm_codont const* codont);

static inline int codon_lprob_idx(struct codon_lprob const* codon_lprob, char const* codon)
{
    return codon_lprob->symbol_idx[(size_t)codon[0]] +
           NSYMBOLS * codon_lprob->symbol_idx[(size_t)codon[1]] +
           NSYMBOLS * NSYMBOLS * codon_lprob->symbol_idx[(size_t)codon[2]];
}

static inline double codon_lprob_value(struct codon_lprob const* codon_lprob,
                                       char const*               codon)
{
    return codon_lprob->lprob[codon_lprob_idx(codon_lprob, codon)];
}

struct nmm_frame_state
{
    struct imm_state*        interface;
    struct nmm_baset const*  baset;
    struct nmm_codont const* codon;
    struct imm_abc const*    abc;
    double                   epsilon;
    double                   leps;
    double                   l1eps;
    double                   zero_lprob;
    char                     any_symbol;
    struct codon_lprob       codon_lprob;
};

#define ARRAY_LENGTH(arr) (sizeof(arr) / sizeof((arr)[0]))

static double frame_state_lprob(struct imm_state const* state, char const* seq, int seq_len);
static int    frame_state_min_seq(struct imm_state const* state);
static int    frame_state_max_seq(struct imm_state const* state);

static double joint_seq_len1(struct nmm_frame_state const* state, char const* seq);
static double joint_seq_len2(struct nmm_frame_state const* state, char const* seq);
static double joint_seq_len3(struct nmm_frame_state const* state, char const* seq);
static double joint_seq_len4(struct nmm_frame_state const* state, char const* seq);
static double joint_seq_len5(struct nmm_frame_state const* state, char const* seq);
static double lprob_frag_given_codon1(struct nmm_frame_state const* state, char const* seq,
                                      struct nmm_codon const* ccode);
static double lprob_frag_given_codon2(struct nmm_frame_state const* state, char const* seq,
                                      struct nmm_codon const* ccode);
static double lprob_frag_given_codon3(struct nmm_frame_state const* state, char const* seq,
                                      struct nmm_codon const* ccode);
static double lprob_frag_given_codon4(struct nmm_frame_state const* state, char const* seq,
                                      struct nmm_codon const* ccode);
static double lprob_frag_given_codon5(struct nmm_frame_state const* state, char const* seq,
                                      struct nmm_codon const* ccode);
static double compute_codon_lprob(char const* bases_comb, struct nmm_codont const* codont,
                                  char const* codon);
static double base_lprob(struct nmm_frame_state const* state, char id);
static inline double logaddexp3(double const a, double const b, double const c)
{
    return logaddexp(logaddexp(a, b), c);
}
static inline double logsumexp(double const* arr, int len) { return imm_lprob_sum(arr, len); }
static inline double ecodon_lprob_value(struct codon_lprob const* codon_lprob,
                                        char const* seq, int const a, int const b,
                                        int const c)
{
    char const codon[3] = {seq[a], seq[b], seq[c]};
    return codon_lprob_value(codon_lprob, codon);
}

struct nmm_frame_state* nmm_frame_state_create(char const*              name,
                                               struct nmm_baset const*  baset,
                                               struct nmm_codont const* codont,
                                               double                   epsilon)
{
    struct nmm_frame_state* state = malloc(sizeof(struct nmm_frame_state));
    struct imm_abc const*   abc = nmm_baset_get_abc(baset);

    if (abc != nmm_codont_get_abc(codont)) {
        free(state);
        imm_error("alphabets from base and codon are different");
        return NULL;
    }

    if (imm_abc_length(abc) != NMM_CODON_NBASES) {
        free(state);
        imm_error("alphabet length is not four");
        return NULL;
    }

    state->abc = abc;
    state->baset = baset;
    state->codon = codont;
    state->epsilon = epsilon;
    state->leps = log(epsilon);
    state->l1eps = log(1 - epsilon);
    state->zero_lprob = imm_lprob_zero();
    state->any_symbol = imm_abc_any_symbol(abc);

    codon_lprob_init(&state->codon_lprob, codont);

    struct imm_state_funcs funcs = {frame_state_lprob, frame_state_min_seq,
                                    frame_state_max_seq};

    state->interface = imm_state_create(name, nmm_baset_get_abc(baset), funcs, state);
    return state;
}

double nmm_frame_state_lposterior(struct nmm_frame_state* state,
                                  struct nmm_codon const* codon, char const* seq, int seq_len)
{
    double lprob = state->zero_lprob;

    if (seq_len == 1)
        lprob = lprob_frag_given_codon1(state, seq, codon);
    else if (seq_len == 2)
        lprob = lprob_frag_given_codon2(state, seq, codon);
    else if (seq_len == 3)
        lprob = lprob_frag_given_codon3(state, seq, codon);
    else if (seq_len == 4)
        lprob = lprob_frag_given_codon4(state, seq, codon);
    else if (seq_len == 5)
        lprob = lprob_frag_given_codon5(state, seq, codon);

    return lprob + nmm_codont_get_lprob(state->codon, codon);
}

double nmm_frame_state_decode(struct nmm_frame_state* state, char const* seq, int seq_len,
                              struct nmm_codon* codon)
{
    struct imm_abc const* abc = nmm_codont_get_abc(state->codon);
    char const*           symbols = imm_abc_symbols(abc);
    int                   n = imm_abc_length(abc);

    double max_lprob = state->zero_lprob;

    for (int i0 = 0; i0 < n; ++i0) {
        for (int i1 = 0; i1 < n; ++i1) {
            for (int i2 = 0; i2 < n; ++i2) {

                struct nmm_codon tmp = {symbols[i0], symbols[i1], symbols[i2]};
                double lprob = nmm_frame_state_lposterior(state, &tmp, seq, seq_len);

                if (lprob >= max_lprob) {
                    max_lprob = lprob;
                    codon->a = tmp.a;
                    codon->b = tmp.b;
                    codon->c = tmp.c;
                }
            }
        }
    }
    return max_lprob;
}

void nmm_frame_state_destroy(struct nmm_frame_state* state)
{
    if (!state)
        return;

    imm_state_destroy(state->interface);
    state->interface = NULL;

    state->baset = NULL;
    state->codon = NULL;
    free(state);
}

static double frame_state_lprob(struct imm_state const* state, char const* seq, int seq_len)
{
    const struct nmm_frame_state* s = imm_state_get_impl_c(state);
    if (seq_len == 1)
        return joint_seq_len1(s, seq);
    else if (seq_len == 2)
        return joint_seq_len2(s, seq);
    else if (seq_len == 3)
        return joint_seq_len3(s, seq);
    else if (seq_len == 4)
        return joint_seq_len4(s, seq);
    else if (seq_len == 5)
        return joint_seq_len5(s, seq);

    return s->zero_lprob;
}

static int frame_state_min_seq(struct imm_state const* state) { return 1; }

static int frame_state_max_seq(struct imm_state const* state) { return 5; }

static double joint_seq_len1(struct nmm_frame_state const* state, char const* seq)
{
    const char _ = state->any_symbol;
    const char c0__[3] = {seq[0], _, _};
    const char c_0_[3] = {_, seq[0], _};
    const char c__0[3] = {_, _, seq[0]};

    double c = 2 * state->leps + 2 * state->l1eps;

    double e0 = codon_lprob_value(&state->codon_lprob, c0__);
    double e1 = codon_lprob_value(&state->codon_lprob, c_0_);
    double e2 = codon_lprob_value(&state->codon_lprob, c__0);
    return c + logaddexp3(e0, e1, e2) - log(3);
}

static double joint_seq_len2(struct nmm_frame_state const* state, char const* seq)
{
#define c_lprob(codon) codon_lprob_value(&state->codon_lprob, codon)
    const char _ = state->any_symbol;

    const char c_01[3] = {_, seq[0], seq[1]};
    const char c0_1[3] = {seq[0], _, seq[1]};
    const char c01_[3] = {seq[0], seq[1], _};

    const char c0__[3] = {seq[0], _, _};
    const char c_0_[3] = {_, seq[0], _};
    const char c__0[3] = {_, _, seq[0]};

    const char c1__[3] = {seq[1], _, _};
    const char c_1_[3] = {_, seq[1], _};
    const char c__1[3] = {_, _, seq[1]};

    const double b_lp0 = base_lprob(state, seq[0]);
    const double b_lp1 = base_lprob(state, seq[1]);

    double v0 = log(2) + state->leps + state->l1eps * 3 - log(3);
    v0 += logaddexp3(c_lprob(c_01), c_lprob(c0_1), c_lprob(c01_));

    double c = 3 * state->leps + state->l1eps - log(3);
    double v1 = c;
    v1 += logaddexp3(c_lprob(c0__), c_lprob(c_0_), c_lprob(c__0)) + b_lp1;

    double v2 = c;
    v2 += logaddexp3(c_lprob(c1__), c_lprob(c_1_), c_lprob(c__1)) + b_lp0;

    return logaddexp3(v0, v1, v2);
#undef c_lprob
}

static double joint_seq_len3(struct nmm_frame_state const* state, char const* seq)
{
#define C(a, b, c) ecodon_lprob_value(&state->codon_lprob, eseq, a, b, c)
    const char eseq[] = {seq[0], seq[1], seq[2], state->any_symbol};
    const char _ = sizeof(eseq) - 1;

    const double B[] = {base_lprob(state, seq[0]), base_lprob(state, seq[1]),
                        base_lprob(state, seq[2])};

    double v0 = 4 * state->l1eps + C(0, 1, 2);

    double c1 = log(4) + 2 * state->leps + 2 * state->l1eps - log(9);
    double v1[] = {logaddexp3(C(_, 1, 2), C(1, _, 2), C(1, 2, _)) + B[0],
                   logaddexp3(C(_, 0, 2), C(0, _, 2), C(0, 2, _)) + B[1],
                   logaddexp3(C(_, 0, 1), C(0, _, 1), C(0, 1, _)) + B[2]};

    double c2 = 4 * state->leps - log(9);
    double v2[] = {logaddexp3(C(2, _, _), C(_, 2, _), C(_, _, 2)) + B[0] + B[1],
                   logaddexp3(C(1, _, _), C(_, 1, _), C(_, _, 1)) + B[0] + B[2],
                   logaddexp3(C(0, _, _), C(_, 0, _), C(_, _, 0)) + B[1] + B[2]};

    return logaddexp3(v0, c1 + logsumexp(v1, 3), c2 + logsumexp(v2, 3));
#undef C
}

static double joint_seq_len4(struct nmm_frame_state const* state, char const* seq)
{
#define C(a, b, c) ecodon_lprob_value(&state->codon_lprob, eseq, a, b, c)
    const char eseq[] = {seq[0], seq[1], seq[2], seq[3], state->any_symbol};
    const char _ = sizeof(eseq) - 1;

    const double B[] = {base_lprob(state, seq[0]), base_lprob(state, seq[1]),
                        base_lprob(state, seq[2]), base_lprob(state, seq[3])};

    double c0 = state->leps + state->l1eps * 3 - log(2);
    double v0[] = {C(1, 2, 3) + B[0], C(0, 2, 3) + B[1], C(0, 1, 3) + B[2],
                   C(0, 1, 2) + B[3]};

    double c1 = 3 * state->leps + state->l1eps - log(9);
    double v1[] = {
        C(_, 2, 3) + B[0] + B[1], C(_, 1, 3) + B[0] + B[2], C(_, 1, 2) + B[0] + B[3],
        C(_, 0, 3) + B[1] + B[2], C(_, 0, 2) + B[1] + B[3], C(_, 0, 1) + B[2] + B[3],
        C(2, _, 3) + B[0] + B[1], C(1, _, 3) + B[0] + B[2], C(1, _, 2) + B[0] + B[3],
        C(0, _, 3) + B[1] + B[2], C(0, _, 2) + B[1] + B[3], C(0, _, 1) + B[2] + B[3],
        C(2, 3, _) + B[0] + B[1], C(1, 3, _) + B[0] + B[2], C(1, 2, _) + B[0] + B[3],
        C(0, 3, _) + B[1] + B[2], C(0, 2, _) + B[1] + B[3], C(0, 1, _) + B[2] + B[3]};

    return logaddexp(c0 + logsumexp(v0, 4), c1 + logsumexp(v1, 18));
#undef C
}

static double joint_seq_len5(struct nmm_frame_state const* state, char const* seq)
{
#define c_lp(codon) codon_lprob_value(&state->codon_lprob, codon)
    const char c012[3] = {seq[0], seq[1], seq[2]};
    const char c013[3] = {seq[0], seq[1], seq[3]};
    const char c014[3] = {seq[0], seq[1], seq[4]};
    const char c023[3] = {seq[0], seq[2], seq[3]};
    const char c024[3] = {seq[0], seq[2], seq[4]};
    const char c034[3] = {seq[0], seq[3], seq[4]};
    const char c123[3] = {seq[1], seq[2], seq[3]};
    const char c124[3] = {seq[1], seq[2], seq[4]};
    const char c134[3] = {seq[1], seq[3], seq[4]};
    const char c234[3] = {seq[2], seq[3], seq[4]};

    const double b_lp0 = base_lprob(state, seq[0]);
    const double b_lp1 = base_lprob(state, seq[1]);
    const double b_lp2 = base_lprob(state, seq[2]);
    const double b_lp3 = base_lprob(state, seq[3]);
    const double b_lp4 = base_lprob(state, seq[4]);

    double v[] = {logaddexp(b_lp0 + b_lp1 + c_lp(c234), b_lp0 + b_lp2 + c_lp(c134)),
                  logaddexp(b_lp0 + b_lp3 + c_lp(c124), b_lp0 + b_lp4 + c_lp(c123)),
                  logaddexp(b_lp1 + b_lp2 + c_lp(c034), b_lp1 + b_lp3 + c_lp(c024)),
                  logaddexp(b_lp1 + b_lp4 + c_lp(c023), b_lp2 + b_lp3 + c_lp(c014)),
                  logaddexp(b_lp2 + b_lp4 + c_lp(c013), b_lp3 + b_lp4 + c_lp(c012))};

    return 2 * state->leps + 2 * state->l1eps - log(10) + logsumexp(v, 5);
#undef c_lprob
}

static double compute_codon_lprob(char const* bases_comb, struct nmm_codont const* codont,
                                  char const* codon)
{
    double lprob = imm_lprob_zero();
    char   bases[3 * 4] = {
        bases_comb[0], bases_comb[1], bases_comb[2], bases_comb[3],
        bases_comb[0], bases_comb[1], bases_comb[2], bases_comb[3],
        bases_comb[0], bases_comb[1], bases_comb[2], bases_comb[3],
    };
    int nbases[3] = {NMM_CODON_NBASES, NMM_CODON_NBASES, NMM_CODON_NBASES};

    char const any_symbol = imm_abc_any_symbol(nmm_codont_get_abc(codont));
    for (int i = 0; i < 3; ++i) {
        if (codon[i] != any_symbol) {
            bases[i * NMM_CODON_NBASES] = codon[i];
            nbases[i] = 1;
        }
    }

    char const* a_id = bases;
    char const* b_id = bases + NMM_CODON_NBASES;
    char const* c_id = bases + 2 * NMM_CODON_NBASES;
    for (int a = 0; a < nbases[0]; ++a) {
        for (int b = 0; b < nbases[1]; ++b) {
            for (int c = 0; c < nbases[2]; ++c) {
                struct nmm_codon const ccode = {a_id[a], b_id[b], c_id[c]};
                double const           t = nmm_codont_get_lprob(codont, &ccode);
                lprob = logaddexp(lprob, t);
            }
        }
    }

    return lprob;
}

static double base_lprob(struct nmm_frame_state const* state, char id)
{
    return nmm_baset_get_lprob(state->baset, id);
}

static double lprob_frag_given_codon1(struct nmm_frame_state const* state, char const* seq,
                                      struct nmm_codon const* codon)
{
    double const loge = state->leps;
    double const log1e = state->l1eps;

    char const x1 = codon->a;
    char const x2 = codon->b;
    char const x3 = codon->c;

    char const z1 = seq[0];

    double const c = 2 * loge + 2 * log1e;

    return c + log((x1 == z1) + (x2 == z1) + (x3 == z1)) - log(3);
}

static double lprob_frag_given_codon2(struct nmm_frame_state const* state, char const* seq,
                                      struct nmm_codon const* codon)
{
    double const loge = state->leps;
    double const log1e = state->l1eps;

    char const x1 = codon->a;
    char const x2 = codon->b;
    char const x3 = codon->c;

    char const z1 = seq[0];
    char const z2 = seq[1];

    double const lprob_z1 = base_lprob(state, z1);
    double const lprob_z2 = base_lprob(state, z2);

    double const c1 = log(2) + loge + log1e * 3 - log(3);
    double const v0 =
        c1 + log((x2 == z1) * (x3 == z2) + (x1 == z1) * (x3 == z2) + (x1 == z1) * (x2 == z2));

    double const c2 = 3 * loge + log1e - log(3);

    double const v1 = c2 + log((x1 == z1) + (x2 == z1) + (x3 == z1)) + lprob_z2;
    double const v2 = c2 + log((x1 == z2) + (x2 == z2) + (x3 == z2)) + lprob_z1;

    return logaddexp3(v0, v1, v2);
}

static double lprob_frag_given_codon3(struct nmm_frame_state const* state, char const* seq,
                                      struct nmm_codon const* codon)
{
    double const loge = state->leps;
    double const log1e = state->l1eps;

    char const x1 = codon->a;
    char const x2 = codon->b;
    char const x3 = codon->c;

    char const z1 = seq[0];
    char const z2 = seq[1];
    char const z3 = seq[2];

    double const lprob_z1 = base_lprob(state, z1);
    double const lprob_z2 = base_lprob(state, z2);
    double const lprob_z3 = base_lprob(state, z3);

    double const v0 = 4 * log1e + log((x1 == z1) * (x2 == z2) * (x3 == z3));

    double const c1 = log(4) + 2 * loge + 2 * log1e - log(9);

    double const v1 =
        c1 +
        log((x2 == z2) * (x3 == z3) + (x1 == z2) * (x3 == z3) + (x1 == z2) * (x2 == z3)) +
        lprob_z1;

    double const v2 =
        c1 +
        log((x2 == z1) * (x3 == z3) + (x1 == z1) * (x3 == z3) + (x1 == z1) * (x2 == z3)) +
        lprob_z2;

    double const v3 =
        c1 +
        log((x2 == z1) * (x3 == z2) + (x1 == z1) * (x3 == z2) + (x1 == z1) * (x2 == z2)) +
        lprob_z3;

    double const c2 = 4 * loge - log(9);

    double const v4 = c2 + log((x1 == z3) + (x2 == z3) + (x3 == z3)) + lprob_z1 + lprob_z2;
    double const v5 = c2 + log((x1 == z2) + (x2 == z2) + (x3 == z2)) + lprob_z1 + lprob_z3;
    double const v6 = c2 + log((x1 == z1) + (x2 == z1) + (x3 == z1)) + lprob_z2 + lprob_z3;

    double v[] = {v0, v1, v2, v3, v4, v5, v6};

    return logsumexp(v, ARRAY_LENGTH(v));
}

static double lprob_frag_given_codon4(struct nmm_frame_state const* state, char const* seq,
                                      struct nmm_codon const* codon)
{
    double const loge = state->leps;
    double const log1e = state->l1eps;

    char const x1 = codon->a;
    char const x2 = codon->b;
    char const x3 = codon->c;

    char const z1 = seq[0];
    char const z2 = seq[1];
    char const z3 = seq[2];
    char const z4 = seq[3];

    double const lprob_z1 = base_lprob(state, z1);
    double const lprob_z2 = base_lprob(state, z2);
    double const lprob_z3 = base_lprob(state, z3);
    double const lprob_z4 = base_lprob(state, z4);

    double const v0[] = {log((x1 == z2) * (x2 == z3) * (x3 == z4)) + lprob_z1,
                         log((x1 == z1) * (x2 == z3) * (x3 == z4)) + lprob_z2,
                         log((x1 == z1) * (x2 == z2) * (x3 == z4)) + lprob_z3,
                         log((x1 == z1) * (x2 == z2) * (x3 == z3)) + lprob_z4};

    double const v1[] = {
        log((x2 == z3) * (x3 == z4)) + lprob_z1 + lprob_z2,
        log((x2 == z2) * (x3 == z4)) + lprob_z1 + lprob_z3,
        log((x2 == z2) * (x3 == z3)) + lprob_z1 + lprob_z4,
        log((x2 == z1) * (x3 == z4)) + lprob_z2 + lprob_z3,
        log((x2 == z1) * (x3 == z3)) + lprob_z2 + lprob_z4,
        log((x2 == z1) * (x3 == z2)) + lprob_z3 + lprob_z4,
        log((x1 == z3) * (x3 == z4)) + lprob_z1 + lprob_z2,
        log((x1 == z2) * (x3 == z4)) + lprob_z1 + lprob_z3,
        log((x1 == z2) * (x3 == z3)) + lprob_z1 + lprob_z4,
        log((x1 == z1) * (x3 == z4)) + lprob_z2 + lprob_z3,
        log((x1 == z1) * (x3 == z3)) + lprob_z2 + lprob_z4,
        log((x1 == z1) * (x3 == z2)) + lprob_z3 + lprob_z4,
        log((x1 == z3) * (x2 == z4)) + lprob_z1 + lprob_z2,
        log((x1 == z2) * (x2 == z4)) + lprob_z1 + lprob_z3,
        log((x1 == z2) * (x2 == z3)) + lprob_z1 + lprob_z4,
        log((x1 == z1) * (x2 == z4)) + lprob_z2 + lprob_z3,
        log((x1 == z1) * (x2 == z3)) + lprob_z2 + lprob_z4,
        log((x1 == z1) * (x2 == z2)) + lprob_z3 + lprob_z4,
    };

    return logaddexp(loge + log1e * 3 - log(2) + logsumexp(v0, ARRAY_LENGTH(v0)),
                     3 * loge + log1e - log(9) + logsumexp(v1, ARRAY_LENGTH(v1)));
}

static double lprob_frag_given_codon5(struct nmm_frame_state const* state, char const* seq,
                                      struct nmm_codon const* codon)
{
    double const loge = state->leps;
    double const log1e = state->l1eps;

    char const x1 = codon->a;
    char const x2 = codon->b;
    char const x3 = codon->c;

    char const z1 = seq[0];
    char const z2 = seq[1];
    char const z3 = seq[2];
    char const z4 = seq[3];
    char const z5 = seq[4];

    double const lprob_z1 = base_lprob(state, z1);
    double const lprob_z2 = base_lprob(state, z2);
    double const lprob_z3 = base_lprob(state, z3);
    double const lprob_z4 = base_lprob(state, z4);
    double const lprob_z5 = base_lprob(state, z5);

    double const v[] = {
        lprob_z1 + lprob_z2 + log((x1 == z3) * (x2 == z4) * (x3 == z5)),
        lprob_z1 + lprob_z3 + log((x1 == z2) * (x2 == z4) * (x3 == z5)),
        lprob_z1 + lprob_z4 + log((x1 == z2) * (x2 == z3) * (x3 == z5)),
        lprob_z1 + lprob_z5 + log((x1 == z2) * (x2 == z3) * (x3 == z4)),
        lprob_z2 + lprob_z3 + log((x1 == z1) * (x2 == z4) * (x3 == z5)),
        lprob_z2 + lprob_z4 + log((x1 == z1) * (x2 == z3) * (x3 == z5)),
        lprob_z2 + lprob_z5 + log((x1 == z1) * (x2 == z3) * (x3 == z4)),
        lprob_z3 + lprob_z4 + log((x1 == z1) * (x2 == z2) * (x3 == z5)),
        lprob_z3 + lprob_z5 + log((x1 == z1) * (x2 == z2) * (x3 == z4)),
        lprob_z4 + lprob_z5 + log((x1 == z1) * (x2 == z2) * (x3 == z3)),
    };

    return 2 * loge + 2 * log1e - log(10) + logsumexp(v, ARRAY_LENGTH(v));
}

static void codon_lprob_init(struct codon_lprob* codon_lprob, struct nmm_codont const* codont)
{
    struct imm_abc const* abc = nmm_codont_get_abc(codont);
    char const            bases_comb[3 * NMM_CODON_NBASES] = {
        imm_abc_symbol_id(abc, 0), imm_abc_symbol_id(abc, 1), imm_abc_symbol_id(abc, 2),
        imm_abc_symbol_id(abc, 3), imm_abc_symbol_id(abc, 0), imm_abc_symbol_id(abc, 1),
        imm_abc_symbol_id(abc, 2), imm_abc_symbol_id(abc, 3), imm_abc_symbol_id(abc, 0),
        imm_abc_symbol_id(abc, 1), imm_abc_symbol_id(abc, 2), imm_abc_symbol_id(abc, 3)};

    char const any_symbol = imm_abc_any_symbol(abc);

    for (int i = 0; i <= ASCII_LAST_STD; ++i)
        codon_lprob->symbol_idx[i] = -1;

    char const* symbols = imm_abc_symbols(abc);
    char const  all_symbols[NSYMBOLS] = {symbols[0], symbols[1], symbols[2], symbols[3],
                                        any_symbol};

    for (int i = 0; i < NMM_CODON_NBASES; ++i)
        codon_lprob->symbol_idx[(size_t)symbols[i]] = imm_abc_symbol_idx(abc, symbols[i]);

    codon_lprob->symbol_idx[(size_t)any_symbol] = NMM_CODON_NBASES;

    double* lprob = codon_lprob->lprob;
    for (int i0 = 0; i0 < NSYMBOLS; ++i0) {
        for (int i1 = 0; i1 < NSYMBOLS; ++i1) {
            for (int i2 = 0; i2 < NSYMBOLS; ++i2) {

                char const codon[3] = {all_symbols[i0], all_symbols[i1], all_symbols[i2]};

                int const i = codon_lprob_idx(codon_lprob, codon);
                lprob[i] = compute_codon_lprob(bases_comb, codont, codon);
            }
        }
    }
}

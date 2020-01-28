#include "codon_lmarg.h"
#include "imm/imm.h"
#include "logaddexp.h"
#include "nmm/nmm.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

struct nmm_frame_state
{
    struct imm_state*        interface;
    struct nmm_baset const*  baset;
    struct nmm_codont const* codont;
    double                   epsilon;
    double                   leps;
    double                   l1eps;
    double                   zero_lprob;
    char                     any_symbol;
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
static inline double base_lprob(struct nmm_frame_state const* state, char id)
{
    return nmm_baset_lprob(state->baset, id);
}
static inline double logaddexp3(double const a, double const b, double const c)
{
    return logaddexp(logaddexp(a, b), c);
}
static inline double logsumexp(double const* arr, int len) { return imm_lprob_sum(arr, len); }

struct nmm_frame_state* nmm_frame_state_create(char const*              name,
                                               struct nmm_baset const*  baset,
                                               struct nmm_codont const* codont,
                                               double                   epsilon)
{
    struct nmm_frame_state* state = malloc(sizeof(struct nmm_frame_state));
    struct imm_abc const*   abc = nmm_base_get_abc(nmm_baset_get_base(baset));

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

    state->baset = baset;
    state->codont = codont;
    state->epsilon = epsilon;
    state->leps = log(epsilon);
    state->l1eps = log(1 - epsilon);
    state->zero_lprob = imm_lprob_zero();
    state->any_symbol = imm_abc_any_symbol(abc);

    struct imm_state_funcs funcs = {frame_state_lprob, frame_state_min_seq,
                                    frame_state_max_seq};

    state->interface =
        imm_state_create(name, nmm_base_get_abc(nmm_baset_get_base(baset)), funcs, state);
    return state;
}

double nmm_frame_state_lposterior(struct nmm_frame_state const* state,
                                  struct nmm_codon const* codon, char const* seq,
                                  int const seq_len)
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

    return lprob + nmm_codont_lprob(state->codont, codon);
}

double nmm_frame_state_decode(struct nmm_frame_state const* state, char const* seq,
                              int const seq_len, struct nmm_codon* codon)
{
    struct imm_abc const* abc = nmm_codont_get_abc(state->codont);
    char const*           symbols = imm_abc_symbols(abc);
    int const             n = imm_abc_length(abc);

    double max_lprob = state->zero_lprob;

    for (int i0 = 0; i0 < n; ++i0) {
        for (int i1 = 0; i1 < n; ++i1) {
            for (int i2 = 0; i2 < n; ++i2) {

                struct nmm_codon const tmp = {symbols[i0], symbols[i1], symbols[i2]};
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
    if (!state) {
        imm_error("state should not be NULL");
        return;
    }

    imm_state_destroy(state->interface);
    state->interface = NULL;

    state->baset = NULL;
    state->codont = NULL;
    free(state);
}

static double frame_state_lprob(struct imm_state const* state, char const* seq,
                                int const seq_len)
{
    struct nmm_frame_state const* s = imm_state_get_impl_c(state);
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
    char const _ = state->any_symbol;

    double const c = 2 * state->leps + 2 * state->l1eps;

    double const e0 = nmm_codont_lprob(state->codont, &NMM_CODON(seq[0], _, _));
    double const e1 = nmm_codont_lprob(state->codont, &NMM_CODON(_, seq[0], _));
    double const e2 = nmm_codont_lprob(state->codont, &NMM_CODON(_, _, seq[0]));

    return c + logaddexp3(e0, e1, e2) - log(3);
}

static double joint_seq_len2(struct nmm_frame_state const* state, char const* seq)
{
#define C(a, b, c) nmm_codont_lprob(state->codont, &NMM_CODON(eseq[a], eseq[b], eseq[c]))
    char const   eseq[] = {seq[0], seq[1], state->any_symbol};
    size_t const _ = sizeof(eseq) - 1;

    double const b_lp0 = base_lprob(state, seq[0]);
    double const b_lp1 = base_lprob(state, seq[1]);

    double const c0 = log(2) + state->leps + state->l1eps * 3 - log(3);
    double const v0 = c0 + logaddexp3(C(_, 0, 1), C(0, _, 1), C(0, 1, _));

    double const c1 = 3 * state->leps + state->l1eps - log(3);
    double const v1 = c1 + logaddexp3(C(0, _, _), C(_, 0, _), C(_, _, 0)) + b_lp1;
    double const v2 = c1 + logaddexp3(C(1, _, _), C(_, 1, _), C(_, _, 1)) + b_lp0;

    return logaddexp3(v0, v1, v2);
#undef C
}

static double joint_seq_len3(struct nmm_frame_state const* state, char const* seq)
{
#define C(a, b, c) nmm_codont_lprob(state->codont, &NMM_CODON(eseq[a], eseq[b], eseq[c]))
    char const   eseq[] = {seq[0], seq[1], seq[2], state->any_symbol};
    size_t const _ = sizeof(eseq) - 1;

    double const B[] = {base_lprob(state, seq[0]), base_lprob(state, seq[1]),
                        base_lprob(state, seq[2])};

    double const v0 = 4 * state->l1eps + C(0, 1, 2);

    double const c1 = log(4) + 2 * state->leps + 2 * state->l1eps - log(9);
    double const v1[] = {logaddexp3(C(_, 1, 2), C(1, _, 2), C(1, 2, _)) + B[0],
                         logaddexp3(C(_, 0, 2), C(0, _, 2), C(0, 2, _)) + B[1],
                         logaddexp3(C(_, 0, 1), C(0, _, 1), C(0, 1, _)) + B[2]};

    double const c2 = 4 * state->leps - log(9);
    double const v2[] = {logaddexp3(C(2, _, _), C(_, 2, _), C(_, _, 2)) + B[0] + B[1],
                         logaddexp3(C(1, _, _), C(_, 1, _), C(_, _, 1)) + B[0] + B[2],
                         logaddexp3(C(0, _, _), C(_, 0, _), C(_, _, 0)) + B[1] + B[2]};

    return logaddexp3(v0, c1 + logsumexp(v1, 3), c2 + logsumexp(v2, 3));
#undef C
}

static double joint_seq_len4(struct nmm_frame_state const* state, char const* seq)
{
#define C(a, b, c) nmm_codont_lprob(state->codont, &NMM_CODON(eseq[a], eseq[b], eseq[c]))
    char const   eseq[] = {seq[0], seq[1], seq[2], seq[3], state->any_symbol};
    size_t const _ = sizeof(eseq) - 1;

    double const B[] = {base_lprob(state, seq[0]), base_lprob(state, seq[1]),
                        base_lprob(state, seq[2]), base_lprob(state, seq[3])};

    double const c0 = state->leps + state->l1eps * 3 - log(2);
    double const v0[] = {C(1, 2, 3) + B[0], C(0, 2, 3) + B[1], C(0, 1, 3) + B[2],
                         C(0, 1, 2) + B[3]};

    double const c1 = 3 * state->leps + state->l1eps - log(9);
    double const v1[] = {
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
#define LPROB(codon) nmm_codont_lprob(state->codont, &codon)
    struct nmm_codon const c012 = NMM_CODON(seq[0], seq[1], seq[2]);
    struct nmm_codon const c013 = NMM_CODON(seq[0], seq[1], seq[3]);
    struct nmm_codon const c014 = NMM_CODON(seq[0], seq[1], seq[4]);
    struct nmm_codon const c023 = NMM_CODON(seq[0], seq[2], seq[3]);
    struct nmm_codon const c024 = NMM_CODON(seq[0], seq[2], seq[4]);
    struct nmm_codon const c034 = NMM_CODON(seq[0], seq[3], seq[4]);
    struct nmm_codon const c123 = NMM_CODON(seq[1], seq[2], seq[3]);
    struct nmm_codon const c124 = NMM_CODON(seq[1], seq[2], seq[4]);
    struct nmm_codon const c134 = NMM_CODON(seq[1], seq[3], seq[4]);
    struct nmm_codon const c234 = NMM_CODON(seq[2], seq[3], seq[4]);

    double const b_lp0 = base_lprob(state, seq[0]);
    double const b_lp1 = base_lprob(state, seq[1]);
    double const b_lp2 = base_lprob(state, seq[2]);
    double const b_lp3 = base_lprob(state, seq[3]);
    double const b_lp4 = base_lprob(state, seq[4]);

    double const v[] = {logaddexp(b_lp0 + b_lp1 + LPROB(c234), b_lp0 + b_lp2 + LPROB(c134)),
                        logaddexp(b_lp0 + b_lp3 + LPROB(c124), b_lp0 + b_lp4 + LPROB(c123)),
                        logaddexp(b_lp1 + b_lp2 + LPROB(c034), b_lp1 + b_lp3 + LPROB(c024)),
                        logaddexp(b_lp1 + b_lp4 + LPROB(c023), b_lp2 + b_lp3 + LPROB(c014)),
                        logaddexp(b_lp2 + b_lp4 + LPROB(c013), b_lp3 + b_lp4 + LPROB(c012))};

    return 2 * state->leps + 2 * state->l1eps - log(10) + logsumexp(v, 5);
#undef c_lprob
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

#include "nmm/frame_state.h"
#include "free.h"
#include "nmm/array_size.h"
#include "nmm/base_abc.h"
#include "nmm/base_table.h"
#include "nmm/codon.h"
#include "nmm/codon_table.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

struct nmm_frame_state
{
    struct imm_state const*       interface;
    struct nmm_base_table const*  baset;
    struct nmm_codon_table const* codont;
    double                        epsilon;
    double                        leps;
    double                        l1eps;
    double                        zero_lprob;
    char                          any_symbol;
};

static double   frame_state_lprob(struct imm_state const* state, struct imm_seq const* seq);
static unsigned frame_state_min_seq(struct imm_state const* state);
static unsigned frame_state_max_seq(struct imm_state const* state);

static double joint_seq_len1(struct nmm_frame_state const* state, struct imm_seq const* seq);
static double joint_seq_len2(struct nmm_frame_state const* state, struct imm_seq const* seq);
static double joint_seq_len3(struct nmm_frame_state const* state, struct imm_seq const* seq);
static double joint_seq_len4(struct nmm_frame_state const* state, struct imm_seq const* seq);
static double joint_seq_len5(struct nmm_frame_state const* state, struct imm_seq const* seq);

static double        lprob_frag_given_codon1(struct nmm_frame_state const* state,
                                             struct imm_seq const*         seq,
                                             struct nmm_codon const*       ccode);
static double        lprob_frag_given_codon2(struct nmm_frame_state const* state,
                                             struct imm_seq const*         seq,
                                             struct nmm_codon const*       ccode);
static double        lprob_frag_given_codon3(struct nmm_frame_state const* state,
                                             struct imm_seq const*         seq,
                                             struct nmm_codon const*       ccode);
static double        lprob_frag_given_codon4(struct nmm_frame_state const* state,
                                             struct imm_seq const*         seq,
                                             struct nmm_codon const*       ccode);
static double        lprob_frag_given_codon5(struct nmm_frame_state const* state,
                                             struct imm_seq const*         seq,
                                             struct nmm_codon const*       ccode);
static inline double base_lprob(struct nmm_frame_state const* state, char id)
{
    return nmm_base_table_lprob(state->baset, id);
}
static inline double logaddexp3(double const a, double const b, double const c)
{
    return imm_lprob_add(imm_lprob_add(a, b), c);
}

static inline struct nmm_codon const* codon_set(struct nmm_codon* codon, char a, char b,
                                                char c)
{
    codon->a = a;
    codon->b = b;
    codon->c = c;
    return codon;
}

struct nmm_frame_state const* nmm_frame_state_create(char const*                   name,
                                                     struct nmm_base_table const*  baset,
                                                     struct nmm_codon_table const* codont,
                                                     double const                  epsilon)
{
    struct nmm_frame_state* state = malloc(sizeof(struct nmm_frame_state));

    if (nmm_base_table_get_base_abc(baset) != nmm_codon_table_get_base(codont)) {
        free_c(state);
        return NULL;
    }

    state->baset = baset;
    state->codont = codont;
    state->epsilon = epsilon;
    state->leps = log(epsilon);
    state->l1eps = log(1 - epsilon);
    state->zero_lprob = imm_lprob_zero();
    state->any_symbol =
        imm_abc_any_symbol(nmm_base_abc_cast(nmm_base_table_get_base_abc(baset)));

    struct imm_state_funcs funcs = {frame_state_lprob, frame_state_min_seq,
                                    frame_state_max_seq};

    state->interface = imm_state_create(
        name, nmm_base_abc_cast(nmm_base_table_get_base_abc(baset)), funcs, state);
    return state;
}

double nmm_frame_state_lposterior(struct nmm_frame_state const* state,
                                  struct nmm_codon const* codon, struct imm_seq const* seq)
{
    double   lprob = state->zero_lprob;
    unsigned length = imm_seq_length(seq);

    if (length == 1)
        lprob = lprob_frag_given_codon1(state, seq, codon);
    else if (length == 2)
        lprob = lprob_frag_given_codon2(state, seq, codon);
    else if (length == 3)
        lprob = lprob_frag_given_codon3(state, seq, codon);
    else if (length == 4)
        lprob = lprob_frag_given_codon4(state, seq, codon);
    else if (length == 5)
        lprob = lprob_frag_given_codon5(state, seq, codon);

    return lprob + nmm_codon_table_lprob(state->codont, codon);
}

double nmm_frame_state_decode(struct nmm_frame_state const* state, struct imm_seq const* seq,
                              struct nmm_codon* codon)
{
    struct imm_abc const* abc = nmm_base_abc_cast(nmm_codon_table_get_base(state->codont));
    char const*           symbols = imm_abc_symbols(abc);
    unsigned const        n = imm_abc_length(abc);

    double            max_lprob = state->zero_lprob;
    struct nmm_codon* tmp = nmm_codon_create(nmm_codon_table_get_base(state->codont));

    for (unsigned i0 = 0; i0 < n; ++i0) {
        for (unsigned i1 = 0; i1 < n; ++i1) {
            for (unsigned i2 = 0; i2 < n; ++i2) {

                struct nmm_triplet triplet = {symbols[i0], symbols[i1], symbols[i2]};
                nmm_codon_set_triplet(tmp, triplet);
                double lprob = nmm_frame_state_lposterior(state, tmp, seq);

                if (lprob >= max_lprob) {
                    max_lprob = lprob;
                    nmm_codon_set_triplet(codon, triplet);
                }
            }
        }
    }
    nmm_codon_destroy(tmp);
    return max_lprob;
}

void nmm_frame_state_destroy(struct nmm_frame_state const* state)
{
    imm_state_destroy(state->interface);
    free_c(state);
}

static double frame_state_lprob(struct imm_state const* state, struct imm_seq const* seq)
{
    struct nmm_frame_state const* s = imm_state_get_impl_c(state);
    unsigned                      length = imm_seq_length(seq);

    if (length == 1)
        return joint_seq_len1(s, seq);
    else if (length == 2)
        return joint_seq_len2(s, seq);
    else if (length == 3)
        return joint_seq_len3(s, seq);
    else if (length == 4)
        return joint_seq_len4(s, seq);
    else if (length == 5)
        return joint_seq_len5(s, seq);

    return s->zero_lprob;
}

static unsigned frame_state_min_seq(struct imm_state const* state) { return 1; }

static unsigned frame_state_max_seq(struct imm_state const* state) { return 5; }

static double joint_seq_len1(struct nmm_frame_state const* state, struct imm_seq const* seq)
{
    char const _ = state->any_symbol;
    NMM_CODON_DECL(codon, nmm_base_table_get_base_abc(state->baset));

    double const c = 2 * state->leps + 2 * state->l1eps;

    char const   nucl = imm_seq_string(seq)[0];
    double const e0 = nmm_codon_table_lprob(state->codont, codon_set(&codon, nucl, _, _));
    double const e1 = nmm_codon_table_lprob(state->codont, codon_set(&codon, _, nucl, _));
    double const e2 = nmm_codon_table_lprob(state->codont, codon_set(&codon, _, _, nucl));

    return c + logaddexp3(e0, e1, e2) - log(3);
}

static double joint_seq_len2(struct nmm_frame_state const* state, struct imm_seq const* seq)
{
#define C(a, b, c)                                                                           \
    nmm_codon_table_lprob(state->codont, codon_set(&codon, eseq[a], eseq[b], eseq[c]))
    char const*  str = imm_seq_string(seq);
    char const   eseq[] = {str[0], str[1], state->any_symbol};
    size_t const _ = sizeof(eseq) - 1;
    NMM_CODON_DECL(codon, nmm_base_table_get_base_abc(state->baset));

    double const b_lp0 = base_lprob(state, str[0]);
    double const b_lp1 = base_lprob(state, str[1]);

    double const c0 = log(2) + state->leps + state->l1eps * 3 - log(3);
    double const v0 = c0 + logaddexp3(C(_, 0, 1), C(0, _, 1), C(0, 1, _));

    double const c1 = 3 * state->leps + state->l1eps - log(3);
    double const v1 = c1 + logaddexp3(C(0, _, _), C(_, 0, _), C(_, _, 0)) + b_lp1;
    double const v2 = c1 + logaddexp3(C(1, _, _), C(_, 1, _), C(_, _, 1)) + b_lp0;

    return logaddexp3(v0, v1, v2);
#undef C
}

static double joint_seq_len3(struct nmm_frame_state const* state, struct imm_seq const* seq)
{
#define C(a, b, c)                                                                           \
    nmm_codon_table_lprob(state->codont, codon_set(&codon, eseq[a], eseq[b], eseq[c]))
    char const*  str = imm_seq_string(seq);
    char const   eseq[] = {str[0], str[1], str[2], state->any_symbol};
    size_t const _ = sizeof(eseq) - 1;
    NMM_CODON_DECL(codon, nmm_base_table_get_base_abc(state->baset));

    double const B[] = {base_lprob(state, str[0]), base_lprob(state, str[1]),
                        base_lprob(state, str[2])};

    double const v0 = 4 * state->l1eps + C(0, 1, 2);

    double const c1 = log(4) + 2 * state->leps + 2 * state->l1eps - log(9);
    double const v1[] = {logaddexp3(C(_, 1, 2), C(1, _, 2), C(1, 2, _)) + B[0],
                         logaddexp3(C(_, 0, 2), C(0, _, 2), C(0, 2, _)) + B[1],
                         logaddexp3(C(_, 0, 1), C(0, _, 1), C(0, 1, _)) + B[2]};

    double const c2 = 4 * state->leps - log(9);
    double const v2[] = {logaddexp3(C(2, _, _), C(_, 2, _), C(_, _, 2)) + B[0] + B[1],
                         logaddexp3(C(1, _, _), C(_, 1, _), C(_, _, 1)) + B[0] + B[2],
                         logaddexp3(C(0, _, _), C(_, 0, _), C(_, _, 0)) + B[1] + B[2]};

    return logaddexp3(v0, c1 + imm_lprob_sum(v1, 3), c2 + imm_lprob_sum(v2, 3));
#undef C
}

static double joint_seq_len4(struct nmm_frame_state const* state, struct imm_seq const* seq)
{
#define C(a, b, c)                                                                           \
    nmm_codon_table_lprob(state->codont, codon_set(&codon, eseq[a], eseq[b], eseq[c]))
    char const*  str = imm_seq_string(seq);
    char const   eseq[] = {str[0], str[1], str[2], str[3], state->any_symbol};
    size_t const _ = sizeof(eseq) - 1;
    NMM_CODON_DECL(codon, nmm_base_table_get_base_abc(state->baset));

    double const B[] = {base_lprob(state, str[0]), base_lprob(state, str[1]),
                        base_lprob(state, str[2]), base_lprob(state, str[3])};

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

    return imm_lprob_add(c0 + imm_lprob_sum(v0, 4), c1 + imm_lprob_sum(v1, 18));
#undef C
}

static double joint_seq_len5(struct nmm_frame_state const* state, struct imm_seq const* seq)
{
#define LPROB(codon) nmm_codon_table_lprob(state->codont, codon)
#define C(a, b, c) codon_set(&codon, str[a], str[b], str[c])
    char const* str = imm_seq_string(seq);
    NMM_CODON_DECL(codon, nmm_base_table_get_base_abc(state->baset));

    double const b_lp0 = base_lprob(state, str[0]);
    double const b_lp1 = base_lprob(state, str[1]);
    double const b_lp2 = base_lprob(state, str[2]);
    double const b_lp3 = base_lprob(state, str[3]);
    double const b_lp4 = base_lprob(state, str[4]);

    double const v[] = {
        imm_lprob_add(b_lp0 + b_lp1 + LPROB(C(2, 3, 4)), b_lp0 + b_lp2 + LPROB(C(1, 3, 4))),
        imm_lprob_add(b_lp0 + b_lp3 + LPROB(C(1, 2, 4)), b_lp0 + b_lp4 + LPROB(C(1, 2, 3))),
        imm_lprob_add(b_lp1 + b_lp2 + LPROB(C(0, 3, 4)), b_lp1 + b_lp3 + LPROB(C(0, 2, 4))),
        imm_lprob_add(b_lp1 + b_lp4 + LPROB(C(0, 2, 3)), b_lp2 + b_lp3 + LPROB(C(0, 1, 4))),
        imm_lprob_add(b_lp2 + b_lp4 + LPROB(C(0, 1, 3)), b_lp3 + b_lp4 + LPROB(C(0, 1, 2)))};

    return 2 * state->leps + 2 * state->l1eps - log(10) + imm_lprob_sum(v, 5);
#undef C
#undef LPROB
}

static double lprob_frag_given_codon1(struct nmm_frame_state const* state,
                                      struct imm_seq const*         seq,
                                      struct nmm_codon const*       codon)
{
    double const loge = state->leps;
    double const log1e = state->l1eps;

    char const x1 = codon->a;
    char const x2 = codon->b;
    char const x3 = codon->c;

    char const z1 = imm_seq_string(seq)[0];

    double const c = 2 * loge + 2 * log1e;

    return c + log((x1 == z1) + (x2 == z1) + (x3 == z1)) - log(3);
}

static double lprob_frag_given_codon2(struct nmm_frame_state const* state,
                                      struct imm_seq const*         seq,
                                      struct nmm_codon const*       codon)
{
    double const loge = state->leps;
    double const log1e = state->l1eps;

    char const x1 = codon->a;
    char const x2 = codon->b;
    char const x3 = codon->c;

    char const z1 = imm_seq_string(seq)[0];
    char const z2 = imm_seq_string(seq)[1];

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

static double lprob_frag_given_codon3(struct nmm_frame_state const* state,
                                      struct imm_seq const*         seq,
                                      struct nmm_codon const*       codon)
{
    double const loge = state->leps;
    double const log1e = state->l1eps;

    char const x1 = codon->a;
    char const x2 = codon->b;
    char const x3 = codon->c;

    char const z1 = imm_seq_string(seq)[0];
    char const z2 = imm_seq_string(seq)[1];
    char const z3 = imm_seq_string(seq)[2];

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

    return imm_lprob_sum(v, NMM_ARRAY_SIZE(v));
}

static double lprob_frag_given_codon4(struct nmm_frame_state const* state,
                                      struct imm_seq const*         seq,
                                      struct nmm_codon const*       codon)
{
    double const loge = state->leps;
    double const log1e = state->l1eps;

    char const x1 = codon->a;
    char const x2 = codon->b;
    char const x3 = codon->c;

    char const z1 = imm_seq_string(seq)[0];
    char const z2 = imm_seq_string(seq)[1];
    char const z3 = imm_seq_string(seq)[2];
    char const z4 = imm_seq_string(seq)[3];

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

    return imm_lprob_add(loge + log1e * 3 - log(2) + imm_lprob_sum(v0, NMM_ARRAY_SIZE(v0)),
                         3 * loge + log1e - log(9) + imm_lprob_sum(v1, NMM_ARRAY_SIZE(v1)));
}

static double lprob_frag_given_codon5(struct nmm_frame_state const* state,
                                      struct imm_seq const*         seq,
                                      struct nmm_codon const*       codon)
{
    double const loge = state->leps;
    double const log1e = state->l1eps;

    char const x1 = codon->a;
    char const x2 = codon->b;
    char const x3 = codon->c;

    char const z1 = imm_seq_string(seq)[0];
    char const z2 = imm_seq_string(seq)[1];
    char const z3 = imm_seq_string(seq)[2];
    char const z4 = imm_seq_string(seq)[3];
    char const z5 = imm_seq_string(seq)[4];

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

    return 2 * loge + 2 * log1e - log(10) + imm_lprob_sum(v, NMM_ARRAY_SIZE(v));
}

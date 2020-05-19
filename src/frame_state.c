#include "nmm/frame_state.h"
#include "free.h"
#include "model.h"
#include "nmm/array_size.h"
#include "nmm/base_abc.h"
#include "nmm/base_table.h"
#include "nmm/codon.h"
#include "nmm/codon_table.h"
#include "nmm/model.h"
#include "nmm/state_types.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

struct nmm_frame_state
{
    struct imm_state const*       super;
    struct nmm_base_table const*  baset;
    struct nmm_codon_table const* codont;
    double                        epsilon;
    double                        leps;
    double                        l1eps;
    double                        zero_lprob;
    char                          any_symbol;
};

static inline double                  base_lprob(struct nmm_frame_state const* state, char id);
static inline struct nmm_codon const* codon_set(struct nmm_codon* codon, char a, char b, char c);
static void                           destroy(struct imm_state const* state);
static double        joint_seq_len1(struct nmm_frame_state const* state, struct imm_seq const* seq);
static double        joint_seq_len2(struct nmm_frame_state const* state, struct imm_seq const* seq);
static double        joint_seq_len3(struct nmm_frame_state const* state, struct imm_seq const* seq);
static double        joint_seq_len4(struct nmm_frame_state const* state, struct imm_seq const* seq);
static double        joint_seq_len5(struct nmm_frame_state const* state, struct imm_seq const* seq);
static inline double logaddexp3(double const a, double const b, double const c);
static double        lprob(struct imm_state const* state, struct imm_seq const* seq);
static double        lprob_frag_given_codon1(struct nmm_frame_state const* state,
                                             struct imm_seq const* seq, struct nmm_codon const* ccode);
static double        lprob_frag_given_codon2(struct nmm_frame_state const* state,
                                             struct imm_seq const* seq, struct nmm_codon const* ccode);
static double        lprob_frag_given_codon3(struct nmm_frame_state const* state,
                                             struct imm_seq const* seq, struct nmm_codon const* ccode);
static double        lprob_frag_given_codon4(struct nmm_frame_state const* state,
                                             struct imm_seq const* seq, struct nmm_codon const* ccode);
static double        lprob_frag_given_codon5(struct nmm_frame_state const* state,
                                             struct imm_seq const* seq, struct nmm_codon const* ccode);
static uint8_t       max_seq(struct imm_state const* state);
static uint8_t       min_seq(struct imm_state const* state);
static uint8_t       type_id(struct imm_state const* state);

static struct imm_state_vtable const __vtable = {destroy, lprob, max_seq, min_seq, type_id};

struct nmm_base_table const* nmm_frame_state_baset(struct nmm_frame_state const* state)
{
    return state->baset;
}

struct nmm_codon_table const* nmm_frame_state_codont(struct nmm_frame_state const* state)
{
    return state->codont;
}

struct nmm_frame_state const* nmm_frame_state_create(char const*                   name,
                                                     struct nmm_base_table const*  baset,
                                                     struct nmm_codon_table const* codont,
                                                     double const                  epsilon)
{
    struct nmm_frame_state* state = malloc(sizeof(*state));

    if (nmm_base_table_abc(baset) != nmm_codon_table_abc(codont)) {
        free_c(state);
        return NULL;
    }

    state->baset = baset;
    state->codont = codont;
    state->epsilon = epsilon;
    state->leps = log(epsilon);
    state->l1eps = log(1 - epsilon);
    state->zero_lprob = imm_lprob_zero();
    state->any_symbol = imm_abc_any_symbol(nmm_base_abc_super(nmm_base_table_abc(baset)));

    struct imm_abc const* abc = nmm_base_abc_super(nmm_base_table_abc(baset));
    state->super = imm_state_create(name, abc, __vtable, state);
    return state;
}

double nmm_frame_state_decode(struct nmm_frame_state const* state, struct imm_seq const* seq,
                              struct nmm_codon* codon)
{
    struct imm_abc const* abc = nmm_base_abc_super(nmm_codon_table_abc(state->codont));
    char const*           symbols = imm_abc_symbols(abc);
    unsigned const        n = imm_abc_length(abc);

    double            max_lprob = state->zero_lprob;
    struct nmm_codon* tmp = nmm_codon_create(nmm_codon_table_abc(state->codont));

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

struct nmm_frame_state const* nmm_frame_state_derived(struct imm_state const* state)
{
    if (imm_state_type_id(state) != NMM_FRAME_STATE_TYPE_ID) {
        imm_error("could not cast to frame_state");
        return NULL;
    }
    return __imm_state_derived(state);
}

void nmm_frame_state_destroy(struct nmm_frame_state const* state)
{
    state->super->vtable.destroy(state->super);
}

double nmm_frame_state_lposterior(struct nmm_frame_state const* state,
                                  struct nmm_codon const* codon, struct imm_seq const* seq)
{
    double   lprob = state->zero_lprob;
    unsigned length = imm_seq_length(seq);

    switch (length) {
    case 1:
        lprob = lprob_frag_given_codon1(state, seq, codon);
        break;
    case 2:
        lprob = lprob_frag_given_codon2(state, seq, codon);
        break;
    case 3:
        lprob = lprob_frag_given_codon3(state, seq, codon);
        break;
    case 4:
        lprob = lprob_frag_given_codon4(state, seq, codon);
        break;
    case 5:
        lprob = lprob_frag_given_codon5(state, seq, codon);
        break;
    default:
        return imm_lprob_invalid();
    }

    return lprob + nmm_codon_table_lprob(state->codont, codon);
}

struct imm_state const* nmm_frame_state_super(struct nmm_frame_state const* state)
{
    return state->super;
}

struct imm_state const* nmm_frame_state_read(FILE* stream, struct nmm_model const* model)
{
    struct imm_abc const* abc = nmm_model_abc(model);
    struct imm_state*     state = __imm_state_read(stream, abc);
    if (!state) {
        imm_error("could not state_read");
        return NULL;
    }

    state->vtable = __vtable;

    struct nmm_frame_state* frame_state = malloc(sizeof(*frame_state));
    frame_state->super = state;
    state->derived = frame_state;

    uint32_t index = 0;
    if (fread(&index, sizeof(index), 1, stream) < 1) {
        imm_error("could not read baset index");
        goto err;
    }

    struct nmm_base_table const* baset = nmm_model_base_table(model, index);
    if (!baset) {
        imm_error("could not get baset");
        goto err;
    }
    frame_state->baset = baset;

    index = 0;
    if (fread(&index, sizeof(index), 1, stream) < 1) {
        imm_error("could not read codont index");
        goto err;
    }

    struct nmm_codon_table const* codont = nmm_model_codon_table(model, index);
    if (!codont) {
        imm_error("could not get codont");
        goto err;
    }
    frame_state->codont = codont;

    if (fread(&frame_state->epsilon, sizeof(frame_state->epsilon), 1, stream) < 1) {
        imm_error("could not read epsilon");
        goto err;
    }

    frame_state->leps = log(frame_state->epsilon);
    frame_state->l1eps = log(1 - frame_state->epsilon);
    frame_state->zero_lprob = imm_lprob_zero();
    frame_state->any_symbol = imm_abc_any_symbol(nmm_base_abc_super(nmm_base_table_abc(baset)));

    return state;

err:
    free_c(frame_state);
    __imm_state_destroy(state);
    return NULL;
}

static inline double base_lprob(struct nmm_frame_state const* state, char id)
{
    return nmm_base_table_lprob(state->baset, id);
}

static inline struct nmm_codon const* codon_set(struct nmm_codon* codon, char a, char b, char c)
{
    codon->a = a;
    codon->b = b;
    codon->c = c;
    return codon;
}

static void destroy(struct imm_state const* state)
{
    free_c(__imm_state_derived(state));
    __imm_state_destroy(state);
}

static double joint_seq_len1(struct nmm_frame_state const* state, struct imm_seq const* seq)
{
    char const _ = state->any_symbol;
    NMM_CODON_DECL(codon, nmm_base_table_abc(state->baset));

    double const c = 2 * state->leps + 2 * state->l1eps;

    char const   nucl = imm_seq_string(seq)[0];
    double const e0 = nmm_codon_table_lprob(state->codont, codon_set(&codon, nucl, _, _));
    double const e1 = nmm_codon_table_lprob(state->codont, codon_set(&codon, _, nucl, _));
    double const e2 = nmm_codon_table_lprob(state->codont, codon_set(&codon, _, _, nucl));

    return c + logaddexp3(e0, e1, e2) - log(3);
}

static double joint_seq_len2(struct nmm_frame_state const* state, struct imm_seq const* seq)
{
#define C(a, b, c)                                                                                 \
    nmm_codon_table_lprob(state->codont, codon_set(&codon, eseq[a], eseq[b], eseq[c]))
    char const*  str = imm_seq_string(seq);
    char const   eseq[] = {str[0], str[1], state->any_symbol};
    size_t const _ = sizeof(eseq) - 1;
    NMM_CODON_DECL(codon, nmm_base_table_abc(state->baset));

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
#define C(a, b, c)                                                                                 \
    nmm_codon_table_lprob(state->codont, codon_set(&codon, eseq[a], eseq[b], eseq[c]))
    char const*  str = imm_seq_string(seq);
    char const   eseq[] = {str[0], str[1], str[2], state->any_symbol};
    size_t const _ = sizeof(eseq) - 1;
    NMM_CODON_DECL(codon, nmm_base_table_abc(state->baset));

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
#define C(a, b, c)                                                                                 \
    nmm_codon_table_lprob(state->codont, codon_set(&codon, eseq[a], eseq[b], eseq[c]))
    char const*  str = imm_seq_string(seq);
    char const   eseq[] = {str[0], str[1], str[2], str[3], state->any_symbol};
    size_t const _ = sizeof(eseq) - 1;
    NMM_CODON_DECL(codon, nmm_base_table_abc(state->baset));

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
    NMM_CODON_DECL(codon, nmm_base_table_abc(state->baset));

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

static inline double logaddexp3(double const a, double const b, double const c)
{
    return imm_lprob_add(imm_lprob_add(a, b), c);
}

static double lprob(struct imm_state const* state, struct imm_seq const* seq)
{
    struct nmm_frame_state const* s = __imm_state_derived(state);
    unsigned                      length = imm_seq_length(seq);

    switch (length) {
    case 1:
        return joint_seq_len1(s, seq);
    case 2:
        return joint_seq_len2(s, seq);
    case 3:
        return joint_seq_len3(s, seq);
    case 4:
        return joint_seq_len4(s, seq);
    case 5:
        return joint_seq_len5(s, seq);
    }

    return s->zero_lprob;
}

static double lprob_frag_given_codon1(struct nmm_frame_state const* state,
                                      struct imm_seq const* seq, struct nmm_codon const* codon)
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
                                      struct imm_seq const* seq, struct nmm_codon const* codon)
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
                                      struct imm_seq const* seq, struct nmm_codon const* codon)
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
        c1 + log((x2 == z2) * (x3 == z3) + (x1 == z2) * (x3 == z3) + (x1 == z2) * (x2 == z3)) +
        lprob_z1;

    double const v2 =
        c1 + log((x2 == z1) * (x3 == z3) + (x1 == z1) * (x3 == z3) + (x1 == z1) * (x2 == z3)) +
        lprob_z2;

    double const v3 =
        c1 + log((x2 == z1) * (x3 == z2) + (x1 == z1) * (x3 == z2) + (x1 == z1) * (x2 == z2)) +
        lprob_z3;

    double const c2 = 4 * loge - log(9);

    double const v4 = c2 + log((x1 == z3) + (x2 == z3) + (x3 == z3)) + lprob_z1 + lprob_z2;
    double const v5 = c2 + log((x1 == z2) + (x2 == z2) + (x3 == z2)) + lprob_z1 + lprob_z3;
    double const v6 = c2 + log((x1 == z1) + (x2 == z1) + (x3 == z1)) + lprob_z2 + lprob_z3;

    double v[] = {v0, v1, v2, v3, v4, v5, v6};

    return imm_lprob_sum(v, NMM_ARRAY_SIZE(v));
}

static double lprob_frag_given_codon4(struct nmm_frame_state const* state,
                                      struct imm_seq const* seq, struct nmm_codon const* codon)
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
                                      struct imm_seq const* seq, struct nmm_codon const* codon)
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

static uint8_t max_seq(struct imm_state const* state) { return 5; }

static uint8_t min_seq(struct imm_state const* state) { return 1; }

static uint8_t type_id(struct imm_state const* state) { return NMM_FRAME_STATE_TYPE_ID; }

int nmm_frame_state_write(struct imm_state const* state, struct nmm_model const* model,
                          FILE* stream)
{
    struct nmm_frame_state const* s = nmm_frame_state_derived(state);
    uint32_t                      baset_idx = model_baset_index(model, s->baset);
    uint32_t                      codont_idx = model_codont_index(model, s->codont);

    if (__imm_state_write(state, stream)) {
        imm_error("could not write super state");
        return 1;
    }

    if (fwrite(&baset_idx, sizeof(baset_idx), 1, stream) < 1) {
        imm_error("could not write baset_idx index");
        return 1;
    }

    if (fwrite(&codont_idx, sizeof(codont_idx), 1, stream) < 1) {
        imm_error("could not write codont_idx index");
        return 1;
    }

    if (fwrite(&s->epsilon, sizeof(s->epsilon), 1, stream) < 1) {
        imm_error("could not write epsilon");
        return 1;
    }

    return 0;
}

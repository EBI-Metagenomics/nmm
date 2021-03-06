#include "free.h"
#include "nmm/nmm.h"
#include "profile.h"
#include <stdlib.h>
#include <string.h>

struct nmm_frame_state
{
    struct imm_state const*      super;
    struct nmm_base_lprob const* basep;
    struct nmm_codon_marg const* codonm;
    imm_float                    epsilon;
    imm_float                    leps;
    imm_float                    l1eps;
    imm_float                    zero_lprob;
    char                         any_symbol;
};

static inline imm_float               base_lprob(struct nmm_frame_state const* state, char id);
static inline struct nmm_codon const* codon_set(struct nmm_codon* codon, char a, char b, char c);
static void                           destroy(struct imm_state const* state);
static imm_float                      joint_seq_len1(struct nmm_frame_state const* state, struct imm_seq const* seq);
static imm_float                      joint_seq_len2(struct nmm_frame_state const* state, struct imm_seq const* seq);
static imm_float                      joint_seq_len3(struct nmm_frame_state const* state, struct imm_seq const* seq);
static imm_float                      joint_seq_len4(struct nmm_frame_state const* state, struct imm_seq const* seq);
static imm_float                      joint_seq_len5(struct nmm_frame_state const* state, struct imm_seq const* seq);
static inline imm_float               logaddexp3(imm_float const a, imm_float const b, imm_float const c);
static imm_float                      lprob(struct imm_state const* state, struct imm_seq const* seq);
static imm_float lprob_frag_given_codon1(struct nmm_frame_state const* state, struct imm_seq const* seq,
                                         struct nmm_codon const* ccode);
static imm_float lprob_frag_given_codon2(struct nmm_frame_state const* state, struct imm_seq const* seq,
                                         struct nmm_codon const* ccode);
static imm_float lprob_frag_given_codon3(struct nmm_frame_state const* state, struct imm_seq const* seq,
                                         struct nmm_codon const* ccode);
static imm_float lprob_frag_given_codon4(struct nmm_frame_state const* state, struct imm_seq const* seq,
                                         struct nmm_codon const* ccode);
static imm_float lprob_frag_given_codon5(struct nmm_frame_state const* state, struct imm_seq const* seq,
                                         struct nmm_codon const* ccode);
static uint8_t   max_seq(struct imm_state const* state);
static uint8_t   min_seq(struct imm_state const* state);
static uint8_t   type_id(struct imm_state const* state);

static struct imm_state_vtable const __vtable = {destroy, lprob, max_seq, min_seq, type_id};

struct nmm_base_lprob const* nmm_frame_state_base_lprob(struct nmm_frame_state const* state) { return state->basep; }

struct nmm_codon_marg const* nmm_frame_state_codon_marg(struct nmm_frame_state const* state) { return state->codonm; }

struct nmm_frame_state const* nmm_frame_state_create(char const* name, struct nmm_base_lprob const* basep,
                                                     struct nmm_codon_marg const* codonm, imm_float const epsilon)
{
    struct nmm_frame_state* state = malloc(sizeof(*state));

    if (nmm_base_lprob_abc(basep) != nmm_codon_marg_abc(codonm)) {
        free_c(state);
        return NULL;
    }

    state->basep = basep;
    state->codonm = codonm;
    state->epsilon = epsilon;
    state->leps = imm_log(epsilon);
    state->l1eps = imm_log(1 - epsilon);
    state->zero_lprob = imm_lprob_zero();
    state->any_symbol = imm_abc_any_symbol(nmm_base_abc_super(nmm_base_lprob_abc(basep)));

    struct imm_abc const* abc = nmm_base_abc_super(nmm_base_lprob_abc(basep));
    state->super = imm_state_create(name, abc, __vtable, state);
    return state;
}

imm_float nmm_frame_state_decode(struct nmm_frame_state const* state, struct imm_seq const* seq,
                                 struct nmm_codon* codon)
{
    struct imm_abc const* abc = nmm_base_abc_super(nmm_codon_marg_abc(state->codonm));
    char const*           symbols = imm_abc_symbols(abc);
    unsigned const        n = imm_abc_length(abc);

    imm_float         max_lprob = state->zero_lprob;
    struct nmm_codon* tmp = nmm_codon_create(nmm_codon_marg_abc(state->codonm));

    for (unsigned i0 = 0; i0 < n; ++i0) {
        for (unsigned i1 = 0; i1 < n; ++i1) {
            for (unsigned i2 = 0; i2 < n; ++i2) {

                struct nmm_triplet triplet = {symbols[i0], symbols[i1], symbols[i2]};
                nmm_codon_set_triplet(tmp, triplet);
                imm_float lprob = nmm_frame_state_lposterior(state, tmp, seq);

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

void nmm_frame_state_destroy(struct nmm_frame_state const* state) { state->super->vtable.destroy(state->super); }

imm_float nmm_frame_state_epsilon(struct nmm_frame_state const* state) { return state->epsilon; }

imm_float nmm_frame_state_lposterior(struct nmm_frame_state const* state, struct nmm_codon const* codon,
                                     struct imm_seq const* seq)
{
    imm_float lprob = state->zero_lprob;
    unsigned  length = imm_seq_length(seq);

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

    return lprob + nmm_codon_marg_lprob(state->codonm, codon);
}

struct imm_state const* nmm_frame_state_super(struct nmm_frame_state const* state) { return state->super; }

struct imm_state const* nmm_frame_state_read(FILE* stream, struct nmm_profile const* prof)
{
    struct imm_abc const* abc = nmm_profile_abc(prof);
    struct imm_state*     state = __imm_state_read(stream, abc);
    if (!state) {
        imm_error("could not state_read");
        return NULL;
    }

    state->vtable = __vtable;

    struct nmm_frame_state* frame_state = malloc(sizeof(*frame_state));
    frame_state->super = state;
    state->derived = frame_state;

    uint16_t index = 0;
    if (fread(&index, sizeof(index), 1, stream) < 1) {
        imm_error("could not read basep index");
        goto err;
    }

    struct nmm_base_lprob const* basep = nmm_profile_base_lprob(prof, index);
    if (!basep) {
        imm_error("could not get basep");
        goto err;
    }
    frame_state->basep = basep;

    index = 0;
    if (fread(&index, sizeof(index), 1, stream) < 1) {
        imm_error("could not read codonm index");
        goto err;
    }

    struct nmm_codon_marg const* codonm = nmm_profile_codon_marg(prof, index);
    if (!codonm) {
        imm_error("could not get codonm");
        goto err;
    }
    frame_state->codonm = codonm;

    if (fread(&frame_state->epsilon, sizeof(frame_state->epsilon), 1, stream) < 1) {
        imm_error("could not read epsilon");
        goto err;
    }

    frame_state->leps = imm_log(frame_state->epsilon);
    frame_state->l1eps = imm_log(1 - frame_state->epsilon);
    frame_state->zero_lprob = imm_lprob_zero();
    frame_state->any_symbol = imm_abc_any_symbol(nmm_base_abc_super(nmm_base_lprob_abc(basep)));

    return state;

err:
    free_c(frame_state);
    __imm_state_destroy(state);
    return NULL;
}

static inline imm_float base_lprob(struct nmm_frame_state const* state, char id)
{
    return nmm_base_lprob_get(state->basep, id);
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

static imm_float joint_seq_len1(struct nmm_frame_state const* state, struct imm_seq const* seq)
{
    char const _ = state->any_symbol;
    NMM_CODON_DECL(codon, nmm_base_lprob_abc(state->basep));

    imm_float const c = 2 * state->leps + 2 * state->l1eps;

    char const      nucl = imm_seq_string(seq)[0];
    imm_float const e0 = nmm_codon_marg_lprob(state->codonm, codon_set(&codon, nucl, _, _));
    imm_float const e1 = nmm_codon_marg_lprob(state->codonm, codon_set(&codon, _, nucl, _));
    imm_float const e2 = nmm_codon_marg_lprob(state->codonm, codon_set(&codon, _, _, nucl));

    return c + logaddexp3(e0, e1, e2) - imm_log(3);
}

static imm_float joint_seq_len2(struct nmm_frame_state const* state, struct imm_seq const* seq)
{
#define C(a, b, c) nmm_codon_marg_lprob(state->codonm, codon_set(&codon, eseq[a], eseq[b], eseq[c]))
    char const*  str = imm_seq_string(seq);
    char const   eseq[] = {str[0], str[1], state->any_symbol};
    size_t const _ = sizeof(eseq) - 1;
    NMM_CODON_DECL(codon, nmm_base_lprob_abc(state->basep));

    imm_float const b_lp0 = base_lprob(state, str[0]);
    imm_float const b_lp1 = base_lprob(state, str[1]);

    imm_float const c0 = imm_log(2) + state->leps + state->l1eps * 3 - imm_log(3);
    imm_float const v0 = c0 + logaddexp3(C(_, 0, 1), C(0, _, 1), C(0, 1, _));

    imm_float const c1 = 3 * state->leps + state->l1eps - imm_log(3);
    imm_float const v1 = c1 + logaddexp3(C(0, _, _), C(_, 0, _), C(_, _, 0)) + b_lp1;
    imm_float const v2 = c1 + logaddexp3(C(1, _, _), C(_, 1, _), C(_, _, 1)) + b_lp0;

    return logaddexp3(v0, v1, v2);
#undef C
}

static imm_float joint_seq_len3(struct nmm_frame_state const* state, struct imm_seq const* seq)
{
#define C(a, b, c) nmm_codon_marg_lprob(state->codonm, codon_set(&codon, eseq[a], eseq[b], eseq[c]))
    char const*  str = imm_seq_string(seq);
    char const   eseq[] = {str[0], str[1], str[2], state->any_symbol};
    size_t const _ = sizeof(eseq) - 1;
    NMM_CODON_DECL(codon, nmm_base_lprob_abc(state->basep));

    imm_float const B[] = {base_lprob(state, str[0]), base_lprob(state, str[1]), base_lprob(state, str[2])};

    imm_float const v0 = 4 * state->l1eps + C(0, 1, 2);

    imm_float const c1 = imm_log(4) + 2 * state->leps + 2 * state->l1eps - imm_log(9);
    imm_float const v1[] = {logaddexp3(C(_, 1, 2), C(1, _, 2), C(1, 2, _)) + B[0],
                            logaddexp3(C(_, 0, 2), C(0, _, 2), C(0, 2, _)) + B[1],
                            logaddexp3(C(_, 0, 1), C(0, _, 1), C(0, 1, _)) + B[2]};

    imm_float const c2 = 4 * state->leps - imm_log(9);
    imm_float const v2[] = {logaddexp3(C(2, _, _), C(_, 2, _), C(_, _, 2)) + B[0] + B[1],
                            logaddexp3(C(1, _, _), C(_, 1, _), C(_, _, 1)) + B[0] + B[2],
                            logaddexp3(C(0, _, _), C(_, 0, _), C(_, _, 0)) + B[1] + B[2]};

    return logaddexp3(v0, c1 + imm_lprob_sum(v1, 3), c2 + imm_lprob_sum(v2, 3));
#undef C
}

static imm_float joint_seq_len4(struct nmm_frame_state const* state, struct imm_seq const* seq)
{
#define C(a, b, c) nmm_codon_marg_lprob(state->codonm, codon_set(&codon, eseq[a], eseq[b], eseq[c]))
    char const*  str = imm_seq_string(seq);
    char const   eseq[] = {str[0], str[1], str[2], str[3], state->any_symbol};
    size_t const _ = sizeof(eseq) - 1;
    NMM_CODON_DECL(codon, nmm_base_lprob_abc(state->basep));

    imm_float const B[] = {base_lprob(state, str[0]), base_lprob(state, str[1]), base_lprob(state, str[2]),
                           base_lprob(state, str[3])};

    imm_float const c0 = state->leps + state->l1eps * 3 - imm_log(2);
    imm_float const v0[] = {C(1, 2, 3) + B[0], C(0, 2, 3) + B[1], C(0, 1, 3) + B[2], C(0, 1, 2) + B[3]};

    imm_float const c1 = 3 * state->leps + state->l1eps - imm_log(9);
    imm_float const v1[] = {C(_, 2, 3) + B[0] + B[1], C(_, 1, 3) + B[0] + B[2], C(_, 1, 2) + B[0] + B[3],
                            C(_, 0, 3) + B[1] + B[2], C(_, 0, 2) + B[1] + B[3], C(_, 0, 1) + B[2] + B[3],
                            C(2, _, 3) + B[0] + B[1], C(1, _, 3) + B[0] + B[2], C(1, _, 2) + B[0] + B[3],
                            C(0, _, 3) + B[1] + B[2], C(0, _, 2) + B[1] + B[3], C(0, _, 1) + B[2] + B[3],
                            C(2, 3, _) + B[0] + B[1], C(1, 3, _) + B[0] + B[2], C(1, 2, _) + B[0] + B[3],
                            C(0, 3, _) + B[1] + B[2], C(0, 2, _) + B[1] + B[3], C(0, 1, _) + B[2] + B[3]};

    return imm_lprob_add(c0 + imm_lprob_sum(v0, 4), c1 + imm_lprob_sum(v1, 18));
#undef C
}

static imm_float joint_seq_len5(struct nmm_frame_state const* state, struct imm_seq const* seq)
{
#define LPROB(codon) nmm_codon_marg_lprob(state->codonm, codon)
#define C(a, b, c) codon_set(&codon, str[a], str[b], str[c])
    char const* str = imm_seq_string(seq);
    NMM_CODON_DECL(codon, nmm_base_lprob_abc(state->basep));

    imm_float const b_lp0 = base_lprob(state, str[0]);
    imm_float const b_lp1 = base_lprob(state, str[1]);
    imm_float const b_lp2 = base_lprob(state, str[2]);
    imm_float const b_lp3 = base_lprob(state, str[3]);
    imm_float const b_lp4 = base_lprob(state, str[4]);

    imm_float const v[] = {imm_lprob_add(b_lp0 + b_lp1 + LPROB(C(2, 3, 4)), b_lp0 + b_lp2 + LPROB(C(1, 3, 4))),
                           imm_lprob_add(b_lp0 + b_lp3 + LPROB(C(1, 2, 4)), b_lp0 + b_lp4 + LPROB(C(1, 2, 3))),
                           imm_lprob_add(b_lp1 + b_lp2 + LPROB(C(0, 3, 4)), b_lp1 + b_lp3 + LPROB(C(0, 2, 4))),
                           imm_lprob_add(b_lp1 + b_lp4 + LPROB(C(0, 2, 3)), b_lp2 + b_lp3 + LPROB(C(0, 1, 4))),
                           imm_lprob_add(b_lp2 + b_lp4 + LPROB(C(0, 1, 3)), b_lp3 + b_lp4 + LPROB(C(0, 1, 2)))};

    return 2 * state->leps + 2 * state->l1eps - imm_log(10) + imm_lprob_sum(v, 5);
#undef C
#undef LPROB
}

static inline imm_float logaddexp3(imm_float const a, imm_float const b, imm_float const c)
{
    return imm_lprob_add(imm_lprob_add(a, b), c);
}

static imm_float lprob(struct imm_state const* state, struct imm_seq const* seq)
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

static imm_float lprob_frag_given_codon1(struct nmm_frame_state const* state, struct imm_seq const* seq,
                                         struct nmm_codon const* codon)
{
    imm_float const loge = state->leps;
    imm_float const log1e = state->l1eps;

    char const x1 = codon->a;
    char const x2 = codon->b;
    char const x3 = codon->c;

    char const z1 = imm_seq_string(seq)[0];

    imm_float const c = 2 * loge + 2 * log1e;

    return c + imm_log((x1 == z1) + (x2 == z1) + (x3 == z1)) - imm_log(3);
}

static imm_float lprob_frag_given_codon2(struct nmm_frame_state const* state, struct imm_seq const* seq,
                                         struct nmm_codon const* codon)
{
    imm_float const loge = state->leps;
    imm_float const log1e = state->l1eps;

    char const x1 = codon->a;
    char const x2 = codon->b;
    char const x3 = codon->c;

    char const z1 = imm_seq_string(seq)[0];
    char const z2 = imm_seq_string(seq)[1];

    imm_float const lprob_z1 = base_lprob(state, z1);
    imm_float const lprob_z2 = base_lprob(state, z2);

    imm_float const c1 = imm_log(2) + loge + log1e * 3 - imm_log(3);
    imm_float const v0 = c1 + imm_log((x2 == z1) * (x3 == z2) + (x1 == z1) * (x3 == z2) + (x1 == z1) * (x2 == z2));

    imm_float const c2 = 3 * loge + log1e - imm_log(3);

    imm_float const v1 = c2 + imm_log((x1 == z1) + (x2 == z1) + (x3 == z1)) + lprob_z2;
    imm_float const v2 = c2 + imm_log((x1 == z2) + (x2 == z2) + (x3 == z2)) + lprob_z1;

    return logaddexp3(v0, v1, v2);
}

static imm_float lprob_frag_given_codon3(struct nmm_frame_state const* state, struct imm_seq const* seq,
                                         struct nmm_codon const* codon)
{
    imm_float const loge = state->leps;
    imm_float const log1e = state->l1eps;

    char const x1 = codon->a;
    char const x2 = codon->b;
    char const x3 = codon->c;

    char const z1 = imm_seq_string(seq)[0];
    char const z2 = imm_seq_string(seq)[1];
    char const z3 = imm_seq_string(seq)[2];

    imm_float const lprob_z1 = base_lprob(state, z1);
    imm_float const lprob_z2 = base_lprob(state, z2);
    imm_float const lprob_z3 = base_lprob(state, z3);

    imm_float const v0 = 4 * log1e + imm_log((x1 == z1) * (x2 == z2) * (x3 == z3));

    imm_float const c1 = imm_log(4) + 2 * loge + 2 * log1e - imm_log(9);

    imm_float const v1 =
        c1 + imm_log((x2 == z2) * (x3 == z3) + (x1 == z2) * (x3 == z3) + (x1 == z2) * (x2 == z3)) + lprob_z1;

    imm_float const v2 =
        c1 + imm_log((x2 == z1) * (x3 == z3) + (x1 == z1) * (x3 == z3) + (x1 == z1) * (x2 == z3)) + lprob_z2;

    imm_float const v3 =
        c1 + imm_log((x2 == z1) * (x3 == z2) + (x1 == z1) * (x3 == z2) + (x1 == z1) * (x2 == z2)) + lprob_z3;

    imm_float const c2 = 4 * loge - imm_log(9);

    imm_float const v4 = c2 + imm_log((x1 == z3) + (x2 == z3) + (x3 == z3)) + lprob_z1 + lprob_z2;
    imm_float const v5 = c2 + imm_log((x1 == z2) + (x2 == z2) + (x3 == z2)) + lprob_z1 + lprob_z3;
    imm_float const v6 = c2 + imm_log((x1 == z1) + (x2 == z1) + (x3 == z1)) + lprob_z2 + lprob_z3;

    imm_float v[] = {v0, v1, v2, v3, v4, v5, v6};

    return imm_lprob_sum(v, IMM_ARRAY_SIZE(v));
}

static imm_float lprob_frag_given_codon4(struct nmm_frame_state const* state, struct imm_seq const* seq,
                                         struct nmm_codon const* codon)
{
    imm_float const loge = state->leps;
    imm_float const log1e = state->l1eps;

    char const x1 = codon->a;
    char const x2 = codon->b;
    char const x3 = codon->c;

    char const z1 = imm_seq_string(seq)[0];
    char const z2 = imm_seq_string(seq)[1];
    char const z3 = imm_seq_string(seq)[2];
    char const z4 = imm_seq_string(seq)[3];

    imm_float const lprob_z1 = base_lprob(state, z1);
    imm_float const lprob_z2 = base_lprob(state, z2);
    imm_float const lprob_z3 = base_lprob(state, z3);
    imm_float const lprob_z4 = base_lprob(state, z4);

    imm_float const v0[] = {imm_log((x1 == z2) * (x2 == z3) * (x3 == z4)) + lprob_z1,
                            imm_log((x1 == z1) * (x2 == z3) * (x3 == z4)) + lprob_z2,
                            imm_log((x1 == z1) * (x2 == z2) * (x3 == z4)) + lprob_z3,
                            imm_log((x1 == z1) * (x2 == z2) * (x3 == z3)) + lprob_z4};

    imm_float const v1[] = {
        imm_log((x2 == z3) * (x3 == z4)) + lprob_z1 + lprob_z2, imm_log((x2 == z2) * (x3 == z4)) + lprob_z1 + lprob_z3,
        imm_log((x2 == z2) * (x3 == z3)) + lprob_z1 + lprob_z4, imm_log((x2 == z1) * (x3 == z4)) + lprob_z2 + lprob_z3,
        imm_log((x2 == z1) * (x3 == z3)) + lprob_z2 + lprob_z4, imm_log((x2 == z1) * (x3 == z2)) + lprob_z3 + lprob_z4,
        imm_log((x1 == z3) * (x3 == z4)) + lprob_z1 + lprob_z2, imm_log((x1 == z2) * (x3 == z4)) + lprob_z1 + lprob_z3,
        imm_log((x1 == z2) * (x3 == z3)) + lprob_z1 + lprob_z4, imm_log((x1 == z1) * (x3 == z4)) + lprob_z2 + lprob_z3,
        imm_log((x1 == z1) * (x3 == z3)) + lprob_z2 + lprob_z4, imm_log((x1 == z1) * (x3 == z2)) + lprob_z3 + lprob_z4,
        imm_log((x1 == z3) * (x2 == z4)) + lprob_z1 + lprob_z2, imm_log((x1 == z2) * (x2 == z4)) + lprob_z1 + lprob_z3,
        imm_log((x1 == z2) * (x2 == z3)) + lprob_z1 + lprob_z4, imm_log((x1 == z1) * (x2 == z4)) + lprob_z2 + lprob_z3,
        imm_log((x1 == z1) * (x2 == z3)) + lprob_z2 + lprob_z4, imm_log((x1 == z1) * (x2 == z2)) + lprob_z3 + lprob_z4,
    };

    return imm_lprob_add(loge + log1e * 3 - imm_log(2) + imm_lprob_sum(v0, IMM_ARRAY_SIZE(v0)),
                         3 * loge + log1e - imm_log(9) + imm_lprob_sum(v1, IMM_ARRAY_SIZE(v1)));
}

static imm_float lprob_frag_given_codon5(struct nmm_frame_state const* state, struct imm_seq const* seq,
                                         struct nmm_codon const* codon)
{
    imm_float const loge = state->leps;
    imm_float const log1e = state->l1eps;

    char const x1 = codon->a;
    char const x2 = codon->b;
    char const x3 = codon->c;

    char const z1 = imm_seq_string(seq)[0];
    char const z2 = imm_seq_string(seq)[1];
    char const z3 = imm_seq_string(seq)[2];
    char const z4 = imm_seq_string(seq)[3];
    char const z5 = imm_seq_string(seq)[4];

    imm_float const lprob_z1 = base_lprob(state, z1);
    imm_float const lprob_z2 = base_lprob(state, z2);
    imm_float const lprob_z3 = base_lprob(state, z3);
    imm_float const lprob_z4 = base_lprob(state, z4);
    imm_float const lprob_z5 = base_lprob(state, z5);

    imm_float const v[] = {
        lprob_z1 + lprob_z2 + imm_log((x1 == z3) * (x2 == z4) * (x3 == z5)),
        lprob_z1 + lprob_z3 + imm_log((x1 == z2) * (x2 == z4) * (x3 == z5)),
        lprob_z1 + lprob_z4 + imm_log((x1 == z2) * (x2 == z3) * (x3 == z5)),
        lprob_z1 + lprob_z5 + imm_log((x1 == z2) * (x2 == z3) * (x3 == z4)),
        lprob_z2 + lprob_z3 + imm_log((x1 == z1) * (x2 == z4) * (x3 == z5)),
        lprob_z2 + lprob_z4 + imm_log((x1 == z1) * (x2 == z3) * (x3 == z5)),
        lprob_z2 + lprob_z5 + imm_log((x1 == z1) * (x2 == z3) * (x3 == z4)),
        lprob_z3 + lprob_z4 + imm_log((x1 == z1) * (x2 == z2) * (x3 == z5)),
        lprob_z3 + lprob_z5 + imm_log((x1 == z1) * (x2 == z2) * (x3 == z4)),
        lprob_z4 + lprob_z5 + imm_log((x1 == z1) * (x2 == z2) * (x3 == z3)),
    };

    return 2 * loge + 2 * log1e - imm_log(10) + imm_lprob_sum(v, IMM_ARRAY_SIZE(v));
}

static uint8_t max_seq(struct imm_state const* state) { return 5; }

static uint8_t min_seq(struct imm_state const* state) { return 1; }

static uint8_t type_id(struct imm_state const* state) { return NMM_FRAME_STATE_TYPE_ID; }

int nmm_frame_state_write(struct imm_state const* state, struct nmm_profile const* prof, FILE* stream)
{
    struct nmm_frame_state const* s = nmm_frame_state_derived(state);
    uint16_t                      basep_idx = profile_base_lprob_index(prof, s->basep);
    uint16_t                      codonm_idx = profile_codon_marg_index(prof, s->codonm);

    if (__imm_state_write(state, stream)) {
        imm_error("could not write super state");
        return 1;
    }

    if (fwrite(&basep_idx, sizeof(basep_idx), 1, stream) < 1) {
        imm_error("could not write basep_idx index");
        return 1;
    }

    if (fwrite(&codonm_idx, sizeof(codonm_idx), 1, stream) < 1) {
        imm_error("could not write codonm_idx index");
        return 1;
    }

    if (fwrite(&s->epsilon, sizeof(s->epsilon), 1, stream) < 1) {
        imm_error("could not write epsilon");
        return 1;
    }

    return 0;
}

#include "imm.h"
#include "nmm.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

struct nmm_frame_state
{
    struct imm_state *interface;
    struct nmm_base const *base;
    const struct nmm_codon *codon;
    double epsilon;
    double leps;
    double l1eps;
};

double frame_state_lprob(const struct imm_state *state, const char *seq, int seq_len);
int frame_state_min_seq(const struct imm_state *state);
int frame_state_max_seq(const struct imm_state *state);

struct nmm_frame_state *nmm_frame_state_create(const char *name,
                                               const struct nmm_base *base,
                                               const struct nmm_codon *codon,
                                               double epsilon)
{
    struct nmm_frame_state *state = malloc(sizeof(struct nmm_frame_state));

    if (nmm_base_get_abc(base) != nmm_codon_get_abc(codon)) {
        free(state);
        imm_error("alphabets from base and codon are different");
        return NULL;
    }

    if (imm_abc_length(nmm_base_get_abc(base)) != 4) {
        free(state);
        imm_error("alphabet length is not four");
        return NULL;
    }

    state->base = base;
    state->codon = codon;
    state->epsilon = epsilon;
    state->leps = log(epsilon);
    state->l1eps = log(1 - epsilon);

    struct imm_state_funcs funcs = {frame_state_lprob, frame_state_min_seq,
                                    frame_state_max_seq};

    state->interface = imm_state_create(name, nmm_base_get_abc(base), funcs, state);
    return state;
}

double joint_seq_len1(const struct nmm_frame_state *state, const char *seq);
double joint_seq_len2(const struct nmm_frame_state *state, const char *seq);
double joint_seq_len3(const struct nmm_frame_state *state, const char *seq);
double joint_seq_len4(const struct nmm_frame_state *state, const char *seq);
double joint_seq_len5(const struct nmm_frame_state *state, const char *seq);
double codon_lprob(const struct nmm_frame_state *state, const char *codon);
double base_lprob(const struct nmm_frame_state *state, char id);
static inline double logaddexp(double a, double b) { return imm_logaddexp(a, b); }
static inline double logaddexp3(double a, double b, double c)
{
    return logaddexp(logaddexp(a, b), c);
}
static inline double logsumexp(double *arr, int len) { return imm_logsumexp(arr, len); }
static inline double ecodon_lprob(const struct nmm_frame_state *state, const char *seq,
                                  int a, int b, int c)
{
    const char codon[3] = {seq[a], seq[b], seq[c]};
    return codon_lprob(state, codon);
}

double frame_state_lprob(const struct imm_state *state, const char *seq, int seq_len)
{
    const struct nmm_frame_state *s = imm_state_get_impl_c(state);
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

    return -INFINITY;
}

int frame_state_min_seq(const struct imm_state *state) { return 1; }

int frame_state_max_seq(const struct imm_state *state) { return 5; }

void nmm_frame_state_destroy(struct nmm_frame_state *state)
{
    if (!state)
        return;

    imm_state_destroy(state->interface);
    state->interface = NULL;

    state->base = NULL;
    state->codon = NULL;
    free(state);
}

double joint_seq_len1(const struct nmm_frame_state *state, const char *seq)
{
    const char _ = IMM_ANY_SYMBOL;
    const char c0__[3] = {seq[0], _, _};
    const char c_0_[3] = {_, seq[0], _};
    const char c__0[3] = {_, _, seq[0]};

    double c = 2 * state->leps + 2 * state->l1eps;

    double e0 = codon_lprob(state, c0__);
    double e1 = codon_lprob(state, c_0_);
    double e2 = codon_lprob(state, c__0);
    return c + logaddexp3(e0, e1, e2) - log(3);
}

double joint_seq_len2(const struct nmm_frame_state *state, const char *seq)
{
#define c_lprob(codon) codon_lprob(state, codon)
    const char _ = IMM_ANY_SYMBOL;

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

double joint_seq_len3(const struct nmm_frame_state *state, const char *seq)
{
#define C(a, b, c) ecodon_lprob(state, eseq, a, b, c)
    const char eseq[] = {seq[0], seq[1], seq[2], IMM_ANY_SYMBOL};
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

double joint_seq_len4(const struct nmm_frame_state *state, const char *seq)
{
#define C(a, b, c) ecodon_lprob(state, eseq, a, b, c)
    const char eseq[] = {seq[0], seq[1], seq[2], seq[3], IMM_ANY_SYMBOL};
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

double joint_seq_len5(const struct nmm_frame_state *state, const char *seq)
{
#define c_lp(codon) codon_lprob(state, codon)
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

double codon_lprob(const struct nmm_frame_state *state, const char *codon)
{
    const struct imm_abc *abc = imm_state_get_abc(imm_state_cast_c(state));
    double lprob = -INFINITY;
    char bases[3 * 4] = {
        imm_abc_symbol_id(abc, 0), imm_abc_symbol_id(abc, 1), imm_abc_symbol_id(abc, 2),
        imm_abc_symbol_id(abc, 3), imm_abc_symbol_id(abc, 0), imm_abc_symbol_id(abc, 1),
        imm_abc_symbol_id(abc, 2), imm_abc_symbol_id(abc, 3), imm_abc_symbol_id(abc, 0),
        imm_abc_symbol_id(abc, 1), imm_abc_symbol_id(abc, 2), imm_abc_symbol_id(abc, 3),
    };
    int nbases[3] = {4, 4, 4};

    for (int i = 0; i < 3; ++i) {
        if (codon[i] != IMM_ANY_SYMBOL) {
            bases[i * 4] = codon[i];
            nbases[i] = 1;
        }
    }

    const char *a_id = bases;
    const char *b_id = bases + 4;
    const char *c_id = bases + 8;
    for (int a = 0; a < nbases[0]; ++a) {
        for (int b = 0; b < nbases[1]; ++b) {
            for (int c = 0; c < nbases[2]; ++c) {
                double t = nmm_codon_get_lprob(state->codon, a_id[a], b_id[b], c_id[c]);
                lprob = logaddexp(lprob, t);
            }
        }
    }

    return lprob;
}

double base_lprob(const struct nmm_frame_state *state, char id)
{
    return nmm_base_get_lprob(state->base, id);
}

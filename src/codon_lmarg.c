#include "codon_lmarg.h"
#include "codon_iter.h"
#include "logaddexp.h"

/**
 * The symbols here are the standard alphabet symbols plus the special any-symbol one.
 */
#define NSYMBOLS (NMM_CODON_NBASES + 1)

static void set_symbol_index(int* symbol_idx, struct imm_abc const* abc);
static int  set_nonmarginal_lprobs(struct codon_lmarg* lmarg, struct nmm_codonp const* lprob);
static void set_marginal_lprobs(struct codon_lmarg* lmarg, struct imm_abc const* abc);

struct codon_lmarg* nmm_codon_lmarg_create(struct imm_abc const*    abc,
                                           struct nmm_codonp const* lprob)
{
    struct codon_lmarg* lmarg = malloc(sizeof(struct codon_lmarg));

    set_symbol_index(lmarg->symbol_idx, abc);

    static int const n = NSYMBOLS;
    lmarg->lprobs = array3d_create(n, n, n);

    if (set_nonmarginal_lprobs(lmarg, lprob)) {
        array3d_destroy(lmarg->lprobs);
        free(lmarg);
        return NULL;
    }

    set_marginal_lprobs(lmarg, abc);

    return lmarg;
}

void nmm_codon_lmarg_destroy(struct codon_lmarg* codon_lprob)
{
    if (!codon_lprob) {
        imm_error("codon_lprob should not be NULL");
        return;
    }

    array3d_destroy(codon_lprob->lprobs);
    free(codon_lprob);
}

static void set_symbol_index(int* symbol_idx, const struct imm_abc* abc)
{
    for (int i = 0; i <= ASCII_LAST_STD; ++i)
        symbol_idx[i] = -1;

    char const* symbols = imm_abc_symbols(abc);

    for (int i = 0; i < NMM_CODON_NBASES; ++i)
        symbol_idx[(size_t)symbols[i]] = imm_abc_symbol_idx(abc, symbols[i]);

    symbol_idx[(size_t)imm_abc_any_symbol(abc)] = NMM_CODON_NBASES;
}

static inline void codon_lmarg_set(struct codon_lmarg*     codon_lprob,
                                   struct nmm_codon const* codon, double lprob)
{

    int const* symbol_idx = codon_lprob->symbol_idx;
    int const  dim[3] = {symbol_idx[(size_t)codon->a], symbol_idx[(size_t)codon->b],
                        symbol_idx[(size_t)codon->c]};
    array3d_set(&codon_lprob->lprobs, dim[0], dim[1], dim[2], lprob);
}

static int set_nonmarginal_lprobs(struct codon_lmarg* lmarg, struct nmm_codonp const* lprob)
{
    struct codon_iter iter = codon_iter_begin(nmm_codonp_get_abc(lprob));
    while (!codon_iter_end(&iter)) {
        struct nmm_codon const codon = codon_iter_next(&iter);
        codon_lmarg_set(lmarg, &codon, nmm_codonp_get(lprob, &codon));
    }
    return 0;
}

static inline int not_marginal(struct nmm_codon const* codon, char const any_symbol)
{
    return codon->a != any_symbol && codon->b != any_symbol && codon->c != any_symbol;
}

static double marginalization(struct codon_lmarg const* lmarg, char const* symbols,
                              struct nmm_codon const* codon);

static void set_marginal_lprobs(struct codon_lmarg* lmarg, struct imm_abc const* abc)
{
    char       any_symbol = imm_abc_any_symbol(abc);
    char const symbols[5] = {imm_abc_symbol_id(abc, 0), imm_abc_symbol_id(abc, 1),
                             imm_abc_symbol_id(abc, 2), imm_abc_symbol_id(abc, 3),
                             any_symbol};

    for (int i0 = 0; i0 < NSYMBOLS; ++i0) {
        for (int i1 = 0; i1 < NSYMBOLS; ++i1) {
            for (int i2 = 0; i2 < NSYMBOLS; ++i2) {

                struct nmm_codon const codon = {symbols[i0], symbols[i1], symbols[i2]};

                if (not_marginal(&codon, any_symbol))
                    continue;

                codon_lmarg_set(lmarg, &codon, marginalization(lmarg, symbols, &codon));
            }
        }
    }
}

static double marginalization(struct codon_lmarg const* lmarg, char const* symbols,
                              struct nmm_codon const* codon)
{
    char const any_symbol = symbols[NSYMBOLS - 1];

    char const seq[3] = {codon->a, codon->b, codon->c};

    char const* arr[3];
    int         shape[3];
    for (int i = 0; i < 3; ++i) {
        if (seq[i] == any_symbol) {
            arr[i] = symbols;
            shape[i] = NMM_CODON_NBASES;
        } else {
            arr[i] = seq + i;
            shape[i] = 1;
        }
    }

    double lprob = imm_lprob_zero();
    for (int a = 0; a < shape[0]; ++a) {
        for (int b = 0; b < shape[1]; ++b) {
            for (int c = 0; c < shape[2]; ++c) {

                struct nmm_codon const t = {arr[0][a], arr[1][b], arr[2][c]};

                lprob = logaddexp(lprob, codon_lmarg_get(lmarg, &t));
            }
        }
    }

    return lprob;
}

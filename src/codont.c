#include "array.h"
#include "codon_iter.h"
#include "free.h"
#include "imm/imm.h"
#include "logaddexp.h"
#include "nmm/base.h"
#include "nmm/nmm.h"
#include <math.h>
#include <stdlib.h>

/**
 * The symbols here are the standard alphabet symbols plus the special any-symbol one.
 */
#define NSYMBOLS (NMM_BASE_SIZE + 1)

#define ASCII_LAST_STD 127

struct nmm_codont
{
    struct nmm_base const* base;
    /**
     * Maps alphabet symbols to the indices 0, 1, 2, and 3 and the any-symbol to 4.
     */
    int symbol_idx[ASCII_LAST_STD + 1];
    /**
     * Pre-computed marginalization forms of
     * p(𝑋₁=𝚡₁,𝑋₂=𝚡₂,𝑋₃=𝚡₃).
     */
    struct array3d lprobs;
};

static void set_symbol_index(int* symbol_idx, struct imm_abc const* abc);
static int set_nonmarginal_lprobs(struct nmm_codont* codont, struct nmm_codonp const* codonp);
static void set_marginal_lprobs(struct nmm_codont* codont, struct nmm_base const* base);

struct nmm_codont const* nmm_codont_create(struct nmm_codonp const* codonp)
{
    struct nmm_codont* codont = malloc(sizeof(struct nmm_codont));
    codont->base = nmm_codonp_get_base(codonp);

    struct imm_abc const* abc = nmm_base_get_abc(codont->base);
    set_symbol_index(codont->symbol_idx, abc);

    codont->lprobs = array3d_create(NSYMBOLS, NSYMBOLS, NSYMBOLS);

    if (set_nonmarginal_lprobs(codont, codonp)) {
        array3d_destroy(codont->lprobs);
        free_c(codont);
        return NULL;
    }

    set_marginal_lprobs(codont, codont->base);

    return codont;
}

/**
 * Calculate any of the marginalization forms of
 * p(𝑋₁=𝚡₁,𝑋₂=𝚡₂,𝑋₃=𝚡₃).
 *
 * The alphabet's any-symbol can be passed to @codon to perform marginalization over
 * the corresponding random variable. Let `"ACGT"` be set of nucleotides and let `'X`'
 * be the any-symbol of the given alphabet. The code
 *
 *     codon_lprob_get(codon_lprob, "AXG")
 *
 * will evaluate the probability p(𝑋₁=𝙰,𝑋₃=𝙶).
 */
static inline double codon_lmarg_get(struct nmm_codont const* codont,
                                     struct nmm_codon const*  codon)
{
    int const* symbol_idx = codont->symbol_idx;
    char       a, b, c;
    nmm_codon_get(codon, &a, &b, &c);

    int const dim[3] = {symbol_idx[(size_t)a], symbol_idx[(size_t)b], symbol_idx[(size_t)c]};

    return array3d_get(&codont->lprobs, (unsigned)dim[0], (unsigned)dim[1], (unsigned)dim[2]);
}

double nmm_codont_lprob(struct nmm_codont const* codont, struct nmm_codon const* codon)
{
    return codon_lmarg_get(codont, codon);
}

void nmm_codont_destroy(struct nmm_codont const* codont)
{
    array3d_destroy(codont->lprobs);
    free_c(codont);
}

struct nmm_base const* nmm_codont_get_base(struct nmm_codont const* codont)
{
    return codont->base;
}

static void set_symbol_index(int* symbol_idx, const struct imm_abc* abc)
{
    for (int i = 0; i <= ASCII_LAST_STD; ++i)
        symbol_idx[i] = -1;

    char const* symbols = imm_abc_symbols(abc);

    for (int i = 0; i < NMM_BASE_SIZE; ++i)
        symbol_idx[(size_t)symbols[i]] = imm_abc_symbol_idx(abc, symbols[i]);

    symbol_idx[(size_t)imm_abc_any_symbol(abc)] = NMM_BASE_SIZE;
}

static inline void codon_lmarg_set(struct nmm_codont* codont, struct nmm_codon const* codon,
                                   double lprob)
{
    int const* symbol_idx = codont->symbol_idx;
    char       a, b, c;
    nmm_codon_get(codon, &a, &b, &c);

    int const dim[3] = {symbol_idx[(size_t)a], symbol_idx[(size_t)b], symbol_idx[(size_t)c]};

    array3d_set(&codont->lprobs, (unsigned)dim[0], (unsigned)dim[1], (unsigned)dim[2], lprob);
}

static int set_nonmarginal_lprobs(struct nmm_codont* codont, struct nmm_codonp const* codonp)
{
    struct codon_iter iter = codon_iter_begin(nmm_codonp_get_base(codonp));
    while (!codon_iter_end(iter)) {
        struct nmm_codon const* codon = codon_iter_next(&iter);
        codon_lmarg_set(codont, codon, nmm_codonp_get(codonp, codon));
    }
    codon_iter_destroy(iter);
    return 0;
}

static double marginalization(struct nmm_codont const* codont, char const* symbols,
                              struct nmm_codon const* codon);

static inline int not_marginal(struct nmm_codon const* codon, char const any_symbol)
{
    char a, b, c;
    nmm_codon_get(codon, &a, &b, &c);
    return a != any_symbol && b != any_symbol && c != any_symbol;
}

static void set_marginal_lprobs(struct nmm_codont* codont, struct nmm_base const* base)
{
    struct imm_abc const * abc = nmm_base_get_abc(base);
    char       any_symbol = imm_abc_any_symbol(abc);
    char const symbols[5] = {imm_abc_symbol_id(abc, 0), imm_abc_symbol_id(abc, 1),
                             imm_abc_symbol_id(abc, 2), imm_abc_symbol_id(abc, 3),
                             any_symbol};

    struct nmm_codon* codon = nmm_codon_create(base);

    for (int i0 = 0; i0 < NSYMBOLS; ++i0) {
        for (int i1 = 0; i1 < NSYMBOLS; ++i1) {
            for (int i2 = 0; i2 < NSYMBOLS; ++i2) {

                nmm_codon_set(codon, symbols[i0], symbols[i1], symbols[i2]);

                if (not_marginal(codon, any_symbol))
                    continue;

                codon_lmarg_set(codont, codon, marginalization(codont, symbols, codon));
            }
        }
    }

    nmm_codon_destroy(codon);
}

static double marginalization(struct nmm_codont const* codont, char const* symbols,
                              struct nmm_codon const* codon)
{
    char const any_symbol = symbols[NSYMBOLS - 1];

    char na, nb, nc;
    nmm_codon_get(codon, &na, &nb, &nc);
    char const seq[3] = {na, nb, nc};

    char const* arr[3];
    int         shape[3];
    for (int i = 0; i < 3; ++i) {
        if (seq[i] == any_symbol) {
            arr[i] = symbols;
            shape[i] = NMM_BASE_SIZE;
        } else {
            arr[i] = seq + i;
            shape[i] = 1;
        }
    }

    struct nmm_codon * tmp = nmm_codon_create(nmm_codon_get_base(codon));
    double lprob = imm_lprob_zero();
    for (int a = 0; a < shape[0]; ++a) {
        for (int b = 0; b < shape[1]; ++b) {
            for (int c = 0; c < shape[2]; ++c) {

                nmm_codon_set(tmp, arr[0][a], arr[1][b], arr[2][c]);

                lprob = logaddexp(lprob, codon_lmarg_get(codont, tmp));
            }
        }
    }
    nmm_codon_destroy(tmp);

    return lprob;
}

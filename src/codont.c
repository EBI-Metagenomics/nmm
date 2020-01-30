#include "nmm/codont.h"
#include "codon_iter.h"
#include "codont_static.h"
#include "free.h"
#include "imm/imm.h"
#include "logaddexp.h"
#include "nmm/base.h"
#include "nmm/codonp.h"
#include <stdlib.h>

/**
 * Number of bases (four) plus the special any-symbol one.
 */
#define NSYMBOLS (NMM_BASE_SIZE + 1)

static void set_symbol_index(struct nmm_codont* codont);
static int set_nonmarginal_lprobs(struct nmm_codont* codont, struct nmm_codonp const* codonp);
static void set_marginal_lprobs(struct nmm_codont* codont, struct nmm_base const* base);

static double marginalization(struct nmm_codont const* codont, char const* symbols,
                              struct nmm_codon const* codon);

static inline int not_marginal(struct nmm_codon const* codon, char const any_symbol)
{
    struct nmm_triplet const t = nmm_codon_get(codon);
    return t.a != any_symbol && t.b != any_symbol && t.c != any_symbol;
}

struct nmm_codont const* nmm_codont_create(struct nmm_codonp const* codonp)
{
    struct nmm_codont* codont = malloc(sizeof(struct nmm_codont));
    codont->base = nmm_codonp_get_base(codonp);

    set_symbol_index(codont);

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
double nmm_codont_lprob(struct nmm_codont const* codont, struct nmm_codon const* codon)
{
    return codont_lprob(codont, codon);
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

static void set_symbol_index(struct nmm_codont* codont)
{
    struct imm_abc const* abc = nmm_base_get_abc(codont->base);

    char const* symbols = imm_abc_symbols(abc);

    for (unsigned i = 0; i < NMM_BASE_SIZE; ++i) {
        size_t j = (size_t)symbols[i];
        codont->symbol_idx[j] = (unsigned)imm_abc_symbol_idx(abc, symbols[i]);
    }

    codont->symbol_idx[(size_t)imm_abc_any_symbol(abc)] = NMM_BASE_SIZE;
}

static inline void set_marginal_lprob(struct nmm_codont*      codont,
                                      struct nmm_codon const* codon, double lprob)
{
    array3d_set(&codont->lprobs, codont_get_array_idx(codont, codon), lprob);
}

static int set_nonmarginal_lprobs(struct nmm_codont* codont, struct nmm_codonp const* codonp)
{
    struct codon_iter iter = codon_iter_begin(codont->base);
    while (!codon_iter_end(iter)) {
        struct nmm_codon const* codon = codon_iter_next(&iter);
        set_marginal_lprob(codont, codon, nmm_codonp_get(codonp, codon));
    }
    codon_iter_destroy(iter);
    return 0;
}

static void set_marginal_lprobs(struct nmm_codont* codont, struct nmm_base const* base)
{
    struct imm_abc const* abc = nmm_base_get_abc(base);
    char const            any_symbol = imm_abc_any_symbol(abc);
    char const            symbols[5] = {imm_abc_symbol_id(abc, 0), imm_abc_symbol_id(abc, 1),
                             imm_abc_symbol_id(abc, 2), imm_abc_symbol_id(abc, 3),
                             any_symbol};

    struct nmm_codon* codon = nmm_codon_create(base);

    for (unsigned i0 = 0; i0 < NSYMBOLS; ++i0) {
        for (unsigned i1 = 0; i1 < NSYMBOLS; ++i1) {
            for (unsigned i2 = 0; i2 < NSYMBOLS; ++i2) {

                struct nmm_triplet const t = {symbols[i0], symbols[i1], symbols[i2]};
                nmm_codon_set(codon, t);

                if (not_marginal(codon, any_symbol))
                    continue;

                set_marginal_lprob(codont, codon, marginalization(codont, symbols, codon));
            }
        }
    }

    nmm_codon_destroy(codon);
}

static double marginalization(struct nmm_codont const* codont, char const* symbols,
                              struct nmm_codon const* codon)
{
    char const any_symbol = symbols[NSYMBOLS - 1];

    struct nmm_triplet t = nmm_codon_get(codon);
    char const         seq[3] = {t.a, t.b, t.c};

    char const* arr[3];
    unsigned    shape[3];
    for (unsigned i = 0; i < 3; ++i) {
        if (seq[i] == any_symbol) {
            arr[i] = symbols;
            shape[i] = NMM_BASE_SIZE;
        } else {
            arr[i] = seq + i;
            shape[i] = 1;
        }
    }

    struct nmm_codon* tmp = nmm_codon_create(nmm_codon_get_base(codon));
    double            lprob = imm_lprob_zero();
    for (unsigned a = 0; a < shape[0]; ++a) {
        for (unsigned b = 0; b < shape[1]; ++b) {
            for (unsigned c = 0; c < shape[2]; ++c) {

                t.a = arr[0][a];
                t.b = arr[1][b];
                t.c = arr[2][c];
                nmm_codon_set(tmp, t);
                lprob = logaddexp(lprob, nmm_codont_lprob(codont, tmp));
            }
        }
    }
    nmm_codon_destroy(tmp);

    return lprob;
}

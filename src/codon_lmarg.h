#ifndef NMM_CODON_LMARG_H
#define NMM_CODON_LMARG_H

#include "array.h"
#include "hide.h"
#include "imm/imm.h"
#include "nmm/nmm.h"

#define ASCII_LAST_STD 127

/** @file codon_lprob.h
 * Compute the marginal probability of emitting a codon.
 *
 * Let p(𝑋₁=𝚡₁,𝑋₂=𝚡₂,𝑋₃=𝚡₃) be the probability of emitting codon (𝚡₁,𝚡₂,𝚡₃), where 𝚡ᵢ𝜖𝒜.
 * This modules implements the computation of p(𝑋₁=𝚡₁,𝑋₂=𝚡₂,𝑋₃=𝚡₃) and any of its
 * marginalization forms (e.g., p(𝑋₁=𝚡₁,𝑋₃=𝚡₃) or p(𝑋₃=𝚡₃)).
 */

struct codon_lmarg
{
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

HIDE struct codon_lmarg* nmm_codon_lmarg_create(struct imm_abc const*    abc,
                                                struct nmm_codonp const* lprob);

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
static inline double codon_lmarg_get(struct codon_lmarg const* lmarg,
                                     struct nmm_codon const*   codon)
{
    int const* symbol_idx = lmarg->symbol_idx;
    int const  dim[3] = {symbol_idx[(size_t)codon->a], symbol_idx[(size_t)codon->b],
                        symbol_idx[(size_t)codon->c]};
    return array3d_get(&lmarg->lprobs, dim[0], dim[1], dim[2]);
}
HIDE void nmm_codon_lmarg_destroy(struct codon_lmarg* lmarg);

#endif

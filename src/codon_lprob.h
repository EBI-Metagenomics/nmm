#ifndef NMM_CODON_LPROB_H
#define NMM_CODON_LPROB_H

#include "hide.h"
#include "imm/imm.h"
#include "nmm/nmm.h"
#include <stddef.h>

#define NSYMBOLS (NMM_CODON_NBASES + 1)
#define ASCII_LAST_STD 127

/**
 * Compute the probability of emitting a codon.
 *
 * Let p(𝑋₁=𝚡₁,𝑋₂=𝚡₂,𝑋₃=𝚡₃) be the probability of emitting codon (𝚡₁,𝚡₂,𝚡₃), where 𝚡ᵢ𝜖𝒜.
 * This modules implements the computation of p(𝑋₁=𝚡₁,𝑋₂=𝚡₂,𝑋₃=𝚡₃ and any of its
 * marginalization forms (e.g., p(𝑋₁=𝚡₁,𝑋₃=𝚡₃) or p(𝑋₃=𝚡₃)).
 */

struct codon_lprob
{
    /**
     * Pre-computed marginalization forms of
     * p(𝑋₁=𝚡₁,𝑋₂=𝚡₂,𝑋₃=𝚡₃).
     */
    double lprob[NSYMBOLS * NSYMBOLS * NSYMBOLS];
    /**
     * Maps alphabet symbols to the indices 0, 1, 2, and 3 and the any-symbol to 4.
     */
    int symbol_idx[ASCII_LAST_STD + 1];
};

HIDE void nmm_codon_lprob_init(struct codon_lprob*      codon_lprob,
                               struct nmm_codont const* codont);

/**
 * Maps codon (containing any-symbol or not) to a position in @ref codon_lprob.symbol_idx.
 */
static inline int codon_lprob_idx(struct codon_lprob const* codon_lprob, char const* codon)
{
    size_t const     i[3] = {(size_t)codon[0], (size_t)codon[1], (size_t)codon[2]};
    int const*       symbol_idx = codon_lprob->symbol_idx;
    int const        dims[3] = {symbol_idx[i[0]], symbol_idx[i[1]], symbol_idx[i[2]]};
    static int const strides[3] = {1, NSYMBOLS, NSYMBOLS * NSYMBOLS};

    return strides[0] * dims[0] + strides[1] * dims[1] + strides[2] * dims[2];
}

/**
 * Calculate any of the marginalization forms of
 * p(𝑋₁=𝚡₁,𝑋₂=𝚡₂,𝑋₃=𝚡₃).
 *
 * The alphabet's any-symbol can be passed to @codon to perform marginalization over
 * the corresponding random variable. Let `"ACGT"` be set of nucleotides and let `'X`'
 * be the any-symbol of the given alphabet. The code
 *
 *     codon_lprob_value(codon_lprob, "AXG")
 *
 * will evaluate the probability p(𝑋₁=𝙰,𝑋₃=𝙶).
 */
static inline double codon_lprob_value(struct codon_lprob const* codon_lprob,
                                       char const*               codon)
{
    return codon_lprob->lprob[codon_lprob_idx(codon_lprob, codon)];
}

#endif

#ifndef NMM_CODONT_H
#define NMM_CODONT_H

#include "nmm/api.h"
#include "nmm/array.h"
#include "nmm/codon.h"

/** @file codont.h
 * Codon table module.
 *
 * A codon table is represented by an (immutable) object of type @ref nmm_codont
 * and is used to compute the marginalization forms of p(𝑋₁=𝚡₁,𝑋₂=𝚡₂,𝑋₃=𝚡₃),
 * the probability of emitting codon (𝚡₁,𝚡₂,𝚡₃).
 */

struct nmm_codonp;

#define NMM_ASCII_LAST_STD 127

struct nmm_codont
{
    struct nmm_base const* base;
    /**
     * Maps alphabet symbols to the indices 0, 1, 2, and 3 and the any-symbol to 4.
     */
    unsigned symbol_idx[NMM_ASCII_LAST_STD + 1];
    /**
     * Pre-computed marginalization forms of
     * p(𝑋₁=𝚡₁,𝑋₂=𝚡₂,𝑋₃=𝚡₃).
     */
    struct nmm_array3d lprobs;
};

NMM_API static inline struct nmm_array3d_idx __nmm_codont_get_array_idx(
    struct nmm_codont const* codont, struct nmm_codon const* codon)
{
    struct nmm_triplet const t = nmm_codon_get_triplet(codon);
    return (struct nmm_array3d_idx){codont->symbol_idx[(size_t)t.a],
                                    codont->symbol_idx[(size_t)t.b],
                                    codont->symbol_idx[(size_t)t.c]};
}

NMM_API struct nmm_codont const* nmm_codont_create(struct nmm_codonp const* codonp);

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
NMM_API static inline double nmm_codont_lprob(struct nmm_codont const* codont,
                                              struct nmm_codon const*  codon)
{
    return nmm_array3d_get(&codont->lprobs, __nmm_codont_get_array_idx(codont, codon));
}

NMM_API void nmm_codont_destroy(struct nmm_codont const* codont);

NMM_API static inline struct nmm_base const* nmm_codont_get_base(
    struct nmm_codont const* codont)
{
    return codont->base;
}

#endif

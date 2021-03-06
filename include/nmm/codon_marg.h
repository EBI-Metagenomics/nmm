#ifndef NMM_CODON_MARG_H
#define NMM_CODON_MARG_H

#include "nmm/array3d.h"
#include "nmm/codon.h"
#include "nmm/export.h"

/** @file codon_marg.h
 * Codon marginalization module.
 *
 * A codon marginalization is represented by an (immutable) object of type @ref nmm_codon_marg
 * and is used to compute the marginalization forms of p(𝑋₁=𝚡₁,𝑋₂=𝚡₂,𝑋₃=𝚡₃),
 * the probability of emitting codon (𝚡₁,𝚡₂,𝚡₃).
 */

struct nmm_codon_lprob;

#define NMM_ASCII_LAST_STD 127

struct nmm_codon_marg
{
    struct nmm_base_abc const* base_abc;
    /**
     * Maps alphabet symbols to the indices 0, 1, 2, and 3 and the any-symbol to 4.
     */
    uint8_t symbol_idx[NMM_ASCII_LAST_STD + 1];
    /**
     * Pre-computed marginalization forms of
     * p(𝑋₁=𝚡₁,𝑋₂=𝚡₂,𝑋₃=𝚡₃).
     */
    struct nmm_array3d lprobs;
};

static inline struct nmm_base_abc const* nmm_codon_marg_abc(struct nmm_codon_marg const* codonm);
NMM_API struct nmm_codon_marg const*     nmm_codon_marg_create(struct nmm_codon_lprob const* codonp);
NMM_API void                             nmm_codon_marg_destroy(struct nmm_codon_marg const* codonm);
static inline imm_float nmm_codon_marg_lprob(struct nmm_codon_marg const* codonm, struct nmm_codon const* codon);

static inline struct nmm_array3d_idx __nmm_codon_marg_array_idx(struct nmm_codon_marg const* tbl,
                                                                struct nmm_codon const*      codon);

static inline struct nmm_base_abc const* nmm_codon_marg_abc(struct nmm_codon_marg const* codonm)
{
    return codonm->base_abc;
}

/**
 * Calculate any of the marginalization forms of
 * p(𝑋₁=𝚡₁,𝑋₂=𝚡₂,𝑋₃=𝚡₃).
 *
 * The alphabet's any-symbol can be passed to @codon to perform marginalization over
 * the corresponding random variable. Let `"ACGT"` be a set of nucleotides and let `'X`'
 * be the any-symbol of the given alphabet. The code
 *
 *     nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'X', 'G'));
 *     nmm_codon_marg_lprob(codonm, codon);
 *
 * will evaluate the probability p(𝑋₁=𝙰,𝑋₃=𝙶).
 */
static inline imm_float nmm_codon_marg_lprob(struct nmm_codon_marg const* codonm, struct nmm_codon const* codon)
{
    return nmm_array3d_get(&codonm->lprobs, __nmm_codon_marg_array_idx(codonm, codon));
}

static inline struct nmm_array3d_idx __nmm_codon_marg_array_idx(struct nmm_codon_marg const* tbl,
                                                                struct nmm_codon const*      codon)
{
    struct nmm_triplet const t = nmm_codon_get_triplet(codon);
    return (struct nmm_array3d_idx){tbl->symbol_idx[(size_t)t.a], tbl->symbol_idx[(size_t)t.b],
                                    tbl->symbol_idx[(size_t)t.c]};
}

#endif

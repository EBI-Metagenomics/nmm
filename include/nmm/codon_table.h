#ifndef NMM_CODON_TABLE_H
#define NMM_CODON_TABLE_H

#include "nmm/array.h"
#include "nmm/codon.h"
#include "nmm/export.h"

/** @file codont.h
 * Codon table module.
 *
 * A codon table is represented by an (immutable) object of type @ref nmm_codon_table
 * and is used to compute the marginalization forms of p(𝑋₁=𝚡₁,𝑋₂=𝚡₂,𝑋₃=𝚡₃),
 * the probability of emitting codon (𝚡₁,𝚡₂,𝚡₃).
 */

struct nmm_codon_lprob;

#define NMM_ASCII_LAST_STD 127

struct nmm_codon_table
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

static inline struct nmm_base_abc const* nmm_codon_table_abc(struct nmm_codon_table const* codont);
NMM_API struct nmm_codon_table const*    nmm_codon_table_create(struct nmm_codon_lprob const* prob);
NMM_API void        nmm_codon_table_destroy(struct nmm_codon_table const* codont);
static inline float nmm_codon_table_lprob(struct nmm_codon_table const* codont,
                                          struct nmm_codon const*       codon);

static inline struct nmm_array3d_idx __nmm_codon_table_array_idx(struct nmm_codon_table const* tbl,
                                                                 struct nmm_codon const* codon);
/**
 * Calculate any of the marginalization forms of
 * p(𝑋₁=𝚡₁,𝑋₂=𝚡₂,𝑋₃=𝚡₃).
 *
 * The alphabet's any-symbol can be passed to @codon to perform marginalization over
 * the corresponding random variable. Let `"ACGT"` be a set of nucleotides and let `'X`'
 * be the any-symbol of the given alphabet. The code
 *
 *     codon_lprob_get(codon_lprob, "AXG")
 *
 * will evaluate the probability p(𝑋₁=𝙰,𝑋₃=𝙶).
 */

static inline struct nmm_base_abc const* nmm_codon_table_abc(struct nmm_codon_table const* codont)
{
    return codont->base_abc;
}

static inline float nmm_codon_table_lprob(struct nmm_codon_table const* codont,
                                          struct nmm_codon const*       codon)
{
    return nmm_array3d_get(&codont->lprobs, __nmm_codon_table_array_idx(codont, codon));
}

static inline struct nmm_array3d_idx __nmm_codon_table_array_idx(struct nmm_codon_table const* tbl,
                                                                 struct nmm_codon const* codon)
{
    struct nmm_triplet const t = nmm_codon_get_triplet(codon);
    return (struct nmm_array3d_idx){tbl->symbol_idx[(size_t)t.a], tbl->symbol_idx[(size_t)t.b],
                                    tbl->symbol_idx[(size_t)t.c]};
}

#endif

#ifndef NMM_CODONT_STATIC_H
#define NMM_CODONT_STATIC_H

#include "array.h"
#include "codon_static.h"

#define ASCII_LAST_STD 127

struct nmm_codont
{
    struct nmm_base const* base;
    /**
     * Maps alphabet symbols to the indices 0, 1, 2, and 3 and the any-symbol to 4.
     */
    unsigned symbol_idx[ASCII_LAST_STD + 1];
    /**
     * Pre-computed marginalization forms of
     * p(𝑋₁=𝚡₁,𝑋₂=𝚡₂,𝑋₃=𝚡₃).
     */
    struct array3d lprobs;
};

static inline struct array3d_idx codont_get_array_idx(struct nmm_codont const* codont,
                                                      struct nmm_codon const*  codon)
{
    struct nmm_triplet const t = codon_get(codon);
    return (struct array3d_idx){codont->symbol_idx[(size_t)t.a],
                                codont->symbol_idx[(size_t)t.b],
                                codont->symbol_idx[(size_t)t.c]};
}

static inline double codont_lprob(struct nmm_codont const* codont,
                                  struct nmm_codon const*  codon)
{
    return array3d_get(&codont->lprobs, codont_get_array_idx(codont, codon));
}

#endif

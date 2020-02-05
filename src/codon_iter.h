#ifndef NMM_CODON_ITER_H
#define NMM_CODON_ITER_H

#include "codon_static.h"
#include "imm/imm.h"
#include "nmm/base.h"
#include "nmm/codon.h"
#include <stdbool.h>

struct codon_iter
{
    struct nmm_base const* base;
    char const*            bases;
    int                    pos;
};

static inline struct codon_iter codon_iter_begin(struct nmm_base const* base)
{
    char const* bases = imm_abc_symbols(nmm_base_get_abc(base));
    return (struct codon_iter){base, bases, 0};
}

static inline struct nmm_codon codon_iter_next(struct codon_iter* iter)
{
    int a = (iter->pos / (NMM_BASE_SIZE * NMM_BASE_SIZE)) % NMM_BASE_SIZE;
    int b = (iter->pos / NMM_BASE_SIZE) % NMM_BASE_SIZE;
    int c = iter->pos % NMM_BASE_SIZE;
    iter->pos++;

    struct nmm_triplet t = {iter->bases[a], iter->bases[b], iter->bases[c]};
    CODON_DECL(codon, iter->base);
    nmm_codon_set(&codon, t);

    return codon;
}

static inline bool codon_iter_end(struct codon_iter const iter)
{
    return iter.pos >= NMM_BASE_SIZE * NMM_BASE_SIZE * NMM_BASE_SIZE;
}

#endif

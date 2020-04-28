#ifndef CODON_ITER_H
#define CODON_ITER_H

#include "imm/imm.h"
#include "nmm/base_abc.h"
#include "nmm/codon.h"
#include <stdbool.h>

struct codon_iter
{
    struct nmm_base_abc const* base_abc;
    char const*                bases;
    int                        pos;
};

static inline struct codon_iter codon_iter_begin(struct nmm_base_abc const* base)
{
    char const* bases = imm_abc_symbols(nmm_base_abc_super(base));
    return (struct codon_iter){base, bases, 0};
}

static inline struct nmm_codon codon_iter_next(struct codon_iter* iter)
{
    int a = (iter->pos / (NMM_BASE_ABC_SIZE * NMM_BASE_ABC_SIZE)) % NMM_BASE_ABC_SIZE;
    int b = (iter->pos / NMM_BASE_ABC_SIZE) % NMM_BASE_ABC_SIZE;
    int c = iter->pos % NMM_BASE_ABC_SIZE;
    iter->pos++;

    struct nmm_triplet t = {iter->bases[a], iter->bases[b], iter->bases[c]};
    NMM_CODON_DECL(codon, iter->base_abc);
    nmm_codon_set_triplet(&codon, t);

    return codon;
}

static inline bool codon_iter_end(struct codon_iter const iter)
{
    return iter.pos >= NMM_BASE_ABC_SIZE * NMM_BASE_ABC_SIZE * NMM_BASE_ABC_SIZE;
}

#endif

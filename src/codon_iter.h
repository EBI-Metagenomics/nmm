#ifndef NMM_CODON_ITER_H
#define NMM_CODON_ITER_H

#include "imm/imm.h"
#include "nmm/nmm.h"

struct codon_iter
{
    char const* bases;
    int         pos;
};

static inline struct codon_iter codon_iter_begin(struct imm_abc const* abc)
{
    return (struct codon_iter){imm_abc_symbols(abc), 0};
}

static inline struct nmm_codon codon_iter_next(struct codon_iter* iter)
{
    int a = (iter->pos / (NMM_BASE_SIZE * NMM_BASE_SIZE)) % NMM_BASE_SIZE;
    int b = (iter->pos / NMM_BASE_SIZE) % NMM_BASE_SIZE;
    int c = iter->pos % NMM_BASE_SIZE;
    iter->pos++;

    return NMM_CODON(iter->bases[a], iter->bases[b], iter->bases[c]);
}

static inline int codon_iter_end(struct codon_iter const* iter)
{
    return iter->pos >= NMM_BASE_SIZE * NMM_BASE_SIZE * NMM_BASE_SIZE;
}

#endif

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
    int a = (iter->pos / (NMM_CODON_NBASES * NMM_CODON_NBASES)) % NMM_CODON_NBASES;
    int b = (iter->pos / NMM_CODON_NBASES) % NMM_CODON_NBASES;
    int c = iter->pos % NMM_CODON_NBASES;
    iter->pos++;

    return NMM_CODON(iter->bases[a], iter->bases[b], iter->bases[c]);
}

static inline int codon_iter_end(struct codon_iter const* iter)
{
    return iter->pos >= NMM_CODON_NBASES * NMM_CODON_NBASES * NMM_CODON_NBASES;
}

#endif

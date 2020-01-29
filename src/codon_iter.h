#ifndef NMM_CODON_ITER_H
#define NMM_CODON_ITER_H

#include "imm/imm.h"
#include "nmm/base.h"
#include "nmm/codon.h"
#include <stdbool.h>

struct codon_iter
{
    struct nmm_codon* codon;
    char const*       bases;
    int               pos;
};

static inline struct codon_iter codon_iter_begin(struct nmm_base const* base)
{
    struct nmm_codon* codon = nmm_codon_create(base);
    char const*       bases = imm_abc_symbols(nmm_base_get_abc(base));
    return (struct codon_iter){codon, bases, 0};
}

static inline struct nmm_codon const* codon_iter_next(struct codon_iter * iter)
{
    int a = (iter->pos / (NMM_BASE_SIZE * NMM_BASE_SIZE)) % NMM_BASE_SIZE;
    int b = (iter->pos / NMM_BASE_SIZE) % NMM_BASE_SIZE;
    int c = iter->pos % NMM_BASE_SIZE;
    iter->pos++;

    nmm_codon_set(iter->codon, iter->bases[a], iter->bases[b], iter->bases[c]);
    return iter->codon;
}

static inline bool codon_iter_end(struct codon_iter const iter)
{
    return iter.pos >= NMM_BASE_SIZE * NMM_BASE_SIZE * NMM_BASE_SIZE;
}

static inline void codon_iter_destroy(struct codon_iter iter)
{
    nmm_codon_destroy(iter.codon);
}

#endif

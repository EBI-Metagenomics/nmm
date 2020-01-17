#ifndef NMM_CODON_LPROB2_H
#define NMM_CODON_LPROB2_H

#include "array.h"
#include "imm/imm.h"
#include "nmm/nmm.h"

#define ASCII_LAST_STD 127

struct codon_lprob2
{
    int            symbol_idx[ASCII_LAST_STD + 1];
    struct array3d lprobs;
};

struct codon_lprob2* codon_lprob2_create(const struct imm_abc*         abc,
                                         struct nmm_codon_lprob const* lprobs,
                                         int                           lprobs_length);
static inline void   codon_lprob2_set(struct codon_lprob2*    codon_lprob,
                                      struct nmm_codon const* codon, double lprob)
{

    int const* symbol_idx = codon_lprob->symbol_idx;
    int const  dim[3] = {symbol_idx[(size_t)codon->a], symbol_idx[(size_t)codon->b],
                        symbol_idx[(size_t)codon->c]};
    array3d_set(&codon_lprob->lprobs, dim[0], dim[1], dim[2], lprob);
}
static inline double codon_lprob2_get(struct codon_lprob2 const* codon_lprob,
                                      struct nmm_codon const*    codon)
{

    int const* symbol_idx = codon_lprob->symbol_idx;
    int const  dim[3] = {symbol_idx[(size_t)codon->a], symbol_idx[(size_t)codon->b],
                        symbol_idx[(size_t)codon->c]};
    return array3d_get(&codon_lprob->lprobs, dim[0], dim[1], dim[2]);
}
void codon_lprob2_destroy(struct codon_lprob2* codon_lprob);

#endif

#ifndef NMM_AMINO_LPROB_H
#define NMM_AMINO_LPROB_H

#include "imm/imm.h"
#include "nmm/amino_abc.h"
#include "nmm/export.h"

struct nmm_amino_abc;

struct nmm_amino_lprob
{
    struct nmm_amino_abc const* amino_abc;
    imm_float                   lprobs[NMM_AMINO_ABC_SIZE];
};

NMM_API struct nmm_amino_lprob const* nmm_amino_lprob_create(struct nmm_amino_abc const* abc,
                                                             imm_float const*            lprobs);
NMM_API void                          nmm_amino_lprob_destroy(struct nmm_amino_lprob const* aminop);
static inline struct nmm_amino_abc const* nmm_amino_lprob_abc(struct nmm_amino_lprob const* aminop);
static inline imm_float nmm_amino_lprob_get(struct nmm_amino_lprob const* aminop, char const amino);

static inline struct nmm_amino_abc const* nmm_amino_lprob_abc(struct nmm_amino_lprob const* aminop)
{
    return aminop->amino_abc;
}

static inline imm_float nmm_amino_lprob_get(struct nmm_amino_lprob const* aminop, char const amino)
{
    uint_fast8_t i = imm_abc_symbol_idx(nmm_amino_abc_super(aminop->amino_abc), amino);
    return aminop->lprobs[i];
}

#endif

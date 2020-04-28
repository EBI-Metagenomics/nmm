#ifndef NMM_AMINO_TABLE_H
#define NMM_AMINO_TABLE_H

#include "imm/imm.h"
#include "nmm/amino_abc.h"
#include "nmm/array_size.h"
#include "nmm/export.h"

struct nmm_amino_abc;

struct nmm_amino_table
{
    struct nmm_amino_abc const* amino_abc;
    double                      lprobs[NMM_AMINO_ABC_SIZE];
};

NMM_EXPORT struct nmm_amino_table const* nmm_amino_table_create(
    struct nmm_amino_abc const* amino_abc, double const* lprobs);

static inline double nmm_amino_table_lprob(struct nmm_amino_table const* aminot, char const amino)
{
    unsigned i = imm_abc_symbol_idx(nmm_amino_abc_super(aminot->amino_abc), amino);
    return aminot->lprobs[i];
}

NMM_EXPORT void nmm_amino_table_destroy(struct nmm_amino_table const* aminot);

static inline struct nmm_amino_abc const* nmm_amino_table_get_amino_abc(
    struct nmm_amino_table const* aminot)
{
    return aminot->amino_abc;
}

#endif

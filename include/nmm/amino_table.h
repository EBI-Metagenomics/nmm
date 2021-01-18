#ifndef NMM_AMINO_TABLE_H
#define NMM_AMINO_TABLE_H

#include "imm/imm.h"
#include "nmm/amino_abc.h"
#include "nmm/export.h"

struct nmm_amino_abc;

struct nmm_amino_table
{
    struct nmm_amino_abc const* amino_abc;
    double                      lprobs[NMM_AMINO_ABC_SIZE];
};

NMM_API struct nmm_amino_table const* nmm_amino_table_create(struct nmm_amino_abc const* abc,
                                                             double const*               lprobs);
NMM_API void                          nmm_amino_table_destroy(struct nmm_amino_table const* tbl);
static inline struct nmm_amino_abc const* nmm_amino_table_abc(struct nmm_amino_table const* tbl);
static inline double nmm_amino_table_lprob(struct nmm_amino_table const* tbl, char const amino);

static inline struct nmm_amino_abc const* nmm_amino_table_abc(struct nmm_amino_table const* tbl)
{
    return tbl->amino_abc;
}

static inline double nmm_amino_table_lprob(struct nmm_amino_table const* tbl, char const amino)
{
    unsigned i = imm_abc_symbol_idx(nmm_amino_abc_super(tbl->amino_abc), amino);
    return tbl->lprobs[i];
}

#endif

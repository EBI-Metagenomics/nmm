#ifndef NMM_AMINO_TABLE_H
#define NMM_AMINO_TABLE_H

#include "imm/imm.h"
#include "nmm/amino_abc.h"
#include "nmm/api.h"
#include "nmm/array_size.h"

struct nmm_amino_abc;

struct nmm_amino_table
{
    struct nmm_amino_abc const* amino_abc;
    double                      lprobs[NMM_AMINO_ABC_SIZE];
};

NMM_API struct nmm_amino_table const* nmm_amino_table_create(
    struct nmm_amino_abc const* amino_abc, double const* lprobs);

NMM_API static inline double nmm_amino_table_lprob(struct nmm_amino_table const* aminot,
                                                   char const                    amino)
{
    int idx = imm_abc_symbol_idx(nmm_amino_abc_cast(aminot->amino_abc), amino);
    if (idx < 0) {
        imm_error("amino not found");
        return imm_lprob_invalid();
    }
    size_t i = (size_t)idx;
    IMM_BUG(i >= NMM_ARRAY_SIZE(aminot->lprobs));
    return aminot->lprobs[i];
}

NMM_API void nmm_amino_table_destroy(struct nmm_amino_table const* aminot);

NMM_API static inline struct nmm_amino_abc const* nmm_amino_table_get_amino_abc(
    struct nmm_amino_table const* aminot)
{
    return aminot->amino_abc;
}

#endif

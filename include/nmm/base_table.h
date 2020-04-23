#ifndef NMM_BASE_TABLE_H
#define NMM_BASE_TABLE_H

#include "imm/imm.h"
#include "nmm/array_size.h"
#include "nmm/base_abc.h"
#include "nmm/export.h"

struct nmm_base_abc;

struct nmm_base_table
{
    struct nmm_base_abc const* base_abc;
    double                     lprobs[NMM_BASE_ABC_SIZE];
};

NMM_EXPORT struct nmm_base_table const* nmm_base_table_create(struct nmm_base_abc const* base_abc,
                                                              double a, double b, double c,
                                                              double d);

static inline double nmm_base_table_lprob(struct nmm_base_table const* baset, char const base)
{
    unsigned i = imm_abc_symbol_idx(nmm_base_abc_parent(baset->base_abc), base);
    return baset->lprobs[i];
}

NMM_EXPORT void nmm_base_table_destroy(struct nmm_base_table const* baset);

static inline struct nmm_base_abc const* nmm_base_table_get_base_abc(
    struct nmm_base_table const* baset)
{
    return baset->base_abc;
}

#endif

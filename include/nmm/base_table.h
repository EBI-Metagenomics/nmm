#ifndef NMM_BASE_TABLE_H
#define NMM_BASE_TABLE_H

#include "imm/imm.h"
#include "nmm/api.h"
#include "nmm/array_size.h"
#include "nmm/base_abc.h"

struct nmm_base_abc;

struct nmm_base_table
{
    struct nmm_base_abc const* base;
    double                     lprobs[4];
};

NMM_API struct nmm_base_table const* nmm_base_table_create(struct nmm_base_abc const* base,
                                                           double a, double b, double c,
                                                           double d);

NMM_API static inline double nmm_base_table_lprob(struct nmm_base_table const* baset,
                                                  char const                   base)
{
    int idx = imm_abc_symbol_idx(nmm_base_abc_cast(baset->base), base);
    if (idx < 0) {
        imm_error("base not found");
        return imm_lprob_invalid();
    }
    size_t i = (size_t)idx;
    IMM_BUG(i >= NMM_ARRAY_SIZE(baset->lprobs));
    return baset->lprobs[i];
}

NMM_API void nmm_base_table_destroy(struct nmm_base_table const* baset);

NMM_API static inline struct nmm_base_abc const* nmm_base_table_get_base(
    struct nmm_base_table const* baset)
{
    return baset->base;
}

#endif

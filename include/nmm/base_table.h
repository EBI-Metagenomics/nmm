#ifndef NMM_BASE_TABLE_H
#define NMM_BASE_TABLE_H

#include "nmm/api.h"

struct nmm_base;

struct nmm_base_table
{
    struct nmm_base const* base;
    double                 lprobs[4];
};

NMM_API struct nmm_base_table const* nmm_base_table_create(struct nmm_base const* base,
                                                           double a, double b, double c,
                                                           double d);
NMM_API double nmm_base_table_lprob(struct nmm_base_table const* baset, char nucleotide);
NMM_API void   nmm_base_table_destroy(struct nmm_base_table const* baset);
NMM_API static inline struct nmm_base const* nmm_base_table_get_base(
    struct nmm_base_table const* baset)
{
    return baset->base;
}

#endif

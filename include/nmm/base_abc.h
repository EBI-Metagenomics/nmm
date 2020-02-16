#ifndef NMM_BASE_ABC_H
#define NMM_BASE_ABC_H

#include "nmm/api.h"

#define NMM_BASE_ABC_SIZE 4

struct imm_abc;

struct nmm_base_abc
{
    struct imm_abc const* abc;
};

NMM_API struct nmm_base_abc const* nmm_base_abc_create(struct imm_abc const* abc);
NMM_API void                       nmm_base_abc_destroy(struct nmm_base_abc const* base_abc);

NMM_API static inline struct imm_abc const* nmm_base_abc_cast(
    struct nmm_base_abc const* base_abc)
{
    return base_abc->abc;
}

#endif
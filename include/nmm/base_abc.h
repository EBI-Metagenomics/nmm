#ifndef NMM_BASE_ABC_H
#define NMM_BASE_ABC_H

#include "nmm/export.h"

#define NMM_BASE_ABC_SIZE 4

struct imm_abc;

struct nmm_base_abc
{
    struct imm_abc const* parent;
};

NMM_EXPORT struct nmm_base_abc const* nmm_base_abc_create(char const* symbols,
                                                          char const  any_symbol);
NMM_EXPORT void                       nmm_base_abc_destroy(struct nmm_base_abc const* base_abc);

static inline struct imm_abc const* nmm_base_abc_parent(struct nmm_base_abc const* base_abc)
{
    return base_abc->parent;
}

NMM_EXPORT struct nmm_base_abc const* nmm_base_abc_child(struct imm_abc const* abc);

#endif

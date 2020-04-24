#ifndef NMM_AMINO_ABC_H
#define NMM_AMINO_ABC_H

#include "nmm/export.h"

#define NMM_AMINO_ABC_SIZE 20

struct imm_abc;

struct nmm_amino_abc
{
    struct imm_abc const* parent;
};

NMM_EXPORT struct nmm_amino_abc const* nmm_amino_abc_create(char const* symbols,
                                                            char const  any_symbol);
NMM_EXPORT void                        nmm_amino_abc_destroy(struct nmm_amino_abc const* amino_abc);

static inline struct imm_abc const* nmm_amino_abc_parent(struct nmm_amino_abc const* amino_abc)
{
    return amino_abc->parent;
}

NMM_EXPORT struct nmm_amino_abc const* nmm_amino_abc_child(struct imm_abc const* abc);

#endif

#ifndef NMM_AMINO_ABC_H
#define NMM_AMINO_ABC_H

#include "nmm/api.h"

#define NMM_AMINO_ABC_SIZE 20

struct imm_abc;

struct nmm_amino_abc
{
    struct imm_abc const* abc;
};

NMM_API struct nmm_amino_abc const* nmm_amino_abc_create(struct imm_abc const* abc);
NMM_API void nmm_amino_abc_destroy(struct nmm_amino_abc const* amino_abc);

NMM_API static inline struct imm_abc const* nmm_amino_abc_cast(
    struct nmm_amino_abc const* amino_abc)
{
    return amino_abc->abc;
}

#endif

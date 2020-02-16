#include "nmm/amino_abc.h"
#include "free.h"
#include "imm/imm.h"
#include <stdlib.h>

struct nmm_amino_abc const* nmm_amino_abc_create(struct imm_abc const* abc)
{
    if (imm_abc_length(abc) != NMM_AMINO_ABC_SIZE) {
        imm_error("alphabet length must be %u", NMM_AMINO_ABC_SIZE);
        return NULL;
    }
    struct nmm_amino_abc* amino_abc = malloc(sizeof(struct nmm_amino_abc));
    amino_abc->abc = abc;
    return amino_abc;
}

void nmm_amino_abc_destroy(struct nmm_amino_abc const* amino_abc) { free_c(amino_abc); }

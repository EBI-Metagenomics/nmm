#include "nmm/amino_table.h"
#include "free.h"
#include <stdlib.h>
#include <string.h>

struct nmm_amino_table const* nmm_amino_table_create(struct nmm_amino_abc const* amino_abc,
                                                     double const*               lprobs)
{
    struct nmm_amino_table* aminot = malloc(sizeof(struct nmm_amino_table));
    aminot->amino_abc = amino_abc;
    memcpy(aminot->lprobs, lprobs, NMM_AMINO_ABC_SIZE * sizeof(double));
    return aminot;
}

void nmm_amino_table_destroy(struct nmm_amino_table const* aminot) { free_c(aminot); }

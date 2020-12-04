#include "nmm/amino_table.h"
#include "free.h"
#include <stdlib.h>
#include <string.h>

struct nmm_amino_table const* nmm_amino_table_create(struct nmm_amino_abc const* abc,
                                                     float const*                lprobs)
{
    struct nmm_amino_table* aminot = malloc(sizeof(*aminot));
    aminot->amino_abc = abc;
    memcpy(aminot->lprobs, lprobs, NMM_AMINO_ABC_SIZE * sizeof(float));
    return aminot;
}

void nmm_amino_table_destroy(struct nmm_amino_table const* tbl) { free_c(tbl); }

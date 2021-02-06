#include "nmm/amino_lprob.h"
#include "free.h"
#include <stdlib.h>
#include <string.h>

struct nmm_amino_lprob const* nmm_amino_lprob_create(struct nmm_amino_abc const* abc,
                                                     imm_float const             lprobs[static NMM_AMINO_ABC_SIZE])
{
    struct nmm_amino_lprob* aminop = malloc(sizeof(*aminop));
    aminop->amino_abc = abc;
    memcpy(aminop->lprobs, lprobs, NMM_AMINO_ABC_SIZE * sizeof(*lprobs));
    return aminop;
}

void nmm_amino_lprob_destroy(struct nmm_amino_lprob const* aminop) { free_c(aminop); }

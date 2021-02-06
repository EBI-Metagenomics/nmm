#include "base_lprob.h"
#include "free.h"
#include "nmm/base_lprob.h"
#include <stdlib.h>

struct nmm_base_lprob const* nmm_base_lprob_create(struct nmm_base_abc const* abc,
                                                   imm_float const a, imm_float const b,
                                                   imm_float const c, imm_float const d)
{
    struct nmm_base_lprob* basep = malloc(sizeof(*basep));
    basep->base_abc = abc;
    basep->lprobs[0] = a;
    basep->lprobs[1] = b;
    basep->lprobs[2] = c;
    basep->lprobs[3] = d;
    return basep;
}

void nmm_base_lprob_destroy(struct nmm_base_lprob const* basep) { free_c(basep); }

struct nmm_base_lprob const* base_lprob_read(FILE* stream, struct nmm_base_abc const* abc)
{
    struct nmm_base_lprob* basep = malloc(sizeof(*basep));
    basep->base_abc = abc;

    if (fread(basep->lprobs, sizeof(*basep->lprobs), NMM_BASE_ABC_SIZE, stream) <
        NMM_BASE_ABC_SIZE) {
        imm_error("could not read lprobs");
        free_c(basep);
        return NULL;
    }

    return basep;
}

int base_lprob_write(struct nmm_base_lprob const* basep, FILE* stream)
{
    if (fwrite(basep->lprobs, sizeof(*basep->lprobs), NMM_BASE_ABC_SIZE, stream) <
        NMM_BASE_ABC_SIZE) {
        imm_error("could not write lprobs");
        return 1;
    }

    return 0;
}

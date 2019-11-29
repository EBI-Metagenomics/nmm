#include "imm/imm.h"
#include "nmm/nmm.h"
#include <math.h>
#include <stdlib.h>

struct nmm_baset
{
    const struct imm_abc* abc;
    double                emiss_lprobs[4];
};

struct nmm_baset* nmm_baset_create(const struct imm_abc* abc)
{
    if (imm_abc_length(abc) != 4) {
        imm_error("alphabet length is not four");
        return NULL;
    }
    struct nmm_baset* baset = malloc(sizeof(struct nmm_baset));
    baset->abc = abc;
    baset->emiss_lprobs[0] = -INFINITY;
    baset->emiss_lprobs[1] = -INFINITY;
    baset->emiss_lprobs[2] = -INFINITY;
    baset->emiss_lprobs[3] = -INFINITY;
    return baset;
}

int nmm_baset_set_lprob(struct nmm_baset* baset, char nucleotide, double lprob)
{
    int idx = imm_abc_symbol_idx(baset->abc, nucleotide);
    if (idx < 0) {
        imm_error("nucleotide not found");
        return 1;
    }
    baset->emiss_lprobs[idx] = lprob;
    return 0;
}

double nmm_baset_get_lprob(const struct nmm_baset* baset, char nucleotide)
{
    int idx = imm_abc_symbol_idx(baset->abc, nucleotide);
    if (idx < 0) {
        imm_error("nucleotide not found");
        return NAN;
    }
    return baset->emiss_lprobs[idx];
}

int nmm_baset_normalize(struct nmm_baset* baset)
{
    return imm_lprob_normalize(baset->emiss_lprobs, 4);
}

void nmm_baset_destroy(struct nmm_baset* baset)
{
    if (!baset)
        return;

    baset->abc = NULL;
    free(baset);
}

const struct imm_abc* nmm_baset_get_abc(const struct nmm_baset* baset) { return baset->abc; }

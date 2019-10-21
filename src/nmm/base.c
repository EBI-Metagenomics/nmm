#include "imm.h"
#include "nmm.h"
#include <math.h>
#include <stdlib.h>

struct nmm_base
{
    const struct imm_abc *abc;
    double emiss_lprobs[4];
};

struct nmm_base *nmm_base_create(const struct imm_abc *abc)
{
    if (imm_abc_length(abc) != 4) {
        imm_error("alphabet length is not four");
        return NULL;
    }
    struct nmm_base *base = malloc(sizeof(struct nmm_base));
    base->abc = abc;
    base->emiss_lprobs[0] = -INFINITY;
    base->emiss_lprobs[1] = -INFINITY;
    base->emiss_lprobs[2] = -INFINITY;
    base->emiss_lprobs[3] = -INFINITY;
    return base;
}

int nmm_base_set_lprob(struct nmm_base *base, char nucleotide, double lprob)
{
    int idx = imm_abc_symbol_idx(base->abc, nucleotide);
    if (idx < 0) {
        imm_error("nucleotide not found");
        return 1;
    }
    base->emiss_lprobs[idx] = lprob;
    return 0;
}

double nmm_base_get_lprob(const struct nmm_base *base, char nucleotide)
{
    int idx = imm_abc_symbol_idx(base->abc, nucleotide);
    if (idx < 0) {
        imm_error("nucleotide not found");
        return NAN;
    }
    return base->emiss_lprobs[idx];
}

int nmm_base_normalize(struct nmm_base *base)
{
    return imm_lognormalize(base->emiss_lprobs, 4);
}

void nmm_base_destroy(struct nmm_base *base)
{
    if (!base)
        return;

    base->abc = NULL;
    free(base);
}

const struct imm_abc *nmm_base_get_abc(const struct nmm_base *base)
{
    return base->abc;
}

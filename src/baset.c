#include "nmm/baset.h"
#include "free.h"
#include "imm/imm.h"
#include "nmm/base.h"
#include <stdlib.h>

struct nmm_baset
{
    struct nmm_base const* base;
    double                 emiss_lprobs[4];
};

struct nmm_baset const* nmm_baset_create(struct nmm_base const* base, double a, double b,
                                         double c, double d)
{
    struct nmm_baset* baset = malloc(sizeof(struct nmm_baset));
    baset->base = base;
    baset->emiss_lprobs[0] = a;
    baset->emiss_lprobs[1] = b;
    baset->emiss_lprobs[2] = c;
    baset->emiss_lprobs[3] = d;
    return baset;
}

double nmm_baset_lprob(const struct nmm_baset* baset, char const nucleotide)
{
    if (baset == NULL)
        return 0.0;
    if (baset->base == NULL)
        return 0.0;
    int idx = imm_abc_symbol_idx(nmm_base_get_abc(baset->base), nucleotide);
    if (idx < 0) {
        imm_error("nucleotide not found");
        return imm_lprob_invalid();
    }
    return baset->emiss_lprobs[idx];
}

void nmm_baset_destroy(struct nmm_baset const* baset) { free_c(baset); }

struct nmm_base const* nmm_baset_get_base(struct nmm_baset const* baset)
{
    return baset->base;
}

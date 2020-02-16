#include "nmm/baset.h"
#include "array_size.h"
#include "baset_static.h"
#include "free.h"
#include "imm/imm.h"
#include "nmm/base.h"
#include <stdlib.h>

struct nmm_baset const* nmm_baset_create(struct nmm_base const* base, double const a,
                                         double const b, double const c, double const d)
{
    struct nmm_baset* baset = malloc(sizeof(struct nmm_baset));
    baset->base = base;
    baset->lprobs[0] = a;
    baset->lprobs[1] = b;
    baset->lprobs[2] = c;
    baset->lprobs[3] = d;
    return baset;
}

double nmm_baset_lprob(struct nmm_baset const* baset, char const nucleotide)
{
    int idx = imm_abc_symbol_idx(nmm_base_get_abc(baset->base), nucleotide);
    if (idx < 0) {
        imm_error("nucleotide not found");
        return imm_lprob_invalid();
    }
    size_t i = (size_t)idx;
    IMM_BUG(i >= ARRAY_SIZE(baset->lprobs));
    return baset->lprobs[i];
}

void nmm_baset_destroy(struct nmm_baset const* baset) { free_c(baset); }

struct nmm_base const* nmm_baset_get_base(struct nmm_baset const* baset)
{
    return baset_get_base(baset);
}

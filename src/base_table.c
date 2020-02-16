#include "array_size.h"
#include "free.h"
#include "imm/imm.h"
#include "nmm/base.h"
#include "nmm/base_table.h"
#include <stdlib.h>

struct nmm_base_table const* nmm_base_table_create(struct nmm_base const* base,
                                                   double const a, double const b,
                                                   double const c, double const d)
{
    struct nmm_base_table* baset = malloc(sizeof(struct nmm_base_table));
    baset->base = base;
    baset->lprobs[0] = a;
    baset->lprobs[1] = b;
    baset->lprobs[2] = c;
    baset->lprobs[3] = d;
    return baset;
}

double nmm_base_table_lprob(struct nmm_base_table const* baset, char const nucleotide)
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

void nmm_base_table_destroy(struct nmm_base_table const* baset) { free_c(baset); }

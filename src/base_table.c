#include "nmm/base_table.h"
#include "free.h"
#include <stdlib.h>

struct nmm_base_table const* nmm_base_table_create(struct nmm_base_abc const* base_abc,
                                                   double const a, double const b,
                                                   double const c, double const d)
{
    struct nmm_base_table* baset = malloc(sizeof(struct nmm_base_table));
    baset->base = base_abc;
    baset->lprobs[0] = a;
    baset->lprobs[1] = b;
    baset->lprobs[2] = c;
    baset->lprobs[3] = d;
    return baset;
}

void nmm_base_table_destroy(struct nmm_base_table const* baset) { free_c(baset); }

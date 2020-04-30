#include "base_table.h"
#include "free.h"
#include "nmm/base_table.h"
#include <stdlib.h>

struct nmm_base_table const* nmm_base_table_create(struct nmm_base_abc const* base_abc,
                                                   double const a, double const b, double const c,
                                                   double const d)
{
    struct nmm_base_table* baset = malloc(sizeof(*baset));
    baset->base_abc = base_abc;
    baset->lprobs[0] = a;
    baset->lprobs[1] = b;
    baset->lprobs[2] = c;
    baset->lprobs[3] = d;
    return baset;
}

void nmm_base_table_destroy(struct nmm_base_table const* baset) { free_c(baset); }

struct nmm_base_table const* base_table_read(FILE* stream, struct nmm_base_abc const* base_abc)
{
    struct nmm_base_table* baset = malloc(sizeof(*baset));
    baset->base_abc = base_abc;

    if (fread(baset->lprobs, sizeof(*baset->lprobs), NMM_BASE_ABC_SIZE, stream) <
        NMM_BASE_ABC_SIZE) {
        imm_error("could not read lprobs");
        free_c(baset);
        return NULL;
    }

    return baset;
}

int base_table_write(struct nmm_base_table const* baset, FILE* stream)
{
    if (fwrite(baset->lprobs, sizeof(*baset->lprobs), NMM_BASE_ABC_SIZE, stream) <
        NMM_BASE_ABC_SIZE) {
        imm_error("could not write lprobs");
        return 1;
    }

    return 0;
}

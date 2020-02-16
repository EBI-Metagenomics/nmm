#include "nmm/base.h"
#include "free.h"
#include "imm/imm.h"
#include <stdlib.h>

struct nmm_base const* nmm_base_create(struct imm_abc const* abc)
{
    if (imm_abc_length(abc) != NMM_BASE_SIZE) {
        imm_error("alphabet length must be %u", NMM_BASE_SIZE);
        return NULL;
    }
    struct nmm_base* base = malloc(sizeof(struct nmm_base));
    base->abc = abc;
    return base;
}

void nmm_base_destroy(struct nmm_base const* base) { free_c(base); }

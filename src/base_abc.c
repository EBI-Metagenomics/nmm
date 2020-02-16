#include "nmm/base_abc.h"
#include "free.h"
#include "imm/imm.h"
#include <stdlib.h>

struct nmm_base_abc const* nmm_base_abc_create(struct imm_abc const* abc)
{
    if (imm_abc_length(abc) != NMM_BASE_ABC_SIZE) {
        imm_error("alphabet length must be %u", NMM_BASE_ABC_SIZE);
        return NULL;
    }
    struct nmm_base_abc* base_abc = malloc(sizeof(struct nmm_base_abc));
    base_abc->abc = abc;
    return base_abc;
}

void nmm_base_abc_destroy(struct nmm_base_abc const* base_abc) { free_c(base_abc); }

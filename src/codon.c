#include "nmm/codon.h"
#include "free.h"
#include "imm/imm.h"
#include "nmm/base_abc.h"
#include <stdbool.h>
#include <stdlib.h>

struct nmm_codon* nmm_codon_create(struct nmm_base_abc const* base_abc)
{
    struct nmm_codon* codon = malloc(sizeof(struct nmm_codon));
    nmm_codon_init(codon, base_abc);
    return codon;
}

void nmm_codon_destroy(struct nmm_codon const* codon) { free_c(codon); }

int nmm_codon_set_triplet(struct nmm_codon* codon, struct nmm_triplet triplet)
{
    struct imm_abc const* abc = nmm_base_abc_parent(codon->base_abc);

    bool ok = (imm_abc_symbol_type(abc, triplet.a) != IMM_SYMBOL_UNKNOWN &&
               imm_abc_symbol_type(abc, triplet.b) != IMM_SYMBOL_UNKNOWN &&
               imm_abc_symbol_type(abc, triplet.c) != IMM_SYMBOL_UNKNOWN);

    if (!ok) {
        imm_error("unknown codon");
        return 1;
    }

    codon->a = triplet.a;
    codon->b = triplet.b;
    codon->c = triplet.c;

    return 0;
}

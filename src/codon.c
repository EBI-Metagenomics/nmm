#include "nmm/codon.h"
#include "codon_static.h"
#include "free.h"
#include "nmm/base.h"
#include <stdbool.h>
#include <stdlib.h>

struct nmm_codon* nmm_codon_create(struct nmm_base const* base)
{
    struct nmm_codon* codon = malloc(sizeof(struct nmm_codon));
    codon_init(codon, base);
    return codon;
}

void nmm_codon_destroy(struct nmm_codon const* codon) { free_c(codon); }

struct nmm_base const* nmm_codon_get_base(struct nmm_codon const* codon)
{
    return codon->base;
}

int nmm_codon_set(struct nmm_codon* codon, struct nmm_triplet triplet)
{
    struct imm_abc const* abc = nmm_base_get_abc(codon->base);

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

struct nmm_triplet nmm_codon_get(struct nmm_codon const* codon)
{
    return (struct nmm_triplet){codon->a, codon->b, codon->c};
}

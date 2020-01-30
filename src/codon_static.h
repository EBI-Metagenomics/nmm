#ifndef NMM_CODON_STATIC_H
#define NMM_CODON_STATIC_H

#include "bug.h"
#include "imm/imm.h"
#include "nmm/base.h"
#include "nmm/codon.h"

struct nmm_codon
{
    struct nmm_base const* base;
    char                   a;
    char                   b;
    char                   c;
};

static inline void codon_init(struct nmm_codon* codon, struct nmm_base const* base)
{
    codon->base = base;

    struct imm_abc const* abc = nmm_base_get_abc(base);
    codon->a = imm_abc_any_symbol(abc);
    codon->b = imm_abc_any_symbol(abc);
    codon->c = imm_abc_any_symbol(abc);
}

#define CODON_DECL(name, base) struct nmm_codon name; codon_init(&(name), (base))

static inline struct nmm_codon const* codon_set(struct nmm_codon* codon, char a, char b,
                                                char c)
{
    BUG(nmm_codon_set(codon, (struct nmm_triplet){a, b, c}) != 0);
    return codon;
}

#endif

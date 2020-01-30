#ifndef NMM_CODON_STATIC_H
#define NMM_CODON_STATIC_H

#include "base_static.h"
#include "imm/imm.h"
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
}

#define CODON_DECL(name, base)                                                               \
    struct nmm_codon name;                                                                   \
    codon_init(&(name), (base))

static inline struct nmm_codon const* codon_set(struct nmm_codon* codon, char a, char b,
                                                char c)
{
    codon->a = a;
    codon->b = b;
    codon->c = c;
    return codon;
}

static inline struct nmm_triplet codon_get(struct nmm_codon const* codon)
{
    return NMM_TRIPLET(codon->a, codon->b, codon->c);
}

#endif

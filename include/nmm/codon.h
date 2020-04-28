#ifndef NMM_CODON_H
#define NMM_CODON_H

#include "imm/imm.h"
#include "nmm/base_abc.h"
#include "nmm/export.h"
#include "nmm/triplet.h"

struct nmm_codon
{
    struct nmm_base_abc const* base_abc;
    char                       a;
    char                       b;
    char                       c;
};

NMM_EXPORT struct nmm_codon* nmm_codon_create(struct nmm_base_abc const* base);
NMM_EXPORT void              nmm_codon_destroy(struct nmm_codon const* codon);

static inline struct nmm_base_abc const* nmm_codon_get_base(struct nmm_codon const* codon)
{
    return codon->base_abc;
}

NMM_EXPORT int nmm_codon_set_triplet(struct nmm_codon* codon, struct nmm_triplet triplet);

static inline struct nmm_triplet nmm_codon_get_triplet(struct nmm_codon const* codon)
{
    return NMM_TRIPLET(codon->a, codon->b, codon->c);
}

static inline void nmm_codon_init(struct nmm_codon* codon, struct nmm_base_abc const* base_abc)
{
    codon->base_abc = base_abc;
    codon->a = imm_abc_any_symbol(nmm_base_abc_super(base_abc));
    codon->b = imm_abc_any_symbol(nmm_base_abc_super(base_abc));
    codon->c = imm_abc_any_symbol(nmm_base_abc_super(base_abc));
}

#define NMM_CODON_DECL(name, base)                                                                 \
    struct nmm_codon name;                                                                         \
    nmm_codon_init(&(name), (base))

#endif

#ifndef NMM_CODON_H
#define NMM_CODON_H

#include "imm/imm.h"
#include "nmm/api.h"
#include "nmm/base.h"
#include "nmm/triplet.h"

struct nmm_codon
{
    struct nmm_base const* base;
    char                   a;
    char                   b;
    char                   c;
};

NMM_API struct nmm_codon* nmm_codon_create(struct nmm_base const* base);
NMM_API void              nmm_codon_destroy(struct nmm_codon const* codon);

NMM_API static inline struct nmm_base const* nmm_codon_get_base(struct nmm_codon const* codon)
{
    return codon->base;
}

NMM_API int nmm_codon_set_triplet(struct nmm_codon* codon, struct nmm_triplet triplet);

NMM_API static inline struct nmm_triplet nmm_codon_get_triplet(struct nmm_codon const* codon)
{
    return NMM_TRIPLET(codon->a, codon->b, codon->c);
}

NMM_API static inline void nmm_codon_init(struct nmm_codon*      codon,
                                          struct nmm_base const* base)
{
    codon->base = base;
    codon->a = imm_abc_any_symbol(nmm_base_get_abc(base));
    codon->b = imm_abc_any_symbol(nmm_base_get_abc(base));
    codon->c = imm_abc_any_symbol(nmm_base_get_abc(base));
}

#define NMM_CODON_DECL(name, base)                                                           \
    struct nmm_codon name;                                                                   \
    nmm_codon_init(&(name), (base))

#endif

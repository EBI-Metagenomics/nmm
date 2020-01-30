#ifndef NMM_CODON_H
#define NMM_CODON_H

#include "nmm/api.h"

struct nmm_base;
struct nmm_codon;

struct nmm_triplet
{
    char a;
    char b;
    char c;
};

NMM_API struct nmm_codon* nmm_codon_create(struct nmm_base const* base);
NMM_API void nmm_codon_destroy(struct nmm_codon const* codon);
NMM_API struct nmm_base const* nmm_codon_get_base(struct nmm_codon const* codon);
NMM_API int                nmm_codon_set(struct nmm_codon* codon, struct nmm_triplet triplet);
NMM_API struct nmm_triplet nmm_codon_get(struct nmm_codon const* codon);

#endif

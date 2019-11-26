#ifndef NMM_CODON_H
#define NMM_CODON_H

#include "nmm/api.h"

struct imm_abc;
struct nmm_codont;

struct nmm_ccode
{
    char a;
    char b;
    char c;
};

#define NMM_CCODE(a, b, c) ((struct nmm_ccode){(a), (b), (c)})

NMM_API struct nmm_codont* nmm_codon_create(struct imm_abc const* abc);
NMM_API int    nmm_codon_set_lprob(struct nmm_codont* codon, struct nmm_ccode const* ccode,
                                   double lprob);
NMM_API double nmm_codon_get_lprob(struct nmm_codont const* codon,
                                   struct nmm_ccode const* ccode);
NMM_API int    nmm_codon_normalize(struct nmm_codont* codon);
NMM_API void   nmm_codon_destroy(struct nmm_codont* codon);
NMM_API struct imm_abc const* nmm_codon_get_abc(struct nmm_codont const* codon);

#endif

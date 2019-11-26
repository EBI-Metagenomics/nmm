#ifndef NMM_CODONT_H
#define NMM_CODONT_H

#include "nmm/api.h"

struct imm_abc;
struct nmm_codont;

struct nmm_codon
{
    char a;
    char b;
    char c;
};

#define NMM_CODON(a, b, c) ((struct nmm_codon){(a), (b), (c)})

NMM_API struct nmm_codont* nmm_codont_create(struct imm_abc const* abc);
NMM_API int    nmm_codont_set_lprob(struct nmm_codont* codont, struct nmm_codon const* codon,
                                    double lprob);
NMM_API double nmm_codont_get_lprob(struct nmm_codont const* codont,
                                    struct nmm_codon const*  codon);
NMM_API int    nmm_codont_normalize(struct nmm_codont* codont);
NMM_API void   nmm_codont_destroy(struct nmm_codont* codont);
NMM_API struct imm_abc const* nmm_codont_get_abc(struct nmm_codont const* codont);

#endif

#ifndef NMM_CODON_H
#define NMM_CODON_H

#include "nmm/api.h"

struct imm_abc;
struct nmm_codon;

NMM_API struct nmm_codon* nmm_codon_create(struct imm_abc const* abc);
NMM_API int               nmm_codon_set_lprob(struct nmm_codon* codon, char a, char b, char c,
                                              double lprob);
NMM_API double nmm_codon_get_lprob(struct nmm_codon const* codon, char a, char b, char c);
NMM_API int    nmm_codon_normalize(struct nmm_codon* codon);
NMM_API void   nmm_codon_destroy(struct nmm_codon* codon);
NMM_API struct imm_abc const* nmm_codon_get_abc(struct nmm_codon const* codon);

#endif

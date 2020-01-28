#ifndef NMM_CODONP_H
#define NMM_CODONP_H

#include "nmm/api.h"

struct imm_abc;
struct nmm_codon;
struct nmm_codonp;

NMM_API struct nmm_codonp* nmm_codonp_create(struct imm_abc const* abc);
NMM_API int    nmm_codonp_set(struct nmm_codonp* codonp, struct nmm_codon const* codon,
                              double lprob);
NMM_API double nmm_codonp_get(struct nmm_codonp const* codonp, struct nmm_codon const* codon);
NMM_API int    nmm_codonp_normalize(struct nmm_codonp* codonp);
NMM_API void   nmm_codonp_destroy(struct nmm_codonp const* codonp);
NMM_API struct imm_abc const* nmm_codonp_get_abc(struct nmm_codonp const* codonp);

#endif

#ifndef NMM_BASET_H
#define NMM_BASET_H

#include "nmm/api.h"

struct imm_abc;
struct nmm_baset;

NMM_API struct nmm_baset* nmm_baset_create(struct imm_abc const* abc);
NMM_API int    nmm_baset_set_lprob(struct nmm_baset* baset, char nucleotide, double lprob);
NMM_API double nmm_baset_get_lprob(struct nmm_baset const* baset, char nucleotide);
NMM_API int    nmm_baset_normalize(struct nmm_baset* baset);
NMM_API void   nmm_baset_destroy(struct nmm_baset* baset);
NMM_API struct imm_abc const* nmm_baset_get_abc(struct nmm_baset const* baset);

#endif

#ifndef NMM_BASET_H
#define NMM_BASET_H

#include "nmm/api.h"

struct nmm_base;
struct nmm_baset;

NMM_API struct nmm_baset const* nmm_baset_create(struct nmm_base const* base, double a,
                                                 double b, double c, double d);
NMM_API double nmm_baset_lprob(struct nmm_baset const* baset, char nucleotide);
NMM_API void   nmm_baset_destroy(struct nmm_baset const* baset);
NMM_API struct nmm_base const* nmm_baset_get_base(struct nmm_baset const* baset);

#endif

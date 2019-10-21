#ifndef NMM_BASE_H_API
#define NMM_BASE_H_API

#include "nmm/api.h"

struct imm_abc;
struct nmm_base;

NMM_API struct nmm_base *nmm_base_create(const struct imm_abc *abc);
NMM_API int nmm_base_set_lprob(struct nmm_base *base, char nucleotide, double lprob);
NMM_API double nmm_base_get_lprob(const struct nmm_base *base, char nucleotide);
NMM_API int nmm_base_normalize(struct nmm_base *base);
NMM_API void nmm_base_destroy(struct nmm_base *base);
NMM_API const struct imm_abc *nmm_base_get_abc(const struct nmm_base *base);

#endif

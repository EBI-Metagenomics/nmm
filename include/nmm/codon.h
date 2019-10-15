#ifndef NMM_CODON_H_API
#define NMM_CODON_H_API

#include "nmm/api.h"

struct nmm_codon;

NMM_API struct nmm_codon *nmm_codon_create(void);
NMM_API struct nmm_codon *nmm_codon_clone(const struct nmm_codon *codon);
NMM_API void nmm_codon_set_lprob(struct nmm_codon *codon, int a, int b, int c,
                                 double lprob);
NMM_API void nmm_codon_set_ninfs(struct nmm_codon *codon);
NMM_API double nmm_codon_get_lprob(const struct nmm_codon *codon, int a, int b, int c);
NMM_API int nmm_codon_normalize(struct nmm_codon *codon);
NMM_API void nmm_codon_destroy(struct nmm_codon *codon);

#endif

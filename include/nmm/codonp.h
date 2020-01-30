#ifndef NMM_CODONP_H
#define NMM_CODONP_H

#include "nmm/api.h"

/** @file codonp.h
 * Codon probability module.
 *
 * An object of type @ref nmm_codonp is used to define the probabilities
 * p(𝑋₁=𝚡₁,𝑋₂=𝚡₂,𝑋₃=𝚡₃) of emitting codon (𝚡₁,𝚡₂,𝚡₃). Its sole purpose is
 * to be used by the constructor of the type @ref nmm_codont.
 */

struct nmm_base;
struct nmm_codon;
struct nmm_codonp;

NMM_API struct nmm_codonp* nmm_codonp_create(struct nmm_base const* base);
NMM_API int    nmm_codonp_set(struct nmm_codonp* codonp, struct nmm_codon const* codon,
                              double lprob);
NMM_API double nmm_codonp_get(struct nmm_codonp const* codonp, struct nmm_codon const* codon);
NMM_API int    nmm_codonp_normalize(struct nmm_codonp* codonp);
NMM_API void   nmm_codonp_destroy(struct nmm_codonp const* codonp);
NMM_API struct nmm_base const* nmm_codonp_get_base(struct nmm_codonp const* codonp);

#endif

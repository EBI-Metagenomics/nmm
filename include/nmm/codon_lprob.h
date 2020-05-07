#ifndef NMM_CODON_LPROB_H
#define NMM_CODON_LPROB_H

#include "nmm/export.h"

/** @file codonp.h
 * Codon probability module.
 *
 * An object of type @ref nmm_codon_lprob is used to define the probabilities
 * p(𝑋₁=𝚡₁,𝑋₂=𝚡₂,𝑋₃=𝚡₃) of emitting codon (𝚡₁,𝚡₂,𝚡₃). Its sole purpose is
 * to be used by the constructor of the type @ref nmm_codont.
 */

struct nmm_base_abc;
struct nmm_codon;
struct nmm_codon_lprob;

NMM_EXPORT struct nmm_base_abc const* nmm_codon_lprob_abc(struct nmm_codon_lprob const* codonp);
NMM_EXPORT struct nmm_codon_lprob*    nmm_codon_lprob_create(struct nmm_base_abc const* abc);
NMM_EXPORT void                       nmm_codon_lprob_destroy(struct nmm_codon_lprob const* codonp);
NMM_EXPORT double                     nmm_codon_lprob_get(struct nmm_codon_lprob const* codonp,
                                                          struct nmm_codon const*       codon);
NMM_EXPORT int                        nmm_codon_lprob_normalize(struct nmm_codon_lprob* codonp);
NMM_EXPORT int nmm_codon_lprob_set(struct nmm_codon_lprob* codonp, struct nmm_codon const* codon,
                                   double lprob);

#endif

#ifndef NMM_CODONT_H
#define NMM_CODONT_H

#include "nmm/api.h"

/** @file codont.h
 * Codon table module.
 *
 * A codon table is represented by an (immutable) object of type @ref nmm_codon
 * and is used to compute the marginalization forms of p(𝑋₁=𝚡₁,𝑋₂=𝚡₂,𝑋₃=𝚡₃),
 * the probability of emitting codon (𝚡₁,𝚡₂,𝚡₃).
 */

struct nmm_codon;
struct nmm_codonp;
struct nmm_codont;

NMM_API struct nmm_codont const* nmm_codont_create(struct nmm_codonp const* codonp);
NMM_API double                   nmm_codont_lprob(struct nmm_codont const* codont,
                                                  struct nmm_codon const*  codon);
NMM_API void                     nmm_codont_destroy(struct nmm_codont const* codont);
NMM_API struct nmm_base const*   nmm_codont_get_base(struct nmm_codont const* codont);

#endif

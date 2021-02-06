#ifndef NMM_BASE_LPROB_H
#define NMM_BASE_LPROB_H

#include "imm/imm.h"
#include "nmm/base_abc.h"
#include "nmm/export.h"
#include <stdio.h>

/** @file base_lprob.h
 * Base probability module.
 *
 * An object of type @ref nmm_base_lprob is used to define the probabilities
 * p(𝑋₁=𝚡₁) of emitting base 𝚡₁.
 */

struct nmm_base_abc;

struct nmm_base_lprob
{
    struct nmm_base_abc const* base_abc;
    imm_float                  lprobs[NMM_BASE_ABC_SIZE];
};

static inline struct nmm_base_abc const* nmm_base_lprob_abc(struct nmm_base_lprob const* basep);
NMM_API struct nmm_base_lprob const*     nmm_base_lprob_create(struct nmm_base_abc const* abc,
                                                               imm_float a, imm_float b, imm_float c,
                                                               imm_float d);
NMM_API void                             nmm_base_lprob_destroy(struct nmm_base_lprob const* basep);
static inline imm_float nmm_base_lprob_get(struct nmm_base_lprob const* basep, char const base);

static inline struct nmm_base_abc const* nmm_base_lprob_abc(struct nmm_base_lprob const* basep)
{
    return basep->base_abc;
}

static inline imm_float nmm_base_lprob_get(struct nmm_base_lprob const* basep, char const base)
{
    uint8_t i = imm_abc_symbol_idx(nmm_base_abc_super(basep->base_abc), base);
    return basep->lprobs[i];
}

#endif

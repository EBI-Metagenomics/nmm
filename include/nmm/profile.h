#ifndef NMM_PROFILE_H
#define NMM_PROFILE_H

#include "nmm/export.h"
#include <inttypes.h>
#include <stdbool.h>
#include <stdio.h>

struct imm_abc;
struct imm_model;
struct nmm_base_lprob;
struct nmm_codon_marg;
struct nmm_profile;

NMM_API struct imm_abc const*         nmm_profile_abc(struct nmm_profile const* prof);
NMM_API struct nmm_base_lprob const*  nmm_profile_base_lprob(struct nmm_profile const* prof, uint16_t index);
NMM_API struct nmm_codon_lprob const* nmm_profile_codon_lprob(struct nmm_profile const* prof, uint16_t index);
NMM_API struct nmm_codon_marg const*  nmm_profile_codon_marg(struct nmm_profile const* prof, uint16_t index);
NMM_API struct nmm_profile*           nmm_profile_create(struct imm_abc const* abc);
NMM_API void                          nmm_profile_append_model(struct nmm_profile* prof, struct imm_model* model);
NMM_API struct imm_model*             nmm_profile_get_model(struct nmm_profile const* prof, uint8_t i);
NMM_API uint8_t                       nmm_profile_nmodels(struct nmm_profile const* prof);
NMM_API void                          nmm_profile_destroy(struct nmm_profile const* prof, bool deep);
NMM_API void                          nmm_profile_free(struct nmm_profile const* prof);
NMM_API uint16_t                      nmm_profile_nbase_lprobs(struct nmm_profile const* prof);
NMM_API uint16_t                      nmm_profile_ncodon_lprobs(struct nmm_profile const* prof);
NMM_API uint16_t                      nmm_profile_ncodon_margs(struct nmm_profile const* prof);
NMM_API struct nmm_profile*           nmm_profile_read(FILE* stream);
NMM_API int                           nmm_profile_write(struct nmm_profile const* prof, FILE* stream);

#endif

#ifndef NMM_CODON_STATE_H
#define NMM_CODON_STATE_H

#include "nmm/export.h"
#include <stdio.h>

struct imm_state;
struct nmm_codon_lprob;
struct nmm_codon_state;
struct nmm_profile;

NMM_API struct nmm_codon_lprob const* nmm_codon_state_codon_lprob(struct nmm_codon_state const* state);
NMM_API struct nmm_codon_state const* nmm_codon_state_create(char const* name, struct nmm_codon_lprob const* codonp);
NMM_API struct nmm_codon_state const* nmm_codon_state_derived(struct imm_state const* state);
NMM_API void                          nmm_codon_state_destroy(struct nmm_codon_state const* state);
NMM_API struct imm_state const*       nmm_codon_state_read(FILE* stream, struct nmm_profile const* prof);
NMM_API struct imm_state const*       nmm_codon_state_super(struct nmm_codon_state const* state);
NMM_API int nmm_codon_state_write(struct imm_state const* state, struct nmm_profile const* prof, FILE* stream);

#endif

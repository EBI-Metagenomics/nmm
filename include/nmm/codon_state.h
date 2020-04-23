#ifndef NMM_CODON_STATE_H
#define NMM_CODON_STATE_H

#include "nmm/export.h"

struct nmm_codon_lprob;
struct nmm_codon_state;

NMM_EXPORT struct nmm_codon_state const* nmm_codon_state_create(
    char const* name, struct nmm_codon_lprob const* codonp);
NMM_EXPORT void                    nmm_codon_state_destroy(struct nmm_codon_state const* state);
NMM_EXPORT struct imm_state const* nmm_codon_state_parent(struct nmm_codon_state const* state);
NMM_EXPORT struct nmm_codon_state const* nmm_codon_state_child(struct imm_state const* state);

#endif

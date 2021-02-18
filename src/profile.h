#ifndef PROFILE_H
#define PROFILE_H

#include <inttypes.h>

struct nmm_base_lprob;
struct nmm_codon_lprob;
struct nmm_codon_marg;
struct nmm_profile;

uint16_t profile_base_lprob_index(struct nmm_profile const* prof, struct nmm_base_lprob const* base_lprob);
uint16_t profile_codon_lprob_index(struct nmm_profile const* prof, struct nmm_codon_lprob const* codon_lprob);
uint16_t profile_codon_marg_index(struct nmm_profile const* prof, struct nmm_codon_marg const* codon_marg);

#endif

#ifndef MODEL_H
#define MODEL_H

#include <inttypes.h>

struct nmm_base_lprob;
struct nmm_codon_lprob;
struct nmm_codon_marg;
struct nmm_model;

uint16_t model_base_lprob_index(struct nmm_model const* io, struct nmm_base_lprob const* base_lprob);
uint16_t model_codon_lprob_index(struct nmm_model const* io, struct nmm_codon_lprob const* codon_lprob);
uint16_t model_codon_marg_index(struct nmm_model const* io, struct nmm_codon_marg const* codon_marg);

#endif

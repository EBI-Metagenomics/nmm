#ifndef MODEL_H
#define MODEL_H

#include <inttypes.h>

struct nmm_base_lprob;
struct nmm_codon_lprob;
struct nmm_codon_marg;
struct nmm_model;

uint16_t model_basep_index(struct nmm_model const* io, struct nmm_base_lprob const* basep);
uint16_t model_codonp_index(struct nmm_model const* io, struct nmm_codon_lprob const* codonp);
uint16_t model_codont_index(struct nmm_model const* io, struct nmm_codon_marg const* codont);

#endif

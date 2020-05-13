#ifndef MODEL_H
#define MODEL_H

#include <inttypes.h>

struct nmm_base_table;
struct nmm_codon_lprob;
struct nmm_codon_table;
struct nmm_model;

uint32_t model_baset_index(struct nmm_model const* io, struct nmm_base_table const* baset);
uint32_t model_codonp_index(struct nmm_model const* io, struct nmm_codon_lprob const* codonp);
uint32_t model_codont_index(struct nmm_model const* io, struct nmm_codon_table const* codont);
struct nmm_model* model_new(void);

#endif

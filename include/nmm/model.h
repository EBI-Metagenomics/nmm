#ifndef NMM_MODEL_H
#define NMM_MODEL_H

#include "nmm/export.h"
#include <inttypes.h>
#include <stdio.h>

struct imm_dp;
struct imm_hmm;
struct imm_model;
struct nmm_base_lprob;
struct nmm_codon_table;
struct nmm_model;

NMM_API struct imm_abc const*         nmm_model_abc(struct nmm_model const* model);
NMM_API struct nmm_base_lprob const*  nmm_model_base_lprob(struct nmm_model const* model, uint16_t index);
NMM_API struct nmm_codon_lprob const* nmm_model_codon_lprob(struct nmm_model const* model, uint16_t index);
NMM_API struct nmm_codon_table const* nmm_model_codon_table(struct nmm_model const* model, uint16_t index);
NMM_API struct nmm_model const*       nmm_model_create(struct imm_hmm* hmm, struct imm_dp const* dp);
NMM_API void                          nmm_model_destroy(struct nmm_model const* model);
NMM_API uint16_t                      nmm_model_nbase_lprobs(struct nmm_model const* model);
NMM_API uint16_t                      nmm_model_ncodon_lprobs(struct nmm_model const* model);
NMM_API uint16_t                      nmm_model_ncodon_tables(struct nmm_model const* model);
NMM_API struct imm_hmm*               nmm_model_hmm(struct nmm_model const* model);
NMM_API struct imm_dp const*          nmm_model_dp(struct nmm_model const* model);
NMM_API struct imm_state const*       nmm_model_state(struct nmm_model const* model, uint16_t i);
NMM_API uint16_t                      nmm_model_nstates(struct nmm_model const* model);
NMM_API struct nmm_model const*       nmm_model_read(FILE* stream);
NMM_API int                           nmm_model_write(struct nmm_model const* io, FILE* stream);

#endif

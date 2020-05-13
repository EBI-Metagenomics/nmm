#ifndef NMM_MODEL_H
#define NMM_MODEL_H

#include "nmm/export.h"
#include <inttypes.h>
#include <stdio.h>

struct imm_dp;
struct imm_hmm;
struct imm_model;
struct nmm_base_table;
struct nmm_codon_table;
struct nmm_model;

NMM_API struct nmm_model const*      nmm_model_create(struct imm_hmm* hmm, struct imm_dp const* dp);
NMM_API struct nmm_model*            nmm_model_derived(struct imm_model* model);
NMM_API struct nmm_model const*      nmm_model_derived_c(struct imm_model const* model);
NMM_API void                         nmm_model_destroy(struct nmm_model const* model);
NMM_API uint32_t                     nmm_model_nbase_tables(struct nmm_model const* model);
NMM_API uint32_t                     nmm_model_ncodon_lprobs(struct nmm_model const* model);
NMM_API uint32_t                     nmm_model_ncodon_tables(struct nmm_model const* model);
NMM_API struct imm_model const*      nmm_model_super(struct nmm_model const* model);
NMM_API struct nmm_base_table const* nmm_model_base_table(struct nmm_model const* model,
                                                          uint32_t                index);
NMM_API struct nmm_codon_lprob const* nmm_model_codon_lprob(struct nmm_model const* model,
                                                            uint32_t                index);
NMM_API struct nmm_codon_table const* nmm_model_codon_table(struct nmm_model const* model,
                                                            uint32_t                index);
NMM_API int                           nmm_model_write(struct nmm_model const* io, FILE* stream);

#endif

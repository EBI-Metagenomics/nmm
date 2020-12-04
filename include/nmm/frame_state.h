#ifndef NMM_FRAME_STATE_H
#define NMM_FRAME_STATE_H

#include "nmm/export.h"
#include <stdio.h>

struct imm_abc;
struct imm_seq;
struct imm_state;
struct nmm_base_table;
struct nmm_codon;
struct nmm_codon_table;
struct nmm_frame_state;
struct nmm_model;

NMM_API struct nmm_base_table const* nmm_frame_state_base_table(
    struct nmm_frame_state const* state);
NMM_API struct nmm_codon_table const* nmm_frame_state_codon_table(
    struct nmm_frame_state const* state);
NMM_API struct nmm_frame_state const* nmm_frame_state_create(char const*                   name,
                                                             struct nmm_base_table const*  baset,
                                                             struct nmm_codon_table const* codont,
                                                             float                         epsilon);
NMM_API float nmm_frame_state_decode(struct nmm_frame_state const* state, struct imm_seq const* seq,
                                     struct nmm_codon* codon);
NMM_API struct nmm_frame_state const* nmm_frame_state_derived(struct imm_state const* state);
NMM_API void                          nmm_frame_state_destroy(struct nmm_frame_state const* state);
NMM_API float                         nmm_frame_state_epsilon(struct nmm_frame_state const* state);
NMM_API float                   nmm_frame_state_lposterior(struct nmm_frame_state const* state,
                                                           struct nmm_codon const* codon, struct imm_seq const* seq);
NMM_API struct imm_state const* nmm_frame_state_read(FILE* stream, struct nmm_model const* model);
NMM_API struct imm_state const* nmm_frame_state_super(struct nmm_frame_state const* state);
NMM_API int nmm_frame_state_write(struct imm_state const* state, struct nmm_model const* model,
                                  FILE* stream);

#endif

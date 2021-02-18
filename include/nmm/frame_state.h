#ifndef NMM_FRAME_STATE_H
#define NMM_FRAME_STATE_H

#include "imm/imm.h"
#include "nmm/export.h"
#include <stdio.h>

struct imm_seq;
struct imm_state;
struct nmm_base_lprob;
struct nmm_codon;
struct nmm_codon_marg;
struct nmm_frame_state;
struct nmm_profile;

NMM_API struct nmm_base_lprob const*  nmm_frame_state_base_lprob(struct nmm_frame_state const* state);
NMM_API struct nmm_codon_marg const*  nmm_frame_state_codon_marg(struct nmm_frame_state const* state);
NMM_API struct nmm_frame_state const* nmm_frame_state_create(char const* name, struct nmm_base_lprob const* basep,
                                                             struct nmm_codon_marg const* codont, imm_float epsilon);
NMM_API imm_float nmm_frame_state_decode(struct nmm_frame_state const* state, struct imm_seq const* seq,
                                         struct nmm_codon* codon);
NMM_API struct nmm_frame_state const* nmm_frame_state_derived(struct imm_state const* state);
NMM_API void                          nmm_frame_state_destroy(struct nmm_frame_state const* state);
NMM_API imm_float                     nmm_frame_state_epsilon(struct nmm_frame_state const* state);
NMM_API imm_float nmm_frame_state_lposterior(struct nmm_frame_state const* state, struct nmm_codon const* codon,
                                             struct imm_seq const* seq);
NMM_API struct imm_state const* nmm_frame_state_read(FILE* stream, struct nmm_profile const* prof);
NMM_API struct imm_state const* nmm_frame_state_super(struct nmm_frame_state const* state);
NMM_API int nmm_frame_state_write(struct imm_state const* state, struct nmm_profile const* prof, FILE* stream);

#endif

#ifndef FRAME_STATE_H
#define FRAME_STATE_H

#include <stdio.h>

struct imm_state;
struct nmm_base_table;
struct nmm_frame_state;
struct nmm_model;

struct nmm_base_table const*  frame_state_baset(struct nmm_frame_state const* state);
struct nmm_codon_table const* frame_state_codont(struct nmm_frame_state const* state);

#endif

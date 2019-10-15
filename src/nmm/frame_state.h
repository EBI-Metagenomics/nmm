#ifndef IMM_STATE_FRAME_H
#define IMM_STATE_FRAME_H

#include "state/state.h"

struct frame_state;

void frame_state_create(struct imm_state *state, const double *base_emiss_lprobs,
                        const struct imm_codon *codon, double epsilon);
double frame_state_emiss_lprob(const struct imm_state *state, const char *seq,
                               int seq_len);
int frame_state_normalize(struct imm_state *state);
int frame_state_min_seq(const struct imm_state *state);
int frame_state_max_seq(const struct imm_state *state);
void frame_state_destroy(struct imm_state *state);

#endif

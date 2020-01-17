#ifndef NMM_FRAME_STATE_H
#define NMM_FRAME_STATE_H

#include "nmm/api.h"

struct imm_abc;
struct nmm_baset;
struct nmm_codon;
struct nmm_codont;
struct nmm_frame_state;

NMM_API struct nmm_frame_state* nmm_frame_state_create(char const*              name,
                                                       struct nmm_baset const*  baset,
                                                       struct nmm_codont const* codont,
                                                       double                   epsilon);
NMM_API double nmm_frame_state_lposterior(struct nmm_frame_state const* state,
                                          struct nmm_codon const* codon, char const* seq,
                                          int seq_len);
NMM_API double nmm_frame_state_decode(struct nmm_frame_state const* state, char const* seq,
                                      int seq_len, struct nmm_codon* codon);
NMM_API void   nmm_frame_state_destroy(struct nmm_frame_state* state);

#endif

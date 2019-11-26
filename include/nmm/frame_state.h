#ifndef NMM_FRAME_STATE_H
#define NMM_FRAME_STATE_H

#include "nmm/api.h"
#include "nmm/codont.h"

struct imm_abc;
struct nmm_base;
struct nmm_frame_state;
struct nmm_codont;

NMM_API struct nmm_frame_state* nmm_frame_state_create(char const*              name,
                                                       struct nmm_base const*   base,
                                                       struct nmm_codont const* codon,
                                                       double                   epsilon);
NMM_API double                  nmm_frame_state_lposterior(struct nmm_frame_state* state,
                                                           struct nmm_codon const* ccode, char const* seq,
                                                           int seq_len);
NMM_API double nmm_frame_state_decode(struct nmm_frame_state* state, char const* seq,
                                      int seq_len, struct nmm_codon* ccode);
NMM_API void   nmm_frame_state_destroy(struct nmm_frame_state* state);

#endif

#ifndef NMM_FRAME_STATE_H
#define NMM_FRAME_STATE_H

#include "nmm/api.h"
#include "nmm/codon.h"

struct imm_abc;
struct nmm_base;
struct nmm_frame_state;
struct nmm_codon;

NMM_API struct nmm_frame_state* nmm_frame_state_create(char const*             name,
                                                       struct nmm_base const*  base,
                                                       struct nmm_codon const* codon,
                                                       double                  epsilon);
NMM_API double                  nmm_frame_state_posterior(struct nmm_frame_state* state,
                                                          struct nmm_ccode const* ccode, char const* seq,
                                                          int seq_len);
NMM_API void                    nmm_frame_state_destroy(struct nmm_frame_state* state);

#endif

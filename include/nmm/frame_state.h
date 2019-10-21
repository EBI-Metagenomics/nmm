#ifndef NMM_FRAME_STATE_H_API
#define NMM_FRAME_STATE_H_API

#include "nmm/api.h"

struct imm_abc;
struct nmm_base;
struct nmm_frame_state;
struct nmm_codon;

NMM_API struct nmm_frame_state *nmm_frame_state_create(const char *name,
                                                       const struct nmm_base *base,
                                                       const struct nmm_codon *codon,
                                                       double epsilon);
NMM_API void nmm_frame_state_destroy(struct nmm_frame_state *state);

#endif

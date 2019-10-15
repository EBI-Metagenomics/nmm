#ifndef NMM_FRAME_STATE_H_API
#define NMM_FRAME_STATE_H_API

#include "nmm/api.h"
#include "nmm/codon.h"

struct imm_abc;
struct nmm_frame_state;

NMM_API struct nmm_frame_state *nmm_frame_state_create(const char *name,
                                                       const struct imm_abc *bases,
                                                       const double *base_lprobs,
                                                       const struct nmm_codon *codon,
                                                       double epsilon);

NMM_API void nmm_frame_state_destroy(struct nmm_frame_state *state);
NMM_API int nmm_frame_state_normalize(struct nmm_frame_state *state);

#endif

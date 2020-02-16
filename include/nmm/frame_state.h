#ifndef NMM_FRAME_STATE_H
#define NMM_FRAME_STATE_H

#include "nmm/api.h"

struct imm_abc;
struct imm_seq;
struct nmm_baset;
struct nmm_codon;
struct nmm_codon_table;
struct nmm_frame_state;

NMM_API struct nmm_frame_state const* nmm_frame_state_create(
    char const* name, struct nmm_baset const* baset, struct nmm_codon_table const* codont,
    double epsilon);
NMM_API double nmm_frame_state_lposterior(struct nmm_frame_state const* state,
                                          struct nmm_codon const*       codon,
                                          struct imm_seq const*         seq);
NMM_API double nmm_frame_state_decode(struct nmm_frame_state const* state,
                                      struct imm_seq const* seq, struct nmm_codon* codon);
NMM_API void   nmm_frame_state_destroy(struct nmm_frame_state const* state);

#endif

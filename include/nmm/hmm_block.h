#ifndef NMM_HMM_BLOCK_H
#define NMM_HMM_BLOCK_H

#include "nmm/export.h"
#include <inttypes.h>

struct imm_dp;
struct imm_hmm;
struct nmm_hmm_block;

NMM_API struct imm_dp const*    nmm_hmm_block_dp(struct nmm_hmm_block const* block);
NMM_API struct imm_hmm*         nmm_hmm_block_hmm(struct nmm_hmm_block const* block);
NMM_API uint16_t                nmm_hmm_block_nstates(struct nmm_hmm_block const* block);
NMM_API struct imm_state const* nmm_hmm_block_state(struct nmm_hmm_block const* block, uint16_t i);

#endif

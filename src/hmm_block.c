#include "nmm/hmm_block.h"
#include "hmm_block.h"
#include "imm/imm.h"

struct imm_dp const* nmm_hmm_block_dp(struct nmm_hmm_block const* block) { return imm_hmm_block_dp(block->super); }

struct imm_hmm* nmm_hmm_block_hmm(struct nmm_hmm_block const* block) { return imm_hmm_block_hmm(block->super); }

uint16_t nmm_hmm_block_nstates(struct nmm_hmm_block const* block) { return imm_hmm_block_nstates(block->super); }

struct imm_state const* nmm_hmm_block_state(struct nmm_hmm_block const* block, uint16_t i)
{
    return imm_hmm_block_state(block->super, i);
}

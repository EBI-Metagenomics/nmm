#ifndef CODON_STATE_H
#define CODON_STATE_H

#include <stdio.h>

struct imm_state;
struct nmm_codon_state;
struct nmm_io;

struct nmm_codon_lprob const* codon_state_codonp(struct nmm_codon_state const* state);
struct imm_state const*       codon_state_read(FILE* stream, struct nmm_io const* io);

#endif

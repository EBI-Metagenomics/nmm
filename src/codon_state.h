#ifndef CODON_STATE_H
#define CODON_STATE_H

#include <stdio.h>

struct imm_state;

struct imm_state const* codon_state_read(FILE* stream);

#endif

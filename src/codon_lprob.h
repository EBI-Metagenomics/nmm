#ifndef CODON_LPROB_H
#define CODON_LPROB_H

#include <stdio.h>

struct nmm_codon_lprob;

int codon_lprob_write(struct nmm_codon_lprob const* codonp, FILE* stream);

#endif

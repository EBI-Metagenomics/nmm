#ifndef CODON_LPROB_H
#define CODON_LPROB_H

#include <stdio.h>

struct nmm_base_abc;
struct nmm_codon_lprob;

struct nmm_codon_lprob const* codon_lprob_read(FILE* stream, struct nmm_base_abc const* base_abc);
int                           codon_lprob_write(struct nmm_codon_lprob const* codonp, FILE* stream);

#endif

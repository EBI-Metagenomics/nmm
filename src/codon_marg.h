#ifndef CODON_MARG_H
#define CODON_MARG_H

#include <stdio.h>

struct nmm_base_abc;
struct nmm_codon_marg;

struct nmm_codon_marg const* codon_marg_read(FILE* stream, struct nmm_base_abc const* base_abc);
int                          codon_marg_write(struct nmm_codon_marg const* codont, FILE* stream);

#endif

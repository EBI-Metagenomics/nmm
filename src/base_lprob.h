#ifndef BASE_LPROB_H
#define BASE_LPROB_H

#include <stdio.h>

struct nmm_base_abc;
struct nmm_base_lprob;

struct nmm_base_lprob const* base_lprob_read(FILE* stream, struct nmm_base_abc const* abc);
int                          base_lprob_write(struct nmm_base_lprob const* basep, FILE* stream);

#endif

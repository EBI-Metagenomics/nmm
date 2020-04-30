#ifndef CODON_TABLE_H
#define CODON_TABLE_H

#include <stdio.h>

struct nmm_base_abc;
struct nmm_codon_table;

struct nmm_codon_table const* codon_table_read(FILE* stream, struct nmm_base_abc const* base_abc);
int                           codon_table_write(struct nmm_codon_table const* codont, FILE* stream);

#endif

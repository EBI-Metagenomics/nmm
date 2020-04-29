#ifndef CODON_TABLE_H
#define CODON_TABLE_H

#include <stdio.h>

struct nmm_codon_table;

int codon_table_write(struct nmm_codon_table const* codont, FILE* stream);

#endif

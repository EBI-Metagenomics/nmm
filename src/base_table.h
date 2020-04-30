#ifndef BASE_TABLE_H
#define BASE_TABLE_H

#include <stdio.h>

struct nmm_base_abc;
struct nmm_base_table;

int base_table_read(FILE* stream, struct nmm_base_abc const* base_abc);
int base_table_write(struct nmm_base_table const* baset, FILE* stream);

#endif

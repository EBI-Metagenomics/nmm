#ifndef BASE_TABLE_H
#define BASE_TABLE_H

#include <stdio.h>

struct nmm_base_table;

int base_table_write(struct nmm_base_table const* baset, FILE* stream);

#endif

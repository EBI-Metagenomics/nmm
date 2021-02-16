#ifndef FIO_H
#define FIO_H

#include <stdint.h>
#include <stdio.h>

int     nmm_fseek(FILE* stream, int64_t offset, int origin);
int64_t nmm_ftell(FILE* stream);

#endif

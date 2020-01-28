#ifndef NMM_FREE_H
#define NMM_FREE_H

#include <stdlib.h>

static inline void free_c(void const* p) { free((void*)p); }

#endif

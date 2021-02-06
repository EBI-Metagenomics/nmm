#ifndef HELPER_H
#define HELPER_H

#include "imm/imm.h"

#define LOG(x) ((imm_float)log((x)))

static inline imm_float zero(void) { return imm_lprob_zero(); }

#endif

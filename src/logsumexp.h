#ifndef NMM_LOGSUMEXP_H
#define NMM_LOGSUMEXP_H

#include "imm/imm.h"
#include "logaddexp.h"

static inline double logsumexp(double const* arr, unsigned len)
{
    double r = IMM_LPROB_ZERO;
    for (size_t i = 0; i < len; ++i)
        r = logaddexp(r, arr[i]);
    return r;
}

#endif

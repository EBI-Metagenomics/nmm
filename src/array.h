#ifndef NMM_ARRAY_H
#define NMM_ARRAY_H

#include "imm/imm.h"
#include <stdlib.h>

struct array3d
{
    double*  values;
    unsigned strides[3];
};

static inline struct array3d array3d_create(unsigned dim0, unsigned dim1, unsigned dim2)
{
    size_t         len = (size_t)(dim0 * dim1 * dim2);
    struct array3d arr = {malloc(sizeof(double) * len), {dim1 * dim2, dim2, 1}};
    return arr;
}

static inline void array3d_set(struct array3d* arr, unsigned i0, unsigned i1, unsigned i2,
                               double val)
{
    unsigned const* s = arr->strides;
    arr->values[i0 * s[0] + i1 * s[1] + i2 * s[2]] = val;
}

static inline void array3d_fill(struct array3d* arr, double val)
{
    size_t const len = arr->strides[0] * arr->strides[1] * arr->strides[2];
    for (size_t i = 0; i < len; ++i)
        arr->values[i] = val;
}

static inline double array3d_get(struct array3d const* arr, unsigned i0, unsigned i1,
                                 unsigned i2)
{
    unsigned const* s = arr->strides;
    return arr->values[i0 * s[0] + i1 * s[1] + i2 * s[2]];
}

static inline void array3d_destroy(struct array3d arr) { free(arr.values); }

static inline int array3d_normalize(struct array3d const* arr)
{
    unsigned len = arr->strides[0] * arr->strides[1] * arr->strides[2];
    return imm_lprob_normalize(arr->values, len);
}

#endif

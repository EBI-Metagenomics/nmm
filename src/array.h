#ifndef NMM_ARRAY_H
#define NMM_ARRAY_H

#include <stdlib.h>

struct array3d
{
    double* values;
    int     strides[3];
};

static inline struct array3d array3d_create(int dim0, int dim1, int dim2)
{
    size_t         len = (size_t)(dim0 * dim1 * dim2);
    struct array3d arr = {malloc(sizeof(double) * len), {dim1 * dim2, dim2, 1}};
    return arr;
}

static inline void array3d_set(struct array3d* arr, int i0, int i1, int i2, double val)
{
    int const* s = arr->strides;
    arr->values[i0 * s[0] + i1 * s[1] + i2 * s[2]] = val;
}

static inline void array3d_fill(struct array3d* arr, double val)
{
    int const len = arr->strides[0] * arr->strides[1] * arr->strides[2];
    for (int i = 0; i < len; ++i)
        arr->values[i] = val;
}

static inline double array3d_get(struct array3d const* arr, int i0, int i1, int i2)
{
    int const* s = arr->strides;
    return arr->values[i0 * s[0] + i1 * s[1] + i2 * s[2]];
}

static inline void array3d_destroy(struct array3d arr) { free(arr.values); }

#endif

#ifndef NMM_ARRAY_H
#define NMM_ARRAY_H

#include "imm/imm.h"
#include "nmm/export.h"
#include <stdio.h>
#include <stdlib.h>

struct nmm_array3d
{
    unsigned strides[3];
    float*   values;
};

struct nmm_array3d_idx
{
    unsigned i0;
    unsigned i1;
    unsigned i2;
};

static inline struct nmm_array3d nmm_array3d_create(unsigned dim0, unsigned dim1, unsigned dim2);
static inline void               nmm_array3d_destroy(struct nmm_array3d arr) { free(arr.values); }
static inline void               nmm_array3d_fill(struct nmm_array3d* arr, float val);
static inline float    nmm_array3d_get(struct nmm_array3d const* arr, struct nmm_array3d_idx idx);
static inline unsigned nmm_array3d_length(struct nmm_array3d const* arr);
static inline int      nmm_array3d_normalize(struct nmm_array3d const* arr);
NMM_API int            nmm_array3d_read(struct nmm_array3d* arr, FILE* stream);
static inline void nmm_array3d_set(struct nmm_array3d* arr, struct nmm_array3d_idx idx, float val);
NMM_API int        nmm_array3d_write(struct nmm_array3d const* arr, FILE* stream);

static inline struct nmm_array3d nmm_array3d_create(unsigned dim0, unsigned dim1, unsigned dim2)
{
    unsigned           len = dim0 * dim1 * dim2;
    struct nmm_array3d arr = {{dim1 * dim2, dim2, 1}, malloc(sizeof(*arr.values) * len)};
    return arr;
}

static inline void nmm_array3d_fill(struct nmm_array3d* arr, float val)
{
    for (unsigned i = 0; i < nmm_array3d_length(arr); ++i)
        arr->values[i] = val;
}

static inline float nmm_array3d_get(struct nmm_array3d const* arr, struct nmm_array3d_idx idx)
{
    unsigned const* s = arr->strides;
    return arr->values[idx.i0 * s[0] + idx.i1 * s[1] + idx.i2 * s[2]];
}

static inline unsigned nmm_array3d_length(struct nmm_array3d const* arr)
{
    return arr->strides[0] * arr->strides[1] * arr->strides[2];
}

static inline int nmm_array3d_normalize(struct nmm_array3d const* arr)
{
    return imm_lprob_normalize(arr->values, nmm_array3d_length(arr));
}

static inline void nmm_array3d_set(struct nmm_array3d* arr, struct nmm_array3d_idx idx, float val)
{
    unsigned const* s = arr->strides;
    arr->values[idx.i0 * s[0] + idx.i1 * s[1] + idx.i2 * s[2]] = val;
}

#endif

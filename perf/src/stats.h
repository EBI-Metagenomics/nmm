#ifndef PERF_STATS_H
#define PERF_STATS_H

#include "imm/imm.h"
#include <math.h>
#include <stdlib.h>

static inline imm_float __mean(imm_float const* arr, unsigned n);
static inline imm_float __median(imm_float const* sorted_arr, unsigned n);
static inline imm_float __sample_std(imm_float const* arr, unsigned n, imm_float mean);
static int              __comp(void const* a, void const* b);

struct stats
{
    imm_float mean;   /**< Arithmetic mean. */
    imm_float median; /**< Arithmetic median. */
    imm_float sem;    /**< Standard error of the mean. */
};

static inline struct stats compute_stats(imm_float* arr, unsigned n)
{
    qsort(arr, n, sizeof(*arr), __comp);
    imm_float mean = __mean(arr, n);
    imm_float median = __median(arr, n);
    imm_float std = __sample_std(arr, n, mean);
    imm_float sem = (imm_float)(std / sqrt(n));
    return (struct stats){mean, median, sem};
}

static inline imm_float __mean(imm_float const* arr, unsigned n)
{
    if (n == 0)
        return NAN;

    imm_float total = arr[0];

    for (unsigned i = 1; i < n; ++i)
        total += arr[i];

    return (imm_float)((double)total / n);
}

static inline imm_float __median(imm_float const* sorted_arr, unsigned n)
{
    if (n == 0)
        return NAN;

    if (n % 2 == 1)
        return sorted_arr[n / 2];

    unsigned i = n / 2;
    return (sorted_arr[i - 1] + sorted_arr[i]) / 2;
}

static inline imm_float __sample_std(imm_float const* arr, unsigned n, imm_float mean)
{
    if (n <= 1)
        return NAN;

    imm_float var = 0.0;

    for (unsigned i = 0; i < n; ++i)
        var += (arr[i] - mean) * (arr[i] - mean);

    return (imm_float)sqrt((double)var / (n - 1));
}

static int __comp(void const* a, void const* b)
{
    imm_float fa = *(imm_float const*)a;
    imm_float fb = *(imm_float const*)b;
    return (fa > fb) - (fa < fb);
}

#endif

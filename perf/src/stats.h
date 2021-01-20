#ifndef PERF_STATS_H
#define PERF_STATS_H

#include <math.h>
#include <stdlib.h>

static inline double __mean(double const* arr, unsigned n);
static inline double __median(double const* sorted_arr, unsigned n);
static inline double __sample_std(double const* arr, unsigned n, double mean);
static int           __comp(void const* a, void const* b);

struct stats
{
    double mean;   /**< Arithmetic mean. */
    double median; /**< Arithmetic median. */
    double sem;    /**< Standard error of the mean. */
};

static inline struct stats compute_stats(double* arr, unsigned n)
{
    qsort(arr, n, sizeof(*arr), __comp);
    double mean = __mean(arr, n);
    double median = __median(arr, n);
    double std = __sample_std(arr, n, mean);
    double sem = std / sqrt(n);
    return (struct stats){mean, median, sem};
}

static inline double __mean(double const* arr, unsigned n)
{
    if (n == 0)
        return NAN;

    double total = arr[0];

    for (unsigned i = 1; i < n; ++i)
        total += arr[i];

    return total / n;
}

static inline double __median(double const* sorted_arr, unsigned n)
{
    if (n == 0)
        return NAN;

    if (n % 2 == 1)
        return sorted_arr[n / 2];

    unsigned i = n / 2;
    return (sorted_arr[i - 1] + sorted_arr[i]) / 2;
}

static inline double __sample_std(double const* arr, unsigned n, double mean)
{
    if (n <= 1)
        return NAN;

    double var = 0.0;

    for (unsigned i = 0; i < n; ++i)
        var += (arr[i] - mean) * (arr[i] - mean);

    return sqrt(var / (n - 1));
}

static int __comp(void const* a, void const* b)
{
    double fa = *(double const*)a;
    double fb = *(double const*)b;
    return (fa > fb) - (fa < fb);
}

#endif

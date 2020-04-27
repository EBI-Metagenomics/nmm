#ifndef LIB_KHASH_PTR_H
#define LIB_KHASH_PTR_H

#include "lib/khash.h"
#include <assert.h>
#include <stdint.h>

/*
 * Source: Linux kernel.
 *
 * This hash multiplies the input by a large odd number and takes the
 * high bits.  Since multiplication propagates changes to the most
 * significant end only, it is essential that the high bits of the
 * product be used for the hash value.
 *
 * Chuck Lever verified the effectiveness of this technique:
 * http://www.citi.umich.edu/techreports/reports/citi-tr-00-1.pdf
 *
 * Although a random odd number will do, it turns out that the golden
 * ratio phi = (sqrt(5)-1)/2, or its negative, has particularly nice
 * properties.  (See Knuth vol 3, section 6.4, exercise 9.)
 *
 * These are the negative, (1 - phi) = phi**2 = (3 - sqrt(5))/2,
 * which is very slightly easier to multiply by and makes no
 * difference to the hash distribution.
 */
#define GOLDEN_RATIO_64 0x61C8864680B583EBull

/* Source: http://graphics.stanford.edu/~seander/bithacks.html#BitReverseObvious */
static kh_inline uint64_t reverse_bits(uint64_t v)
{
    unsigned long r = v;                        // r will be reversed bits of v; first get LSB of v
    int           s = sizeof(v) * CHAR_BIT - 1; // extra shift needed at end

    for (v >>= 1; v; v >>= 1) {
        r <<= 1;
        r |= v & 1;
        s--;
    }
    r <<= s; // shift when v's highest bits are zero
    return r;
}

static kh_inline khint_t ptr_hash_func(void const* ptr)
{
    uint64_t val = (uint64_t)ptr;
    return (khint_t)reverse_bits(val * GOLDEN_RATIO_64);
}

#define ptr_hash_equal(a, b) ((a) == (b))

#define KHASH_MAP_INIT_PTR(name, khval_t)                                                          \
    KHASH_INIT(name, void const*, khval_t, 1, ptr_hash_func, ptr_hash_equal)

#endif

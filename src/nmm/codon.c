#include "nmm.h"
#include "imm.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

struct nmm_codon
{
    double emiss_lprobs[4 * 4 * 4];
};

int codon_check_range(int a, int b, int c);

struct nmm_codon *nmm_codon_create(void)
{
    return malloc(sizeof(struct nmm_codon));
}

struct nmm_codon *nmm_codon_clone(const struct nmm_codon *codon)
{
    struct nmm_codon *c = malloc(sizeof(struct nmm_codon));
    memcpy(c->emiss_lprobs, codon->emiss_lprobs, sizeof(double) * 4 * 4 * 4);
    return c;
}

void nmm_codon_set_lprob(struct nmm_codon *codon, int a, int b, int c,
                                   double lprob)
{
    codon_check_range(a, b, c);

    if (a < 0 || a > 3 || b < 0 || b > 3 || c < 0 || c > 3)
        imm_error("base index outside the range [0, 3]");

    codon->emiss_lprobs[4 * 4 * a + 4 * b + c] = lprob;
}

void nmm_codon_set_ninfs(struct nmm_codon *codon)
{
    for (int i = 0; i < 4 * 4 * 4; ++i)
        codon->emiss_lprobs[i] = -INFINITY;
}

double nmm_codon_get_lprob(const struct nmm_codon *codon, int a, int b,
                                     int c)
{
    if (codon_check_range(a, b, c))
        return NAN;

    return codon->emiss_lprobs[4 * 4 * a + 4 * b + c];
}

int nmm_codon_normalize(struct nmm_codon *codon)
{
    return imm_lognormalize(codon->emiss_lprobs, 4 * 4 * 4);
}

void nmm_codon_destroy(struct nmm_codon *codon)
{
    if (codon)
        free(codon);
}

int codon_check_range(int a, int b, int c)
{
    if (a < 0 || a > 3 || b < 0 || b > 3 || c < 0 || c > 3) {
        imm_error("base index outside the range [0, 3]");
        return -1;
    }
    return 0;
}

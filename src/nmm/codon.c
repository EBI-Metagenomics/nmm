#include "imm.h"
#include "nmm.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

struct nmm_codon
{
    const struct imm_abc *abc;
    double emiss_lprobs[4 * 4 * 4];
};

void codon_set_ninfs(struct nmm_codon *codon);

struct nmm_codon *nmm_codon_create(const struct imm_abc *abc)
{
    if (imm_abc_length(abc) != 4) {
        imm_error("alphabet length is not four");
        return NULL;
    }
    struct nmm_codon *codon = malloc(sizeof(struct nmm_codon));
    codon->abc = abc;
    codon_set_ninfs(codon);
    return codon;
}

struct nmm_codon *nmm_codon_clone(const struct nmm_codon *codon)
{
    struct nmm_codon *c = malloc(sizeof(struct nmm_codon));
    c->abc = codon->abc;
    memcpy(c->emiss_lprobs, codon->emiss_lprobs, sizeof(double) * 4 * 4 * 4);
    return c;
}

int nmm_codon_set_lprob(struct nmm_codon *codon, char a, char b, char c, double lprob)
{
    int idx[3] = {imm_abc_symbol_idx(codon->abc, a), imm_abc_symbol_idx(codon->abc, b),
                  imm_abc_symbol_idx(codon->abc, c)};

    if (idx[0] < 0 || idx[1] < 0 || idx[2] < 0) {
        imm_error("nucleotide not found");
        return 1;
    }

    codon->emiss_lprobs[4 * 4 * idx[0] + 4 * idx[1] + idx[2]] = lprob;
    return 0;
}

void codon_set_ninfs(struct nmm_codon *codon)
{
    for (int i = 0; i < 4 * 4 * 4; ++i)
        codon->emiss_lprobs[i] = -INFINITY;
}

double nmm_codon_get_lprob(const struct nmm_codon *codon, char a, char b, char c)
{
    int idx[3] = {imm_abc_symbol_idx(codon->abc, a), imm_abc_symbol_idx(codon->abc, b),
                  imm_abc_symbol_idx(codon->abc, c)};

    if (idx[0] < 0 || idx[1] < 0 || idx[2] < 0) {
        imm_error("nucleotide not found");
        return NAN;
    }

    return codon->emiss_lprobs[4 * 4 * idx[0] + 4 * idx[1] + idx[2]];
}

int nmm_codon_normalize(struct nmm_codon *codon)
{
    return imm_lognormalize(codon->emiss_lprobs, 4 * 4 * 4);
}

void nmm_codon_destroy(struct nmm_codon *codon)
{
    if (!codon)
        return;

    codon->abc = NULL;
    free(codon);
}

#include "imm/imm.h"
#include "nmm/nmm.h"
#include <math.h>
#include <stdlib.h>

struct nmm_codont
{
    const struct imm_abc* abc;
    double                emiss_lprobs[4 * 4 * 4];
};

static void codon_set_ninfs(struct nmm_codont* codon);

struct nmm_codont* nmm_codont_create(const struct imm_abc* abc)
{
    if (imm_abc_length(abc) != 4) {
        imm_error("alphabet length is not four");
        return NULL;
    }
    struct nmm_codont* codon = malloc(sizeof(struct nmm_codont));
    codon->abc = abc;
    codon_set_ninfs(codon);
    return codon;
}

int nmm_codont_set_lprob(struct nmm_codont* codon, struct nmm_codon const* ccode,
                         double lprob)
{
    int idx[3] = {imm_abc_symbol_idx(codon->abc, ccode->a),
                  imm_abc_symbol_idx(codon->abc, ccode->b),
                  imm_abc_symbol_idx(codon->abc, ccode->c)};

    if (idx[0] < 0 || idx[1] < 0 || idx[2] < 0) {
        imm_error("nucleotide not found");
        return 1;
    }

    codon->emiss_lprobs[4 * 4 * idx[0] + 4 * idx[1] + idx[2]] = lprob;
    return 0;
}

double nmm_codont_get_lprob(const struct nmm_codont* codon, struct nmm_codon const* ccode)
{
    int idx[3] = {imm_abc_symbol_idx(codon->abc, ccode->a),
                  imm_abc_symbol_idx(codon->abc, ccode->b),
                  imm_abc_symbol_idx(codon->abc, ccode->c)};

    if (idx[0] < 0 || idx[1] < 0 || idx[2] < 0) {
        imm_error("nucleotide not found");
        return NAN;
    }

    return codon->emiss_lprobs[4 * 4 * idx[0] + 4 * idx[1] + idx[2]];
}

int nmm_codont_normalize(struct nmm_codont* codon)
{
    return imm_lprob_normalize(codon->emiss_lprobs, 4 * 4 * 4);
}

void nmm_codont_destroy(struct nmm_codont* codon)
{
    if (!codon)
        return;

    codon->abc = NULL;
    free(codon);
}

const struct imm_abc* nmm_codont_get_abc(const struct nmm_codont* codon)
{
    return codon->abc;
}

static void codon_set_ninfs(struct nmm_codont* codon)
{
    for (int i = 0; i < 4 * 4 * 4; ++i)
        codon->emiss_lprobs[i] = -INFINITY;
}

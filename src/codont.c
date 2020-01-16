#include "imm/imm.h"
#include "nmm/nmm.h"
#include <math.h>
#include <stdlib.h>

struct nmm_codont
{
    const struct imm_abc* abc;
    double                emiss_lprobs[4 * 4 * 4];
};

static void set_zero_lprobs(double* lprobs);

struct nmm_codont* nmm_codont_create(const struct imm_abc*         abc,
                                     struct nmm_codon_lprob const* lprobs, int lprobs_length)
{
    if (imm_abc_length(abc) != 4) {
        imm_error("alphabet length is not four");
        return NULL;
    }

    struct nmm_codont* codont = malloc(sizeof(struct nmm_codont));
    codont->abc = abc;

    set_zero_lprobs(codont->emiss_lprobs);
    for (int i = 0; i < lprobs_length; ++i) {
        int idx[3] = {imm_abc_symbol_idx(abc, lprobs[i].codon.a),
                      imm_abc_symbol_idx(abc, lprobs[i].codon.b),
                      imm_abc_symbol_idx(abc, lprobs[i].codon.c)};

        if (idx[0] < 0 || idx[1] < 0 || idx[2] < 0) {
            free(codont);
            imm_error("nucleotide not found");
            return NULL;
        }

        codont->emiss_lprobs[4 * 4 * idx[0] + 4 * idx[1] + idx[2]] = lprobs[i].lprob;
    }

    return codont;
}

/* int nmm_codont_set_lprob(struct nmm_codont* codont, struct nmm_codon const* codon, */
/*                          double lprob) */
/* { */
/*     int idx[3] = {imm_abc_symbol_idx(codont->abc, codon->a), */
/*                   imm_abc_symbol_idx(codont->abc, codon->b), */
/*                   imm_abc_symbol_idx(codont->abc, codon->c)}; */

/*     if (idx[0] < 0 || idx[1] < 0 || idx[2] < 0) { */
/*         imm_error("nucleotide not found"); */
/*         return 1; */
/*     } */

/*     codont->emiss_lprobs[4 * 4 * idx[0] + 4 * idx[1] + idx[2]] = lprob; */
/*     return 0; */
/* } */

double nmm_codont_lprob(struct nmm_codont const* codont, struct nmm_codon const* codon)
{
    int idx[3] = {imm_abc_symbol_idx(codont->abc, codon->a),
                  imm_abc_symbol_idx(codont->abc, codon->b),
                  imm_abc_symbol_idx(codont->abc, codon->c)};

    if (idx[0] < 0 || idx[1] < 0 || idx[2] < 0) {
        imm_error("nucleotide not found");
        return imm_lprob_invalid();
    }

    return codont->emiss_lprobs[4 * 4 * idx[0] + 4 * idx[1] + idx[2]];
}

int nmm_codont_normalize(struct nmm_codont* codont)
{
    return imm_lprob_normalize(codont->emiss_lprobs, 4 * 4 * 4);
}

void nmm_codont_destroy(struct nmm_codont* codont)
{
    if (!codont)
        return;

    codont->abc = NULL;
    free(codont);
}

struct imm_abc const* nmm_codont_get_abc(const struct nmm_codont* codont)
{
    return codont->abc;
}

static void set_zero_lprobs(double* lprobs)
{
    double const zero_lprob = imm_lprob_zero();
    for (int i = 0; i < 4 * 4 * 4; ++i)
        lprobs[i] = zero_lprob;
}

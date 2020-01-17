#include "array.h"
#include "codon_lprob.h"
#include "imm/imm.h"
#include "nmm/nmm.h"
#include <math.h>
#include <stdlib.h>

#define ASCII_LAST_STD 127
#define NBASES NMM_CODON_NBASES

struct nmm_codont
{
    const struct imm_abc* abc;
    struct codon_lprob*   lprob;
};

struct nmm_codont* nmm_codont_create(const struct imm_abc*         abc,
                                     struct nmm_codon_lprob const* lprobs, int lprobs_length)
{
    if (imm_abc_length(abc) != NMM_CODON_NBASES) {
        imm_error("alphabet length is not four");
        return NULL;
    }

    struct nmm_codont* codont = malloc(sizeof(struct nmm_codont));
    codont->abc = abc;
    codont->lprob = nmm_codon_lprob_create(abc, lprobs, lprobs_length);

    if (!codont->lprob) {
        nmm_codon_lprob_destroy(codont->lprob);
        free(codont);
        return NULL;
    }

    return codont;
}

double nmm_codont_lprob(struct nmm_codont const* codont, struct nmm_codon const* codon)
{
    return codon_lprob_get(codont->lprob, codon);
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

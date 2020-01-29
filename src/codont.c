#include "nmm/base.h"
#include "array.h"
#include "codon_lmarg.h"
#include "imm/imm.h"
#include "nmm/nmm.h"
#include <math.h>
#include <stdlib.h>

struct nmm_codont
{
    struct imm_abc const* abc;
    struct codon_lmarg*   lmarg;
};

struct nmm_codont* nmm_codont_create(struct imm_abc const*    abc,
                                     struct nmm_codonp const* lprobs)
{
    if (imm_abc_length(abc) != NMM_BASE_SIZE) {
        imm_error("alphabet length is not four");
        return NULL;
    }

    struct nmm_codont* codont = malloc(sizeof(struct nmm_codont));
    codont->abc = abc;
    codont->lmarg = nmm_codon_lmarg_create(abc, lprobs);

    if (!codont->lmarg) {
        nmm_codon_lmarg_destroy(codont->lmarg);
        free(codont);
        return NULL;
    }

    return codont;
}

double nmm_codont_lprob(struct nmm_codont const* codont, struct nmm_codon const* codon)
{
    return codon_lmarg_get(codont->lmarg, codon);
}

void nmm_codont_destroy(struct nmm_codont* codont)
{
    if (!codont)
        return;

    codont->abc = NULL;
    nmm_codon_lmarg_destroy(codont->lmarg);
    free(codont);
}

struct imm_abc const* nmm_codont_get_abc(struct nmm_codont const* codont)
{
    return codont->abc;
}

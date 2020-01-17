#include "array.h"
#include "imm/imm.h"
#include "nmm/nmm.h"
#include <stdlib.h>

struct nmm_codonp
{
    struct imm_abc const* abc;
    struct array3d        lprobs;
};

struct nmm_codonp* nmm_codonp_create(struct imm_abc const* abc)
{
    if (imm_abc_length(abc) != NMM_CODON_NBASES) {
        imm_error("alphabet length is not four");
        return NULL;
    }

    struct nmm_codonp* codonp = malloc(sizeof(struct nmm_codonp));
    codonp->abc = NULL;

    codonp->lprobs = array3d_create(NMM_CODON_NBASES, NMM_CODON_NBASES, NMM_CODON_NBASES);
    array3d_fill(&codonp->lprobs, imm_lprob_zero());

    return codonp;
}

int nmm_codonp_set(struct nmm_codonp* codonp, struct nmm_codon const* codon, double lprob)
{

    int a = imm_abc_symbol_idx(codonp->abc, codon->a);
    int b = imm_abc_symbol_idx(codonp->abc, codon->b);
    int c = imm_abc_symbol_idx(codonp->abc, codon->c);

    if (a < 0 || b < 0 || c < 0) {
        imm_error("codon not found");
        return 1;
    }

    array3d_set(&codonp->lprobs, a, b, c, lprob);

    return 0;
}

double nmm_codonp_get(struct nmm_codonp const* codonp, struct nmm_codon const* codon)
{
    int a = imm_abc_symbol_idx(codonp->abc, codon->a);
    int b = imm_abc_symbol_idx(codonp->abc, codon->b);
    int c = imm_abc_symbol_idx(codonp->abc, codon->c);

    if (a < 0 || b < 0 || c < 0)
        return imm_lprob_invalid();

    return array3d_get(&codonp->lprobs, a, b, c);
}

int nmm_codonp_normalize(struct nmm_codonp* codonp) { return array3d_normalize(&codonp->lprobs); }

void nmm_codonp_destroy(struct nmm_codonp* codonp)
{
    if (!codonp) {
        imm_error("codonp should not be NULL");
        return;
    }

    codonp->abc = NULL;
    array3d_destroy(codonp->lprobs);
    free(codonp);
}

struct imm_abc const* nmm_codonp_get_abc(struct nmm_codonp const* codonp)
{
    return codonp->abc;
}

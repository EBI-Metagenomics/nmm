#include "nmm/codonp.h"
#include "array.h"
#include "free.h"
#include "imm/imm.h"
#include "nmm/codon.h"
#include <stdlib.h>

/** @file codonp.c
 * Compute the probability of emitting a codon.
 *
 * Let p(𝑋₁=𝚡₁,𝑋₂=𝚡₂,𝑋₃=𝚡₃) be the probability of emitting codon (𝚡₁,𝚡₂,𝚡₃), where 𝚡ᵢ𝜖𝒜.
 * This modules implements the computation of p(𝑋₁=𝚡₁,𝑋₂=𝚡₂,𝑋₃=𝚡₃).
 */

struct nmm_codonp
{
    struct imm_abc const* abc;
    /**
     * Pre-computed probability p(𝑋₁=𝚡₁,𝑋₂=𝚡₂,𝑋₃=𝚡₃).
     */
    struct array3d lprobs;
};

struct nmm_codonp* nmm_codonp_create(struct imm_abc const* abc)
{
    if (imm_abc_length(abc) != NMM_CODON_NBASES) {
        imm_error("alphabet length is not four");
        return NULL;
    }

    struct nmm_codonp* codonp = malloc(sizeof(struct nmm_codonp));
    codonp->abc = abc;

    codonp->lprobs = array3d_create(NMM_CODON_NBASES, NMM_CODON_NBASES, NMM_CODON_NBASES);
    array3d_fill(&codonp->lprobs, imm_lprob_zero());

    return codonp;
}

int nmm_codonp_set(struct nmm_codonp* codonp, struct nmm_codon const* codon, double lprob)
{
    int const a = imm_abc_symbol_idx(codonp->abc, codon->a);
    int const b = imm_abc_symbol_idx(codonp->abc, codon->b);
    int const c = imm_abc_symbol_idx(codonp->abc, codon->c);

    if (a < 0 || b < 0 || c < 0) {
        imm_error("codon not found");
        return 1;
    }

    array3d_set(&codonp->lprobs, (unsigned)a, (unsigned)b, (unsigned)c, lprob);

    return 0;
}

double nmm_codonp_get(struct nmm_codonp const* codonp, struct nmm_codon const* codon)
{
    int const a = imm_abc_symbol_idx(codonp->abc, codon->a);
    int const b = imm_abc_symbol_idx(codonp->abc, codon->b);
    int const c = imm_abc_symbol_idx(codonp->abc, codon->c);

    if (a < 0 || b < 0 || c < 0)
        return imm_lprob_invalid();

    return array3d_get(&codonp->lprobs, (unsigned)a, (unsigned)b, (unsigned)c);
}

int nmm_codonp_normalize(struct nmm_codonp* codonp)
{
    return array3d_normalize(&codonp->lprobs);
}

void nmm_codonp_destroy(struct nmm_codonp const* codonp)
{
    array3d_destroy(codonp->lprobs);
    free_c(codonp);
}

struct imm_abc const* nmm_codonp_get_abc(struct nmm_codonp const* codonp)
{
    return codonp->abc;
}

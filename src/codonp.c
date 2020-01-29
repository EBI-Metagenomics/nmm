#include "nmm/codonp.h"
#include "array.h"
#include "free.h"
#include "imm/imm.h"
#include "nmm/base.h"
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
    struct nmm_base const* base;
    /**
     * Pre-computed probability p(𝑋₁=𝚡₁,𝑋₂=𝚡₂,𝑋₃=𝚡₃).
     */
    struct array3d lprobs;
};

struct nmm_codonp* nmm_codonp_create(struct nmm_base const* base)
{
    struct nmm_codonp* codonp = malloc(sizeof(struct nmm_codonp));
    codonp->base = base;

    codonp->lprobs = array3d_create(NMM_BASE_SIZE, NMM_BASE_SIZE, NMM_BASE_SIZE);
    array3d_fill(&codonp->lprobs, imm_lprob_zero());

    return codonp;
}

int nmm_codonp_set(struct nmm_codonp* codonp, struct nmm_codon const* codon, double lprob)
{
    struct imm_abc const* abc = nmm_base_get_abc(codonp->base);
    int const             a = imm_abc_symbol_idx(abc, codon->a);
    int const             b = imm_abc_symbol_idx(abc, codon->b);
    int const             c = imm_abc_symbol_idx(abc, codon->c);

    if (a < 0 || b < 0 || c < 0) {
        imm_error("codon not found");
        return 1;
    }

    array3d_set(&codonp->lprobs, (unsigned)a, (unsigned)b, (unsigned)c, lprob);

    return 0;
}

double nmm_codonp_get(struct nmm_codonp const* codonp, struct nmm_codon const* codon)
{
    struct imm_abc const* abc = nmm_base_get_abc(codonp->base);
    int const             a = imm_abc_symbol_idx(abc, codon->a);
    int const             b = imm_abc_symbol_idx(abc, codon->b);
    int const             c = imm_abc_symbol_idx(abc, codon->c);

    if (a < 0 || b < 0 || c < 0) {
        imm_error("codon not found");
        return imm_lprob_invalid();
    }

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

struct nmm_base const* nmm_codonp_get_base(struct nmm_codonp const* codonp)
{
    return codonp->base;
}

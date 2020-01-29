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
    if (codonp->base != nmm_codon_get_base(codon)) {
        imm_error("bases must be the same");
        return 1;
    }

    struct nmm_triplet t = nmm_codon_get(codon);
    char const         any_symbol = imm_abc_any_symbol(nmm_base_get_abc(codonp->base));
    if (any_symbol == t.a || any_symbol == t.b || any_symbol == t.c) {
        imm_error("any-symbol is not allowed");
        return 1;
    }

    struct imm_abc const* abc = nmm_base_get_abc(codonp->base);
    struct array3d_idx    idx = {(unsigned)imm_abc_symbol_idx(abc, t.a),
                              (unsigned)imm_abc_symbol_idx(abc, t.b),
                              (unsigned)imm_abc_symbol_idx(abc, t.c)};

    array3d_set(&codonp->lprobs, idx, lprob);

    return 0;
}

double nmm_codonp_get(struct nmm_codonp const* codonp, struct nmm_codon const* codon)
{
    if (codonp->base != nmm_codon_get_base(codon)) {
        imm_error("bases must be the same");
        return imm_lprob_invalid();
    }

    struct nmm_triplet t = nmm_codon_get(codon);
    char const         any_symbol = imm_abc_any_symbol(nmm_base_get_abc(codonp->base));
    if (any_symbol == t.a || any_symbol == t.b || any_symbol == t.c) {
        imm_error("any-symbol is not allowed");
        return 1;
    }

    struct imm_abc const* abc = nmm_base_get_abc(codonp->base);
    struct array3d_idx    idx = {(unsigned)imm_abc_symbol_idx(abc, t.a),
                              (unsigned)imm_abc_symbol_idx(abc, t.b),
                              (unsigned)imm_abc_symbol_idx(abc, t.c)};

    return array3d_get(&codonp->lprobs, idx);
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

#include "codon_lprob.h"
#include "free.h"
#include "imm/imm.h"
#include "nmm/array.h"
#include "nmm/base_abc.h"
#include "nmm/codon.h"
#include "nmm/codon_lprob.h"
#include <stdbool.h>
#include <stdlib.h>

/** @file codonp.c
 * Compute the probability of emitting a codon.
 *
 * Let p(𝑋₁=𝚡₁,𝑋₂=𝚡₂,𝑋₃=𝚡₃) be the probability of emitting codon (𝚡₁,𝚡₂,𝚡₃), where 𝚡ᵢ𝜖𝒜.
 * This modules implements the computation of p(𝑋₁=𝚡₁,𝑋₂=𝚡₂,𝑋₃=𝚡₃).
 */

struct nmm_codon_lprob
{
    struct nmm_base_abc const* base_abc;
    /**
     * Pre-computed probability p(𝑋₁=𝚡₁,𝑋₂=𝚡₂,𝑋₃=𝚡₃).
     */
    struct nmm_array3d lprobs;
};

static inline bool triplet_has(struct nmm_triplet const triplet, char const symbol);
static inline struct nmm_array3d_idx triplet_index(struct imm_abc const*    abc,
                                                   struct nmm_triplet const triplet);

struct nmm_base_abc const* nmm_codon_lprob_abc(struct nmm_codon_lprob const* codonp)
{
    return codonp->base_abc;
}

struct nmm_codon_lprob* nmm_codon_lprob_create(struct nmm_base_abc const* abc)
{
    struct nmm_codon_lprob* codonp = malloc(sizeof(*codonp));
    codonp->base_abc = abc;

    codonp->lprobs = nmm_array3d_create(NMM_BASE_ABC_SIZE, NMM_BASE_ABC_SIZE, NMM_BASE_ABC_SIZE);
    nmm_array3d_fill(&codonp->lprobs, imm_lprob_zero());

    return codonp;
}

void nmm_codon_lprob_destroy(struct nmm_codon_lprob const* codonp)
{
    nmm_array3d_destroy(codonp->lprobs);
    free_c(codonp);
}

float nmm_codon_lprob_get(struct nmm_codon_lprob const* codonp, struct nmm_codon const* codon)
{
    if (codonp->base_abc != nmm_codon_abc(codon)) {
        imm_error("bases must be the same");
        return imm_lprob_invalid();
    }

    struct nmm_triplet const triplet = nmm_codon_get_triplet(codon);
    char const               any_symbol = imm_abc_any_symbol(nmm_base_abc_super(codonp->base_abc));
    if (triplet_has(triplet, any_symbol)) {
        imm_error("any-symbol is not allowed");
        return imm_lprob_invalid();
    }

    struct nmm_array3d_idx const idx = triplet_index(nmm_base_abc_super(codonp->base_abc), triplet);

    return nmm_array3d_get(&codonp->lprobs, idx);
}

int nmm_codon_lprob_normalize(struct nmm_codon_lprob* codonp)
{
    return nmm_array3d_normalize(&codonp->lprobs);
}

int nmm_codon_lprob_set(struct nmm_codon_lprob* codonp, struct nmm_codon const* codon,
                        float const lprob)
{
    if (codonp->base_abc != nmm_codon_abc(codon)) {
        imm_error("bases must be the same");
        return 1;
    }

    struct nmm_triplet const triplet = nmm_codon_get_triplet(codon);
    char const               any_symbol = imm_abc_any_symbol(nmm_base_abc_super(codonp->base_abc));
    if (triplet_has(triplet, any_symbol)) {
        imm_error("any-symbol is not allowed");
        return 1;
    }

    struct nmm_array3d_idx const idx = triplet_index(nmm_base_abc_super(codonp->base_abc), triplet);

    nmm_array3d_set(&codonp->lprobs, idx, lprob);

    return 0;
}

struct nmm_codon_lprob const* codon_lprob_read(FILE* stream, struct nmm_base_abc const* abc)
{
    struct nmm_codon_lprob* codonp = malloc(sizeof(*codonp));
    codonp->base_abc = abc;

    if (nmm_array3d_read(&codonp->lprobs, stream)) {
        imm_error("could not read lprobs");
        free_c(codonp);
        return NULL;
    }

    return codonp;
}

int codon_lprob_write(struct nmm_codon_lprob const* codonp, FILE* stream)
{
    if (nmm_array3d_write(&codonp->lprobs, stream)) {
        imm_error("could not write lprobs");
        return 1;
    }

    return 0;
}

static inline bool triplet_has(struct nmm_triplet const triplet, char const symbol)
{
    return symbol == triplet.a || symbol == triplet.b || symbol == triplet.c;
}

static inline struct nmm_array3d_idx triplet_index(struct imm_abc const*    abc,
                                                   struct nmm_triplet const triplet)
{
    return (struct nmm_array3d_idx){imm_abc_symbol_idx(abc, triplet.a),
                                    imm_abc_symbol_idx(abc, triplet.b),
                                    imm_abc_symbol_idx(abc, triplet.c)};
}

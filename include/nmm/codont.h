#ifndef NMM_CODONT_H
#define NMM_CODONT_H

#include "nmm/api.h"

#define NMM_CODON_NBASES 4

/** @file codont.h
 * Codon table module.
 *
 * A codon table is represented by an (immutable) object of type @ref nmm_codon.
 */

struct imm_abc;
struct nmm_codont;

struct nmm_codon
{
    char a;
    char b;
    char c;
};

#define NMM_CODON(a, b, c) ((struct nmm_codon){(a), (b), (c)})

struct nmm_codon_lprob
{
    struct nmm_codon codon;
    double           lprob;
};

NMM_API struct nmm_codont*    nmm_codont_create(struct imm_abc const*         abc,
                                                struct nmm_codon_lprob const* lprobs,
                                                int                           lprobs_length);
NMM_API double                nmm_codont_lprob(struct nmm_codont const* codont,
                                               struct nmm_codon const*  codon);
NMM_API void                  nmm_codont_destroy(struct nmm_codont* codont);
NMM_API struct imm_abc const* nmm_codont_get_abc(struct nmm_codont const* codont);

#endif

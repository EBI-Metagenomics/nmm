#ifndef NMM_CODON_LPROB_H
#define NMM_CODON_LPROB_H

#include "hide.h"
#include "imm/imm.h"
#include "nmm/nmm.h"
#include <stddef.h>

#define NSYMBOLS (NMM_CODON_NBASES + 1)
#define ASCII_LAST_STD 127

struct codon_lprob
{
    double lprob[NSYMBOLS * NSYMBOLS * NSYMBOLS];
    int    symbol_idx[ASCII_LAST_STD + 1];
};

HIDE void nmm_codon_lprob_init(struct codon_lprob*      codon_lprob,
                               struct nmm_codont const* codont);

static inline int codon_lprob_idx(struct codon_lprob const* codon_lprob,
                                      char const*               codon)
{
    return codon_lprob->symbol_idx[(size_t)codon[0]] +
           NSYMBOLS * codon_lprob->symbol_idx[(size_t)codon[1]] +
           NSYMBOLS * NSYMBOLS * codon_lprob->symbol_idx[(size_t)codon[2]];
}

static inline double codon_lprob_value(struct codon_lprob const* codon_lprob,
                                           char const*               codon)
{
    return codon_lprob->lprob[codon_lprob_idx(codon_lprob, codon)];
}

#endif

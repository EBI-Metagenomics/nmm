#ifndef NMM_CODON_LPROB_H
#define NMM_CODON_LPROB_H

#include "nmm/api.h"
#include "nmm/codon.h"

struct nmm_codon_lprob
{
    struct nmm_codon codon;
    double           lprob;
};

#endif

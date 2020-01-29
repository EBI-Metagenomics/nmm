#ifndef NMM_CODON_H
#define NMM_CODON_H

struct nmm_codon
{
    char a;
    char b;
    char c;
};

#define NMM_CODON(a, b, c) ((struct nmm_codon){(a), (b), (c)})

#endif

#include "cass.h"
#include "imm/imm.h"
#include "nmm/nmm.h"

#define CODON(a, b, c) NMM_CODON(a, b, c)

void test_codont_works(void);
void test_codont_fails(void);

int main(void)
{
    test_codont_works();
    test_codont_fails();
    return cass_status();
}

void test_codont_works(void)
{
    struct imm_abc* abc = imm_abc_create("ACGT", 'X');

    struct nmm_codon_lprob const lprobs[2] = {{CODON('A', 'T', 'G'), log(0.8)},
                                              {CODON('A', 'T', 'T'), log(0.1)}};
    struct nmm_codont*           codont = nmm_codont_create(abc, lprobs, 2);

    cass_cond(codont != NULL);

    cass_close(nmm_codont_lprob(codont, &CODON('A', 'T', 'G')), log(0.8));
    cass_close(nmm_codont_lprob(codont, &CODON('A', 'T', 'T')), log(0.1));
    cass_cond(imm_lprob_is_zero(nmm_codont_lprob(codont, &CODON('T', 'T', 'T'))));
    cass_cond(!imm_lprob_is_valid(nmm_codont_lprob(codont, &CODON('H', 'T', 'T'))));

    imm_abc_destroy(abc);
    nmm_codont_destroy(codont);
}

void test_codont_fails(void)
{
    struct imm_abc* abc = imm_abc_create("ACGT", 'X');

    struct nmm_codon_lprob const lprobs[3] = {{CODON('A', 'T', 'G'), log(0.8)},
                                              {CODON('A', 'T', 'T'), log(0.1)},
                                              {CODON('A', 'H', 'T'), log(0.1)}};
    struct nmm_codont*           codont = nmm_codont_create(abc, lprobs, 3);

    cass_cond(codont == NULL);

    imm_abc_destroy(abc);
    nmm_codont_destroy(codont);
}

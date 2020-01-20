#include "cass.h"
#include "imm/imm.h"
#include "nmm/nmm.h"

#define CODON(a, b, c) NMM_CODON(a, b, c)

void test_codont_nonmarginal(void);
void test_codont_marginal(void);

int main(void)
{
    test_codont_nonmarginal();
    test_codont_marginal();
    return cass_status();
}

void test_codont_nonmarginal(void)
{
    struct imm_abc* abc = imm_abc_create("ACGT", 'X');

    struct nmm_codonp* codonp = nmm_codonp_create(abc);
    cass_cond(nmm_codonp_set(codonp, &CODON('A', 'T', 'G'), log(0.8)) == 0);
    cass_cond(nmm_codonp_set(codonp, &CODON('A', 'T', 'T'), log(0.1)) == 0);
    cass_cond(nmm_codonp_set(codonp, &CODON('A', 'H', 'T'), log(0.1)) == 1);
    struct nmm_codont* codont = nmm_codont_create(abc, codonp);
    cass_cond(codont != NULL);
    nmm_codonp_destroy(codonp);

    cass_close(nmm_codont_lprob(codont, &CODON('A', 'T', 'G')), log(0.8));
    cass_close(nmm_codont_lprob(codont, &CODON('A', 'T', 'T')), log(0.1));
    cass_cond(imm_lprob_is_zero(nmm_codont_lprob(codont, &CODON('T', 'T', 'T'))));
    cass_cond(!imm_lprob_is_valid(nmm_codont_lprob(codont, &CODON('H', 'T', 'T'))));

    imm_abc_destroy(abc);
    nmm_codont_destroy(codont);
}

void test_codont_marginal(void)
{
    struct imm_abc* abc = imm_abc_create("ACGT", 'X');

    struct nmm_codonp* codonp = nmm_codonp_create(abc);
    nmm_codonp_set(codonp, &NMM_CODON('A', 'T', 'G'), log(0.8));
    nmm_codonp_set(codonp, &NMM_CODON('A', 'T', 'T'), log(0.1));
    struct nmm_codont* codont = nmm_codont_create(abc, codonp);
    cass_cond(codont != NULL);
    nmm_codonp_destroy(codonp);

    cass_close(nmm_codont_lprob(codont, &CODON('A', 'T', 'G')), log(0.8));
    cass_close(nmm_codont_lprob(codont, &CODON('A', 'T', 'T')), log(0.1));
    cass_close(nmm_codont_lprob(codont, &CODON('A', 'T', 'X')), log(0.9));
    cass_close(nmm_codont_lprob(codont, &CODON('A', 'X', 'X')), log(0.9));
    cass_close(nmm_codont_lprob(codont, &CODON('X', 'X', 'X')), log(0.9));
    cass_close(nmm_codont_lprob(codont, &CODON('X', 'T', 'X')), log(0.9));
    cass_close(nmm_codont_lprob(codont, &CODON('X', 'X', 'G')), log(0.8));
    cass_close(nmm_codont_lprob(codont, &CODON('X', 'X', 'T')), log(0.1));

    imm_abc_destroy(abc);
    nmm_codont_destroy(codont);
}

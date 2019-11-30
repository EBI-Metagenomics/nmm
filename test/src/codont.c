#include "cass.h"
#include "imm/imm.h"
#include "nmm/nmm.h"

void test_codont(void);

int main(void)
{
    test_codont();
    return cass_status();
}

void test_codont(void)
{
#define CODON(a, b, c) &NMM_CODON(a, b, c)
    struct imm_abc*    abc = imm_abc_create("ACGT", 'X');
    struct nmm_codont* codont = nmm_codont_create(abc);

    cass_cond(nmm_codont_set_lprob(codont, CODON('A', 'T', 'G'), log(0.8)) == 0);
    cass_cond(nmm_codont_set_lprob(codont, CODON('A', 'T', 'T'), log(0.1)) == 0);
    cass_cond(nmm_codont_set_lprob(codont, CODON('A', 'H', 'T'), log(0.1)) == 1);

    cass_close(nmm_codont_get_lprob(codont, CODON('A', 'T', 'G')), log(0.8));
    cass_close(nmm_codont_get_lprob(codont, CODON('A', 'T', 'T')), log(0.1));
    cass_cond(imm_lprob_is_zero(nmm_codont_get_lprob(codont, CODON('T', 'T', 'T'))));
    cass_cond(!imm_lprob_is_valid(nmm_codont_get_lprob(codont, CODON('H', 'T', 'T'))));

    nmm_codont_normalize(codont);

    cass_close(nmm_codont_get_lprob(codont, CODON('A', 'T', 'G')), log(0.8) - log(0.9));
    cass_close(nmm_codont_get_lprob(codont, CODON('A', 'T', 'T')), log(0.1) - log(0.9));
    cass_cond(imm_lprob_is_zero(nmm_codont_get_lprob(codont, CODON('T', 'T', 'T'))));

    imm_abc_destroy(abc);
    nmm_codont_destroy(codont);
#undef CODON
}

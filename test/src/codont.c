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
#define CCODE(a, b, c) &NMM_CCODE(a, b, c)
    struct imm_abc*    abc = imm_abc_create("ACGT");
    struct nmm_codont* codon = nmm_codont_create(abc);

    cass_cond(nmm_codont_set_lprob(codon, CCODE('A', 'T', 'G'), log(0.8)) == 0);
    cass_cond(nmm_codont_set_lprob(codon, CCODE('A', 'T', 'T'), log(0.1)) == 0);
    cass_cond(nmm_codont_set_lprob(codon, CCODE('A', 'H', 'T'), log(0.1)) == 1);

    cass_close(nmm_codont_get_lprob(codon, CCODE('A', 'T', 'G')), log(0.8));
    cass_close(nmm_codont_get_lprob(codon, CCODE('A', 'T', 'T')), log(0.1));
    cass_cond(imm_lprob_is_zero(nmm_codont_get_lprob(codon, CCODE('T', 'T', 'T'))));
    cass_cond(!imm_lprob_is_valid(nmm_codont_get_lprob(codon, CCODE('H', 'T', 'T'))));

    nmm_codont_normalize(codon);

    cass_close(nmm_codont_get_lprob(codon, CCODE('A', 'T', 'G')), log(0.8) - log(0.9));
    cass_close(nmm_codont_get_lprob(codon, CCODE('A', 'T', 'T')), log(0.1) - log(0.9));
    cass_cond(imm_lprob_is_zero(nmm_codont_get_lprob(codon, CCODE('T', 'T', 'T'))));

    imm_abc_destroy(abc);
    nmm_codont_destroy(codon);
#undef CCODE
}

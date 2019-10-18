#include "cass/cass.h"
#include "imm.h"
#include "nmm.h"

void test_codon(void);

int main(void)
{
    test_codon();
    return cass_status();
}

void test_codon(void)
{
    struct imm_abc *abc = imm_abc_create("ACGT");
    struct nmm_codon *codon = nmm_codon_create(abc);

    nmm_codon_set_lprob(codon, 'A', 'T', 'G', log(0.8));
    nmm_codon_set_lprob(codon, 'A', 'T', 'T', log(0.1));

    cass_close(nmm_codon_get_lprob(codon, 'A', 'T', 'G'), log(0.8));
    cass_close(nmm_codon_get_lprob(codon, 'A', 'T', 'T'), log(0.1));
    cass_condition(imm_isninf(nmm_codon_get_lprob(codon, 'T', 'T', 'T')));

    nmm_codon_normalize(codon);

    cass_close(nmm_codon_get_lprob(codon, 'A', 'T', 'G'), log(0.8) - log(0.9));
    cass_close(nmm_codon_get_lprob(codon, 'A', 'T', 'T'), log(0.1) - log(0.9));
    cass_condition(imm_isninf(nmm_codon_get_lprob(codon, 'T', 'T', 'T')));

    imm_abc_destroy(abc);
    nmm_codon_destroy(codon);
}

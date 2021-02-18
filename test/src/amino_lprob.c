#include "cass/cass.h"
#include "nmm/nmm.h"

void test_amino_lprob(void);

int main(void)
{
    test_amino_lprob();
    return cass_status();
}

void test_amino_lprob(void)
{
    struct nmm_amino_abc const*   amino_abc = nmm_amino_abc_create("ACDEFGHIKLMNPQRSTVWY", 'X');
    imm_float const               zero = imm_lprob_zero();
    imm_float const               lprobs[NMM_AMINO_ABC_SIZE] = {zero, imm_log(1), [19] = imm_log(19)};
    struct nmm_amino_lprob const* aminop = nmm_amino_lprob_create(amino_abc, lprobs);

    cass_cond(imm_lprob_is_zero(nmm_amino_lprob_get(aminop, 'A')));
    cass_close(nmm_amino_lprob_get(aminop, 'C'), log(1));
    cass_close(nmm_amino_lprob_get(aminop, 'Y'), log(19));

    nmm_amino_abc_destroy(amino_abc);
    nmm_amino_lprob_destroy(aminop);
}

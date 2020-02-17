#include "cass.h"
#include "imm/imm.h"
#include "nmm/nmm.h"

void test_amino_table(void);

int main(void)
{
    test_amino_table();
    return cass_status();
}

void test_amino_table(void)
{
    struct imm_abc const*         abc = imm_abc_create("ACDEFGHIKLMNPQRSTVWY", 'X');
    struct nmm_amino_abc const*   amino_abc = nmm_amino_abc_create(abc);
    double const                  zero = imm_lprob_zero();
    double const                  lprobs[NMM_AMINO_ABC_SIZE] = {zero, log(1), [19] = log(19)};
    struct nmm_amino_table const* aminot = nmm_amino_table_create(amino_abc, lprobs);

    cass_cond(imm_lprob_is_zero(nmm_amino_table_lprob(aminot, 'A')));
    cass_close(nmm_amino_table_lprob(aminot, 'C'), log(1));
    cass_close(nmm_amino_table_lprob(aminot, 'Y'), log(19));
    cass_cond(!imm_lprob_is_valid(nmm_amino_table_lprob(aminot, 'B')));
    cass_cond(!imm_lprob_is_valid(nmm_amino_table_lprob(aminot, 'X')));

    imm_abc_destroy(abc);
    nmm_amino_abc_destroy(amino_abc);
    nmm_amino_table_destroy(aminot);
}

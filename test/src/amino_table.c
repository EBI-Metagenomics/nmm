#include "cass.h"
#include "imm/imm.h"
#include "nmm/amino_abc.h"
#include "nmm/nmm.h"

void test_amino_table(void);

int main(void)
{
    test_amino_table();
    return cass_status();
}

void test_amino_table(void)
{
    struct nmm_amino_abc const*   amino_abc = nmm_amino_abc_create("ACDEFGHIKLMNPQRSTVWY", 'X');
    float const                   zero = imm_lprob_zero();
    float const                   lprobs[NMM_AMINO_ABC_SIZE] = {zero, logf(1), [19] = logf(19)};
    struct nmm_amino_table const* aminot = nmm_amino_table_create(amino_abc, lprobs);

    cass_cond(imm_lprob_is_zero(nmm_amino_table_lprob(aminot, 'A')));
    cass_close(nmm_amino_table_lprob(aminot, 'C'), logf(1));
    cass_close(nmm_amino_table_lprob(aminot, 'Y'), logf(19));

    nmm_amino_abc_destroy(amino_abc);
    nmm_amino_table_destroy(aminot);
}

#include "cass.h"
#include "imm/imm.h"
#include "nmm/nmm.h"

void test_base_table(void);

int main(void)
{
    test_base_table();
    return cass_status();
}

void test_base_table(void)
{
    struct imm_abc const*        abc = imm_abc_create("ACGT", 'X');
    struct nmm_base_abc const*   base_abc = nmm_base_abc_create(abc);
    double const                 zero = imm_lprob_zero();
    struct nmm_base_table const* baset =
        nmm_base_table_create(base_abc, log(8), log(1), zero, log(3));

    cass_close(nmm_base_table_lprob(baset, 'A'), log(8));
    cass_close(nmm_base_table_lprob(baset, 'C'), log(1));
    cass_cond(imm_lprob_is_zero(nmm_base_table_lprob(baset, 'G')));
    cass_close(nmm_base_table_lprob(baset, 'T'), log(3));
    cass_cond(!imm_lprob_is_valid(nmm_base_table_lprob(baset, 'X')));

    imm_abc_destroy(abc);
    nmm_base_abc_destroy(base_abc);
    nmm_base_table_destroy(baset);
}

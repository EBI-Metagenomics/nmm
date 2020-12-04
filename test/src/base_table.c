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
    struct nmm_base_abc const*   base_abc = nmm_base_abc_create("ACGT", 'X');
    float const                  zero = imm_lprob_zero();
    struct nmm_base_table const* baset =
        nmm_base_table_create(base_abc, logf(8), logf(1), zero, logf(3));

    cass_close(nmm_base_table_lprob(baset, 'A'), logf(8));
    cass_close(nmm_base_table_lprob(baset, 'C'), logf(1));
    cass_cond(imm_lprob_is_zero(nmm_base_table_lprob(baset, 'G')));
    cass_close(nmm_base_table_lprob(baset, 'T'), logf(3));

    nmm_base_abc_destroy(base_abc);
    nmm_base_table_destroy(baset);
}

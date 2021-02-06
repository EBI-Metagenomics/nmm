#include "cass/cass.h"
#include "helper.h"
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
    struct nmm_base_table const* baset =
        nmm_base_table_create(base_abc, LOG(8), LOG(1), zero(), LOG(3));

    cass_close(nmm_base_table_lprob(baset, 'A'), LOG(8));
    cass_close(nmm_base_table_lprob(baset, 'C'), LOG(1));
    cass_cond(imm_lprob_is_zero(nmm_base_table_lprob(baset, 'G')));
    cass_close(nmm_base_table_lprob(baset, 'T'), LOG(3));

    nmm_base_abc_destroy(base_abc);
    nmm_base_table_destroy(baset);
}

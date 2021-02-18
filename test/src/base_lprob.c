#include "cass/cass.h"
#include "helper.h"
#include "nmm/nmm.h"

void test_base_lprob(void);

int main(void)
{
    test_base_lprob();
    return cass_status();
}

void test_base_lprob(void)
{
    struct nmm_base_abc const*   base_abc = nmm_base_abc_create("ACGT", 'X');
    struct nmm_base_lprob const* basep = nmm_base_lprob_create(base_abc, imm_log(8), imm_log(1), zero(), imm_log(3));

    cass_close(nmm_base_lprob_get(basep, 'A'), imm_log(8));
    cass_close(nmm_base_lprob_get(basep, 'C'), imm_log(1));
    cass_cond(imm_lprob_is_zero(nmm_base_lprob_get(basep, 'G')));
    cass_close(nmm_base_lprob_get(basep, 'T'), imm_log(3));

    nmm_base_abc_destroy(base_abc);
    nmm_base_lprob_destroy(basep);
}

#include "cass.h"
#include "imm/imm.h"
#include "nmm/nmm.h"

void test_base(void);

int main(void)
{
    test_base();
    return cass_status();
}

void test_base(void)
{
    struct imm_abc*  abc = imm_abc_create("ACGT");
    struct nmm_base* base = nmm_base_create(abc);

    cass_cond(nmm_base_set_lprob(base, 'A', log(0.8)) == 0);
    cass_cond(nmm_base_set_lprob(base, 'C', log(0.1)) == 0);
    cass_cond(nmm_base_set_lprob(base, 'T', log(0.3)) == 0);
    cass_cond(nmm_base_set_lprob(base, 'H', log(0.3)) == 1);

    cass_close(nmm_base_get_lprob(base, 'A'), log(0.8));
    cass_close(nmm_base_get_lprob(base, 'C'), log(0.1));
    cass_cond(imm_lprob_is_zero(nmm_base_get_lprob(base, 'G')));
    cass_cond(!imm_lprob_is_valid(nmm_base_get_lprob(base, 'H')));

    nmm_base_normalize(base);

    cass_close(nmm_base_get_lprob(base, 'A'), log(0.8) - log(1.2));
    cass_close(nmm_base_get_lprob(base, 'C'), log(0.1) - log(1.2));
    cass_close(nmm_base_get_lprob(base, 'T'), log(0.3) - log(1.2));
    cass_cond(imm_lprob_is_zero(nmm_base_get_lprob(base, 'G')));
    cass_cond(!imm_lprob_is_valid(nmm_base_get_lprob(base, 'H')));

    imm_abc_destroy(abc);
    nmm_base_destroy(base);
}

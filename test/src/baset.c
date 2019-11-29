#include "cass.h"
#include "imm/imm.h"
#include "nmm/nmm.h"

void test_baset(void);

int main(void)
{
    test_baset();
    return cass_status();
}

void test_baset(void)
{
    struct imm_abc*   abc = imm_abc_create("ACGT");
    struct nmm_baset* baset = nmm_baset_create(abc);

    cass_cond(nmm_baset_set_lprob(baset, 'A', log(0.8)) == 0);
    cass_cond(nmm_baset_set_lprob(baset, 'C', log(0.1)) == 0);
    cass_cond(nmm_baset_set_lprob(baset, 'T', log(0.3)) == 0);
    cass_cond(nmm_baset_set_lprob(baset, 'H', log(0.3)) == 1);

    cass_close(nmm_baset_get_lprob(baset, 'A'), log(0.8));
    cass_close(nmm_baset_get_lprob(baset, 'C'), log(0.1));
    cass_cond(imm_lprob_is_zero(nmm_baset_get_lprob(baset, 'G')));
    cass_cond(!imm_lprob_is_valid(nmm_baset_get_lprob(baset, 'H')));

    nmm_baset_normalize(baset);

    cass_close(nmm_baset_get_lprob(baset, 'A'), log(0.8) - log(1.2));
    cass_close(nmm_baset_get_lprob(baset, 'C'), log(0.1) - log(1.2));
    cass_close(nmm_baset_get_lprob(baset, 'T'), log(0.3) - log(1.2));
    cass_cond(imm_lprob_is_zero(nmm_baset_get_lprob(baset, 'G')));
    cass_cond(!imm_lprob_is_valid(nmm_baset_get_lprob(baset, 'H')));

    imm_abc_destroy(abc);
    nmm_baset_destroy(baset);
}

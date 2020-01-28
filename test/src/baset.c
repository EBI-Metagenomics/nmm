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
    struct imm_abc const*   abc = imm_abc_create("ACGT", 'X');
    struct nmm_base const*  base = nmm_base_create(abc);
    double const            zero = imm_lprob_zero();
    struct nmm_baset const* baset =
        nmm_baset_create(base, log(0.8), log(0.1), zero, log(0.3));

    cass_close(nmm_baset_lprob(baset, 'A'), log(0.8));
    /* cass_close(nmm_baset_lprob(baset, 'C'), log(0.1)); */
    /* cass_cond(imm_lprob_is_zero(nmm_baset_lprob(baset, 'G'))); */
    /* cass_close(nmm_baset_lprob(baset, 'T'), log(0.3)); */
    /* cass_cond(!imm_lprob_is_valid(nmm_baset_lprob(baset, 'X'))); */

    imm_abc_destroy(abc);
    nmm_baset_destroy(baset);
}

#include "cass.h"
#include "imm/imm.h"
#include "nmm/nmm.h"

void test_codonp(void);

int main(void)
{
    test_codonp();
    return cass_status();
}

void test_codonp(void)
{
    struct imm_abc const*  abc = imm_abc_create("ACGT", 'X');
    struct nmm_base const* base = nmm_base_create(abc);
    struct nmm_codonp*     codonp = nmm_codonp_create(base);
    cass_cond(codonp != NULL);

    cass_cond(nmm_codonp_normalize(codonp) == 1);

    struct nmm_codon* codon = nmm_codon_create(base);

    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'C', 'C')) == 0);

    cass_cond(imm_lprob_is_zero(nmm_codonp_get_lprob(codonp, codon)));
    cass_cond(nmm_codonp_set_lprob(codonp, codon, log(0.5)) == 0);
    cass_close(nmm_codonp_get_lprob(codonp, codon), log(0.5));

    cass_cond(nmm_codonp_normalize(codonp) == 0);
    cass_close(nmm_codonp_get_lprob(codonp, codon), log(1.0));

    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'C', 'X')) == 0);
    cass_cond(!imm_lprob_is_valid(nmm_codonp_get_lprob(codonp, codon)));

    nmm_codon_destroy(codon);

    nmm_codonp_destroy(codonp);
    nmm_base_destroy(base);
    imm_abc_destroy(abc);
}

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
    struct nmm_base_abc const* base = nmm_base_abc_create("ACGT", 'X');
    struct nmm_codon_lprob*    codonp = nmm_codon_lprob_create(base);
    cass_cond(codonp != NULL);

    cass_cond(nmm_codon_lprob_normalize(codonp) == 1);

    struct nmm_codon* codon = nmm_codon_create(base);

    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'C', 'C')) == 0);

    cass_cond(imm_lprob_is_zero(nmm_codon_lprob_get(codonp, codon)));
    cass_cond(nmm_codon_lprob_set(codonp, codon, log(0.5)) == 0);
    cass_close(nmm_codon_lprob_get(codonp, codon), log(0.5));

    cass_cond(nmm_codon_lprob_normalize(codonp) == 0);
    cass_close(nmm_codon_lprob_get(codonp, codon), log(1.0));

    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'C', 'X')) == 0);
    cass_cond(!imm_lprob_is_valid(nmm_codon_lprob_get(codonp, codon)));

    nmm_codon_destroy(codon);

    nmm_codon_lprob_destroy(codonp);
    nmm_base_abc_destroy(base);
}

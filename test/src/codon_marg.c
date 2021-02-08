#include "cass/cass.h"
#include "imm/imm.h"
#include "nmm/nmm.h"

void test_codonm_nonmarginal(void);
void test_codonm_marginal(void);

int main(void)
{
    test_codonm_nonmarginal();
    test_codonm_marginal();
    return cass_status();
}

void test_codonm_nonmarginal(void)
{
    struct nmm_base_abc const* base = nmm_base_abc_create("ACGT", 'X');
    struct nmm_codon_lprob*    codonp = nmm_codon_lprob_create(base);

    struct nmm_codon* codon = nmm_codon_create(base);
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'T', 'G')) == 0);
    cass_cond(nmm_codon_lprob_set(codonp, codon, imm_log(0.8)) == 0);
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'T', 'T')) == 0);
    cass_cond(nmm_codon_lprob_set(codonp, codon, imm_log(0.1)) == 0);

    struct nmm_codon_marg const* codonm = nmm_codon_marg_create(codonp);
    cass_cond(codonm != NULL);
    nmm_codon_lprob_destroy(codonp);

    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'T', 'G')) == 0);
    cass_close(nmm_codon_marg_lprob(codonm, codon), imm_log(0.8));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'T', 'T')) == 0);
    cass_close(nmm_codon_marg_lprob(codonm, codon), imm_log(0.1));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('T', 'T', 'T')) == 0);
    cass_cond(imm_lprob_is_zero(nmm_codon_marg_lprob(codonm, codon)));
    nmm_codon_destroy(codon);

    nmm_codon_marg_destroy(codonm);
    nmm_base_abc_destroy(base);
}

void test_codonm_marginal(void)
{
    struct nmm_base_abc const* base = nmm_base_abc_create("ACGT", 'X');
    struct nmm_codon_lprob*    codonp = nmm_codon_lprob_create(base);

    struct nmm_codon* codon = nmm_codon_create(base);
    nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'T', 'G'));
    nmm_codon_lprob_set(codonp, codon, imm_log(0.8));
    nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'T', 'T'));
    nmm_codon_lprob_set(codonp, codon, imm_log(0.1));

    struct nmm_codon_marg const* codonm = nmm_codon_marg_create(codonp);
    cass_cond(codonm != NULL);
    nmm_codon_lprob_destroy(codonp);

    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'T', 'G')) == 0);
    cass_close(nmm_codon_marg_lprob(codonm, codon), imm_log(0.8));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'T', 'T')) == 0);
    cass_close(nmm_codon_marg_lprob(codonm, codon), imm_log(0.1));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'T', 'X')) == 0);
    cass_close(nmm_codon_marg_lprob(codonm, codon), imm_log(0.9));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'X', 'X')) == 0);
    cass_close(nmm_codon_marg_lprob(codonm, codon), imm_log(0.9));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('X', 'X', 'X')) == 0);
    cass_close(nmm_codon_marg_lprob(codonm, codon), imm_log(0.9));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('X', 'T', 'X')) == 0);
    cass_close(nmm_codon_marg_lprob(codonm, codon), imm_log(0.9));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('X', 'X', 'G')) == 0);
    cass_close(nmm_codon_marg_lprob(codonm, codon), imm_log(0.8));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('X', 'X', 'T')) == 0);
    cass_close(nmm_codon_marg_lprob(codonm, codon), imm_log(0.1));
    nmm_codon_destroy(codon);

    nmm_codon_marg_destroy(codonm);
    nmm_base_abc_destroy(base);
}

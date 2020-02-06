#include "cass.h"
#include "imm/imm.h"
#include "nmm/nmm.h"

void test_codont_nonmarginal(void);
void test_codont_marginal(void);

int main(void)
{
    test_codont_nonmarginal();
    test_codont_marginal();
    return cass_status();
}

void test_codont_nonmarginal(void)
{
    struct imm_abc const*  abc = imm_abc_create("ACGT", 'X');
    struct nmm_base const* base = nmm_base_create(abc);
    struct nmm_codonp*     codonp = nmm_codonp_create(base);

    struct nmm_codon* codon = nmm_codon_create(base);
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'T', 'G')) == 0);
    cass_cond(nmm_codonp_set_lprob(codonp, codon, log(0.8)) == 0);
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'T', 'T')) == 0);
    cass_cond(nmm_codonp_set_lprob(codonp, codon, log(0.1)) == 0);

    struct nmm_codont const* codont = nmm_codont_create(codonp);
    cass_cond(codont != NULL);
    nmm_codonp_destroy(codonp);

    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'T', 'G')) == 0);
    cass_close(nmm_codont_lprob(codont, codon), log(0.8));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'T', 'T')) == 0);
    cass_close(nmm_codont_lprob(codont, codon), log(0.1));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('T', 'T', 'T')) == 0);
    cass_cond(imm_lprob_is_zero(nmm_codont_lprob(codont, codon)));
    nmm_codon_destroy(codon);

    nmm_codont_destroy(codont);
    nmm_base_destroy(base);
    imm_abc_destroy(abc);
}

void test_codont_marginal(void)
{
    struct imm_abc const*  abc = imm_abc_create("ACGT", 'X');
    struct nmm_base const* base = nmm_base_create(abc);
    struct nmm_codonp*     codonp = nmm_codonp_create(base);

    struct nmm_codon* codon = nmm_codon_create(base);
    nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'T', 'G'));
    nmm_codonp_set_lprob(codonp, codon, log(0.8));
    nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'T', 'T'));
    nmm_codonp_set_lprob(codonp, codon, log(0.1));

    struct nmm_codont const* codont = nmm_codont_create(codonp);
    cass_cond(codont != NULL);
    nmm_codonp_destroy(codonp);

    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'T', 'G')) == 0);
    cass_close(nmm_codont_lprob(codont, codon), log(0.8));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'T', 'T')) == 0);
    cass_close(nmm_codont_lprob(codont, codon), log(0.1));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'T', 'X')) == 0);
    cass_close(nmm_codont_lprob(codont, codon), log(0.9));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'X', 'X')) == 0);
    cass_close(nmm_codont_lprob(codont, codon), log(0.9));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('X', 'X', 'X')) == 0);
    cass_close(nmm_codont_lprob(codont, codon), log(0.9));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('X', 'T', 'X')) == 0);
    cass_close(nmm_codont_lprob(codont, codon), log(0.9));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('X', 'X', 'G')) == 0);
    cass_close(nmm_codont_lprob(codont, codon), log(0.8));
    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('X', 'X', 'T')) == 0);
    cass_close(nmm_codont_lprob(codont, codon), log(0.1));
    nmm_codon_destroy(codon);

    nmm_codont_destroy(codont);
    nmm_base_destroy(base);
    imm_abc_destroy(abc);
}

#include "cass.h"
#include "imm/imm.h"
#include "nmm/nmm.h"

void test_codonp_wrong_alphabet_length(void);
void test_codonp_set_invalid_codon(void);
void test_codonp_get_invalid_codon(void);
void test_codonp_normalize(void);

int main(void)
{
    test_codonp_wrong_alphabet_length();
    test_codonp_set_invalid_codon();
    test_codonp_get_invalid_codon();
    test_codonp_normalize();
    return cass_status();
}

void test_codonp_wrong_alphabet_length(void)
{
    struct imm_abc const* abc = imm_abc_create("ACGTE", 'X');
    struct nmm_codonp*    codonp = nmm_codonp_create(abc);
    cass_cond(codonp == NULL);
    imm_abc_destroy(abc);
}

void test_codonp_set_invalid_codon(void)
{
    struct imm_abc const* abc = imm_abc_create("ACGT", 'X');
    struct nmm_codonp*    codonp = nmm_codonp_create(abc);
    cass_cond(codonp != NULL);
    cass_cond(nmm_codonp_set(codonp, &NMM_CODON('A', 'C', 'U'), log(0.5)) == 1);
    imm_abc_destroy(abc);
    nmm_codonp_destroy(codonp);
}

void test_codonp_get_invalid_codon(void)
{
    struct imm_abc const* abc = imm_abc_create("ACGT", 'X');
    struct nmm_codonp*    codonp = nmm_codonp_create(abc);
    cass_cond(codonp != NULL);
    cass_cond(nmm_codonp_set(codonp, &NMM_CODON('A', 'C', 'C'), log(0.5)) == 0);
    cass_cond(!imm_lprob_is_valid(nmm_codonp_get(codonp, &NMM_CODON('A', 'C', 'U'))));
    imm_abc_destroy(abc);
    nmm_codonp_destroy(codonp);
}

void test_codonp_normalize(void)
{
    struct imm_abc const* abc = imm_abc_create("ACGT", 'X');
    struct nmm_codonp*    codonp = nmm_codonp_create(abc);
    cass_cond(codonp != NULL);
    cass_cond(nmm_codonp_normalize(codonp) == 1);
    cass_cond(nmm_codonp_set(codonp, &NMM_CODON('A', 'C', 'C'), log(0.5)) == 0);
    cass_cond(nmm_codonp_normalize(codonp) == 0);
    cass_close(nmm_codonp_get(codonp, &NMM_CODON('A', 'C', 'C')), log(1.0));
    imm_abc_destroy(abc);
    nmm_codonp_destroy(codonp);
}

#include "cass.h"
#include "imm/imm.h"
#include "nmm/nmm.h"

void test_amino_abc_success1(void);
void test_amino_abc_success2(void);
void test_amino_abc_failure(void);

int main(void)
{
    test_amino_abc_success1();
    test_amino_abc_success2();
    test_amino_abc_failure();
    return cass_status();
}

void test_amino_abc_success1(void)
{
    struct nmm_amino_abc const* amino_abc = nmm_amino_abc_create("ACDEFGHIKLMNPQRSTVWY", 'X');
    cass_cond(amino_abc != NULL);
    nmm_amino_abc_destroy(amino_abc);
}

void test_amino_abc_success2(void)
{
    struct nmm_amino_abc const* amino_abc = nmm_amino_abc_create("ACDEFGHIKLMNPQRSTVWY", 'X');
    cass_cond(amino_abc != NULL);
    struct imm_abc const* abc = nmm_amino_abc_super(amino_abc);
    cass_cond(abc != NULL);
    imm_abc_destroy(abc);
}

void test_amino_abc_failure(void)
{
    struct nmm_amino_abc const* amino_abc = nmm_amino_abc_create("ACT", 'X');
    cass_cond(amino_abc == NULL);
}

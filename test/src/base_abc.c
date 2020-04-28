#include "cass.h"
#include "imm/imm.h"
#include "nmm/nmm.h"

void test_base_abc_success1(void);
void test_base_abc_success2(void);
void test_base_abc_failure(void);

int main(void)
{
    test_base_abc_success1();
    test_base_abc_success2();
    test_base_abc_failure();
    return cass_status();
}

void test_base_abc_success1(void)
{
    struct nmm_base_abc const* base_abc = nmm_base_abc_create("ACGT", 'X');
    cass_cond(base_abc != NULL);
    nmm_base_abc_destroy(base_abc);
}

void test_base_abc_success2(void)
{
    struct nmm_base_abc const* base_abc = nmm_base_abc_create("ACGT", 'X');
    cass_cond(base_abc != NULL);
    struct imm_abc const* abc = nmm_base_abc_super(base_abc);
    cass_cond(abc != NULL);
    imm_abc_destroy(abc);
}

void test_base_abc_failure(void)
{
    struct nmm_base_abc const* base_abc = nmm_base_abc_create("ACT", 'X');
    cass_cond(base_abc == NULL);
}

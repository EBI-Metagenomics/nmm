#include "cass.h"
#include "imm/imm.h"
#include "nmm/nmm.h"

void test_base_abc_success(void);
void test_base_abc_failure(void);

int main(void)
{
    test_base_abc_success();
    test_base_abc_failure();
    return cass_status();
}

void test_base_abc_success(void)
{
    struct imm_abc const*      abc = imm_abc_create("ACGT", 'X');
    struct nmm_base_abc const* base = nmm_base_abc_create(abc);
    cass_cond(base != NULL);
    nmm_base_abc_destroy(base);
    imm_abc_destroy(abc);
}

void test_base_abc_failure(void)
{
    struct imm_abc const*      abc = imm_abc_create("ACT", 'X');
    struct nmm_base_abc const* base = nmm_base_abc_create(abc);
    cass_cond(base == NULL);
    imm_abc_destroy(abc);
}

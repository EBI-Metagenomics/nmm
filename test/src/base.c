#include "cass.h"
#include "imm/imm.h"
#include "nmm/nmm.h"

void test_base_success(void);
void test_base_failure(void);

int main(void)
{
    test_base_success();
    test_base_failure();
    return cass_status();
}

void test_base_success(void)
{
    struct imm_abc const*   abc = imm_abc_create("ACGT", 'X');
    struct nmm_base const* base = nmm_base_create(abc);
    cass_cond(base != NULL);
    nmm_base_destroy(base);
    imm_abc_destroy(abc);
}

void test_base_failure(void)
{
    struct imm_abc const*   abc = imm_abc_create("ACT", 'X');
    struct nmm_base const* base = nmm_base_create(abc);
    cass_cond(base == NULL);
    imm_abc_destroy(abc);
}

#include "cass.h"
#include "imm/imm.h"
#include "nmm/nmm.h"

void test_amino_abc_success(void);
void test_amino_abc_failure(void);

int main(void)
{
    test_amino_abc_success();
    test_amino_abc_failure();
    return cass_status();
}

void test_amino_abc_success(void)
{
    struct imm_abc const*       abc = imm_abc_create("ACDEFGHIKLMNPQRSTVWY", 'X');
    struct nmm_amino_abc const* amino_abc = nmm_amino_abc_create(abc);
    cass_cond(amino_abc != NULL);
    nmm_amino_abc_destroy(amino_abc);
    imm_abc_destroy(abc);
}

void test_amino_abc_failure(void)
{
    struct imm_abc const*       abc = imm_abc_create("ACT", 'X');
    struct nmm_amino_abc const* amino_abc = nmm_amino_abc_create(abc);
    cass_cond(amino_abc == NULL);
    imm_abc_destroy(abc);
}

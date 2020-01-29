#include "cass.h"
#include "imm/imm.h"
#include "nmm/nmm.h"

void test_codon(void);

int main(void)
{
    test_codon();
    return cass_status();
}

void test_codon(void)
{
    struct imm_abc const*  abc = imm_abc_create("ACGT", 'X');
    struct nmm_base const* base = nmm_base_create(abc);
    char a, b, c;

    struct nmm_codon* codon = nmm_codon_create(base);

    cass_cond(nmm_codon_set(codon, 'A', 'T', 'G') == 0);
    nmm_codon_get(codon, &a, &b, &c);
    cass_cond(a == 'A' && b == 'T' && c == 'G');

    cass_cond(nmm_codon_set(codon, 'X', 'T', 'X') == 0);
    nmm_codon_get(codon, &a, &b, &c);
    cass_cond(a == 'X' && b == 'T' && c == 'X');

    cass_cond(nmm_codon_set(codon, 'A', 'T', 'L') == 1);

    nmm_codon_destroy(codon);
    nmm_base_destroy(base);
    imm_abc_destroy(abc);
}

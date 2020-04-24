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
    struct nmm_base_abc const* base = nmm_base_abc_create("ACGT", 'X');

    struct nmm_codon* codon = nmm_codon_create(base);

    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'T', 'G')) == 0);
    struct nmm_triplet t = nmm_codon_get_triplet(codon);
    cass_cond(t.a == 'A' && t.b == 'T' && t.c == 'G');

    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('X', 'T', 'X')) == 0);
    t = nmm_codon_get_triplet(codon);
    cass_cond(t.a == 'X' && t.b == 'T' && t.c == 'X');

    cass_cond(nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'T', 'L')) == 1);

    nmm_codon_destroy(codon);
    nmm_base_abc_destroy(base);
}

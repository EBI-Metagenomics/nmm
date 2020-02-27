#include "cass.h"
#include "imm/imm.h"
#include "nmm/nmm.h"

void test_codon_state1(void);

int main(void)
{
    test_codon_state1();
    return cass_status();
}

void test_codon_state1(void)
{
    struct imm_abc const*        abc = imm_abc_create("ACGT", 'X');
    struct imm_abc const*        another_abc = imm_abc_create("AUT", 'X');
    struct nmm_base_abc const*   base = nmm_base_abc_create(abc);

    struct nmm_codon_lprob* codonp = nmm_codon_lprob_create(base);
    struct nmm_codon*       codon = nmm_codon_create(base);

    nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'T', 'G'));
    nmm_codon_lprob_set(codonp, codon, log(0.8 / 0.9));
    nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'T', 'T'));
    nmm_codon_lprob_set(codonp, codon, log(0.1 / 0.9));

    struct nmm_codon_state const* state = nmm_codon_state_create("State", codonp);
    const struct imm_state*       s = imm_state_cast_c(state);

    struct imm_seq const* seq = imm_seq_create("ATG", abc);
    cass_close(imm_state_lprob(s, seq), log(0.8 / 0.9));
    imm_seq_destroy(seq);

    seq = imm_seq_create("AG", abc);
    cass_cond(!imm_lprob_is_valid(imm_state_lprob(s, seq)));
    imm_seq_destroy(seq);

    seq = imm_seq_create("UUU", another_abc);
    cass_cond(!imm_lprob_is_valid(imm_state_lprob(s, seq)));
    imm_seq_destroy(seq);

    nmm_codon_destroy(codon);
    nmm_codon_state_destroy(state);
    nmm_codon_lprob_destroy(codonp);
    nmm_base_abc_destroy(base);
    imm_abc_destroy(abc);
    imm_abc_destroy(another_abc);
}

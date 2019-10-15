#include "cass/cass.h"
#include "imm.h"
#include "nmm.h"

void test_hmm1(void);

int main(void)
{
    test_hmm1();
    return cass_status();
}

void test_hmm1(void)
{
    struct imm_abc *bases = imm_abc_create("ACGT");
    struct imm_hmm *hmm = imm_hmm_create(bases);
    const int A = imm_abc_symbol_idx(bases, 'A');
    const int C = imm_abc_symbol_idx(bases, 'C');
    const int G = imm_abc_symbol_idx(bases, 'G');
    const int T = imm_abc_symbol_idx(bases, 'T');

    struct nmm_codon *codon = nmm_codon_create();
    nmm_codon_set_ninfs(codon);

    nmm_codon_set_lprob(codon, A, T, G, log(0.8));
    nmm_codon_set_lprob(codon, A, T, T, log(0.1));
    nmm_codon_set_lprob(codon, C, C, C, log(0.1));

    double base_lprobs[] = {log(0.25), log(0.25), log(0.5), LOG0};
    struct nmm_frame_state *state =
        nmm_frame_state_create("M", bases, base_lprobs, codon, 0.0);

    imm_hmm_add_state(hmm, imm_state_cast_c(state), log(1.0));

    struct imm_path *path = imm_path_create();
    imm_path_add(path, imm_state_cast_c(state), 3);
    cass_close(imm_hmm_likelihood(hmm, "ATT", path), 0);
    imm_path_destroy(path);

    imm_hmm_destroy(hmm);
    nmm_frame_state_destroy(state);
    nmm_codon_destroy(codon);
    imm_abc_destroy(bases);
}

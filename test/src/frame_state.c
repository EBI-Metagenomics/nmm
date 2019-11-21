#include "cass.h"
#include "imm/imm.h"
#include "nmm/nmm.h"

void test_frame_state1(void);
void test_frame_state2(void);
void test_frame_state3(void);
void test_frame_state_posterior(void);

int main(void)
{
    test_frame_state1();
    test_frame_state2();
    test_frame_state3();
    test_frame_state_posterior();
    return cass_status();
}

void test_frame_state1(void)
{
    struct imm_abc* abc = imm_abc_create("ACGT");

    struct nmm_base* base = nmm_base_create(abc);
    nmm_base_set_lprob(base, 'A', log(0.2));
    nmm_base_set_lprob(base, 'C', log(0.2));
    nmm_base_set_lprob(base, 'G', log(0.2));
    nmm_base_set_lprob(base, 'T', log(0.2));
    cass_cond(nmm_base_normalize(base) == 0);

    struct nmm_codon* codon = nmm_codon_create(abc);
    nmm_codon_set_lprob(codon, &NMM_CCODE('A', 'T', 'G'), log(0.8 / 0.9));
    nmm_codon_set_lprob(codon, &NMM_CCODE('A', 'T', 'T'), log(0.1 / 0.9));

    struct nmm_frame_state* state = nmm_frame_state_create("State", base, codon, 0.1);

    const struct imm_state* s = imm_state_cast_c(state);
    cass_close(imm_state_lprob(s, "A", 1), -5.914503505971854);
    cass_close(imm_state_lprob(s, "AT", 2), -2.915843423869834);
    cass_close(imm_state_lprob(s, "ATA", 3), -6.905597115665666);
    cass_close(imm_state_lprob(s, "ATG", 3), -0.534773288204706);
    cass_close(imm_state_lprob(s, "ATT", 3), -2.590237330499946);
    cass_close(imm_state_lprob(s, "ATTA", 4), -6.881032208841384);
    cass_close(imm_state_lprob(s, "ATTAA", 5), -12.08828960987379);
    cass_cond(imm_lprob_is_zero(imm_state_lprob(s, "ATTAAT", 6)));

    nmm_frame_state_destroy(state);
    nmm_base_destroy(base);
    nmm_codon_destroy(codon);
    imm_abc_destroy(abc);
}

void test_frame_state2(void)
{
    struct imm_abc* abc = imm_abc_create("ACGT");

    struct nmm_base* base = nmm_base_create(abc);
    nmm_base_set_lprob(base, 'A', log(0.1));
    nmm_base_set_lprob(base, 'C', log(0.2));
    nmm_base_set_lprob(base, 'G', log(0.3));
    nmm_base_set_lprob(base, 'T', log(0.4));
    cass_cond(nmm_base_normalize(base) == 0);

    struct nmm_codon* codon = nmm_codon_create(abc);

    nmm_codon_set_lprob(codon, &NMM_CCODE('A', 'T', 'G'), log(0.8 / 0.9));
    nmm_codon_set_lprob(codon, &NMM_CCODE('A', 'T', 'T'), log(0.1 / 0.9));

    struct nmm_frame_state* state = nmm_frame_state_create("State", base, codon, 0.1);

    const struct imm_state* s = imm_state_cast_c(state);
    cass_close(imm_state_lprob(s, "A", 1), -5.914503505971854);
    cass_cond(imm_lprob_is_zero(imm_state_lprob(s, "C", 1)));
    cass_close(imm_state_lprob(s, "G", 1), -6.032286541628237);
    cass_close(imm_state_lprob(s, "T", 1), -5.809142990314027);

    cass_close(imm_state_lprob(s, "AT", 2), -2.9159357500274385);
    cass_close(imm_state_lprob(s, "ATA", 3), -7.821518343902165);
    cass_close(imm_state_lprob(s, "ATG", 3), -0.5344319079005616);
    cass_close(imm_state_lprob(s, "ATC", 3), -7.129480084106424);
    cass_close(imm_state_lprob(s, "ATT", 3), -2.57514520832882);

    cass_close(imm_state_lprob(s, "ATTA", 4), -7.789644584138959);
    cass_close(imm_state_lprob(s, "ACTG", 4), -5.036637096635257);

    cass_close(imm_state_lprob(s, "ATTAA", 5), -13.920871073622099);

    cass_cond(imm_lprob_is_zero(imm_state_lprob(s, "ATTAAT", 6)));

    nmm_frame_state_destroy(state);
    nmm_base_destroy(base);
    nmm_codon_destroy(codon);
    imm_abc_destroy(abc);
}

void test_frame_state3(void)
{
    struct imm_abc* abc = imm_abc_create("ACGT");

    struct nmm_base* base = nmm_base_create(abc);
    nmm_base_set_lprob(base, 'A', log(0.1));
    nmm_base_set_lprob(base, 'C', log(0.2));
    nmm_base_set_lprob(base, 'G', log(0.3));
    nmm_base_set_lprob(base, 'T', log(0.4));
    cass_cond(nmm_base_normalize(base) == 0);

    struct nmm_codon* codon = nmm_codon_create(abc);
    nmm_codon_set_lprob(codon, &NMM_CCODE('A', 'T', 'G'), log(0.8));
    nmm_codon_set_lprob(codon, &NMM_CCODE('A', 'T', 'T'), log(0.1));
    nmm_codon_set_lprob(codon, &NMM_CCODE('G', 'T', 'C'), log(0.4));
    cass_cond(nmm_codon_normalize(codon) == 0);

    struct nmm_frame_state* state = nmm_frame_state_create("State", base, codon, 0.1);

    const struct imm_state* s = imm_state_cast_c(state);
    cass_close(imm_state_lprob(s, "A", 1), -6.282228286097171);
    cass_close(imm_state_lprob(s, "C", 1), -7.0931585023135);
    cass_close(imm_state_lprob(s, "G", 1), -5.99454621364539);
    cass_close(imm_state_lprob(s, "T", 1), -5.840395533818132);
    cass_close(imm_state_lprob(s, "AT", 2), -3.283414346005771);
    cass_close(imm_state_lprob(s, "CG", 2), -9.395743595307545);
    cass_close(imm_state_lprob(s, "ATA", 3), -8.18911998648269);
    cass_close(imm_state_lprob(s, "ATG", 3), -0.9021560981322401);
    cass_close(imm_state_lprob(s, "ATT", 3), -2.9428648000333952);
    cass_close(imm_state_lprob(s, "ATC", 3), -7.314811395663229);
    cass_close(imm_state_lprob(s, "GTC", 3), -1.5951613351178675);
    cass_close(imm_state_lprob(s, "ATTA", 4), -8.157369364264277);
    cass_close(imm_state_lprob(s, "GTTC", 4), -4.711642430498609);
    cass_close(imm_state_lprob(s, "ACTG", 4), -5.404361876760574);
    cass_close(imm_state_lprob(s, "ATTAA", 5), -14.288595853747417);
    cass_close(imm_state_lprob(s, "GTCAA", 5), -12.902301492627526);

    nmm_frame_state_destroy(state);
    nmm_base_destroy(base);
    nmm_codon_destroy(codon);
    imm_abc_destroy(abc);
}

void test_frame_state_posterior(void)
{

    struct imm_abc* abc = imm_abc_create("ACGT");

    struct nmm_base* base = nmm_base_create(abc);
    nmm_base_set_lprob(base, 'A', log(0.1));
    nmm_base_set_lprob(base, 'C', log(0.2));
    nmm_base_set_lprob(base, 'G', log(0.3));
    nmm_base_set_lprob(base, 'T', log(0.4));
    cass_cond(nmm_base_normalize(base) == 0);

    struct nmm_codon* codon = nmm_codon_create(abc);
    nmm_codon_set_lprob(codon, &NMM_CCODE('A', 'T', 'G'), log(0.8));
    nmm_codon_set_lprob(codon, &NMM_CCODE('A', 'T', 'T'), log(0.1));
    nmm_codon_set_lprob(codon, &NMM_CCODE('G', 'T', 'C'), log(0.4));
    cass_cond(nmm_codon_normalize(codon) == 0);

    struct nmm_frame_state* state = nmm_frame_state_create("State", base, codon, 0.1);

    /* printf("%c: %16f\n", seq[0], nmm_frame_state_posterior(state, seq, 1, &ccode)); */
    struct nmm_ccode ccode = {'A', 'G', 'A'};
    char             seq[5] = "01234";
    double           total = imm_lprob_zero();
    double           total_f[5] = {
        imm_lprob_zero(), imm_lprob_zero(), imm_lprob_zero(),
        imm_lprob_zero(), imm_lprob_zero(),
    };
    for (int i0 = 0; i0 < 4; ++i0) {
        seq[0] = imm_abc_symbol_id(abc, i0);
        total = imm_lprob_add(total, nmm_frame_state_posterior(state, seq, 1, &ccode));
        total_f[0] =
            imm_lprob_add(total_f[0], nmm_frame_state_posterior(state, seq, 1, &ccode));
    }
    printf("Total F=1: %.16f\n", exp(total_f[0]));

    for (int i0 = 0; i0 < 4; ++i0) {
        for (int i1 = 0; i1 < 4; ++i1) {
            seq[0] = imm_abc_symbol_id(abc, i0);
            seq[1] = imm_abc_symbol_id(abc, i1);
            total = imm_lprob_add(total, nmm_frame_state_posterior(state, seq, 2, &ccode));
            total_f[1] =
                imm_lprob_add(total_f[1], nmm_frame_state_posterior(state, seq, 2, &ccode));
        }
    }
    printf("Total F=2: %.16f\n", exp(total_f[1]));

    for (int i0 = 0; i0 < 4; ++i0) {
        for (int i1 = 0; i1 < 4; ++i1) {
            for (int i2 = 0; i2 < 4; ++i2) {
                seq[0] = imm_abc_symbol_id(abc, i0);
                seq[1] = imm_abc_symbol_id(abc, i1);
                seq[2] = imm_abc_symbol_id(abc, i2);
                total =
                    imm_lprob_add(total, nmm_frame_state_posterior(state, seq, 3, &ccode));
                total_f[2] = imm_lprob_add(total_f[2],
                                           nmm_frame_state_posterior(state, seq, 3, &ccode));
            }
        }
    }
    printf("Total F=3: %.16f\n", exp(total_f[2]));

    for (int i0 = 0; i0 < 4; ++i0) {
        for (int i1 = 0; i1 < 4; ++i1) {
            for (int i2 = 0; i2 < 4; ++i2) {
                for (int i3 = 0; i3 < 4; ++i3) {
                    seq[0] = imm_abc_symbol_id(abc, i0);
                    seq[1] = imm_abc_symbol_id(abc, i1);
                    seq[2] = imm_abc_symbol_id(abc, i2);
                    seq[3] = imm_abc_symbol_id(abc, i3);
                    total = imm_lprob_add(total,
                                          nmm_frame_state_posterior(state, seq, 4, &ccode));
                    total_f[3] = imm_lprob_add(
                        total_f[3], nmm_frame_state_posterior(state, seq, 4, &ccode));
                }
            }
        }
    }
    printf("Total F=4: %.16f\n", exp(total_f[3]));

    for (int i0 = 0; i0 < 4; ++i0) {
        for (int i1 = 0; i1 < 4; ++i1) {
            for (int i2 = 0; i2 < 4; ++i2) {
                for (int i3 = 0; i3 < 4; ++i3) {
                    for (int i4 = 0; i4 < 4; ++i4) {
                        seq[0] = imm_abc_symbol_id(abc, i0);
                        seq[1] = imm_abc_symbol_id(abc, i1);
                        seq[2] = imm_abc_symbol_id(abc, i2);
                        seq[3] = imm_abc_symbol_id(abc, i3);
                        seq[4] = imm_abc_symbol_id(abc, i4);
                        total = imm_lprob_add(
                            total, nmm_frame_state_posterior(state, seq, 5, &ccode));
                        total_f[4] = imm_lprob_add(
                            total_f[4], nmm_frame_state_posterior(state, seq, 5, &ccode));
                    }
                }
            }
        }
    }
    printf("Total F=5: %.16f\n", exp(total_f[4]));

    printf("Total: %.16f\n", exp(total));

    nmm_frame_state_destroy(state);
    nmm_base_destroy(base);
    nmm_codon_destroy(codon);
    imm_abc_destroy(abc);
}

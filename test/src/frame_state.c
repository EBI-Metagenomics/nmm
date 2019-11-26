#include "cartes.h"
#include "cass.h"
#include "imm/imm.h"
#include "nmm/nmm.h"

void test_frame_state1(void);
void test_frame_state2(void);
void test_frame_state3(void);
void test_frame_state_lposterior(void);
void test_frame_state_decode(void);

int main(void)
{
    test_frame_state1();
    test_frame_state2();
    test_frame_state3();
    test_frame_state_lposterior();
    test_frame_state_decode();
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

    struct nmm_codont* codon = nmm_codont_create(abc);
    nmm_codont_set_lprob(codon, &NMM_CCODE('A', 'T', 'G'), log(0.8 / 0.9));
    nmm_codont_set_lprob(codon, &NMM_CCODE('A', 'T', 'T'), log(0.1 / 0.9));

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
    nmm_codont_destroy(codon);
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

    struct nmm_codont* codon = nmm_codont_create(abc);

    nmm_codont_set_lprob(codon, &NMM_CCODE('A', 'T', 'G'), log(0.8 / 0.9));
    nmm_codont_set_lprob(codon, &NMM_CCODE('A', 'T', 'T'), log(0.1 / 0.9));

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
    nmm_codont_destroy(codon);
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

    struct nmm_codont* codon = nmm_codont_create(abc);
    nmm_codont_set_lprob(codon, &NMM_CCODE('A', 'T', 'G'), log(0.8));
    nmm_codont_set_lprob(codon, &NMM_CCODE('A', 'T', 'T'), log(0.1));
    nmm_codont_set_lprob(codon, &NMM_CCODE('G', 'T', 'C'), log(0.4));
    cass_cond(nmm_codont_normalize(codon) == 0);

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
    nmm_codont_destroy(codon);
    imm_abc_destroy(abc);
}

void test_frame_state_lposterior(void)
{
    struct imm_abc* abc = imm_abc_create("ACGT");

    struct nmm_base* base = nmm_base_create(abc);
    nmm_base_set_lprob(base, 'A', log(0.1));
    nmm_base_set_lprob(base, 'C', log(0.2));
    nmm_base_set_lprob(base, 'G', log(0.3));
    nmm_base_set_lprob(base, 'T', log(0.4));
    cass_cond(nmm_base_normalize(base) == 0);

    char const* symbols = imm_abc_symbols(abc);
    int         length = imm_abc_length(abc);

    struct nmm_codont* codon = nmm_codont_create(abc);
    struct cartes*     ccode_iter = cartes_create(symbols, length, 3);
    char const*        ccode_item = NULL;

    while ((ccode_item = cartes_next(ccode_iter)) != NULL) {
        struct nmm_ccode ccode = {ccode_item[0], ccode_item[1], ccode_item[2]};
        nmm_codont_set_lprob(codon, &ccode, log(0.001));
    }
    cartes_destroy(ccode_iter);

    nmm_codont_set_lprob(codon, &NMM_CCODE('A', 'T', 'G'), log(0.8));
    nmm_codont_set_lprob(codon, &NMM_CCODE('A', 'T', 'T'), log(0.1));
    nmm_codont_set_lprob(codon, &NMM_CCODE('G', 'T', 'C'), log(0.4));
    cass_cond(nmm_codont_normalize(codon) == 0);

    struct nmm_frame_state* state = nmm_frame_state_create("State", base, codon, 0.1);

    ccode_iter = cartes_create(symbols, length, 3);

    while ((ccode_item = cartes_next(ccode_iter)) != NULL) {
        struct nmm_ccode ccode = {ccode_item[0], ccode_item[1], ccode_item[2]};

        double total = imm_lprob_zero();
        for (int times = 1; times < 6; ++times) {

            struct cartes* seq_iter = cartes_create(symbols, length, times);
            char const*    seq = NULL;

            while ((seq = cartes_next(seq_iter)) != NULL) {
                double lprob = nmm_frame_state_lposterior(state, &ccode, seq, times);
                lprob -= nmm_codont_get_lprob(codon, &ccode);
                total = imm_lprob_add(total, lprob);
            }
            cartes_destroy(seq_iter);
        }
        cass_close(exp(total), 1.0);
    }
    cartes_destroy(ccode_iter);

    nmm_frame_state_destroy(state);
    nmm_base_destroy(base);
    nmm_codont_destroy(codon);
    imm_abc_destroy(abc);
}

void test_frame_state_decode(void)
{
    struct imm_abc* abc = imm_abc_create("ACGT");

    struct nmm_base* base = nmm_base_create(abc);
    nmm_base_set_lprob(base, 'A', log(0.1));
    nmm_base_set_lprob(base, 'C', log(0.2));
    nmm_base_set_lprob(base, 'G', log(0.3));
    nmm_base_set_lprob(base, 'T', log(0.4));
    cass_cond(nmm_base_normalize(base) == 0);

    char const* symbols = imm_abc_symbols(abc);
    int         length = imm_abc_length(abc);

    struct nmm_codont* codon = nmm_codont_create(abc);
    nmm_codont_set_lprob(codon, &NMM_CCODE('A', 'T', 'G'), log(0.8));
    nmm_codont_set_lprob(codon, &NMM_CCODE('A', 'T', 'T'), log(0.1));
    nmm_codont_set_lprob(codon, &NMM_CCODE('G', 'T', 'C'), log(0.4));
    cass_cond(nmm_codont_normalize(codon) == 0);

    struct nmm_frame_state* state = nmm_frame_state_create("State", base, codon, 0.1);

    struct nmm_ccode ccode = NMM_CCODE('X', 'X', 'X');

    cass_close(nmm_frame_state_decode(state, "ATG", 3, &ccode), -0.902566706136);
    cass_cond(ccode.a == 'A' && ccode.b == 'T' && ccode.c == 'G');
    cass_close(nmm_frame_state_decode(state, "ATGT", 4, &ccode), -4.710599080052);
    cass_cond(ccode.a == 'A' && ccode.b == 'T' && ccode.c == 'G');
    cass_close(nmm_frame_state_decode(state, "ATGA", 4, &ccode), -6.097714346951);
    cass_cond(ccode.a == 'A' && ccode.b == 'T' && ccode.c == 'G');
    cass_close(nmm_frame_state_decode(state, "ATGGT", 5, &ccode), -9.031100481720);
    cass_cond(ccode.a == 'A' && ccode.b == 'T' && ccode.c == 'G');
    cass_close(nmm_frame_state_decode(state, "ATT", 3, &ccode), -2.977101440300);
    cass_cond(ccode.a == 'A' && ccode.b == 'T' && ccode.c == 'T');
    cass_close(nmm_frame_state_decode(state, "ATC", 3, &ccode), -7.720225141384);
    cass_cond(ccode.a == 'A' && ccode.b == 'T' && ccode.c == 'G');
    cass_close(nmm_frame_state_decode(state, "TC", 2, &ccode), -4.199089882536);
    cass_cond(ccode.a == 'G' && ccode.b == 'T' && ccode.c == 'C');
    cass_close(nmm_frame_state_decode(state, "A", 1, &ccode), -6.400011321754);
    cass_cond(ccode.a == 'A' && ccode.b == 'T' && ccode.c == 'G');
    cass_close(nmm_frame_state_decode(state, "AG", 2, &ccode), -3.507173471362);
    cass_cond(ccode.a == 'A' && ccode.b == 'T' && ccode.c == 'G');
    cass_close(nmm_frame_state_decode(state, "GC", 2, &ccode), -4.199705077880);
    cass_cond(ccode.a == 'G' && ccode.b == 'T' && ccode.c == 'C');

    nmm_frame_state_destroy(state);
    nmm_base_destroy(base);
    nmm_codont_destroy(codon);
    imm_abc_destroy(abc);
}

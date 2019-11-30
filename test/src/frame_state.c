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
    struct imm_abc* abc = imm_abc_create("ACGT", 'X');

    struct nmm_baset* baset = nmm_baset_create(abc);
    nmm_baset_set_lprob(baset, 'A', log(0.2));
    nmm_baset_set_lprob(baset, 'C', log(0.2));
    nmm_baset_set_lprob(baset, 'G', log(0.2));
    nmm_baset_set_lprob(baset, 'T', log(0.2));
    cass_cond(nmm_baset_normalize(baset) == 0);

    struct nmm_codont* codont = nmm_codont_create(abc);
    nmm_codont_set_lprob(codont, &NMM_CODON('A', 'T', 'G'), log(0.8 / 0.9));
    nmm_codont_set_lprob(codont, &NMM_CODON('A', 'T', 'T'), log(0.1 / 0.9));

    struct nmm_frame_state* state = nmm_frame_state_create("State", baset, codont, 0.1);

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
    nmm_baset_destroy(baset);
    nmm_codont_destroy(codont);
    imm_abc_destroy(abc);
}

void test_frame_state2(void)
{
    struct imm_abc* abc = imm_abc_create("ACGT", 'X');

    struct nmm_baset* baset = nmm_baset_create(abc);
    nmm_baset_set_lprob(baset, 'A', log(0.1));
    nmm_baset_set_lprob(baset, 'C', log(0.2));
    nmm_baset_set_lprob(baset, 'G', log(0.3));
    nmm_baset_set_lprob(baset, 'T', log(0.4));
    cass_cond(nmm_baset_normalize(baset) == 0);

    struct nmm_codont* codont = nmm_codont_create(abc);

    nmm_codont_set_lprob(codont, &NMM_CODON('A', 'T', 'G'), log(0.8 / 0.9));
    nmm_codont_set_lprob(codont, &NMM_CODON('A', 'T', 'T'), log(0.1 / 0.9));

    struct nmm_frame_state* state = nmm_frame_state_create("State", baset, codont, 0.1);

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
    nmm_baset_destroy(baset);
    nmm_codont_destroy(codont);
    imm_abc_destroy(abc);
}

void test_frame_state3(void)
{
    struct imm_abc* abc = imm_abc_create("ACGT", 'X');

    struct nmm_baset* baset = nmm_baset_create(abc);
    nmm_baset_set_lprob(baset, 'A', log(0.1));
    nmm_baset_set_lprob(baset, 'C', log(0.2));
    nmm_baset_set_lprob(baset, 'G', log(0.3));
    nmm_baset_set_lprob(baset, 'T', log(0.4));
    cass_cond(nmm_baset_normalize(baset) == 0);

    struct nmm_codont* codont = nmm_codont_create(abc);
    nmm_codont_set_lprob(codont, &NMM_CODON('A', 'T', 'G'), log(0.8));
    nmm_codont_set_lprob(codont, &NMM_CODON('A', 'T', 'T'), log(0.1));
    nmm_codont_set_lprob(codont, &NMM_CODON('G', 'T', 'C'), log(0.4));
    cass_cond(nmm_codont_normalize(codont) == 0);

    struct nmm_frame_state* state = nmm_frame_state_create("State", baset, codont, 0.1);

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
    nmm_baset_destroy(baset);
    nmm_codont_destroy(codont);
    imm_abc_destroy(abc);
}

void test_frame_state_lposterior(void)
{
    struct imm_abc* abc = imm_abc_create("ACGT", 'X');

    struct nmm_baset* baset = nmm_baset_create(abc);
    nmm_baset_set_lprob(baset, 'A', log(0.1));
    nmm_baset_set_lprob(baset, 'C', log(0.2));
    nmm_baset_set_lprob(baset, 'G', log(0.3));
    nmm_baset_set_lprob(baset, 'T', log(0.4));
    cass_cond(nmm_baset_normalize(baset) == 0);

    char const* symbols = imm_abc_symbols(abc);
    int         length = imm_abc_length(abc);

    struct nmm_codont* codont = nmm_codont_create(abc);
    struct cartes*     codon_iter = cartes_create(symbols, length, 3);
    char const*        codon_item = NULL;

    while ((codon_item = cartes_next(codon_iter)) != NULL) {
        struct nmm_codon codon = NMM_CODON(codon_item[0], codon_item[1], codon_item[2]);
        nmm_codont_set_lprob(codont, &codon, log(0.001));
    }
    cartes_destroy(codon_iter);

    nmm_codont_set_lprob(codont, &NMM_CODON('A', 'T', 'G'), log(0.8));
    nmm_codont_set_lprob(codont, &NMM_CODON('A', 'T', 'T'), log(0.1));
    nmm_codont_set_lprob(codont, &NMM_CODON('G', 'T', 'C'), log(0.4));
    cass_cond(nmm_codont_normalize(codont) == 0);

    struct nmm_frame_state* state = nmm_frame_state_create("State", baset, codont, 0.1);

    codon_iter = cartes_create(symbols, length, 3);

    while ((codon_item = cartes_next(codon_iter)) != NULL) {
        struct nmm_codon codon = {codon_item[0], codon_item[1], codon_item[2]};

        double total = imm_lprob_zero();
        for (int times = 1; times < 6; ++times) {

            struct cartes* seq_iter = cartes_create(symbols, length, times);
            char const*    seq = NULL;

            while ((seq = cartes_next(seq_iter)) != NULL) {
                double lprob = nmm_frame_state_lposterior(state, &codon, seq, times);
                lprob -= nmm_codont_get_lprob(codont, &codon);
                total = imm_lprob_add(total, lprob);
            }
            cartes_destroy(seq_iter);
        }
        cass_close(exp(total), 1.0);
    }
    cartes_destroy(codon_iter);

    nmm_frame_state_destroy(state);
    nmm_baset_destroy(baset);
    nmm_codont_destroy(codont);
    imm_abc_destroy(abc);
}

void test_frame_state_decode(void)
{
    struct imm_abc* abc = imm_abc_create("ACGT", 'X');

    struct nmm_baset* baset = nmm_baset_create(abc);
    nmm_baset_set_lprob(baset, 'A', log(0.1));
    nmm_baset_set_lprob(baset, 'C', log(0.2));
    nmm_baset_set_lprob(baset, 'G', log(0.3));
    nmm_baset_set_lprob(baset, 'T', log(0.4));
    cass_cond(nmm_baset_normalize(baset) == 0);

    char const* symbols = imm_abc_symbols(abc);
    int         length = imm_abc_length(abc);

    struct nmm_codont* codont = nmm_codont_create(abc);
    nmm_codont_set_lprob(codont, &NMM_CODON('A', 'T', 'G'), log(0.8));
    nmm_codont_set_lprob(codont, &NMM_CODON('A', 'T', 'T'), log(0.1));
    nmm_codont_set_lprob(codont, &NMM_CODON('G', 'T', 'C'), log(0.4));
    cass_cond(nmm_codont_normalize(codont) == 0);

    struct nmm_frame_state* state = nmm_frame_state_create("State", baset, codont, 0.1);

    struct nmm_codon codon = NMM_CODON('X', 'X', 'X');

    cass_close(nmm_frame_state_decode(state, "ATG", 3, &codon), -0.902566706136);
    cass_cond(codon.a == 'A' && codon.b == 'T' && codon.c == 'G');
    cass_close(nmm_frame_state_decode(state, "ATGT", 4, &codon), -4.710599080052);
    cass_cond(codon.a == 'A' && codon.b == 'T' && codon.c == 'G');
    cass_close(nmm_frame_state_decode(state, "ATGA", 4, &codon), -6.097714346951);
    cass_cond(codon.a == 'A' && codon.b == 'T' && codon.c == 'G');
    cass_close(nmm_frame_state_decode(state, "ATGGT", 5, &codon), -9.031100481720);
    cass_cond(codon.a == 'A' && codon.b == 'T' && codon.c == 'G');
    cass_close(nmm_frame_state_decode(state, "ATT", 3, &codon), -2.977101440300);
    cass_cond(codon.a == 'A' && codon.b == 'T' && codon.c == 'T');
    cass_close(nmm_frame_state_decode(state, "ATC", 3, &codon), -7.720225141384);
    cass_cond(codon.a == 'A' && codon.b == 'T' && codon.c == 'G');
    cass_close(nmm_frame_state_decode(state, "TC", 2, &codon), -4.199089882536);
    cass_cond(codon.a == 'G' && codon.b == 'T' && codon.c == 'C');
    cass_close(nmm_frame_state_decode(state, "A", 1, &codon), -6.400011321754);
    cass_cond(codon.a == 'A' && codon.b == 'T' && codon.c == 'G');
    cass_close(nmm_frame_state_decode(state, "AG", 2, &codon), -3.507173471362);
    cass_cond(codon.a == 'A' && codon.b == 'T' && codon.c == 'G');
    cass_close(nmm_frame_state_decode(state, "GC", 2, &codon), -4.199705077880);
    cass_cond(codon.a == 'G' && codon.b == 'T' && codon.c == 'C');

    nmm_frame_state_destroy(state);
    nmm_baset_destroy(baset);
    nmm_codont_destroy(codont);
    imm_abc_destroy(abc);
}

#include "cass/cass.h"
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
    struct nmm_base_abc const*   base = nmm_base_abc_create("ACGT", 'X');
    struct imm_abc const*        abc = nmm_base_abc_super(base);
    struct nmm_base_lprob const* basep =
        nmm_base_lprob_create(base, imm_log(0.25), imm_log(0.25), imm_log(0.25), imm_log(0.25));

    struct nmm_codon_lprob* codonp = nmm_codon_lprob_create(base);
    struct nmm_codon*       codon = nmm_codon_create(base);
    nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'T', 'G'));
    nmm_codon_lprob_set(codonp, codon, imm_log(0.8 / 0.9));
    nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'T', 'T'));
    nmm_codon_lprob_set(codonp, codon, imm_log(0.1 / 0.9));
    struct nmm_codon_marg const* codont = nmm_codon_marg_create(codonp);
    nmm_codon_lprob_destroy(codonp);

    struct nmm_frame_state const* state = nmm_frame_state_create("State", basep, codont, (imm_float)0.1);
    const struct imm_state*       s = nmm_frame_state_super(state);

    struct imm_seq const* seq = imm_seq_create("A", abc);
    cass_close(imm_state_lprob(s, seq), -5.914503505971854);
    imm_seq_destroy(seq);
    seq = imm_seq_create("AT", abc);
    cass_close(imm_state_lprob(s, seq), -2.915843423869834);
    imm_seq_destroy(seq);
    seq = imm_seq_create("ATA", abc);
    cass_close(imm_state_lprob(s, seq), -6.905597115665666);
    imm_seq_destroy(seq);
    seq = imm_seq_create("ATG", abc);
    cass_close(imm_state_lprob(s, seq), -0.534773288204706);
    imm_seq_destroy(seq);
    seq = imm_seq_create("ATT", abc);
    cass_close(imm_state_lprob(s, seq), -2.590237330499946);
    imm_seq_destroy(seq);
    seq = imm_seq_create("ATTA", abc);
    cass_close(imm_state_lprob(s, seq), -6.881032208841384);
    imm_seq_destroy(seq);
    seq = imm_seq_create("ATTAA", abc);
    cass_close(imm_state_lprob(s, seq), -12.08828960987379);
    imm_seq_destroy(seq);
    seq = imm_seq_create("ATTAAT", abc);
    cass_cond(imm_lprob_is_zero(imm_state_lprob(s, seq)));
    imm_seq_destroy(seq);

    struct imm_state const* super = nmm_frame_state_super(state);
    cass_cond(nmm_codon_state_derived(super) == NULL);
    cass_cond((state = nmm_frame_state_derived(super)) != NULL);

    nmm_codon_destroy(codon);
    nmm_frame_state_destroy(state);
    nmm_base_lprob_destroy(basep);
    nmm_codon_marg_destroy(codont);
    nmm_base_abc_destroy(base);
}

void test_frame_state2(void)
{
    struct nmm_base_abc const*   base = nmm_base_abc_create("ACGT", 'X');
    struct imm_abc const*        abc = nmm_base_abc_super(base);
    struct nmm_base_lprob const* basep =
        nmm_base_lprob_create(base, imm_log(0.1), imm_log(0.2), imm_log(0.3), imm_log(0.4));

    struct nmm_codon_lprob* codonp = nmm_codon_lprob_create(base);
    struct nmm_codon*       codon = nmm_codon_create(base);
    nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'T', 'G'));
    nmm_codon_lprob_set(codonp, codon, imm_log(0.8 / 0.9));
    nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'T', 'T'));
    nmm_codon_lprob_set(codonp, codon, imm_log(0.1 / 0.9));
    struct nmm_codon_marg const* codont = nmm_codon_marg_create(codonp);
    nmm_codon_lprob_destroy(codonp);

    struct nmm_frame_state const* state = nmm_frame_state_create("State", basep, codont, (imm_float)0.1);
    const struct imm_state*       s = nmm_frame_state_super(state);

    struct imm_seq const* seq = imm_seq_create("A", abc);
    cass_close(imm_state_lprob(s, seq), -5.914503505971854);
    imm_seq_destroy(seq);
    seq = imm_seq_create("C", abc);
    cass_cond(imm_lprob_is_zero(imm_state_lprob(s, seq)));
    imm_seq_destroy(seq);
    seq = imm_seq_create("G", abc);
    cass_close(imm_state_lprob(s, seq), -6.032286541628237);
    imm_seq_destroy(seq);
    seq = imm_seq_create("T", abc);
    cass_close(imm_state_lprob(s, seq), -5.809142990314027);
    imm_seq_destroy(seq);

    seq = imm_seq_create("AT", abc);
    cass_close(imm_state_lprob(s, seq), -2.9159357500274385);
    imm_seq_destroy(seq);
    seq = imm_seq_create("ATA", abc);
    cass_close(imm_state_lprob(s, seq), -7.821518343902165);
    imm_seq_destroy(seq);
    seq = imm_seq_create("ATG", abc);
    cass_close(imm_state_lprob(s, seq), -0.5344319079005616);
    imm_seq_destroy(seq);
    seq = imm_seq_create("ATC", abc);
    cass_close(imm_state_lprob(s, seq), -7.129480084106424);
    imm_seq_destroy(seq);
    seq = imm_seq_create("ATT", abc);
    cass_close(imm_state_lprob(s, seq), -2.57514520832882);
    imm_seq_destroy(seq);

    seq = imm_seq_create("ATTA", abc);
    cass_close(imm_state_lprob(s, seq), -7.789644584138959);
    imm_seq_destroy(seq);
    seq = imm_seq_create("ACTG", abc);
    cass_close(imm_state_lprob(s, seq), -5.036637096635257);
    imm_seq_destroy(seq);

    seq = imm_seq_create("ATTAA", abc);
    cass_close(imm_state_lprob(s, seq), -13.920871073622099);
    imm_seq_destroy(seq);

    seq = imm_seq_create("ATTAAT", abc);
    cass_cond(imm_lprob_is_zero(imm_state_lprob(s, seq)));
    imm_seq_destroy(seq);

    nmm_codon_destroy(codon);
    nmm_frame_state_destroy(state);
    nmm_base_lprob_destroy(basep);
    nmm_codon_marg_destroy(codont);
    nmm_base_abc_destroy(base);
}

void test_frame_state3(void)
{
    struct nmm_base_abc const*   base = nmm_base_abc_create("ACGT", 'X');
    struct imm_abc const*        abc = nmm_base_abc_super(base);
    struct nmm_base_lprob const* basep =
        nmm_base_lprob_create(base, imm_log(0.1), imm_log(0.2), imm_log(0.3), imm_log(0.4));

    struct nmm_codon_lprob* codonp = nmm_codon_lprob_create(base);
    struct nmm_codon*       codon = nmm_codon_create(base);
    nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'T', 'G'));
    nmm_codon_lprob_set(codonp, codon, imm_log(0.8) - imm_log(1.3));
    nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'T', 'T'));
    nmm_codon_lprob_set(codonp, codon, imm_log(0.1) - imm_log(1.3));
    nmm_codon_set_triplet(codon, NMM_TRIPLET('G', 'T', 'C'));
    nmm_codon_lprob_set(codonp, codon, imm_log(0.4) - imm_log(1.3));
    struct nmm_codon_marg const* codont = nmm_codon_marg_create(codonp);
    nmm_codon_lprob_destroy(codonp);

    struct nmm_frame_state const* state = nmm_frame_state_create("State", basep, codont, (imm_float)0.1);
    const struct imm_state*       s = nmm_frame_state_super(state);

    struct imm_seq const* seq = imm_seq_create("A", abc);
    cass_close(imm_state_lprob(s, seq), -6.282228286097171);
    imm_seq_destroy(seq);
    seq = imm_seq_create("C", abc);
    cass_close(imm_state_lprob(s, seq), -7.0931585023135);
    imm_seq_destroy(seq);
    seq = imm_seq_create("G", abc);
    cass_close(imm_state_lprob(s, seq), -5.99454621364539);
    imm_seq_destroy(seq);
    seq = imm_seq_create("T", abc);
    cass_close(imm_state_lprob(s, seq), -5.840395533818132);
    imm_seq_destroy(seq);
    seq = imm_seq_create("AT", abc);
    cass_close(imm_state_lprob(s, seq), -3.283414346005771);
    imm_seq_destroy(seq);
    seq = imm_seq_create("CG", abc);
    cass_close(imm_state_lprob(s, seq), -9.395743595307545);
    imm_seq_destroy(seq);
    seq = imm_seq_create("ATA", abc);
    cass_close(imm_state_lprob(s, seq), -8.18911998648269);
    imm_seq_destroy(seq);
    seq = imm_seq_create("ATG", abc);
    cass_close(imm_state_lprob(s, seq), -0.9021560981322401);
    imm_seq_destroy(seq);
    seq = imm_seq_create("ATT", abc);
    cass_close(imm_state_lprob(s, seq), -2.9428648000333952);
    imm_seq_destroy(seq);
    seq = imm_seq_create("ATC", abc);
    cass_close(imm_state_lprob(s, seq), -7.314811395663229);
    imm_seq_destroy(seq);
    seq = imm_seq_create("GTC", abc);
    cass_close(imm_state_lprob(s, seq), -1.5951613351178675);
    imm_seq_destroy(seq);
    seq = imm_seq_create("ATTA", abc);
    cass_close(imm_state_lprob(s, seq), -8.157369364264277);
    imm_seq_destroy(seq);
    seq = imm_seq_create("GTTC", abc);
    cass_close(imm_state_lprob(s, seq), -4.711642430498609);
    imm_seq_destroy(seq);
    seq = imm_seq_create("ACTG", abc);
    cass_close(imm_state_lprob(s, seq), -5.404361876760574);
    imm_seq_destroy(seq);
    seq = imm_seq_create("ATTAA", abc);
    cass_close(imm_state_lprob(s, seq), -14.288595853747417);
    imm_seq_destroy(seq);
    seq = imm_seq_create("GTCAA", abc);
    cass_close(imm_state_lprob(s, seq), -12.902301492627526);
    imm_seq_destroy(seq);

    nmm_codon_destroy(codon);
    nmm_frame_state_destroy(state);
    nmm_base_lprob_destroy(basep);
    nmm_codon_marg_destroy(codont);
    nmm_base_abc_destroy(base);
}

void test_frame_state_lposterior(void)
{
    struct nmm_base_abc const*   base = nmm_base_abc_create("ACGT", 'X');
    struct imm_abc const*        abc = nmm_base_abc_super(base);
    struct nmm_base_lprob const* basep =
        nmm_base_lprob_create(base, imm_log(0.1), imm_log(0.2), imm_log(0.3), imm_log(0.4));

    char const* symbols = imm_abc_symbols(abc);
    uint16_t    length = imm_abc_length(abc);

    struct imm_cartes* codon_iter = imm_cartes_create(symbols, length, 3);
    char const*        codon_item = NULL;

    struct nmm_codon_lprob* codonp = nmm_codon_lprob_create(base);

    struct nmm_codon* codon = nmm_codon_create(base);
    while ((codon_item = imm_cartes_next(codon_iter)) != NULL) {
        nmm_codon_set_triplet(codon, (struct nmm_triplet){codon_item[0], codon_item[1], codon_item[2]});
        nmm_codon_lprob_set(codonp, codon, imm_log(0.001));
    }
    imm_cartes_destroy(codon_iter);

    nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'T', 'G'));
    nmm_codon_lprob_set(codonp, codon, imm_log(0.8));
    nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'T', 'T'));
    nmm_codon_lprob_set(codonp, codon, imm_log(0.1));
    nmm_codon_set_triplet(codon, NMM_TRIPLET('G', 'T', 'C'));
    nmm_codon_lprob_set(codonp, codon, imm_log(0.4));

    cass_cond(nmm_codon_lprob_normalize(codonp) == 0);
    struct nmm_codon_marg const* codont = nmm_codon_marg_create(codonp);
    nmm_codon_lprob_destroy(codonp);

    struct nmm_frame_state const* state = nmm_frame_state_create("State", basep, codont, (imm_float)0.1);

    codon_iter = imm_cartes_create(symbols, length, 3);

    while ((codon_item = imm_cartes_next(codon_iter)) != NULL) {
        nmm_codon_set_triplet(codon, (struct nmm_triplet){codon_item[0], codon_item[1], codon_item[2]});

        imm_float total = imm_lprob_zero();
        for (uint16_t times = 1; times < 6; ++times) {

            struct imm_cartes* seq_iter = imm_cartes_create(symbols, length, times);
            char const*        seq = NULL;

            while ((seq = imm_cartes_next(seq_iter)) != NULL) {
                struct imm_seq const* tmp = imm_seq_create(seq, abc);
                cass_cond(tmp != NULL);
                imm_float lprob = nmm_frame_state_lposterior(state, codon, tmp);
                imm_seq_destroy(tmp);
                lprob -= nmm_codon_marg_lprob(codont, codon);
                total = imm_lprob_add(total, lprob);
            }
            imm_cartes_destroy(seq_iter);
        }
        cass_close((imm_float)exp(total), 1.0);
    }
    imm_cartes_destroy(codon_iter);

    nmm_base_abc_destroy(base);
    nmm_frame_state_destroy(state);
    nmm_base_lprob_destroy(basep);
    nmm_codon_marg_destroy(codont);
    nmm_codon_destroy(codon);
}

void test_frame_state_decode(void)
{
    struct nmm_base_abc const*   base = nmm_base_abc_create("ACGT", 'X');
    struct imm_abc const*        abc = nmm_base_abc_super(base);
    struct nmm_base_lprob const* basep =
        nmm_base_lprob_create(base, imm_log(0.1), imm_log(0.2), imm_log(0.3), imm_log(0.4));

    struct nmm_codon_lprob* codonp = nmm_codon_lprob_create(base);
    struct nmm_codon*       codon = nmm_codon_create(base);
    nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'T', 'G'));
    nmm_codon_lprob_set(codonp, codon, imm_log(0.8) - imm_log(1.3));
    nmm_codon_set_triplet(codon, NMM_TRIPLET('A', 'T', 'T'));
    nmm_codon_lprob_set(codonp, codon, imm_log(0.1) - imm_log(1.3));
    nmm_codon_set_triplet(codon, NMM_TRIPLET('G', 'T', 'C'));
    nmm_codon_lprob_set(codonp, codon, imm_log(0.4) - imm_log(1.3));
    struct nmm_codon_marg const* codont = nmm_codon_marg_create(codonp);
    nmm_codon_lprob_destroy(codonp);

    struct nmm_frame_state const* state = nmm_frame_state_create("State", basep, codont, (imm_float)0.1);

    struct imm_seq const* seq = imm_seq_create("ATG", abc);
    cass_close(nmm_frame_state_decode(state, seq, codon), -0.902566706136);
    imm_seq_destroy(seq);
    struct nmm_triplet t = nmm_codon_get_triplet(codon);
    cass_cond(t.a == 'A' && t.b == 'T' && t.c == 'G');

    seq = imm_seq_create("ATGT", abc);
    cass_close(nmm_frame_state_decode(state, seq, codon), -4.710599080052);
    imm_seq_destroy(seq);
    t = nmm_codon_get_triplet(codon);
    cass_cond(t.a == 'A' && t.b == 'T' && t.c == 'G');

    seq = imm_seq_create("ATGA", abc);
    cass_close(nmm_frame_state_decode(state, seq, codon), -6.097714346951);
    imm_seq_destroy(seq);
    t = nmm_codon_get_triplet(codon);
    cass_cond(t.a == 'A' && t.b == 'T' && t.c == 'G');

    seq = imm_seq_create("ATGGT", abc);
    cass_close(nmm_frame_state_decode(state, seq, codon), -9.031100481720);
    imm_seq_destroy(seq);
    t = nmm_codon_get_triplet(codon);
    cass_cond(t.a == 'A' && t.b == 'T' && t.c == 'G');

    seq = imm_seq_create("ATT", abc);
    cass_close(nmm_frame_state_decode(state, seq, codon), -2.977101440300);
    imm_seq_destroy(seq);
    t = nmm_codon_get_triplet(codon);
    cass_cond(t.a == 'A' && t.b == 'T' && t.c == 'T');

    seq = imm_seq_create("ATC", abc);
    cass_close(nmm_frame_state_decode(state, seq, codon), -7.720225141384);
    imm_seq_destroy(seq);
    t = nmm_codon_get_triplet(codon);
    cass_cond(t.a == 'A' && t.b == 'T' && t.c == 'G');

    seq = imm_seq_create("TC", abc);
    cass_close(nmm_frame_state_decode(state, seq, codon), -4.199089882536);
    imm_seq_destroy(seq);
    t = nmm_codon_get_triplet(codon);
    cass_cond(t.a == 'G' && t.b == 'T' && t.c == 'C');

    seq = imm_seq_create("A", abc);
    cass_close(nmm_frame_state_decode(state, seq, codon), -6.400011321754);
    imm_seq_destroy(seq);
    t = nmm_codon_get_triplet(codon);
    cass_cond(t.a == 'A' && t.b == 'T' && t.c == 'G');

    seq = imm_seq_create("AG", abc);
    cass_close(nmm_frame_state_decode(state, seq, codon), -3.507173471362);
    imm_seq_destroy(seq);
    t = nmm_codon_get_triplet(codon);
    cass_cond(t.a == 'A' && t.b == 'T' && t.c == 'G');

    seq = imm_seq_create("GC", abc);
    cass_close(nmm_frame_state_decode(state, seq, codon), -4.199705077880);
    imm_seq_destroy(seq);
    t = nmm_codon_get_triplet(codon);
    cass_cond(t.a == 'G' && t.b == 'T' && t.c == 'C');

    nmm_base_abc_destroy(base);
    nmm_frame_state_destroy(state);
    nmm_base_lprob_destroy(basep);
    nmm_codon_marg_destroy(codont);
    nmm_codon_destroy(codon);
}

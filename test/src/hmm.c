#include "cass.h"
#include "imm/imm.h"
#include "nmm/nmm.h"

void test_hmm_frame_state_0eps(void);
void test_hmm_frame_state_len1(void);
void test_hmm_frame_state_len2(void);
void test_hmm_frame_state_len3(void);
void test_hmm_frame_state_len4(void);
void test_hmm_frame_state_len5(void);

int main(void)
{
    test_hmm_frame_state_0eps();
    test_hmm_frame_state_len1();
    test_hmm_frame_state_len2();
    test_hmm_frame_state_len3();
    test_hmm_frame_state_len4();
    test_hmm_frame_state_len5();
    return cass_status();
}

static inline struct imm_state const* cast_c(void const* s) { return imm_state_cast_c(s); }

void test_hmm_frame_state_0eps(void)
{
    struct imm_abc const*   abc = imm_abc_create("ACGT", 'X');
    struct nmm_base const*  base = nmm_base_create(abc);
    double const            zero = imm_lprob_zero();
    struct nmm_baset const* baset =
        nmm_baset_create(base, log(0.25), log(0.25), log(0.5), zero);

    struct nmm_codonp* codonp = nmm_codonp_create(base);
    struct nmm_codon*  codon = nmm_codon_create(base);
    cass_cond(nmm_codon_set(codon, NMM_TRIPLET('A', 'T', 'G')) == 0);
    nmm_codonp_set(codonp, codon, log(0.8));
    cass_cond(nmm_codon_set(codon, NMM_TRIPLET('A', 'T', 'T')) == 0);
    nmm_codonp_set(codonp, codon, log(0.1));
    cass_cond(nmm_codon_set(codon, NMM_TRIPLET('C', 'C', 'C')) == 0);
    nmm_codonp_set(codonp, codon, log(0.1));
    struct nmm_codont const* codont = nmm_codont_create(codonp);
    nmm_codonp_destroy(codonp);

    struct imm_hmm* hmm = imm_hmm_create(abc);

    struct nmm_frame_state const* state = nmm_frame_state_create("M", baset, codont, 0.0);

    imm_hmm_add_state(hmm, cast_c(state), log(1.0));

    struct imm_path* path = imm_path_create();
    imm_path_append(path, cast_c(state), 3);
    struct imm_seq const* seq = imm_seq_create("ATT", abc);
    cass_close(imm_hmm_likelihood(hmm, seq, path), -2.3025850930);
    imm_seq_destroy(seq);
    seq = imm_seq_create("ATG", abc);
    cass_close(imm_hmm_likelihood(hmm, seq, path), -0.2231435513142097);
    imm_seq_destroy(seq);
    seq = imm_seq_create("AT", abc);
    cass_cond(!imm_lprob_is_valid(imm_hmm_likelihood(hmm, seq, path)));
    imm_seq_destroy(seq);
    imm_path_destroy(path);

    path = imm_path_create();
    imm_path_append(path, cast_c(state), 2);
    seq = imm_seq_create("AT", abc);
    cass_cond(imm_lprob_is_zero(imm_hmm_likelihood(hmm, seq, path)));
    imm_seq_destroy(seq);
    imm_path_destroy(path);

    nmm_base_destroy(base);
    nmm_codon_destroy(codon);
    imm_hmm_destroy(hmm);
    nmm_frame_state_destroy(state);
    nmm_codont_destroy(codont);
    nmm_baset_destroy(baset);
    imm_abc_destroy(abc);
}

void test_hmm_frame_state_len1(void)
{
    struct imm_abc const*   abc = imm_abc_create("ACGT", 'X');
    struct nmm_base const*  base = nmm_base_create(abc);
    double const            zero = imm_lprob_zero();
    struct nmm_baset const* baset =
        nmm_baset_create(base, log(0.25), log(0.25), log(0.5), zero);

    struct nmm_codonp* codonp = nmm_codonp_create(base);
    struct nmm_codon*  codon = nmm_codon_create(base);
    cass_cond(nmm_codon_set(codon, NMM_TRIPLET('A', 'T', 'G')) == 0);
    nmm_codonp_set(codonp, codon, log(0.8));
    cass_cond(nmm_codon_set(codon, NMM_TRIPLET('A', 'T', 'T')) == 0);
    nmm_codonp_set(codonp, codon, log(0.1));
    cass_cond(nmm_codon_set(codon, NMM_TRIPLET('C', 'C', 'C')) == 0);
    nmm_codonp_set(codonp, codon, log(0.1));
    struct nmm_codont const* codont = nmm_codont_create(codonp);
    nmm_codonp_destroy(codonp);

    struct imm_hmm* hmm = imm_hmm_create(abc);

    struct nmm_frame_state* state = nmm_frame_state_create("M", baset, codont, 0.1);

    imm_hmm_add_state(hmm, cast_c(state), log(1.0));

    struct imm_path* path = imm_path_create();
    imm_path_append(path, cast_c(state), 1);
    struct imm_seq const* seq = imm_seq_create("A", abc);
    cass_close(imm_hmm_likelihood(hmm, seq, path), -6.0198640216);
    imm_seq_destroy(seq);
    imm_path_destroy(path);

    path = imm_path_create();
    imm_path_append(path, cast_c(state), 1);
    seq = imm_seq_create("C", abc);
    cass_close(imm_hmm_likelihood(hmm, seq, path), -7.118476310297789);
    imm_seq_destroy(seq);
    imm_path_destroy(path);

    path = imm_path_create();
    seq = imm_seq_create("A", abc);
    cass_close(imm_hmm_viterbi(hmm, seq, cast_c(state), path), -6.0198640216);
    cass_close(imm_hmm_likelihood(hmm, seq, path), -6.0198640216);
    imm_seq_destroy(seq);
    imm_path_destroy(path);

    path = imm_path_create();
    seq = imm_seq_create("C", abc);
    cass_close(imm_hmm_viterbi(hmm, seq, cast_c(state), path), -7.1184763103);
    cass_close(imm_hmm_likelihood(hmm, seq, path), -7.1184763103);
    imm_seq_destroy(seq);
    imm_path_destroy(path);

    nmm_base_destroy(base);
    nmm_codon_destroy(codon);
    imm_hmm_destroy(hmm);
    nmm_frame_state_destroy(state);
    nmm_codont_destroy(codont);
    nmm_baset_destroy(baset);
    imm_abc_destroy(abc);
}

void test_hmm_frame_state_len2(void)
{
    struct imm_abc const*   abc = imm_abc_create("ACGT", 'X');
    struct nmm_base const*  base = nmm_base_create(abc);
    double const            zero = imm_lprob_zero();
    struct nmm_baset const* baset =
        nmm_baset_create(base, log(0.25), log(0.25), log(0.5), zero);

    struct nmm_codonp* codonp = nmm_codonp_create(base);
    struct nmm_codon*  codon = nmm_codon_create(base);
    cass_cond(nmm_codon_set(codon, NMM_TRIPLET('A', 'T', 'G')) == 0);
    nmm_codonp_set(codonp, codon, log(0.8));
    cass_cond(nmm_codon_set(codon, NMM_TRIPLET('A', 'T', 'T')) == 0);
    nmm_codonp_set(codonp, codon, log(0.1));
    cass_cond(nmm_codon_set(codon, NMM_TRIPLET('C', 'C', 'C')) == 0);
    nmm_codonp_set(codonp, codon, log(0.1));
    struct nmm_codont const* codont = nmm_codont_create(codonp);
    nmm_codonp_destroy(codonp);

    struct imm_hmm* hmm = imm_hmm_create(abc);

    struct nmm_frame_state* state = nmm_frame_state_create("M", baset, codont, 0.1);

    imm_hmm_add_state(hmm, cast_c(state), log(1.0));

    struct imm_path* path = imm_path_create();
    imm_path_append(path, cast_c(state), 2);
    struct imm_seq const* seq = imm_seq_create("AA", abc);
    cass_close(imm_hmm_likelihood(hmm, seq, path), -8.910235779525845);
    imm_seq_destroy(seq);
    imm_path_destroy(path);

    path = imm_path_create();
    imm_path_append(path, cast_c(state), 2);
    seq = imm_seq_create("TG", abc);
    cass_close(imm_hmm_likelihood(hmm, seq, path), -3.2434246977896133);
    imm_seq_destroy(seq);
    imm_path_destroy(path);

    path = imm_path_create();
    imm_path_append(path, cast_c(state), 2);
    seq = imm_seq_create("CC", abc);
    cass_close(imm_hmm_likelihood(hmm, seq, path), -4.225022885864217);
    imm_seq_destroy(seq);
    imm_path_destroy(path);

    path = imm_path_create();
    imm_path_append(path, cast_c(state), 2);
    seq = imm_seq_create("TT", abc);
    cass_close(imm_hmm_likelihood(hmm, seq, path), -5.326716841069734);
    imm_seq_destroy(seq);
    imm_path_destroy(path);

    path = imm_path_create();
    seq = imm_seq_create("AA", abc);
    cass_close(imm_hmm_viterbi(hmm, seq, cast_c(state), path), -8.910235779525845);
    cass_close(imm_hmm_likelihood(hmm, seq, path), -8.910235779525845);
    imm_seq_destroy(seq);
    imm_path_destroy(path);

    path = imm_path_create();
    seq = imm_seq_create("TG", abc);
    cass_close(imm_hmm_viterbi(hmm, seq, cast_c(state), path), -3.2434246977896133);
    cass_close(imm_hmm_likelihood(hmm, seq, path), -3.2434246977896133);
    imm_seq_destroy(seq);
    imm_path_destroy(path);

    path = imm_path_create();
    seq = imm_seq_create("CC", abc);
    cass_close(imm_hmm_viterbi(hmm, seq, cast_c(state), path), -4.225022885864217);
    cass_close(imm_hmm_likelihood(hmm, seq, path), -4.225022885864217);
    imm_seq_destroy(seq);
    imm_path_destroy(path);

    path = imm_path_create();
    seq = imm_seq_create("TT", abc);
    cass_close(imm_hmm_viterbi(hmm, seq, cast_c(state), path), -5.326716841069734);
    cass_close(imm_hmm_likelihood(hmm, seq, path), -5.326716841069734);
    imm_seq_destroy(seq);
    imm_path_destroy(path);

    nmm_base_destroy(base);
    nmm_codon_destroy(codon);
    imm_hmm_destroy(hmm);
    nmm_frame_state_destroy(state);
    nmm_codont_destroy(codont);
    nmm_baset_destroy(baset);
    imm_abc_destroy(abc);
}

void test_hmm_frame_state_len3(void)
{
    struct imm_abc const*   abc = imm_abc_create("ACGT", 'X');
    struct nmm_base const*  base = nmm_base_create(abc);
    double const            zero = imm_lprob_zero();
    struct nmm_baset const* baset =
        nmm_baset_create(base, log(0.25), log(0.25), log(0.5), zero);

    struct nmm_codonp* codonp = nmm_codonp_create(base);
    struct nmm_codon*  codon = nmm_codon_create(base);
    cass_cond(nmm_codon_set(codon, NMM_TRIPLET('A', 'T', 'G')) == 0);
    nmm_codonp_set(codonp, codon, log(0.8));
    cass_cond(nmm_codon_set(codon, NMM_TRIPLET('A', 'T', 'T')) == 0);
    nmm_codonp_set(codonp, codon, log(0.1));
    cass_cond(nmm_codon_set(codon, NMM_TRIPLET('C', 'C', 'C')) == 0);
    nmm_codonp_set(codonp, codon, log(0.1));
    struct nmm_codont const* codont = nmm_codont_create(codonp);
    nmm_codonp_destroy(codonp);

    struct imm_hmm* hmm = imm_hmm_create(abc);

    struct nmm_frame_state* state = nmm_frame_state_create("M", baset, codont, 0.1);

    imm_hmm_add_state(hmm, cast_c(state), log(1.0));

    struct imm_path* path = imm_path_create();
    imm_path_append(path, cast_c(state), 3);
    struct imm_seq const* seq = imm_seq_create("ATC", abc);
    cass_close(imm_hmm_likelihood(hmm, seq, path), -7.012344487235739);
    imm_seq_destroy(seq);
    imm_path_destroy(path);

    path = imm_path_create();
    imm_path_append(path, cast_c(state), 3);
    seq = imm_seq_create("ATG", abc);
    cass_close(imm_hmm_likelihood(hmm, seq, path), -0.639793371602465);
    imm_seq_destroy(seq);
    imm_path_destroy(path);

    path = imm_path_create();
    seq = imm_seq_create("ATC", abc);
    cass_close(imm_hmm_viterbi(hmm, seq, cast_c(state), path), -7.012344487235739);
    cass_close(imm_hmm_likelihood(hmm, seq, path), -7.012344487235739);
    imm_seq_destroy(seq);
    imm_path_destroy(path);

    path = imm_path_create();
    seq = imm_seq_create("ATG", abc);
    cass_close(imm_hmm_viterbi(hmm, seq, cast_c(state), path), -0.639793371602465);
    cass_close(imm_hmm_likelihood(hmm, seq, path), -0.639793371602465);
    imm_seq_destroy(seq);
    imm_path_destroy(path);

    nmm_base_destroy(base);
    nmm_codon_destroy(codon);
    imm_hmm_destroy(hmm);
    nmm_frame_state_destroy(state);
    nmm_codont_destroy(codont);
    nmm_baset_destroy(baset);
    imm_abc_destroy(abc);
}

void test_hmm_frame_state_len4(void)
{
    struct imm_abc const*   abc = imm_abc_create("ACGT", 'X');
    struct nmm_base const*  base = nmm_base_create(abc);
    double const            zero = imm_lprob_zero();
    struct nmm_baset const* baset =
        nmm_baset_create(base, log(0.25), log(0.25), log(0.5), zero);

    struct nmm_codonp* codonp = nmm_codonp_create(base);
    struct nmm_codon*  codon = nmm_codon_create(base);
    cass_cond(nmm_codon_set(codon, NMM_TRIPLET('A', 'T', 'G')) == 0);
    nmm_codonp_set(codonp, codon, log(0.8));
    cass_cond(nmm_codon_set(codon, NMM_TRIPLET('A', 'T', 'T')) == 0);
    nmm_codonp_set(codonp, codon, log(0.1));
    cass_cond(nmm_codon_set(codon, NMM_TRIPLET('C', 'C', 'C')) == 0);
    nmm_codonp_set(codonp, codon, log(0.1));
    struct nmm_codont const* codont = nmm_codont_create(codonp);
    nmm_codonp_destroy(codonp);

    struct imm_hmm* hmm = imm_hmm_create(abc);

    struct nmm_frame_state* state = nmm_frame_state_create("M", baset, codont, 0.1);

    imm_hmm_add_state(hmm, cast_c(state), log(1.0));

    struct imm_path* path = imm_path_create();
    imm_path_append(path, cast_c(state), 4);
    struct imm_seq const* seq = imm_seq_create("ATCC", abc);
    cass_close(imm_hmm_likelihood(hmm, seq, path), -11.982929094215963);
    imm_seq_destroy(seq);
    imm_path_destroy(path);

    path = imm_path_create();
    seq = imm_seq_create("ATCC", abc);
    cass_close(imm_hmm_viterbi(hmm, seq, cast_c(state), path), -11.982929094215963);
    cass_close(imm_hmm_likelihood(hmm, seq, path), -11.982929094215963);
    imm_seq_destroy(seq);
    imm_path_destroy(path);

    nmm_base_destroy(base);
    nmm_codon_destroy(codon);
    imm_hmm_destroy(hmm);
    nmm_frame_state_destroy(state);
    nmm_codont_destroy(codont);
    nmm_baset_destroy(baset);
    imm_abc_destroy(abc);
}

void test_hmm_frame_state_len5(void)
{
    struct imm_abc const*   abc = imm_abc_create("ACGT", 'X');
    struct nmm_base const*  base = nmm_base_create(abc);
    double const            zero = imm_lprob_zero();
    struct nmm_baset const* baset =
        nmm_baset_create(base, log(0.25), log(0.25), log(0.5), zero);

    struct nmm_codonp* codonp = nmm_codonp_create(base);
    struct nmm_codon*  codon = nmm_codon_create(base);
    cass_cond(nmm_codon_set(codon, NMM_TRIPLET('A', 'T', 'G')) == 0);
    nmm_codonp_set(codonp, codon, log(0.8));
    cass_cond(nmm_codon_set(codon, NMM_TRIPLET('A', 'T', 'T')) == 0);
    nmm_codonp_set(codonp, codon, log(0.1));
    cass_cond(nmm_codon_set(codon, NMM_TRIPLET('C', 'C', 'C')) == 0);
    nmm_codonp_set(codonp, codon, log(0.1));
    struct nmm_codont const* codont = nmm_codont_create(codonp);
    nmm_codonp_destroy(codonp);

    struct imm_hmm* hmm = imm_hmm_create(abc);

    struct nmm_frame_state* state = nmm_frame_state_create("M", baset, codont, 0.1);

    imm_hmm_add_state(hmm, cast_c(state), log(1.0));

    struct imm_path* path = imm_path_create();
    imm_path_append(path, cast_c(state), 5);
    struct imm_seq const* seq = imm_seq_create("ACGTA", abc);
    cass_cond(imm_lprob_is_zero(imm_hmm_likelihood(hmm, seq, path)));
    imm_seq_destroy(seq);
    imm_path_destroy(path);

    path = imm_path_create();
    imm_path_append(path, cast_c(state), 5);
    seq = imm_seq_create("ACTAG", abc);
    cass_close(imm_hmm_likelihood(hmm, seq, path), -10.11420858385178);
    imm_seq_destroy(seq);
    imm_path_destroy(path);

    path = imm_path_create();
    seq = imm_seq_create("ACGTA", abc);
    cass_cond(!imm_lprob_is_valid(imm_hmm_viterbi(hmm, seq, cast_c(state), path)));
    cass_cond(!imm_lprob_is_valid(imm_hmm_likelihood(hmm, seq, path)));
    imm_seq_destroy(seq);
    imm_path_destroy(path);

    path = imm_path_create();
    seq = imm_seq_create("ACTAG", abc);
    cass_close(imm_hmm_viterbi(hmm, seq, cast_c(state), path), -10.11420858385178);
    cass_close(imm_hmm_likelihood(hmm, seq, path), -10.11420858385178);
    imm_seq_destroy(seq);
    imm_path_destroy(path);

    nmm_base_destroy(base);
    nmm_codon_destroy(codon);
    imm_hmm_destroy(hmm);
    nmm_frame_state_destroy(state);
    nmm_codont_destroy(codont);
    nmm_baset_destroy(baset);
    imm_abc_destroy(abc);
}

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
    struct imm_abc* abc = imm_abc_create("ACGT");

    struct nmm_base* base = nmm_base_create(abc);
    nmm_base_set_lprob(base, 'A', log(0.25));
    nmm_base_set_lprob(base, 'C', log(0.25));
    nmm_base_set_lprob(base, 'G', log(0.5));
    nmm_base_set_lprob(base, 'T', imm_lprob_zero());

    struct nmm_codon* codon = nmm_codon_create(abc);
    nmm_codon_set_lprob(codon, &NMM_CCODE('A', 'T', 'G'), log(0.8));
    nmm_codon_set_lprob(codon, &NMM_CCODE('A', 'T', 'T'), log(0.1));
    nmm_codon_set_lprob(codon, &NMM_CCODE('C', 'C', 'C'), log(0.1));

    struct imm_hmm* hmm = imm_hmm_create(abc);

    struct nmm_frame_state* state = nmm_frame_state_create("M", base, codon, 0.0);

    imm_hmm_add_state(hmm, cast_c(state), log(1.0));

    struct imm_path* path = imm_path_create();
    imm_path_append(path, cast_c(state), 3);
    cass_close(imm_hmm_likelihood(hmm, "ATT", path), -2.3025850930);
    cass_close(imm_hmm_likelihood(hmm, "ATG", path), -0.2231435513142097);
    cass_cond(!imm_lprob_is_valid(imm_hmm_likelihood(hmm, "AT", path)));
    imm_path_destroy(path);

    path = imm_path_create();
    imm_path_append(path, cast_c(state), 2);
    cass_cond(imm_lprob_is_zero(imm_hmm_likelihood(hmm, "AT", path)));
    imm_path_destroy(path);

    imm_hmm_destroy(hmm);
    nmm_frame_state_destroy(state);
    nmm_codon_destroy(codon);
    nmm_base_destroy(base);
    imm_abc_destroy(abc);
}

void test_hmm_frame_state_len1(void)
{
    struct imm_abc* abc = imm_abc_create("ACGT");

    struct nmm_base* base = nmm_base_create(abc);
    nmm_base_set_lprob(base, 'A', log(0.25));
    nmm_base_set_lprob(base, 'C', log(0.25));
    nmm_base_set_lprob(base, 'G', log(0.5));
    nmm_base_set_lprob(base, 'T', imm_lprob_zero());

    struct nmm_codon* codon = nmm_codon_create(abc);
    nmm_codon_set_lprob(codon, &NMM_CCODE('A', 'T', 'G'), log(0.8));
    nmm_codon_set_lprob(codon, &NMM_CCODE('A', 'T', 'T'), log(0.1));
    nmm_codon_set_lprob(codon, &NMM_CCODE('C', 'C', 'C'), log(0.1));

    struct imm_hmm* hmm = imm_hmm_create(abc);

    struct nmm_frame_state* state = nmm_frame_state_create("M", base, codon, 0.1);

    imm_hmm_add_state(hmm, cast_c(state), log(1.0));

    struct imm_path* path = imm_path_create();
    imm_path_append(path, cast_c(state), 1);
    cass_close(imm_hmm_likelihood(hmm, "A", path), -6.0198640216);
    imm_path_destroy(path);

    path = imm_path_create();
    imm_path_append(path, cast_c(state), 1);
    cass_close(imm_hmm_likelihood(hmm, "C", path), -7.118476310297789);
    imm_path_destroy(path);

    path = imm_path_create();
    cass_close(imm_hmm_viterbi(hmm, "A", cast_c(state), path), -6.0198640216);
    cass_close(imm_hmm_likelihood(hmm, "A", path), -6.0198640216);
    imm_path_destroy(path);

    path = imm_path_create();
    cass_close(imm_hmm_viterbi(hmm, "C", cast_c(state), path), -7.1184763103);
    cass_close(imm_hmm_likelihood(hmm, "C", path), -7.1184763103);
    imm_path_destroy(path);

    imm_hmm_destroy(hmm);
    nmm_frame_state_destroy(state);
    nmm_codon_destroy(codon);
    nmm_base_destroy(base);
    imm_abc_destroy(abc);
}

void test_hmm_frame_state_len2(void)
{
    struct imm_abc* abc = imm_abc_create("ACGT");

    struct nmm_base* base = nmm_base_create(abc);
    nmm_base_set_lprob(base, 'A', log(0.25));
    nmm_base_set_lprob(base, 'C', log(0.25));
    nmm_base_set_lprob(base, 'G', log(0.5));
    nmm_base_set_lprob(base, 'T', imm_lprob_zero());

    struct nmm_codon* codon = nmm_codon_create(abc);
    nmm_codon_set_lprob(codon, &NMM_CCODE('A', 'T', 'G'), log(0.8));
    nmm_codon_set_lprob(codon, &NMM_CCODE('A', 'T', 'T'), log(0.1));
    nmm_codon_set_lprob(codon, &NMM_CCODE('C', 'C', 'C'), log(0.1));

    struct imm_hmm* hmm = imm_hmm_create(abc);

    struct nmm_frame_state* state = nmm_frame_state_create("M", base, codon, 0.1);

    imm_hmm_add_state(hmm, cast_c(state), log(1.0));

    struct imm_path* path = imm_path_create();
    imm_path_append(path, cast_c(state), 2);
    cass_close(imm_hmm_likelihood(hmm, "AA", path), -8.910235779525845);
    imm_path_destroy(path);

    path = imm_path_create();
    imm_path_append(path, cast_c(state), 2);
    cass_close(imm_hmm_likelihood(hmm, "TG", path), -3.2434246977896133);
    imm_path_destroy(path);

    path = imm_path_create();
    imm_path_append(path, cast_c(state), 2);
    cass_close(imm_hmm_likelihood(hmm, "CC", path), -4.225022885864217);
    imm_path_destroy(path);

    path = imm_path_create();
    imm_path_append(path, cast_c(state), 2);
    cass_close(imm_hmm_likelihood(hmm, "TT", path), -5.326716841069734);
    imm_path_destroy(path);

    path = imm_path_create();
    cass_close(imm_hmm_viterbi(hmm, "AA", cast_c(state), path), -8.910235779525845);
    cass_close(imm_hmm_likelihood(hmm, "AA", path), -8.910235779525845);
    imm_path_destroy(path);

    path = imm_path_create();
    cass_close(imm_hmm_viterbi(hmm, "TG", cast_c(state), path), -3.2434246977896133);
    cass_close(imm_hmm_likelihood(hmm, "TG", path), -3.2434246977896133);
    imm_path_destroy(path);

    path = imm_path_create();
    cass_close(imm_hmm_viterbi(hmm, "CC", cast_c(state), path), -4.225022885864217);
    cass_close(imm_hmm_likelihood(hmm, "CC", path), -4.225022885864217);
    imm_path_destroy(path);

    path = imm_path_create();
    cass_close(imm_hmm_viterbi(hmm, "TT", cast_c(state), path), -5.326716841069734);
    cass_close(imm_hmm_likelihood(hmm, "TT", path), -5.326716841069734);
    imm_path_destroy(path);

    imm_hmm_destroy(hmm);
    nmm_frame_state_destroy(state);
    nmm_codon_destroy(codon);
    nmm_base_destroy(base);
    imm_abc_destroy(abc);
}

void test_hmm_frame_state_len3(void)
{
    struct imm_abc* abc = imm_abc_create("ACGT");

    struct nmm_base* base = nmm_base_create(abc);
    nmm_base_set_lprob(base, 'A', log(0.25));
    nmm_base_set_lprob(base, 'C', log(0.25));
    nmm_base_set_lprob(base, 'G', log(0.5));
    nmm_base_set_lprob(base, 'T', imm_lprob_zero());

    struct nmm_codon* codon = nmm_codon_create(abc);
    nmm_codon_set_lprob(codon, &NMM_CCODE('A', 'T', 'G'), log(0.8));
    nmm_codon_set_lprob(codon, &NMM_CCODE('A', 'T', 'T'), log(0.1));
    nmm_codon_set_lprob(codon, &NMM_CCODE('C', 'C', 'C'), log(0.1));

    struct imm_hmm* hmm = imm_hmm_create(abc);

    struct nmm_frame_state* state = nmm_frame_state_create("M", base, codon, 0.1);

    imm_hmm_add_state(hmm, cast_c(state), log(1.0));

    struct imm_path* path = imm_path_create();
    imm_path_append(path, cast_c(state), 3);
    cass_close(imm_hmm_likelihood(hmm, "ATC", path), -7.012344487235739);
    imm_path_destroy(path);

    path = imm_path_create();
    imm_path_append(path, cast_c(state), 3);
    cass_close(imm_hmm_likelihood(hmm, "ATG", path), -0.639793371602465);
    imm_path_destroy(path);

    path = imm_path_create();
    cass_close(imm_hmm_viterbi(hmm, "ATC", cast_c(state), path), -7.012344487235739);
    cass_close(imm_hmm_likelihood(hmm, "ATC", path), -7.012344487235739);
    imm_path_destroy(path);

    path = imm_path_create();
    cass_close(imm_hmm_viterbi(hmm, "ATG", cast_c(state), path), -0.639793371602465);
    cass_close(imm_hmm_likelihood(hmm, "ATG", path), -0.639793371602465);
    imm_path_destroy(path);

    imm_hmm_destroy(hmm);
    nmm_frame_state_destroy(state);
    nmm_codon_destroy(codon);
    nmm_base_destroy(base);
    imm_abc_destroy(abc);
}

void test_hmm_frame_state_len4(void)
{
    struct imm_abc* abc = imm_abc_create("ACGT");

    struct nmm_base* base = nmm_base_create(abc);
    nmm_base_set_lprob(base, 'A', log(0.25));
    nmm_base_set_lprob(base, 'C', log(0.25));
    nmm_base_set_lprob(base, 'G', log(0.5));
    nmm_base_set_lprob(base, 'T', imm_lprob_zero());

    struct nmm_codon* codon = nmm_codon_create(abc);
    nmm_codon_set_lprob(codon, &NMM_CCODE('A', 'T', 'G'), log(0.8));
    nmm_codon_set_lprob(codon, &NMM_CCODE('A', 'T', 'T'), log(0.1));
    nmm_codon_set_lprob(codon, &NMM_CCODE('C', 'C', 'C'), log(0.1));

    struct imm_hmm* hmm = imm_hmm_create(abc);

    struct nmm_frame_state* state = nmm_frame_state_create("M", base, codon, 0.1);

    imm_hmm_add_state(hmm, cast_c(state), log(1.0));

    struct imm_path* path = imm_path_create();
    imm_path_append(path, cast_c(state), 4);
    cass_close(imm_hmm_likelihood(hmm, "ATCC", path), -11.982929094215963);
    imm_path_destroy(path);

    path = imm_path_create();
    cass_close(imm_hmm_viterbi(hmm, "ATCC", cast_c(state), path), -11.982929094215963);
    cass_close(imm_hmm_likelihood(hmm, "ATCC", path), -11.982929094215963);
    imm_path_destroy(path);

    imm_hmm_destroy(hmm);
    nmm_frame_state_destroy(state);
    nmm_codon_destroy(codon);
    nmm_base_destroy(base);
    imm_abc_destroy(abc);
}

void test_hmm_frame_state_len5(void)
{

    struct imm_abc* abc = imm_abc_create("ACGT");

    struct nmm_base* base = nmm_base_create(abc);
    nmm_base_set_lprob(base, 'A', log(0.25));
    nmm_base_set_lprob(base, 'C', log(0.25));
    nmm_base_set_lprob(base, 'G', log(0.5));
    nmm_base_set_lprob(base, 'T', imm_lprob_zero());

    struct nmm_codon* codon = nmm_codon_create(abc);
    nmm_codon_set_lprob(codon, &NMM_CCODE('A', 'T', 'G'), log(0.8));
    nmm_codon_set_lprob(codon, &NMM_CCODE('A', 'T', 'T'), log(0.1));
    nmm_codon_set_lprob(codon, &NMM_CCODE('C', 'C', 'C'), log(0.1));

    struct imm_hmm* hmm = imm_hmm_create(abc);

    struct nmm_frame_state* state = nmm_frame_state_create("M", base, codon, 0.1);

    imm_hmm_add_state(hmm, cast_c(state), log(1.0));

    struct imm_path* path = imm_path_create();
    imm_path_append(path, cast_c(state), 5);
    cass_cond(imm_lprob_is_zero(imm_hmm_likelihood(hmm, "ACGTA", path)));
    imm_path_destroy(path);

    path = imm_path_create();
    imm_path_append(path, cast_c(state), 5);
    cass_close(imm_hmm_likelihood(hmm, "ACTAG", path), -10.11420858385178);
    imm_path_destroy(path);

    path = imm_path_create();
    cass_cond(!imm_lprob_is_valid(imm_hmm_viterbi(hmm, "ACGTA", cast_c(state), path)));
    cass_cond(!imm_lprob_is_valid(imm_hmm_likelihood(hmm, "ACGTA", path)));
    imm_path_destroy(path);

    path = imm_path_create();
    cass_close(imm_hmm_viterbi(hmm, "ACTAG", cast_c(state), path), -10.11420858385178);
    cass_close(imm_hmm_likelihood(hmm, "ACTAG", path), -10.11420858385178);
    imm_path_destroy(path);

    imm_hmm_destroy(hmm);
    nmm_frame_state_destroy(state);
    nmm_codon_destroy(codon);
    nmm_base_destroy(base);
    imm_abc_destroy(abc);
}

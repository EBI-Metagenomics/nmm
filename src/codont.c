#include "imm/imm.h"
#include "nmm/nmm.h"
#include <math.h>
#include <stdlib.h>

#define ASCII_LAST_STD 127
#define NBASES NMM_CODON_NBASES

struct nmm_codont
{
    const struct imm_abc* abc;
    int                   symbol_idx[ASCII_LAST_STD + 1];
    double                emiss_lprobs[(NBASES + 1) * (NBASES + 1) * (NBASES + 1)];
};

static int  set_lprobs(double* dst_lprobs, int const* symbol_idx,
                       struct nmm_codon_lprob const* lprobs, int lprobs_length);
static void set_symbol_idx(int* symbol_idx, const struct imm_abc* abc);

struct nmm_codont* nmm_codont_create(const struct imm_abc*         abc,
                                     struct nmm_codon_lprob const* lprobs, int lprobs_length)
{
    if (imm_abc_length(abc) != NMM_CODON_NBASES) {
        imm_error("alphabet length is not four");
        return NULL;
    }

    struct nmm_codont* codont = malloc(sizeof(struct nmm_codont));
    codont->abc = abc;
    set_symbol_idx(codont->symbol_idx, abc);
    if (set_lprobs(codont->emiss_lprobs, codont->symbol_idx, lprobs, lprobs_length)) {
        free(codont);
        return NULL;
    }

    return codont;
}

double nmm_codont_lprob(struct nmm_codont const* codont, struct nmm_codon const* codon)
{
    int idx[3] = {imm_abc_symbol_idx(codont->abc, codon->a),
                  imm_abc_symbol_idx(codont->abc, codon->b),
                  imm_abc_symbol_idx(codont->abc, codon->c)};

    if (idx[0] < 0 || idx[1] < 0 || idx[2] < 0) {
        imm_error("nucleotide not found");
        return imm_lprob_invalid();
    }

    int const strides[3] = {(NBASES + 1) * (NBASES + 1), NBASES + 1, 1};
    return codont
        ->emiss_lprobs[strides[0] * idx[0] + strides[1] * idx[1] + strides[2] * idx[2]];
}

void nmm_codont_destroy(struct nmm_codont* codont)
{
    if (!codont)
        return;

    codont->abc = NULL;
    free(codont);
}

struct imm_abc const* nmm_codont_get_abc(const struct nmm_codont* codont)
{
    return codont->abc;
}

static inline void set_zero_lprobs(double* lprobs)
{
    double const zero_lprob = imm_lprob_zero();
    for (int i = 0; i < (NBASES + 1) * (NBASES + 1) * (NBASES + 1); ++i)
        lprobs[i] = zero_lprob;
}

static int set_lprobs(double* dst_lprobs, int const* symbol_idx,
                      struct nmm_codon_lprob const* lprobs, int const lprobs_length)
{
    set_zero_lprobs(dst_lprobs);
    int const strides[3] = {(NBASES + 1) * (NBASES + 1), NBASES + 1, 1};

    for (int i = 0; i < lprobs_length; ++i) {
        int idx[3] = {symbol_idx[(size_t)lprobs[i].codon.a],
                      symbol_idx[(size_t)lprobs[i].codon.b],
                      symbol_idx[(size_t)lprobs[i].codon.c]};

        if (idx[0] < 0 || idx[1] < 0 || idx[2] < 0) {
            imm_error("nucleotide not found");
            return 1;
        }

        dst_lprobs[strides[0] * idx[0] + strides[1] * idx[1] + strides[2] * idx[2]] =
            lprobs[i].lprob;
    }

    return 0;
}

static void set_symbol_idx(int* symbol_idx, const struct imm_abc* abc)
{

    for (int i = 0; i <= ASCII_LAST_STD; ++i)
        symbol_idx[i] = -1;

    char const* symbols = imm_abc_symbols(abc);

    for (int i = 0; i < NMM_CODON_NBASES; ++i)
        symbol_idx[(size_t)symbols[i]] = imm_abc_symbol_idx(abc, symbols[i]);

    symbol_idx[(size_t)imm_abc_any_symbol(abc)] = NMM_CODON_NBASES;
}

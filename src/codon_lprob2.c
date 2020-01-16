#include "codon_lprob2.h"

static void set_symbol_idx(int* symbol_idx, const struct imm_abc* abc);
static int  set_lprobs(struct codon_lprob2* lprob, struct nmm_codon_lprob const* lprobs,
                       int lprobs_length);

struct codon_lprob2* codon_lprob2_create(const struct imm_abc*         abc,
                                         struct nmm_codon_lprob const* lprobs,
                                         int                           lprobs_length)
{

    struct codon_lprob2* codon_lprob = malloc(sizeof(struct codon_lprob2));

    set_symbol_idx(codon_lprob->symbol_idx, abc);

    static int const n = NMM_CODON_NBASES + 1;
    codon_lprob->lprobs = array3d_create(n, n, n);

    if (set_lprobs(codon_lprob, lprobs, lprobs_length)) {
        array3d_destroy(codon_lprob->lprobs);
        free(codon_lprob);
        return NULL;
    }

    return codon_lprob;
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

static inline int is_codon_valid(struct codon_lprob2 const* lprob,
                                 struct nmm_codon const*    codon)
{
    if (lprob->symbol_idx[(size_t)codon->a] < 0)
        return 0;

    if (lprob->symbol_idx[(size_t)codon->b] < 0)
        return 0;

    if (lprob->symbol_idx[(size_t)codon->c] < 0)
        return 0;

    return 1;
}

static int set_lprobs(struct codon_lprob2* lprob, struct nmm_codon_lprob const* lprobs,
                      int const lprobs_length)
{
    array3d_fill(&lprob->lprobs, imm_lprob_zero());

    for (int i = 0; i < lprobs_length; ++i) {

        if (!is_codon_valid(lprob, &lprobs[i].codon)) {
            imm_error("nucleotide not found");
            return 1;
        }
        codon_lprob2_set(lprob, &lprobs[i].codon, lprobs[i].lprob);
    }

    return 0;
}

void codon_lprob2_destroy(struct codon_lprob2* codon_lprob)
{
    if (!codon_lprob) {
        imm_error("codon_lprob should not be NULL");
        return;
    }

    array3d_destroy(codon_lprob->lprobs);
    free(codon_lprob);
}

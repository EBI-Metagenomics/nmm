#include "codon_marg.h"
#include "codon_iter.h"
#include "free.h"
#include "nmm/nmm.h"
#include <stdlib.h>
#include <string.h>

/**
 * Number of bases (four) plus the special any-symbol one.
 */
#define NSYMBOLS (NMM_BASE_ABC_SIZE + 1)

static imm_float              marginalization(struct nmm_codon_marg const* codonm, char const* symbols,
                                              struct nmm_codon const* codon);
static struct nmm_codon_marg* new_codon_marg(struct nmm_base_abc const* base_abc);
static inline int             not_marginal(struct nmm_codon const* codon, char const any_symbol);
static inline void set_marginal_lprob(struct nmm_codon_marg* codonm, struct nmm_codon const* codon, imm_float lprob);
static void        set_marginal_lprobs(struct nmm_codon_marg* codonm, struct nmm_base_abc const* base_abc);
static int         set_nonmarginal_lprobs(struct nmm_codon_marg* codonm, struct nmm_codon_lprob const* codonp);
static void        set_symbol_index(struct nmm_codon_marg* codonm);

struct nmm_codon_marg const* nmm_codon_marg_create(struct nmm_codon_lprob const* codonp)
{
    struct nmm_codon_marg* codonm = new_codon_marg(nmm_codon_lprob_abc(codonp));

    set_symbol_index(codonm);

    codonm->lprobs = nmm_array3d_create(NSYMBOLS, NSYMBOLS, NSYMBOLS);

    if (set_nonmarginal_lprobs(codonm, codonp)) {
        nmm_array3d_destroy(codonm->lprobs);
        free_c(codonm);
        return NULL;
    }

    set_marginal_lprobs(codonm, codonm->base_abc);

    return codonm;
}

void nmm_codon_marg_destroy(struct nmm_codon_marg const* codonm)
{
    nmm_array3d_destroy(codonm->lprobs);
    free_c(codonm);
}

struct nmm_codon_marg const* codon_marg_read(FILE* stream, struct nmm_base_abc const* base_abc)
{
    struct nmm_codon_marg* codonm = new_codon_marg(base_abc);

    if (fread(codonm->symbol_idx, sizeof(*codonm->symbol_idx), NMM_ASCII_LAST_STD + 1, stream) <
        NMM_ASCII_LAST_STD + 1) {
        imm_error("could not read symbol_idx");
        free_c(codonm);
        return NULL;
    }

    if (nmm_array3d_read(&codonm->lprobs, stream)) {
        imm_error("could not read lprobs");
        free_c(codonm);
        return NULL;
    }

    return codonm;
}

int codon_marg_write(struct nmm_codon_marg const* codonm, FILE* stream)
{
    if (fwrite(codonm->symbol_idx, sizeof(*codonm->symbol_idx), NMM_ASCII_LAST_STD + 1, stream) <
        NMM_ASCII_LAST_STD + 1) {
        imm_error("could not write symbol_idx");
        return 1;
    }

    if (nmm_array3d_write(&codonm->lprobs, stream)) {
        imm_error("could not write lprobs");
        return 1;
    }

    return 0;
}

static imm_float marginalization(struct nmm_codon_marg const* codonm, char const* symbols,
                                 struct nmm_codon const* codon)
{
    char const any_symbol = symbols[NSYMBOLS - 1];

    struct nmm_triplet t = nmm_codon_get_triplet(codon);
    char const         seq[3] = {t.a, t.b, t.c};

    char const* arr[3];
    unsigned    shape[3];
    for (unsigned i = 0; i < 3; ++i) {
        if (seq[i] == any_symbol) {
            arr[i] = symbols;
            shape[i] = NMM_BASE_ABC_SIZE;
        } else {
            arr[i] = seq + i;
            shape[i] = 1;
        }
    }

    struct nmm_codon* tmp = nmm_codon_create(nmm_codon_abc(codon));
    imm_float         lprob = imm_lprob_zero();
    for (unsigned a = 0; a < shape[0]; ++a) {
        for (unsigned b = 0; b < shape[1]; ++b) {
            for (unsigned c = 0; c < shape[2]; ++c) {

                t.a = arr[0][a];
                t.b = arr[1][b];
                t.c = arr[2][c];
                nmm_codon_set_triplet(tmp, t);
                lprob = imm_lprob_add(lprob, nmm_codon_marg_lprob(codonm, tmp));
            }
        }
    }
    nmm_codon_destroy(tmp);

    return lprob;
}

static struct nmm_codon_marg* new_codon_marg(struct nmm_base_abc const* base_abc)
{
    struct nmm_codon_marg* codonm = malloc(sizeof(*codonm));
    codonm->base_abc = base_abc;
    memset(codonm->symbol_idx, 0, NMM_ASCII_LAST_STD + 1);
    return codonm;
}

static inline int not_marginal(struct nmm_codon const* codon, char const any_symbol)
{
    struct nmm_triplet const t = nmm_codon_get_triplet(codon);
    return t.a != any_symbol && t.b != any_symbol && t.c != any_symbol;
}

static inline void set_marginal_lprob(struct nmm_codon_marg* codonm, struct nmm_codon const* codon, imm_float lprob)
{
    nmm_array3d_set(&codonm->lprobs, __nmm_codon_marg_array_idx(codonm, codon), lprob);
}

static void set_marginal_lprobs(struct nmm_codon_marg* codonm, struct nmm_base_abc const* base_abc)
{
    struct imm_abc const* abc = nmm_base_abc_super(base_abc);
    char const            any_symbol = imm_abc_any_symbol(abc);
    char const            symbols[5] = {imm_abc_symbol_id(abc, 0), imm_abc_symbol_id(abc, 1), imm_abc_symbol_id(abc, 2),
                             imm_abc_symbol_id(abc, 3), any_symbol};

    struct nmm_codon* codon = nmm_codon_create(base_abc);

    for (unsigned i0 = 0; i0 < NSYMBOLS; ++i0) {
        for (unsigned i1 = 0; i1 < NSYMBOLS; ++i1) {
            for (unsigned i2 = 0; i2 < NSYMBOLS; ++i2) {

                struct nmm_triplet const t = {symbols[i0], symbols[i1], symbols[i2]};
                nmm_codon_set_triplet(codon, t);

                if (not_marginal(codon, any_symbol))
                    continue;

                set_marginal_lprob(codonm, codon, marginalization(codonm, symbols, codon));
            }
        }
    }

    nmm_codon_destroy(codon);
}

static int set_nonmarginal_lprobs(struct nmm_codon_marg* codonm, struct nmm_codon_lprob const* codonp)
{
    struct codon_iter iter = codon_iter_begin(codonm->base_abc);
    while (!codon_iter_end(iter)) {
        struct nmm_codon const codon = codon_iter_next(&iter);
        set_marginal_lprob(codonm, &codon, nmm_codon_lprob_get(codonp, &codon));
    }
    return 0;
}

static void set_symbol_index(struct nmm_codon_marg* codonm)
{
    struct imm_abc const* abc = nmm_base_abc_super(codonm->base_abc);

    char const* symbols = imm_abc_symbols(abc);

    for (unsigned i = 0; i < NMM_BASE_ABC_SIZE; ++i) {
        size_t j = (size_t)symbols[i];
        codonm->symbol_idx[j] = imm_abc_symbol_idx(abc, symbols[i]);
    }

    codonm->symbol_idx[(size_t)imm_abc_any_symbol(abc)] = NMM_BASE_ABC_SIZE;
}

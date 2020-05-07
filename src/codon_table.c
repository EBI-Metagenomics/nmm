#include "codon_table.h"
#include "codon_iter.h"
#include "free.h"
#include "imm/imm.h"
#include "nmm/base_abc.h"
#include "nmm/codon_lprob.h"
#include "nmm/codon_table.h"
#include <stdlib.h>

/**
 * Number of bases (four) plus the special any-symbol one.
 */
#define NSYMBOLS (NMM_BASE_ABC_SIZE + 1)

static void set_symbol_index(struct nmm_codon_table* codont);
static int  set_nonmarginal_lprobs(struct nmm_codon_table*       codont,
                                   struct nmm_codon_lprob const* codonp);
static void set_marginal_lprobs(struct nmm_codon_table*    codont,
                                struct nmm_base_abc const* base_abc);

static double marginalization(struct nmm_codon_table const* codont, char const* symbols,
                              struct nmm_codon const* codon);

static inline int not_marginal(struct nmm_codon const* codon, char const any_symbol)
{
    struct nmm_triplet const t = nmm_codon_get_triplet(codon);
    return t.a != any_symbol && t.b != any_symbol && t.c != any_symbol;
}

struct nmm_codon_table const* nmm_codon_table_create(struct nmm_codon_lprob const* codonp)
{
    struct nmm_codon_table* codont = malloc(sizeof(struct nmm_codon_table));
    codont->base_abc = nmm_codon_lprob_abc(codonp);

    set_symbol_index(codont);

    codont->lprobs = nmm_array3d_create(NSYMBOLS, NSYMBOLS, NSYMBOLS);

    if (set_nonmarginal_lprobs(codont, codonp)) {
        nmm_array3d_destroy(codont->lprobs);
        free_c(codont);
        return NULL;
    }

    set_marginal_lprobs(codont, codont->base_abc);

    return codont;
}

void nmm_codon_table_destroy(struct nmm_codon_table const* codont)
{
    nmm_array3d_destroy(codont->lprobs);
    free_c(codont);
}

struct nmm_codon_table const* codon_table_read(FILE* stream, struct nmm_base_abc const* base_abc)
{
    struct nmm_codon_table* codont = malloc(sizeof(*codont));
    codont->base_abc = base_abc;

    if (fread(codont->symbol_idx, sizeof(*codont->symbol_idx), NMM_ASCII_LAST_STD + 1, stream) <
        NMM_ASCII_LAST_STD + 1) {
        imm_error("could not read symbol_idx");
        free_c(codont);
        return NULL;
    }

    if (nmm_array3d_read(&codont->lprobs, stream)) {
        imm_error("could not read lprobs");
        free_c(codont);
        return NULL;
    }

    return codont;
}

int codon_table_write(struct nmm_codon_table const* codont, FILE* stream)
{
    if (fwrite(codont->symbol_idx, sizeof(*codont->symbol_idx), NMM_ASCII_LAST_STD + 1, stream) <
        NMM_ASCII_LAST_STD + 1) {
        imm_error("could not write symbol_idx");
        return 1;
    }

    if (nmm_array3d_write(&codont->lprobs, stream)) {
        imm_error("could not write lprobs");
        return 1;
    }

    return 0;
}

static void set_symbol_index(struct nmm_codon_table* codont)
{
    struct imm_abc const* abc = nmm_base_abc_super(codont->base_abc);

    char const* symbols = imm_abc_symbols(abc);

    for (unsigned i = 0; i < NMM_BASE_ABC_SIZE; ++i) {
        size_t j = (size_t)symbols[i];
        codont->symbol_idx[j] = imm_abc_symbol_idx(abc, symbols[i]);
    }

    codont->symbol_idx[(size_t)imm_abc_any_symbol(abc)] = NMM_BASE_ABC_SIZE;
}

static inline void set_marginal_lprob(struct nmm_codon_table* codont, struct nmm_codon const* codon,
                                      double lprob)
{
    nmm_array3d_set(&codont->lprobs, __nmm_codon_table_get_array_idx(codont, codon), lprob);
}

static int set_nonmarginal_lprobs(struct nmm_codon_table*       codont,
                                  struct nmm_codon_lprob const* codonp)
{
    struct codon_iter iter = codon_iter_begin(codont->base_abc);
    while (!codon_iter_end(iter)) {
        struct nmm_codon const codon = codon_iter_next(&iter);
        set_marginal_lprob(codont, &codon, nmm_codon_lprob_get(codonp, &codon));
    }
    return 0;
}

static void set_marginal_lprobs(struct nmm_codon_table* codont, struct nmm_base_abc const* base_abc)
{
    struct imm_abc const* abc = nmm_base_abc_super(base_abc);
    char const            any_symbol = imm_abc_any_symbol(abc);
    char const            symbols[5] = {imm_abc_symbol_id(abc, 0), imm_abc_symbol_id(abc, 1),
                             imm_abc_symbol_id(abc, 2), imm_abc_symbol_id(abc, 3), any_symbol};

    struct nmm_codon* codon = nmm_codon_create(base_abc);

    for (unsigned i0 = 0; i0 < NSYMBOLS; ++i0) {
        for (unsigned i1 = 0; i1 < NSYMBOLS; ++i1) {
            for (unsigned i2 = 0; i2 < NSYMBOLS; ++i2) {

                struct nmm_triplet const t = {symbols[i0], symbols[i1], symbols[i2]};
                nmm_codon_set_triplet(codon, t);

                if (not_marginal(codon, any_symbol))
                    continue;

                set_marginal_lprob(codont, codon, marginalization(codont, symbols, codon));
            }
        }
    }

    nmm_codon_destroy(codon);
}

static double marginalization(struct nmm_codon_table const* codont, char const* symbols,
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
    double            lprob = imm_lprob_zero();
    for (unsigned a = 0; a < shape[0]; ++a) {
        for (unsigned b = 0; b < shape[1]; ++b) {
            for (unsigned c = 0; c < shape[2]; ++c) {

                t.a = arr[0][a];
                t.b = arr[1][b];
                t.c = arr[2][c];
                nmm_codon_set_triplet(tmp, t);
                lprob = imm_lprob_add(lprob, nmm_codon_table_lprob(codont, tmp));
            }
        }
    }
    nmm_codon_destroy(tmp);

    return lprob;
}

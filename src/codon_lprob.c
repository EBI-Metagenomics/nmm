#include "codon_lprob.h"
#include "logaddexp.h"

static double compute_codon_lprob(char const* bases_comb, struct nmm_codont const* codont,
                                  char const* codon);

void nmm_codon_lprob_init(struct codon_lprob* codon_lprob, struct nmm_codont const* codont)
{
    struct imm_abc const* abc = nmm_codont_get_abc(codont);
    char const            bases_comb[3 * NMM_CODON_NBASES] = {
        imm_abc_symbol_id(abc, 0), imm_abc_symbol_id(abc, 1), imm_abc_symbol_id(abc, 2),
        imm_abc_symbol_id(abc, 3), imm_abc_symbol_id(abc, 0), imm_abc_symbol_id(abc, 1),
        imm_abc_symbol_id(abc, 2), imm_abc_symbol_id(abc, 3), imm_abc_symbol_id(abc, 0),
        imm_abc_symbol_id(abc, 1), imm_abc_symbol_id(abc, 2), imm_abc_symbol_id(abc, 3)};

    char const any_symbol = imm_abc_any_symbol(abc);

    for (int i = 0; i <= ASCII_LAST_STD; ++i)
        codon_lprob->symbol_idx[i] = -1;

    char const* symbols = imm_abc_symbols(abc);
    char const  all_symbols[NSYMBOLS] = {symbols[0], symbols[1], symbols[2], symbols[3],
                                        any_symbol};

    for (int i = 0; i < NMM_CODON_NBASES; ++i)
        codon_lprob->symbol_idx[(size_t)symbols[i]] = imm_abc_symbol_idx(abc, symbols[i]);

    codon_lprob->symbol_idx[(size_t)any_symbol] = NMM_CODON_NBASES;

    double* lprob = codon_lprob->lprob;
    for (int i0 = 0; i0 < NSYMBOLS; ++i0) {
        for (int i1 = 0; i1 < NSYMBOLS; ++i1) {
            for (int i2 = 0; i2 < NSYMBOLS; ++i2) {

                char const codon[3] = {all_symbols[i0], all_symbols[i1], all_symbols[i2]};

                int const i = codon_lprob_idx(codon_lprob, codon);
                lprob[i] = compute_codon_lprob(bases_comb, codont, codon);
            }
        }
    }
}

static double compute_codon_lprob(char const* bases_comb, struct nmm_codont const* codont,
                                  char const* codon)
{
    double lprob = imm_lprob_zero();
    char   bases[3 * 4] = {
        bases_comb[0], bases_comb[1], bases_comb[2], bases_comb[3],
        bases_comb[0], bases_comb[1], bases_comb[2], bases_comb[3],
        bases_comb[0], bases_comb[1], bases_comb[2], bases_comb[3],
    };
    int nbases[3] = {NMM_CODON_NBASES, NMM_CODON_NBASES, NMM_CODON_NBASES};

    char const any_symbol = imm_abc_any_symbol(nmm_codont_get_abc(codont));
    for (int i = 0; i < 3; ++i) {
        if (codon[i] != any_symbol) {
            bases[i * NMM_CODON_NBASES] = codon[i];
            nbases[i] = 1;
        }
    }

    char const* a_id = bases;
    char const* b_id = bases + NMM_CODON_NBASES;
    char const* c_id = bases + 2 * NMM_CODON_NBASES;
    for (int a = 0; a < nbases[0]; ++a) {
        for (int b = 0; b < nbases[1]; ++b) {
            for (int c = 0; c < nbases[2]; ++c) {
                struct nmm_codon const ccode = {a_id[a], b_id[b], c_id[c]};
                double const           t = nmm_codont_lprob(codont, &ccode);
                lprob = logaddexp(lprob, t);
            }
        }
    }

    return lprob;
}

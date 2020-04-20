#include "nmm/codon_state.h"
#include "free.h"
#include "imm/imm.h"
#include "nmm/base_abc.h"
#include "nmm/codon.h"
#include "nmm/codon_lprob.h"
#include "nmm/state_types.h"
#include <stdlib.h>

struct nmm_codon_state
{
    struct imm_state const*       interface;
    struct nmm_codon_lprob const* codonp;
    struct nmm_base_abc const*    base_abc;
};

static double   codon_state_lprob(struct imm_state const* state, struct imm_seq const* seq);
static unsigned codon_state_min_seq(struct imm_state const* state);
static unsigned codon_state_max_seq(struct imm_state const* state);
static int      codon_state_write(struct imm_state const* state, FILE* stream);

struct nmm_codon_state const* nmm_codon_state_create(char const*                   name,
                                                     struct nmm_codon_lprob const* codonp)
{
    struct nmm_codon_state* state = malloc(sizeof(struct nmm_codon_state));

    struct imm_state_funcs funcs = {codon_state_lprob, codon_state_min_seq,
                                    codon_state_max_seq, codon_state_write};

    state->interface =
        imm_state_create(name, nmm_base_abc_cast(nmm_codon_lprob_get_base_abc(codonp)), funcs,
                         NMM_CODON_STATE_TYPE_ID, state);

    state->codonp = codonp;
    state->base_abc = nmm_codon_lprob_get_base_abc(codonp);

    return state;
}

void nmm_codon_state_destroy(struct nmm_codon_state const* state)
{
    imm_state_destroy(state->interface);
    free_c(state);
}

static double codon_state_lprob(struct imm_state const* state, struct imm_seq const* seq)
{
    struct nmm_codon_state const* s = imm_state_get_impl(state);
    unsigned                      length = imm_seq_length(seq);

    if (length != 3)
        return imm_lprob_invalid();

    NMM_CODON_DECL(codon, s->base_abc);

    struct nmm_triplet const triplet =
        NMM_TRIPLET(seq->string[0], seq->string[1], seq->string[2]);

    if (nmm_codon_set_triplet(&codon, triplet))
        return imm_lprob_invalid();

    return nmm_codon_lprob_get(s->codonp, &codon);
}

static unsigned codon_state_min_seq(struct imm_state const* state) { return 3; }

static unsigned codon_state_max_seq(struct imm_state const* state) { return 3; }

static int codon_state_write(struct imm_state const* state, FILE* stream) { return 0; }

#include "codon_state.h"
#include "free.h"
#include "imm/imm.h"
#include "nmm/base_abc.h"
#include "nmm/codon.h"
#include "nmm/codon_lprob.h"
#include "nmm/codon_state.h"
#include "nmm/state_types.h"
#include <stdlib.h>

struct nmm_codon_state
{
    struct imm_state const*       parent;
    struct nmm_codon_lprob const* codonp;
    struct nmm_base_abc const*    base_abc;
};

static uint8_t  codon_state_type_id(struct imm_state const* state);
static double   codon_state_lprob(struct imm_state const* state, struct imm_seq const* seq);
static unsigned codon_state_min_seq(struct imm_state const* state);
static unsigned codon_state_max_seq(struct imm_state const* state);
static int      codon_state_write(struct imm_state const* state, FILE* stream);
static void     codon_state_destroy(struct imm_state const* state);

static struct imm_state_vtable const vtable = {codon_state_type_id, codon_state_lprob,
                                               codon_state_min_seq, codon_state_max_seq,
                                               codon_state_write,   codon_state_destroy};

struct nmm_codon_state const* nmm_codon_state_create(char const*                   name,
                                                     struct nmm_codon_lprob const* codonp)
{
    struct nmm_codon_state* state = malloc(sizeof(struct nmm_codon_state));

    state->codonp = codonp;
    state->base_abc = nmm_codon_lprob_get_base_abc(codonp);

    struct imm_abc const* abc = nmm_base_abc_cast(nmm_codon_lprob_get_base_abc(codonp));
    state->parent = imm_state_create(name, abc, vtable, state);
    return state;
}

void nmm_codon_state_destroy(struct nmm_codon_state const* state)
{
    struct imm_state const* parent = state->parent;
    codon_state_destroy(parent);
    __imm_state_destroy_parent(parent);
}

struct imm_state const* nmm_codon_state_parent(struct nmm_codon_state const* state)
{
    return state->parent;
}

struct nmm_codon_state const* nmm_codon_state_child(struct imm_state const* state)
{
    if (imm_state_type_id(state) != NMM_CODON_STATE_TYPE_ID) {
        imm_error("could not cast to codon_state");
        return NULL;
    }
    return __imm_state_child(state);
}

static uint8_t codon_state_type_id(struct imm_state const* state)
{
    return NMM_CODON_STATE_TYPE_ID;
}

static double codon_state_lprob(struct imm_state const* state, struct imm_seq const* seq)
{
    struct nmm_codon_state const* s = __imm_state_child(state);
    unsigned                      length = imm_seq_length(seq);

    if (length != 3)
        return imm_lprob_invalid();

    NMM_CODON_DECL(codon, s->base_abc);

    struct nmm_triplet const triplet = NMM_TRIPLET(seq->string[0], seq->string[1], seq->string[2]);

    if (nmm_codon_set_triplet(&codon, triplet))
        return imm_lprob_invalid();

    return nmm_codon_lprob_get(s->codonp, &codon);
}

static unsigned codon_state_min_seq(struct imm_state const* state) { return 3; }

static unsigned codon_state_max_seq(struct imm_state const* state) { return 3; }

static int codon_state_write(struct imm_state const* state, FILE* stream) { return 0; }

static void codon_state_destroy(struct imm_state const* state) { free_c(__imm_state_child(state)); }

struct imm_state const* codon_state_read(FILE* stream) { return NULL; }

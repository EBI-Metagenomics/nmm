#include "free.h"
#include "nmm/nmm.h"
#include "profile.h"
#include <stdlib.h>

struct nmm_codon_state
{
    struct imm_state const*       super;
    struct nmm_codon_lprob const* codonp;
    struct nmm_base_abc const*    base_abc;
};

static void      destroy(struct imm_state const* state);
static imm_float lprob(struct imm_state const* state, struct imm_seq const* seq);
static uint8_t   max_seq(struct imm_state const* state);
static uint8_t   min_seq(struct imm_state const* state);
static uint8_t   type_id(struct imm_state const* state);

static struct imm_state_vtable const __vtable = {destroy, lprob, max_seq, min_seq, type_id};

struct nmm_codon_lprob const* nmm_codon_state_codon_lprob(struct nmm_codon_state const* state) { return state->codonp; }

struct nmm_codon_state const* nmm_codon_state_create(char const* name, struct nmm_codon_lprob const* codonp)
{
    struct nmm_codon_state* state = malloc(sizeof(*state));

    state->codonp = codonp;
    state->base_abc = nmm_codon_lprob_abc(codonp);

    struct imm_abc const* abc = nmm_base_abc_super(nmm_codon_lprob_abc(codonp));
    state->super = imm_state_create(name, abc, __vtable, state);
    return state;
}

struct nmm_codon_state const* nmm_codon_state_derived(struct imm_state const* state)
{
    if (imm_state_type_id(state) != NMM_CODON_STATE_TYPE_ID) {
        imm_error("could not cast to codon_state");
        return NULL;
    }
    return __imm_state_derived(state);
}

void nmm_codon_state_destroy(struct nmm_codon_state const* state) { state->super->vtable.destroy(state->super); }

struct imm_state const* nmm_codon_state_super(struct nmm_codon_state const* state) { return state->super; }

struct imm_state const* nmm_codon_state_read(FILE* stream, struct nmm_profile const* prof)
{
    struct imm_abc const* abc = nmm_profile_abc(prof);
    struct imm_state*     state = __imm_state_read(stream, abc);
    if (!state) {
        imm_error("could not state_read");
        return NULL;
    }

    struct nmm_codon_state* codon_state = malloc(sizeof(*codon_state));
    codon_state->super = state;
    state->vtable = __vtable;
    state->derived = codon_state;

    if (!(codon_state->base_abc = nmm_base_abc_derived(abc))) {
        imm_error("expected base_abc");
        goto err;
    }

    uint16_t index = 0;
    if (fread(&index, sizeof(index), 1, stream) < 1) {
        imm_error("could not read codonp index");
        goto err;
    }

    struct nmm_codon_lprob const* codonp = nmm_profile_codon_lprob(prof, index);
    if (!codonp) {
        imm_error("could not get codonp");
        goto err;
    }
    codon_state->codonp = codonp;

    return state;

err:
    free_c(codon_state);
    __imm_state_destroy(state);
    return NULL;
}

static void destroy(struct imm_state const* state)
{
    free_c(__imm_state_derived(state));
    __imm_state_destroy(state);
}

static imm_float lprob(struct imm_state const* state, struct imm_seq const* seq)
{
    struct nmm_codon_state const* s = __imm_state_derived(state);
    uint32_t                      length = imm_seq_length(seq);

    if (length != 3)
        return imm_lprob_invalid();

    NMM_CODON_DECL(codon, s->base_abc);

    struct nmm_triplet const triplet = NMM_TRIPLET(seq->string[0], seq->string[1], seq->string[2]);

    if (nmm_codon_set_triplet(&codon, triplet))
        return imm_lprob_invalid();

    return nmm_codon_lprob_get(s->codonp, &codon);
}

static uint8_t max_seq(struct imm_state const* state) { return 3; }

static uint8_t min_seq(struct imm_state const* state) { return 3; }

static uint8_t type_id(struct imm_state const* state) { return NMM_CODON_STATE_TYPE_ID; }

int nmm_codon_state_write(struct imm_state const* state, struct nmm_profile const* prof, FILE* stream)
{
    struct nmm_codon_state const* s = nmm_codon_state_derived(state);
    uint16_t                      index = profile_codon_lprob_index(prof, s->codonp);

    if (__imm_state_write(state, stream)) {
        imm_error("could not write super state");
        return 1;
    }

    if (fwrite(&index, sizeof(index), 1, stream) < 1) {
        imm_error("could not write codonp index");
        return 1;
    }

    return 0;
}

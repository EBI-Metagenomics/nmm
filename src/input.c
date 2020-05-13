#include "nmm/input.h"
#include "free.h"
#include "imm/imm.h"
#include "model.h"
#include "nmm/io.h"
#include <stdlib.h>
#include <string.h>

struct nmm_input
{
    struct imm_input* super;

    uint32_t                       nbaset;
    struct nmm_base_table const**  baset_ptrs;
    uint32_t                       ncodonp;
    struct nmm_codon_lprob const** codonp_ptrs;
    uint32_t                       ncodont;
    struct nmm_codon_table const** codont_ptrs;
};

/* model->nbaset = 0; */
/* model->baset_ptrs = NULL; */
/* model->ncodonp = 0; */
/* model->codonp_ptrs = NULL; */
/* model->ncodont = 0; */
/* model->codont_ptrs = NULL; */

struct imm_model const* read_block(struct imm_input* input, uint8_t block_type);

static struct imm_input_vtable __vtable = {read_block};

struct nmm_input* nmm_input_create(char const* filepath)
{
    struct imm_input* super = imm_input_create(filepath);
    if (!super)
        return NULL;

    struct nmm_input* input = malloc(sizeof(*input));
    input->super = super;

    return input;
}

int nmm_input_destroy(struct nmm_input const* input)
{
    int err = imm_input_destroy(input->super);
    free_c(input);
    return err;
}

bool nmm_input_eof(struct nmm_input const* input) { return imm_input_eof(input->super); }

struct nmm_model const* nmm_input_read(struct nmm_input* input)
{
    struct imm_model const* model = imm_input_read(input->super);
    if (!model) {
        imm_error("could not read model");
        return NULL;
    }

    return __imm_model_derived(model);
}

#if 0
struct nmm_model const* nmm_model_create_from_file(FILE* stream)
{
    struct nmm_model* io = malloc(sizeof(*io));
    io->baset_ptrs = NULL;
    io->codonp_ptrs = NULL;
    io->codont_ptrs = NULL;
    io->baset_map = NULL;
    io->codonp_map = NULL;
    io->codont_map = NULL;
    io->super = __imm_io_new(io);
    *__imm_io_vtable(io->super) = __vtable;

    if (read_abc(io, stream)) {
        imm_error("could not read abc");
        goto err;
    }

    if (read_baset(io, stream)) {
        imm_error("could not read baset");
        goto err;
    }

    if (read_codonp(io, stream)) {
        imm_error("could not read codonp");
        goto err;
    }

    if (read_codont(io, stream)) {
        imm_error("could not read codont");
        goto err;
    }

    if (__imm_io_read_hmm(io->super, stream)) {
        imm_error("could not read hmm");
        goto err;
    }

    if (__imm_io_read_dp(io->super, stream)) {
        imm_error("could not read dp");
        goto err;
    }

    __imm_dp_create_from_io(io->super);
    return io;

err:
    __imm_io_vtable(io->super)->destroy_on_read_failure(io->super);
    return NULL;
}
#endif

struct imm_model const* read_block(struct imm_input* input, uint8_t block_type) {}

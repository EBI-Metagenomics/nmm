#include "nmm/io.h"
#include "amino_abc.h"
#include "base_abc.h"
#include "free.h"
#include "imm/imm.h"
#include "nmm/abc_types.h"

struct nmm_io
{
    struct imm_io const* parent;
};

static struct imm_abc const* read_abc(FILE* stream, uint8_t type_id);
static void                  destroy(struct imm_io const* io);

struct imm_io_vtable const __vtable = {read_abc, destroy};

struct nmm_io const* nmm_io_create(struct imm_hmm* hmm, struct imm_dp const* dp)
{

    struct nmm_io* io = malloc(sizeof(*io));

    io->parent = __imm_io_create_parent(hmm, dp, __vtable, io);
    if (!io->parent) {
        imm_error("could not io_create_parent");
        free_c(io);
        return NULL;
    }

    return io;
}

void nmm_io_destroy(struct nmm_io const* io)
{
    struct imm_io const* parent = io->parent;
    destroy(parent);
    __imm_io_destroy_parent(parent);
}

int nmm_io_write(struct nmm_io const* io, FILE* stream) { return imm_io_write(io->parent, stream); }

static struct imm_abc const* read_abc(FILE* stream, uint8_t type_id)
{
    struct imm_abc const* abc = NULL;

    switch (type_id) {
    case IMM_ABC_TYPE_ID:
        if (!(abc = imm_abc_read(stream)))
            imm_error("could not read abc");
        break;
    case NMM_AMINO_ABC_TYPE_ID:
        if (!(abc = amino_abc_read(stream)))
            imm_error("could not read amino_abc");
        break;
    case NMM_BASE_ABC_TYPE_ID:
        if (!(abc = base_abc_read(stream)))
            imm_error("could not read base_abc");
        break;
    default:
        imm_error("unknown abc type_id");
        return NULL;
    }

    return abc;
}

static void destroy(struct imm_io const* io) { free_c(__imm_io_child(io)); }

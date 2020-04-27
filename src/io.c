#include "nmm/io.h"
#include "amino_abc.h"
#include "base_abc.h"
#include "codon_state.h"
#include "free.h"
#include "imm/imm.h"
#include "lib/khash_ptr.h"
#include "nmm/abc_types.h"
#include "nmm/codon_state.h"
#include "nmm/state_types.h"

struct codonp_node
{
    uint32_t                      index;
    struct nmm_codon_lprob const* codonp;
};

KHASH_MAP_INIT_PTR(codonp, struct codonp_node*)

struct nmm_io
{
    struct imm_io const* parent;

    khash_t(codonp) * codonp_map;
};

static struct imm_abc const* read_abc(FILE* stream, uint8_t type_id);
static void                  destroy(struct imm_io const* io);
void                         create_codon_lprobs(struct nmm_io* io);

struct imm_io_vtable const __vtable = {read_abc, destroy};

struct nmm_io const* nmm_io_create(struct imm_hmm* hmm, struct imm_dp const* dp)
{

    struct nmm_io* io = malloc(sizeof(*io));
    io->codonp_map = NULL;

    io->parent = __imm_io_create_parent(hmm, dp, __vtable, io);
    if (!io->parent) {
        imm_error("could not io_create_parent");
        free_c(io);
        return NULL;
    }

    create_codon_lprobs(io);

    return io;
}

void create_codon_lprobs(struct nmm_io* io)
{
    khash_t(codonp)* map = io->codonp_map = kh_init(codonp);

    uint32_t idx = 0;
    for (uint32_t i = 0; i < imm_io_nstates(io->parent); ++i) {
        struct imm_state const* state = imm_io_state(io->parent, i);
        if (imm_state_type_id(state) == NMM_CODON_STATE_TYPE_ID) {

            struct nmm_codon_state const* s = nmm_codon_state_child(state);
            struct nmm_codon_lprob const* codonp = codon_state_codonp(s);
            khint_t                       k = kh_get(codonp, map, codonp);
            if (k != kh_end(map))
                continue;

            struct codonp_node* node = malloc(sizeof(*node));
            node->index = idx++;
            node->codonp = codonp;
            int      ret = 0;
            khiter_t iter = kh_put(codonp, map, node->codonp, &ret);
            IMM_BUG(ret == -1 || ret == 0);
            kh_key(map, iter) = node->codonp;
            kh_val(map, iter) = node;
        }
    }
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

static void destroy(struct imm_io const* io)
{
    struct nmm_io* this = __imm_io_child(io);
    if (this->codonp_map) {
        for (khint_t k = kh_begin(this->codonp_map); k < kh_end(this->codonp_map); ++k)
            if (kh_exist(this->codonp_map, k)) {
                free_c(kh_val(this->codonp_map, k));
                free_c(this);
            }
        kh_destroy(codonp, this->codonp_map);
    }
}

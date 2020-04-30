#include "io.h"
#include "amino_abc.h"
#include "base_abc.h"
#include "base_table.h"
#include "codon_lprob.h"
#include "codon_state.h"
#include "codon_table.h"
#include "frame_state.h"
#include "free.h"
#include "imm/imm.h"
#include "lib/khash_ptr.h"
#include "nmm/abc_types.h"
#include "nmm/codon_state.h"
#include "nmm/frame_state.h"
#include "nmm/io.h"
#include "nmm/state_types.h"

struct baset_node
{
    uint32_t                     index;
    struct nmm_base_table const* baset;
};
KHASH_MAP_INIT_PTR(baset, struct baset_node*)

struct codonp_node
{
    uint32_t                      index;
    struct nmm_codon_lprob const* codonp;
};
KHASH_MAP_INIT_PTR(codonp, struct codonp_node*)

struct codont_node
{
    uint32_t                      index;
    struct nmm_codon_table const* codont;
};
KHASH_MAP_INIT_PTR(codont, struct codont_node*)

struct nmm_io
{
    struct imm_io* super;

    struct nmm_base_table**  baset_ptrs;
    struct nmm_codon_lprob** codonp_ptrs;
    struct nmm_codon_table** codont_ptrs;

    khash_t(baset) * baset_map;
    khash_t(codonp) * codonp_map;
    khash_t(codont) * codont_map;
};

static void                  create_baset_map(struct nmm_io* io);
static void                  create_codonp_map(struct nmm_io* io);
static void                  create_codont_map(struct nmm_io* io);
static struct imm_abc const* read_abc(FILE* stream, uint8_t type_id);
static void                  destroy(struct imm_io const* io);
static void                  destroy_on_read_failure(struct imm_io const* io);
static int                   write(struct imm_io const* io, FILE* stream);
static int                   write_baset(struct nmm_io const* io, FILE* stream);
static int                   write_codonp(struct nmm_io const* io, FILE* stream);
static int                   write_codont(struct nmm_io const* io, FILE* stream);

struct imm_io_vtable const __vtable = {destroy, write, destroy_on_read_failure};

struct nmm_io const* nmm_io_create(struct imm_hmm* hmm, struct imm_dp const* dp)
{
    struct nmm_io* io = malloc(sizeof(*io));

    io->baset_map = NULL;
    io->codonp_map = NULL;
    io->codont_map = NULL;

    io->baset_ptrs = NULL;
    io->codonp_ptrs = NULL;
    io->codont_ptrs = NULL;

    io->super = __imm_io_create(hmm, dp, io);
    if (!io->super) {
        imm_error("could not __imm_io_create");
        free_c(io);
        return NULL;
    }

    create_baset_map(io);
    create_codonp_map(io);
    create_codont_map(io);

    *__imm_io_vtable(io->super) = __vtable;

    return io;
}

struct nmm_io const* nmm_io_create_from_file(FILE* stream)
{
    struct nmm_io* io = malloc(sizeof(*io));
    io->baset_ptrs = NULL;
    io->codonp_ptrs = NULL;
    io->codont_ptrs = NULL;
    io->super = __imm_io_new(io);
    /* TODO: finish this */

    uint32_t nbaset = 0;
    if (fread(&nbaset, sizeof(nbaset), 1, stream) < 1) {
        imm_error("could not read nbaset");
        goto err;
    }
    io->baset_ptrs = malloc(sizeof(*io->baset_ptrs) * nbaset);
    for (uint32_t i = 0; i < nbaset; ++i)
        io->baset_ptrs[i] = NULL;

    for (uint32_t i = 0; i < nbaset; ++i)
        io->baset_ptrs[i] = base_table_read(stream, NULL);

    __imm_io_read(io->super, stream);
    return io;

err:

    return NULL;
}

/* struct imm_io* __imm_io_new(void* derived) */

void nmm_io_destroy(struct nmm_io const* io) { __imm_io_vtable(io->super)->destroy(io->super); }

int nmm_io_write(struct nmm_io const* io, FILE* stream)
{
    return __imm_io_vtable(io->super)->write(io->super, stream);
}

struct nmm_io const* nmm_io_derived(struct imm_io const* io)
{
    /* TODO: add type info to check its validity? */
    return __imm_io_derived(io);
}

uint32_t io_baset_index(struct nmm_io const* io, struct nmm_base_table const* baset)
{
    khiter_t i = kh_get(baset, io->baset_map, baset);
    IMM_BUG(i == kh_end(io->baset_map));

    struct baset_node* node = kh_val(io->baset_map, i);

    return node->index;
}

uint32_t io_codonp_index(struct nmm_io const* io, struct nmm_codon_lprob const* codonp)
{
    khiter_t i = kh_get(codonp, io->codonp_map, codonp);
    IMM_BUG(i == kh_end(io->codonp_map));

    struct codonp_node* node = kh_val(io->codonp_map, i);

    return node->index;
}

uint32_t io_codont_index(struct nmm_io const* io, struct nmm_codon_table const* codont)
{
    khiter_t i = kh_get(codont, io->codont_map, codont);
    IMM_BUG(i == kh_end(io->codont_map));

    struct codont_node* node = kh_val(io->codont_map, i);

    return node->index;
}

static void create_baset_map(struct nmm_io* io)
{
    khash_t(baset)* map = io->baset_map = kh_init(baset);

    uint32_t idx = 0;
    for (uint32_t i = 0; i < imm_io_nstates(io->super); ++i) {
        struct imm_state const* state = imm_io_state(io->super, i);
        if (imm_state_type_id(state) != NMM_FRAME_STATE_TYPE_ID)
            continue;

        struct nmm_frame_state const* s = nmm_frame_state_derived(state);
        struct nmm_base_table const*  baset = frame_state_baset(s);
        khiter_t                      k = kh_get(baset, map, baset);
        if (k != kh_end(map))
            continue;

        struct baset_node* node = malloc(sizeof(*node));
        node->index = idx++;
        node->baset = baset;
        int      ret = 0;
        khiter_t iter = kh_put(baset, map, node->baset, &ret);
        IMM_BUG(ret == -1 || ret == 0);
        kh_key(map, iter) = node->baset;
        kh_val(map, iter) = node;
    }
}

static void create_codonp_map(struct nmm_io* io)
{
    khash_t(codonp)* map = io->codonp_map = kh_init(codonp);

    uint32_t idx = 0;
    for (uint32_t i = 0; i < imm_io_nstates(io->super); ++i) {
        struct imm_state const* state = imm_io_state(io->super, i);
        if (imm_state_type_id(state) != NMM_CODON_STATE_TYPE_ID)
            continue;

        struct nmm_codon_state const* s = nmm_codon_state_derived(state);
        struct nmm_codon_lprob const* codonp = codon_state_codonp(s);
        khiter_t                      k = kh_get(codonp, map, codonp);
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

static void create_codont_map(struct nmm_io* io)
{
    khash_t(codont)* map = io->codont_map = kh_init(codont);

    uint32_t idx = 0;
    for (uint32_t i = 0; i < imm_io_nstates(io->super); ++i) {
        struct imm_state const* state = imm_io_state(io->super, i);
        if (imm_state_type_id(state) != NMM_FRAME_STATE_TYPE_ID)
            continue;

        struct nmm_frame_state const* s = nmm_frame_state_derived(state);
        struct nmm_codon_table const* codont = frame_state_codont(s);
        khiter_t                      k = kh_get(codont, map, codont);
        if (k != kh_end(map))
            continue;

        struct codont_node* node = malloc(sizeof(*node));
        node->index = idx++;
        node->codont = codont;
        int      ret = 0;
        khiter_t iter = kh_put(codont, map, node->codont, &ret);
        IMM_BUG(ret == -1 || ret == 0);
        kh_key(map, iter) = node->codont;
        kh_val(map, iter) = node;
    }
}

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
    struct nmm_io* this = __imm_io_derived(io);
    if (this->codonp_map) {
        for (khint_t k = kh_begin(this->codonp_map); k < kh_end(this->codonp_map); ++k)
            if (kh_exist(this->codonp_map, k)) {
                free_c(kh_val(this->codonp_map, k));
                free_c(this);
            }
        kh_destroy(codonp, this->codonp_map);
    }

    __imm_io_destroy(io);
}

static void destroy_on_read_failure(struct imm_io const* io)
{
    __imm_io_destroy_on_read_failure(io);
}

static int write(struct imm_io const* io, FILE* stream)
{
    if (__imm_io_write_abc(io, stream)) {
        imm_error("could not write abc");
        return 1;
    }

    if (write_baset(nmm_io_derived(io), stream)) {
        imm_error("could not write_baset");
        return 1;
    }

    if (write_codonp(nmm_io_derived(io), stream)) {
        imm_error("could not write_codonp");
        return 1;
    }

    if (write_codont(nmm_io_derived(io), stream)) {
        imm_error("could not write_codont");
        return 1;
    }

    if (__imm_io_write_hmm(io, stream)) {
        imm_error("could not write hmm");
        return 1;
    }

    if (__imm_io_write_dp(io, stream)) {
        imm_error("could not write dp");
        return 1;
    }

    return 0;
}

static int write_baset(struct nmm_io const* io, FILE* stream)
{
    khash_t(baset)* map = io->baset_map;
    IMM_BUG(kh_size(map) > UINT32_MAX);
    uint32_t n = (uint32_t)kh_size(map);

    if (fwrite(&n, sizeof(n), 1, stream) < 1) {
        imm_error("could not write nbaset");
        return 1;
    }

    for (khiter_t i = kh_begin(map); i < kh_end(map); ++i) {
        if (!kh_exist(map, i))
            continue;

        struct baset_node const* node = kh_val(map, i);

        if (base_table_write(node->baset, stream))
            return 1;
    }

    return 0;
}

static int write_codonp(struct nmm_io const* io, FILE* stream)
{
    khash_t(codonp)* map = io->codonp_map;
    IMM_BUG(kh_size(map) > UINT32_MAX);
    uint32_t n = (uint32_t)kh_size(map);

    if (fwrite(&n, sizeof(n), 1, stream) < 1) {
        imm_error("could not write ncodonp");
        return 1;
    }

    for (khiter_t i = kh_begin(map); i < kh_end(map); ++i) {
        if (!kh_exist(map, i))
            continue;

        struct codonp_node const* node = kh_val(map, i);

        if (codon_lprob_write(node->codonp, stream))
            return 1;
    }

    return 0;
}

static int write_codont(struct nmm_io const* io, FILE* stream)
{
    khash_t(codont)* map = io->codont_map;
    IMM_BUG(kh_size(map) > UINT32_MAX);
    uint32_t n = (uint32_t)kh_size(map);

    if (fwrite(&n, sizeof(n), 1, stream) < 1) {
        imm_error("could not write ncodont");
        return 1;
    }

    for (khiter_t i = kh_begin(map); i < kh_end(map); ++i) {
        if (!kh_exist(map, i))
            continue;

        struct codont_node const* node = kh_val(map, i);

        if (codon_table_write(node->codont, stream))
            return 1;
    }

    return 0;
}

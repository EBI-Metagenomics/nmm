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
#include "nmm/base_abc.h"
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

    uint32_t                       nbaset;
    struct nmm_base_table const**  baset_ptrs;
    uint32_t                       ncodonp;
    struct nmm_codon_lprob const** codonp_ptrs;
    uint32_t                       ncodont;
    struct nmm_codon_table const** codont_ptrs;

    khash_t(baset) * baset_map;
    khash_t(codonp) * codonp_map;
    khash_t(codont) * codont_map;
};

static void                    create_baset_map(struct nmm_io* io);
static void                    create_codonp_map(struct nmm_io* io);
static void                    create_codont_map(struct nmm_io* io);
static void                    destroy(struct imm_io const* io);
static void                    destroy_baset_map(khash_t(baset) * baset_map);
static void                    destroy_codonp_map(khash_t(codonp) * codonp_map);
static void                    destroy_codont_map(khash_t(codont) * codont_map);
static void                    destroy_on_read_failure(struct imm_io const* io);
static int                     read_abc(struct nmm_io* io, FILE* stream);
static int                     read_baset(struct nmm_io* io, FILE* stream);
static int                     read_codonp(struct nmm_io* io, FILE* stream);
static int                     read_codont(struct nmm_io* io, FILE* stream);
static struct imm_state const* read_state(struct imm_io const* io, FILE* stream, uint8_t type_id);
static int                     write(struct imm_io const* io, FILE* stream);
static int                     write_baset(struct nmm_io const* io, FILE* stream);
static int                     write_codonp(struct nmm_io const* io, FILE* stream);
static int                     write_codont(struct nmm_io const* io, FILE* stream);

struct imm_io_vtable const __vtable = {destroy, read_state, write, destroy_on_read_failure};

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
    /* TODO: fix memory leaks */
    struct nmm_io* io = malloc(sizeof(*io));
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
    return NULL;
}

static int read_baset(struct nmm_io* io, FILE* stream)
{
    io->nbaset = 0;
    if (fread(&io->nbaset, sizeof(io->nbaset), 1, stream) < 1) {
        imm_error("could not read nbaset");
        goto err;
    }
    if (io->nbaset == 0) {
        io->baset_ptrs = NULL;
        return 0;
    }

    io->baset_ptrs = malloc(sizeof(*io->baset_ptrs) * io->nbaset);
    for (uint32_t i = 0; i < io->nbaset; ++i)
        io->baset_ptrs[i] = NULL;

    for (uint32_t i = 0; i < io->nbaset; ++i) {
        struct nmm_base_abc const* base_abc = nmm_base_abc_derived(imm_io_abc(io->super));
        if (!base_abc)
            goto err;

        io->baset_ptrs[i] = base_table_read(stream, base_abc);
    }

    return 0;

err:
    /* TODO: fix memory leak on error */
    return 1;
}

static int read_codonp(struct nmm_io* io, FILE* stream)
{
    io->ncodonp = 0;
    if (fread(&io->ncodonp, sizeof(io->ncodonp), 1, stream) < 1) {
        imm_error("could not read ncodonp");
        goto err;
    }
    if (io->ncodonp == 0) {
        io->codonp_ptrs = NULL;
        return 0;
    }

    io->codonp_ptrs = malloc(sizeof(*io->codonp_ptrs) * io->ncodonp);
    for (uint32_t i = 0; i < io->ncodonp; ++i)
        io->codonp_ptrs[i] = NULL;

    for (uint32_t i = 0; i < io->ncodonp; ++i) {
        struct nmm_base_abc const* base_abc = nmm_base_abc_derived(imm_io_abc(io->super));
        if (!base_abc)
            goto err;

        io->codonp_ptrs[i] = codon_lprob_read(stream, base_abc);
    }

    return 0;

err:
    /* TODO: fix memory leak */
    return 1;
}

static int read_codont(struct nmm_io* io, FILE* stream)
{
    io->ncodont = 0;
    if (fread(&io->ncodont, sizeof(io->ncodont), 1, stream) < 1) {
        imm_error("could not read ncodont");
        goto err;
    }
    if (io->ncodont == 0) {
        io->codont_ptrs = NULL;
        return 0;
    }

    io->codont_ptrs = malloc(sizeof(*io->codont_ptrs) * io->ncodont);
    for (uint32_t i = 0; i < io->ncodont; ++i)
        io->codont_ptrs[i] = NULL;

    for (uint32_t i = 0; i < io->ncodont; ++i) {
        struct nmm_base_abc const* base_abc = nmm_base_abc_derived(imm_io_abc(io->super));
        if (!base_abc)
            goto err;

        io->codont_ptrs[i] = codon_table_read(stream, base_abc);
    }

    return 0;

err:
    /* TODO: fix memory leak */
    return 1;
}

static struct imm_state const* read_state(struct imm_io const* io, FILE* stream, uint8_t type_id)
{
    if (type_id == IMM_MUTE_STATE_TYPE_ID || type_id == IMM_NORMAL_STATE_TYPE_ID ||
        type_id == IMM_TABLE_STATE_TYPE_ID)
        return __imm_io_read_state(io, stream, type_id);

    struct imm_state const* state = NULL;

    switch (type_id) {
    case NMM_CODON_STATE_TYPE_ID:
        if (!(state = codon_state_read(stream, nmm_io_derived(io))))
            imm_error("could not read codon_state");
        break;
    case NMM_FRAME_STATE_TYPE_ID:
        if (!(state = frame_state_read(stream, nmm_io_derived(io))))
            imm_error("could not read frame_state");
        break;
    default:
        imm_error("unknown state type_id");
    }

    return state;
}

void nmm_io_destroy(struct nmm_io const* io) { __imm_io_vtable(io->super)->destroy(io->super); }

struct imm_io const* nmm_io_super(struct nmm_io const* io) { return io->super; }

int nmm_io_write(struct nmm_io const* io, FILE* stream)
{
    return __imm_io_vtable(io->super)->write(io->super, stream);
}

struct nmm_io const* nmm_io_derived(struct imm_io const* io)
{
    /* TODO: add type info to check its validity? */
    return __imm_io_derived(io);
}

uint32_t nmm_io_nbase_tables(struct nmm_io const* io) { return io->nbaset; }

struct nmm_base_table const* nmm_io_base_table(struct nmm_io const* io, uint32_t index)
{
    return io->baset_ptrs[index];
}

uint32_t nmm_io_ncodon_tables(struct nmm_io const* io) { return io->ncodont; }

struct nmm_codon_table const* nmm_io_codon_table(struct nmm_io const* io, uint32_t index)
{
    return io->codont_ptrs[index];
}

uint32_t nmm_io_ncodon_lprobs(struct nmm_io const* io) { return io->ncodonp; }

struct nmm_codon_lprob const* nmm_io_codon_lprob(struct nmm_io const* io, uint32_t index)
{
    return io->codonp_ptrs[index];
}

uint32_t io_baset_index(struct nmm_io const* io, struct nmm_base_table const* baset)
{
    khiter_t i = kh_get(baset, io->baset_map, baset);
    IMM_BUG(i == kh_end(io->baset_map));

    struct baset_node* node = kh_val(io->baset_map, i);

    return node->index;
}

struct nmm_base_table const* io_get_baset(struct nmm_io const* io, uint32_t index)
{
    if (index > io->nbaset) {
        imm_error("baset index overflow");
        return NULL;
    }

    return io->baset_ptrs[index];
}

uint32_t io_codonp_index(struct nmm_io const* io, struct nmm_codon_lprob const* codonp)
{
    khiter_t i = kh_get(codonp, io->codonp_map, codonp);
    IMM_BUG(i == kh_end(io->codonp_map));

    struct codonp_node* node = kh_val(io->codonp_map, i);

    return node->index;
}

struct nmm_codon_lprob const* io_get_codonp(struct nmm_io const* io, uint32_t index)
{
    if (index > io->ncodonp) {
        imm_error("codonp index overflow");
        return NULL;
    }

    return io->codonp_ptrs[index];
}

uint32_t io_codont_index(struct nmm_io const* io, struct nmm_codon_table const* codont)
{
    khiter_t i = kh_get(codont, io->codont_map, codont);
    IMM_BUG(i == kh_end(io->codont_map));

    struct codont_node* node = kh_val(io->codont_map, i);

    return node->index;
}

struct nmm_codon_table const* io_get_codont(struct nmm_io const* io, uint32_t index)
{
    if (index > io->ncodont) {
        imm_error("codont index overflow");
        return NULL;
    }

    return io->codont_ptrs[index];
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

static int read_abc(struct nmm_io* io, FILE* stream)
{
    uint8_t type_id = 0;

    if (fread(&type_id, sizeof(type_id), 1, stream) < 1) {
        imm_error("could not read type_id");
        return 1;
    }

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
        return 1;
    }

    __imm_io_set_abc(io->super, abc);

    return 0;
}

static void destroy(struct imm_io const* io)
{
    struct nmm_io* this = __imm_io_derived(io);
    if (this->baset_map)
        destroy_baset_map(this->baset_map);

    if (this->codonp_map)
        destroy_codonp_map(this->codonp_map);

    if (this->codont_map)
        destroy_codont_map(this->codont_map);

    if (this->baset_ptrs)
        free_c(this->baset_ptrs);

    if (this->codont_ptrs)
        free_c(this->codont_ptrs);

    if (this->codonp_ptrs)
        free_c(this->codonp_ptrs);

    free_c(this);
    __imm_io_destroy(io);
}

static void destroy_baset_map(khash_t(baset) * baset_map)
{
    for (khint_t k = kh_begin(baset_map); k < kh_end(baset_map); ++k)
        if (kh_exist(baset_map, k)) {
            free_c(kh_val(baset_map, k));
        }
    kh_destroy(baset, baset_map);
}

static void destroy_codonp_map(khash_t(codonp) * codonp_map)
{
    for (khint_t k = kh_begin(codonp_map); k < kh_end(codonp_map); ++k)
        if (kh_exist(codonp_map, k)) {
            free_c(kh_val(codonp_map, k));
        }
    kh_destroy(codonp, codonp_map);
}

static void destroy_codont_map(khash_t(codont) * codont_map)
{
    for (khint_t k = kh_begin(codont_map); k < kh_end(codont_map); ++k)
        if (kh_exist(codont_map, k)) {
            free_c(kh_val(codont_map, k));
        }
    kh_destroy(codont, codont_map);
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

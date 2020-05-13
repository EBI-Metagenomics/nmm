#include "model.h"
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
#include "nmm/base_table.h"
#include "nmm/codon_lprob.h"
#include "nmm/codon_state.h"
#include "nmm/codon_table.h"
#include "nmm/frame_state.h"
#include "nmm/model.h"
#include "nmm/state_types.h"
#include <imm/model.h>

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

struct nmm_model
{
    struct imm_model* super;

    khash_t(baset) * baset_map;
    khash_t(codonp) * codonp_map;
    khash_t(codont) * codont_map;
};

static void                    create_baset_map(struct nmm_model* model);
static void                    create_codonp_map(struct nmm_model* model);
static void                    create_codont_map(struct nmm_model* model);
static void                    destroy(struct imm_model const* model);
static void                    destroy_baset_map(khash_t(baset) * baset_map);
static void                    destroy_codonp_map(khash_t(codonp) * codonp_map);
static void                    destroy_codont_map(khash_t(codont) * codont_map);
static void                    deep_destroy(struct imm_model const* model);
static int                     read(struct imm_model* model, FILE* stream);
static struct imm_abc const*   read_abc(struct imm_model* model, FILE* stream, uint8_t type_id);
static int                     read_baset(struct nmm_model* model, FILE* stream);
static int                     read_codonp(struct nmm_model* model, FILE* stream);
static int                     read_codont(struct nmm_model* model, FILE* stream);
static struct imm_state const* read_state(struct imm_model const* model, FILE* stream,
                                          uint8_t type_id);
static int                     write(struct imm_model const* model, FILE* stream);
static int                     write_abc(struct imm_model const* model, FILE* stream);
static int                     write_baset(struct nmm_model const* model, FILE* stream);
static int                     write_codonp(struct nmm_model const* model, FILE* stream);
static int                     write_codont(struct nmm_model const* model, FILE* stream);

static struct imm_model_vtable const __vtable = {deep_destroy, destroy, read,     read_abc,
                                                 read_state,   write,   write_abc};

struct nmm_model const* nmm_model_create(struct imm_hmm* hmm, struct imm_dp const* dp)
{
    struct nmm_model* model = model_new();

    if (!(model->super = __imm_model_create(hmm, dp, __vtable, model))) {
        imm_error("could not create model");
        destroy_baset_map(model->baset_map);
        destroy_codonp_map(model->codonp_map);
        destroy_codont_map(model->codont_map);
        free_c(model);
        return NULL;
    }

    create_baset_map(model);
    create_codonp_map(model);
    create_codont_map(model);

    return model;
}

struct nmm_model* nmm_model_derived(struct imm_model* model) { return __imm_model_derived(model); }

struct nmm_model const* nmm_model_derived_c(struct imm_model const* model)
{
    return __imm_model_derived_c(model);
}

void nmm_model_destroy(struct nmm_model const* model) { imm_model_destroy(model->super); }

uint32_t nmm_model_nbase_tables(struct nmm_model const* model) { return kh_size(model->baset_map); }

uint32_t nmm_model_ncodon_lprobs(struct nmm_model const* model)
{
    return kh_size(model->codonp_map);
}

uint32_t nmm_model_ncodon_tables(struct nmm_model const* model)
{
    return kh_size(model->codont_map);
}

struct imm_model const* nmm_model_super(struct nmm_model const* model) { return model->super; }

int nmm_model_write(struct nmm_model const* model, FILE* stream)
{
    return imm_model_write(model->super, stream);
}

uint32_t model_baset_index(struct nmm_model const* model, struct nmm_base_table const* baset)
{
    khiter_t i = kh_get(baset, model->baset_map, baset);
    IMM_BUG(i == kh_end(model->baset_map));

    struct baset_node* node = kh_val(model->baset_map, i);

    return node->index;
}

uint32_t model_codonp_index(struct nmm_model const* model, struct nmm_codon_lprob const* codonp)
{
    khiter_t i = kh_get(codonp, model->codonp_map, codonp);
    IMM_BUG(i == kh_end(model->codonp_map));

    struct codonp_node* node = kh_val(model->codonp_map, i);

    return node->index;
}

uint32_t model_codont_index(struct nmm_model const* model, struct nmm_codon_table const* codont)
{
    khiter_t i = kh_get(codont, model->codont_map, codont);
    IMM_BUG(i == kh_end(model->codont_map));

    struct codont_node* node = kh_val(model->codont_map, i);

    return node->index;
}

struct nmm_model* model_new(void)
{
    struct nmm_model* model = malloc(sizeof(*model));
    model->super = NULL;
    model->baset_map = kh_init(baset);
    model->codonp_map = kh_init(codonp);
    model->codont_map = kh_init(codont);
    return model;
}

struct nmm_base_table const* nmm_model_base_table(struct nmm_model const* model, uint32_t index)
{
    for (khint_t k = kh_begin(model->baset_map); k < kh_end(model->baset_map); ++k)
        if (kh_exist(model->baset_map, k)) {
            struct baset_node const* node = kh_val(model->baset_map, k);
            if (node->index == index)
                return node->baset;
        }

    imm_error("baset index not found");
    return NULL;
}

struct nmm_codon_lprob const* nmm_model_codon_lprob(struct nmm_model const* model, uint32_t index)
{
    for (khint_t k = kh_begin(model->codonp_map); k < kh_end(model->codonp_map); ++k)
        if (kh_exist(model->codonp_map, k)) {
            struct codonp_node const* node = kh_val(model->codonp_map, k);
            if (node->index == index)
                return node->codonp;
        }

    imm_error("codonp index not found");
    return NULL;
}

struct nmm_codon_table const* nmm_model_codon_table(struct nmm_model const* model, uint32_t index)
{
    for (khint_t k = kh_begin(model->codont_map); k < kh_end(model->codont_map); ++k)
        if (kh_exist(model->codont_map, k)) {
            struct codont_node const* node = kh_val(model->codont_map, k);
            if (node->index == index)
                return node->codont;
        }

    imm_error("codont index not found");
    return NULL;
}

static void create_baset_map(struct nmm_model* model)
{
    khash_t(baset)* map = model->baset_map;

    uint32_t idx = 0;
    for (uint32_t i = 0; i < imm_model_nstates(model->super); ++i) {
        struct imm_state const* state = imm_model_state(model->super, i);
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

static void create_codonp_map(struct nmm_model* model)
{
    khash_t(codonp)* map = model->codonp_map;

    uint32_t idx = 0;
    for (uint32_t i = 0; i < imm_model_nstates(model->super); ++i) {
        struct imm_state const* state = imm_model_state(model->super, i);
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

static void create_codont_map(struct nmm_model* model)
{
    khash_t(codont)* map = model->codont_map;

    uint32_t idx = 0;
    for (uint32_t i = 0; i < imm_model_nstates(model->super); ++i) {
        struct imm_state const* state = imm_model_state(model->super, i);
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

static void destroy(struct imm_model const* model)
{
    struct nmm_model const* this = nmm_model_derived_c(model);
    imm_model_vtable.destroy(model);
    destroy_baset_map(this->baset_map);
    destroy_codonp_map(this->codonp_map);
    destroy_codont_map(this->codont_map);
    free_c(this);
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

static void deep_destroy(struct imm_model const* model)
{
    struct nmm_model const* this = __imm_model_derived_c(model);
    imm_model_vtable.deep_destroy(model);

    for (khint_t k = kh_begin(this->baset_map); k < kh_end(this->baset_map); ++k) {
        if (kh_exist(this->baset_map, k)) {
            struct baset_node const* node = kh_val(this->baset_map, k);
            nmm_base_table_destroy(node->baset);
            free_c(node);
        }
    }
    kh_destroy(baset, this->baset_map);

    for (khint_t k = kh_begin(this->codonp_map); k < kh_end(this->codonp_map); ++k) {
        if (kh_exist(this->codonp_map, k)) {
            struct codonp_node const* node = kh_val(this->codonp_map, k);
            nmm_codon_lprob_destroy(node->codonp);
            free_c(node);
        }
    }
    kh_destroy(codonp, this->codonp_map);

    for (khint_t k = kh_begin(this->codont_map); k < kh_end(this->codont_map); ++k) {
        if (kh_exist(this->codont_map, k)) {
            struct codont_node const* node = kh_val(this->codont_map, k);
            nmm_codon_table_destroy(node->codont);
            free_c(node);
        }
    }
    kh_destroy(codont, this->codont_map);

    free_c(this);
}

static int read(struct imm_model* model, FILE* stream)
{
    uint8_t abc_type_id = 0;
    if (fread(&abc_type_id, sizeof(abc_type_id), 1, stream) < 1) {
        imm_error("could not read abc type id");
        goto err;
    }

    struct imm_abc const* abc = __imm_model_vtable(model)->read_abc(model, stream, abc_type_id);
    if (!abc) {
        imm_error("could not read abc");
        goto err;
    }
    __imm_model_set_abc(model, abc);

    if (read_baset(nmm_model_derived(model), stream)) {
        imm_error("could not read baset");
        goto err;
    }

    if (read_codonp(nmm_model_derived(model), stream)) {
        imm_error("could not read codonp");
        goto err;
    }

    if (read_codont(nmm_model_derived(model), stream)) {
        imm_error("could not read codont");
        goto err;
    }

err:
    __imm_model_vtable(model)->deep_destroy(model);
    return 1;
}

static struct imm_abc const* read_abc(struct imm_model* model, FILE* stream, uint8_t type_id)
{
    struct imm_abc const* abc = NULL;

    switch (type_id) {
    case NMM_AMINO_ABC_TYPE_ID:
        if (!(abc = amino_abc_read(stream)))
            imm_error("could not read amino_abc");
        break;
    case NMM_BASE_ABC_TYPE_ID:
        if (!(abc = base_abc_read(stream)))
            imm_error("could not read base_abc");
        break;
    default:
        abc = imm_model_vtable.read_abc(model, stream, type_id);
    }

    return abc;
}

static int read_baset(struct nmm_model* model, FILE* stream)
{
    uint32_t nbaset = 0;
    if (fread(&nbaset, sizeof(nbaset), 1, stream) < 1) {
        imm_error("could not read the number of base tables");
        return 1;
    }

    for (uint32_t i = 0; i < nbaset; ++i) {
        struct nmm_base_abc const* base_abc = nmm_base_abc_derived(imm_model_abc(model->super));
        if (!base_abc)
            return 1;

        struct nmm_base_table const* baset = base_table_read(stream, base_abc);
        if (!baset)
            return 1;

        struct baset_node* node = malloc(sizeof(*node));
        node->index = i;
        node->baset = baset;

        int      ret = 0;
        khiter_t iter = kh_put(baset, model->baset_map, node->baset, &ret);
        IMM_BUG(ret == -1 || ret == 0);
        kh_key(model->baset_map, iter) = node->baset;
        kh_val(model->baset_map, iter) = node;
    }

    return 0;
}

static int read_codonp(struct nmm_model* model, FILE* stream)
{
    uint32_t ncodonp = 0;
    if (fread(&ncodonp, sizeof(ncodonp), 1, stream) < 1) {
        imm_error("could not read the number of codon lprobs");
        return 1;
    }

    for (uint32_t i = 0; i < ncodonp; ++i) {
        struct nmm_base_abc const* base_abc = nmm_base_abc_derived(imm_model_abc(model->super));
        if (!base_abc)
            return 1;

        struct nmm_codon_lprob const* codonp = codon_lprob_read(stream, base_abc);
        if (!codonp)
            return 1;

        struct codonp_node* node = malloc(sizeof(*node));
        node->index = i;
        node->codonp = codonp;

        int      ret = 0;
        khiter_t iter = kh_put(codonp, model->codonp_map, node->codonp, &ret);
        IMM_BUG(ret == -1 || ret == 0);
        kh_key(model->codonp_map, iter) = node->codonp;
        kh_val(model->codonp_map, iter) = node;
    }

    return 0;
}

static int read_codont(struct nmm_model* model, FILE* stream)
{
    uint32_t ncodont = 0;
    if (fread(&ncodont, sizeof(ncodont), 1, stream) < 1) {
        imm_error("could not read the number of codon tables");
        return 1;
    }

    for (uint32_t i = 0; i < ncodont; ++i) {
        struct nmm_base_abc const* base_abc = nmm_base_abc_derived(imm_model_abc(model->super));
        if (!base_abc)
            return 1;

        struct nmm_codon_table const* codont = codon_table_read(stream, base_abc);
        if (!codont)
            return 1;

        struct codont_node* node = malloc(sizeof(*node));
        node->index = i;
        node->codont = codont;

        int      ret = 0;
        khiter_t iter = kh_put(codont, model->codont_map, node->codont, &ret);
        IMM_BUG(ret == -1 || ret == 0);
        kh_key(model->codont_map, iter) = node->codont;
        kh_val(model->codont_map, iter) = node;
    }

    return 0;
}

static struct imm_state const* read_state(struct imm_model const* model, FILE* stream,
                                          uint8_t type_id)
{
    struct imm_state const* state = NULL;

    switch (type_id) {
    case NMM_CODON_STATE_TYPE_ID:
        if (!(state = codon_state_read(stream, nmm_model_derived_c(model))))
            imm_error("could not read codon_state");
        break;
    case NMM_FRAME_STATE_TYPE_ID:
        if (!(state = frame_state_read(stream, nmm_model_derived_c(model))))
            imm_error("could not read frame_state");
        break;
    default:
        state = imm_model_vtable.read_state(model, stream, type_id);
    }

    return state;
}

static int write(struct imm_model const* model, FILE* stream)
{
    if (__imm_model_vtable(model)->write_abc(model, stream)) {
        imm_error("could not write abc");
        return 1;
    }

    if (write_baset(nmm_model_derived_c(model), stream)) {
        imm_error("could not write_baset");
        return 1;
    }

    if (write_codonp(nmm_model_derived_c(model), stream)) {
        imm_error("could not write_codonp");
        return 1;
    }

    if (write_codont(nmm_model_derived_c(model), stream)) {
        imm_error("could not write_codont");
        return 1;
    }

    if (__imm_model_write_hmm(model, stream)) {
        imm_error("could not write hmm");
        return 1;
    }

    if (__imm_model_write_dp(model, stream)) {
        imm_error("could not write dp");
        return 1;
    }

    return 0;
}

static int write_abc(struct imm_model const* model, FILE* stream)
{
    return imm_model_vtable.write_abc(model, stream);
}

static int write_baset(struct nmm_model const* model, FILE* stream)
{
    khash_t(baset)* map = model->baset_map;
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

static int write_codonp(struct nmm_model const* model, FILE* stream)
{
    khash_t(codonp)* map = model->codonp_map;
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

static int write_codont(struct nmm_model const* model, FILE* stream)
{
    khash_t(codont)* map = model->codont_map;
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

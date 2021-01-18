#include "model.h"
#include "amino_abc.h"
#include "base_abc.h"
#include "base_table.h"
#include "codon_lprob.h"
#include "codon_state.h"
#include "codon_table.h"
#include "free.h"
#include "imm/imm.h"
#include "lib/khash.h"
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
KHASH_MAP_INIT_INT64(baset_idx, struct baset_node*)

struct codonp_node
{
    uint32_t                      index;
    struct nmm_codon_lprob const* codonp;
};
KHASH_MAP_INIT_PTR(codonp, struct codonp_node*)
KHASH_MAP_INIT_INT64(codonp_idx, struct codonp_node*)

struct codont_node
{
    uint32_t                      index;
    struct nmm_codon_table const* codont;
};
KHASH_MAP_INIT_PTR(codont, struct codont_node*)
KHASH_MAP_INIT_INT64(codont_idx, struct codont_node*)

struct nmm_model
{
    struct imm_model* super;

    khash_t(baset) * baset_map;
    khash_t(codonp) * codonp_map;
    khash_t(codont) * codont_map;

    khash_t(baset_idx) * baset_idx;
    khash_t(codonp_idx) * codonp_idx;
    khash_t(codont_idx) * codont_idx;
};

static void                    create_baset_map(struct nmm_model* model);
static void                    create_codonp_map(struct nmm_model* model);
static void                    create_codont_map(struct nmm_model* model);
static void                    deep_destroy(struct nmm_model const* model);
static void                    destroy_baset_map(khash_t(baset) * baset_map);
static void                    destroy_codonp_map(khash_t(codonp) * codonp_map);
static void                    destroy_codont_map(khash_t(codont) * codont_map);
struct nmm_model*              model_new(void);
static int                     read(struct imm_model* model, FILE* stream);
static struct imm_abc const*   read_abc(FILE* stream, uint8_t type_id);
static int                     read_baset(struct nmm_model* model, FILE* stream);
static int                     read_codonp(struct nmm_model* model, FILE* stream);
static int                     read_codont(struct nmm_model* model, FILE* stream);
static struct imm_state const* read_state(struct imm_model const* model, FILE* stream, void* args);
static int                     write(struct imm_model const* model, FILE* stream);
static int                     write_abc(struct nmm_model const* model, FILE* stream);
static int                     write_baset(struct nmm_model const* model, FILE* stream);
static int                     write_codonp(struct nmm_model const* model, FILE* stream);
static int                     write_codont(struct nmm_model const* model, FILE* stream);
static int write_state(struct imm_model const* model, FILE* stream, struct imm_state const* state,
                       void* args);

struct imm_abc const* nmm_model_abc(struct nmm_model const* model)
{
    return imm_model_abc(model->super);
}

struct nmm_model const* nmm_model_create(struct imm_hmm* hmm, struct imm_dp const* dp)
{
    struct nmm_model* model = malloc(sizeof(*model));

    model->super = __imm_model_create(hmm, dp, read_state, model, write_state, model);
    if (!model->super) {
        free_c(model);
        return NULL;
    }

    model->baset_map = kh_init(baset);
    model->codonp_map = kh_init(codonp);
    model->codont_map = kh_init(codont);

    model->baset_idx = kh_init(baset_idx);
    model->codonp_idx = kh_init(codonp_idx);
    model->codont_idx = kh_init(codont_idx);

    create_baset_map(model);
    create_codonp_map(model);
    create_codont_map(model);

    return model;
}

void nmm_model_destroy(struct nmm_model const* model)
{
    imm_model_destroy(model->super);
    destroy_baset_map(model->baset_map);
    destroy_codonp_map(model->codonp_map);
    destroy_codont_map(model->codont_map);

    kh_destroy(baset_idx, model->baset_idx);
    kh_destroy(codonp_idx, model->codonp_idx);
    kh_destroy(codont_idx, model->codont_idx);
    free_c(model);
}

uint32_t nmm_model_nbase_tables(struct nmm_model const* model) { return kh_size(model->baset_map); }

uint32_t nmm_model_ncodon_lprobs(struct nmm_model const* model)
{
    return kh_size(model->codonp_map);
}

uint32_t nmm_model_ncodon_tables(struct nmm_model const* model)
{
    return kh_size(model->codont_map);
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
    model->super = __imm_model_new(read_state, model, write_state, model);

    model->baset_map = kh_init(baset);
    model->codonp_map = kh_init(codonp);
    model->codont_map = kh_init(codont);

    model->baset_idx = kh_init(baset_idx);
    model->codonp_idx = kh_init(codonp_idx);
    model->codont_idx = kh_init(codont_idx);

    return model;
}

struct nmm_base_table const* nmm_model_base_table(struct nmm_model const* model, uint32_t index)
{
    khiter_t k = kh_get(baset_idx, model->baset_idx, index);
    if (k == kh_end(model->baset_idx)) {
        imm_error("baset index not found");
        return NULL;
    }
    return kh_val(model->baset_idx, k)->baset;
}

struct nmm_codon_lprob const* nmm_model_codon_lprob(struct nmm_model const* model, uint32_t index)
{
    khiter_t k = kh_get(codonp_idx, model->codonp_idx, index);
    if (k == kh_end(model->codonp_idx)) {
        imm_error("codonp index not found");
        return NULL;
    }
    return kh_val(model->codonp_idx, k)->codonp;
}

struct nmm_codon_table const* nmm_model_codon_table(struct nmm_model const* model, uint32_t index)
{
    khiter_t k = kh_get(codont_idx, model->codont_idx, index);
    if (k == kh_end(model->codont_idx)) {
        imm_error("codont index not found");
        return NULL;
    }
    return kh_val(model->codont_idx, k)->codont;
}

static void create_baset_map(struct nmm_model* model)
{
    khash_t(baset)* map = model->baset_map;
    khash_t(baset_idx)* idx = model->baset_idx;

    uint16_t j = 0;
    for (uint16_t i = 0; i < imm_model_nstates(model->super); ++i) {
        struct imm_state const* state = imm_model_state(model->super, i);
        if (imm_state_type_id(state) != NMM_FRAME_STATE_TYPE_ID)
            continue;

        struct nmm_frame_state const* s = nmm_frame_state_derived(state);
        struct nmm_base_table const*  baset = nmm_frame_state_base_table(s);
        khiter_t                      k = kh_get(baset, map, baset);
        if (k != kh_end(map))
            continue;

        struct baset_node* node = malloc(sizeof(*node));
        node->index = j++;
        node->baset = baset;
        int      ret = 0;
        khiter_t iter = kh_put(baset, map, node->baset, &ret);
        IMM_BUG(ret == -1 || ret == 0);
        kh_key(map, iter) = node->baset;
        kh_val(map, iter) = node;

        iter = kh_put(baset_idx, model->baset_idx, node->index, &ret);
        IMM_BUG(ret == -1 || ret == 0);
        kh_key(idx, iter) = node->index;
        kh_val(idx, iter) = node;
    }
}

static void create_codonp_map(struct nmm_model* model)
{
    khash_t(codonp)* map = model->codonp_map;
    khash_t(codonp_idx)* idx = model->codonp_idx;

    uint16_t j = 0;
    for (uint16_t i = 0; i < imm_model_nstates(model->super); ++i) {
        struct imm_state const* state = imm_model_state(model->super, i);
        if (imm_state_type_id(state) != NMM_CODON_STATE_TYPE_ID)
            continue;

        struct nmm_codon_state const* s = nmm_codon_state_derived(state);
        struct nmm_codon_lprob const* codonp = codon_state_codonp(s);
        khiter_t                      k = kh_get(codonp, map, codonp);
        if (k != kh_end(map))
            continue;

        struct codonp_node* node = malloc(sizeof(*node));
        node->index = j++;
        node->codonp = codonp;
        int      ret = 0;
        khiter_t iter = kh_put(codonp, map, node->codonp, &ret);
        IMM_BUG(ret == -1 || ret == 0);
        kh_key(map, iter) = node->codonp;
        kh_val(map, iter) = node;

        iter = kh_put(codonp_idx, model->codonp_idx, node->index, &ret);
        IMM_BUG(ret == -1 || ret == 0);
        kh_key(idx, iter) = node->index;
        kh_val(idx, iter) = node;
    }
}

static void create_codont_map(struct nmm_model* model)
{
    khash_t(codont)* map = model->codont_map;
    khash_t(codont_idx)* idx = model->codont_idx;

    uint16_t j = 0;
    for (uint16_t i = 0; i < imm_model_nstates(model->super); ++i) {
        struct imm_state const* state = imm_model_state(model->super, i);
        if (imm_state_type_id(state) != NMM_FRAME_STATE_TYPE_ID)
            continue;

        struct nmm_frame_state const* s = nmm_frame_state_derived(state);
        struct nmm_codon_table const* codont = nmm_frame_state_codon_table(s);
        khiter_t                      k = kh_get(codont, map, codont);
        if (k != kh_end(map))
            continue;

        struct codont_node* node = malloc(sizeof(*node));
        node->index = j++;
        node->codont = codont;
        int      ret = 0;
        khiter_t iter = kh_put(codont, map, node->codont, &ret);
        IMM_BUG(ret == -1 || ret == 0);
        kh_key(map, iter) = node->codont;
        kh_val(map, iter) = node;

        iter = kh_put(codont_idx, model->codont_idx, node->index, &ret);
        IMM_BUG(ret == -1 || ret == 0);
        kh_key(idx, iter) = node->index;
        kh_val(idx, iter) = node;
    }
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

static void deep_destroy(struct nmm_model const* model)
{
    __imm_model_deep_destroy(model->super);

    for (khint_t k = kh_begin(this->baset_map); k < kh_end(model->baset_map); ++k) {
        if (kh_exist(model->baset_map, k)) {
            struct baset_node const* node = kh_val(model->baset_map, k);
            nmm_base_table_destroy(node->baset);
            free_c(node);
        }
    }
    kh_destroy(baset, model->baset_map);

    for (khint_t k = kh_begin(model->codonp_map); k < kh_end(model->codonp_map); ++k) {
        if (kh_exist(model->codonp_map, k)) {
            struct codonp_node const* node = kh_val(model->codonp_map, k);
            nmm_codon_lprob_destroy(node->codonp);
            free_c(node);
        }
    }
    kh_destroy(codonp, model->codonp_map);

    for (khint_t k = kh_begin(model->codont_map); k < kh_end(model->codont_map); ++k) {
        if (kh_exist(model->codont_map, k)) {
            struct codont_node const* node = kh_val(model->codont_map, k);
            nmm_codon_table_destroy(node->codont);
            free_c(node);
        }
    }
    kh_destroy(codont, model->codont_map);

    free_c(model);
}

struct imm_hmm* nmm_model_hmm(struct nmm_model const* model) { return imm_model_hmm(model->super); }

struct imm_dp const* nmm_model_dp(struct nmm_model const* model)
{
    return imm_model_dp(model->super);
}

struct imm_state const* nmm_model_state(struct nmm_model const* model, uint16_t i)
{
    return imm_model_state(model->super, i);
}

uint16_t nmm_model_nstates(struct nmm_model const* model)
{
    return imm_model_nstates(model->super);
}

struct nmm_model const* nmm_model_read(FILE* stream)
{
    struct nmm_model* model = model_new();

    uint8_t abc_type_id = 0;
    if (fread(&abc_type_id, sizeof(abc_type_id), 1, stream) < 1) {
        imm_error("could not read abc type id");
        goto err;
    }

    struct imm_abc const* abc = read_abc(stream, abc_type_id);
    if (!abc) {
        imm_error("could not read abc");
        goto err;
    }
    __imm_model_set_abc(model->super, abc);

    if (read_baset(model, stream)) {
        imm_error("could not read baset");
        goto err;
    }

    if (read_codonp(model, stream)) {
        imm_error("could not read codonp");
        goto err;
    }

    if (read_codont(model, stream)) {
        imm_error("could not read codont");
        goto err;
    }

    if (__imm_model_read_hmm(model->super, stream)) {
        imm_error("could not read hmm");
        goto err;
    }

    if (__imm_model_read_dp(model->super, stream)) {
        imm_error("could not read dp");
        goto err;
    }

    return model;
err:
    deep_destroy(model);
    return NULL;
}

static struct imm_abc const* read_abc(FILE* stream, uint8_t type_id)
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
        abc = imm_abc_read(stream);
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

        iter = kh_put(baset_idx, model->baset_idx, node->index, &ret);
        IMM_BUG(ret == -1 || ret == 0);
        kh_key(model->baset_idx, iter) = node->index;
        kh_val(model->baset_idx, iter) = node;
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

        iter = kh_put(codonp_idx, model->codonp_idx, node->index, &ret);
        IMM_BUG(ret == -1 || ret == 0);
        kh_key(model->codonp_idx, iter) = node->index;
        kh_val(model->codonp_idx, iter) = node;
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

        iter = kh_put(codont_idx, model->codont_idx, node->index, &ret);
        IMM_BUG(ret == -1 || ret == 0);
        kh_key(model->codont_idx, iter) = node->index;
        kh_val(model->codont_idx, iter) = node;
    }

    return 0;
}

static struct imm_state const* read_state(struct imm_model const* model, FILE* stream, void* args)
{
    struct imm_state const* state = NULL;
    uint8_t                 type_id = 0;

    if (fread(&type_id, sizeof(type_id), 1, stream) < 1) {
        imm_error("could not read type id");
        return NULL;
    }

    switch (type_id) {
    case IMM_MUTE_STATE_TYPE_ID:
        if (!(state = imm_mute_state_read(stream, imm_model_abc(model))))
            imm_error("could not read mute state");
        break;
    case IMM_NORMAL_STATE_TYPE_ID:
        if (!(state = imm_normal_state_read(stream, imm_model_abc(model))))
            imm_error("could not read normal state");
        break;
    case IMM_TABLE_STATE_TYPE_ID:
        if (!(state = imm_table_state_read(stream, imm_model_abc(model))))
            imm_error("could not read table state");
        break;
    case NMM_CODON_STATE_TYPE_ID:
        if (!(state = nmm_codon_state_read(stream, (struct nmm_model const*)args)))
            imm_error("could not read codon_state");
        break;
    case NMM_FRAME_STATE_TYPE_ID:
        if (!(state = nmm_frame_state_read(stream, (struct nmm_model const*)args)))
            imm_error("could not read frame_state");
        break;
    default:
        imm_error("unknown state type id");
    }

    return state;
}

int nmm_model_write(struct nmm_model const* model, FILE* stream)
{
    if (write_abc(model, stream)) {
        imm_error("could not write abc");
        return 1;
    }

    if (write_baset(model, stream)) {
        imm_error("could not write_baset");
        return 1;
    }

    if (write_codonp(model, stream)) {
        imm_error("could not write_codonp");
        return 1;
    }

    if (write_codont(model, stream)) {
        imm_error("could not write_codont");
        return 1;
    }

    if (__imm_model_write_hmm(model->super, stream)) {
        imm_error("could not write hmm");
        return 1;
    }

    if (__imm_model_write_dp(model->super, stream)) {
        imm_error("could not write dp");
        return 1;
    }

    return 0;
}

static int write_abc(struct nmm_model const* model, FILE* stream)
{
    uint8_t type_id = imm_abc_type_id(imm_model_abc(model->super));
    if (fwrite(&type_id, sizeof(type_id), 1, stream) < 1) {
        imm_error("could not write abc type id");
        return 1;
    }

    if (imm_abc_write(imm_model_abc(model->super), stream)) {
        imm_error("could not write abc");
        return 1;
    }

    return 0;
}

static int write_baset(struct nmm_model const* model, FILE* stream)
{
    khash_t(baset_idx)* idx = model->baset_idx;
    IMM_BUG(kh_size(idx) > UINT32_MAX);

    uint32_t n = (uint32_t)kh_size(idx);

    if (fwrite(&n, sizeof(n), 1, stream) < 1) {
        imm_error("could not write nbaset");
        return 1;
    }

    for (uint32_t i = 0; i < n; ++i) {
        khiter_t iter = kh_get(baset_idx, idx, i);
        IMM_BUG(iter == kh_end(idx));

        struct baset_node const* node = kh_val(idx, iter);

        if (base_table_write(node->baset, stream))
            return 1;
    }

    return 0;
}

static int write_codonp(struct nmm_model const* model, FILE* stream)
{
    khash_t(codonp_idx)* idx = model->codonp_idx;
    IMM_BUG(kh_size(idx) > UINT32_MAX);

    uint32_t n = (uint32_t)kh_size(idx);

    if (fwrite(&n, sizeof(n), 1, stream) < 1) {
        imm_error("could not write ncodonp");
        return 1;
    }

    for (uint32_t i = 0; i < n; ++i) {
        khiter_t iter = kh_get(codonp_idx, idx, i);
        IMM_BUG(iter == kh_end(idx));

        struct codonp_node const* node = kh_val(idx, iter);

        if (codon_lprob_write(node->codonp, stream))
            return 1;
    }

    return 0;
}

static int write_codont(struct nmm_model const* model, FILE* stream)
{
    khash_t(codont_idx)* idx = model->codont_idx;
    IMM_BUG(kh_size(idx) > UINT32_MAX);

    uint32_t n = (uint32_t)kh_size(idx);

    if (fwrite(&n, sizeof(n), 1, stream) < 1) {
        imm_error("could not write ncodont");
        return 1;
    }

    for (uint32_t i = 0; i < n; ++i) {
        khiter_t iter = kh_get(codont_idx, idx, i);
        IMM_BUG(iter == kh_end(idx));

        struct codont_node const* node = kh_val(idx, iter);

        if (codon_table_write(node->codont, stream))
            return 1;
    }

    return 0;
}

static int write_state(struct imm_model const* model, FILE* stream, struct imm_state const* state,
                       void* args)
{
    uint8_t type_id = imm_state_type_id(state);
    if (fwrite(&type_id, sizeof(type_id), 1, stream) < 1) {
        imm_error("could not write state type id");
        return 1;
    }

    int errno = 0;
    switch (type_id) {
    case IMM_MUTE_STATE_TYPE_ID:
        errno = imm_mute_state_write(state, model, stream);
        break;
    case IMM_NORMAL_STATE_TYPE_ID:
        errno = imm_normal_state_write(state, model, stream);
        break;
    case IMM_TABLE_STATE_TYPE_ID:
        errno = imm_table_state_write(state, model, stream);
        break;
    case NMM_CODON_STATE_TYPE_ID:
        errno = nmm_codon_state_write(state, (struct nmm_model const*)args, stream);
        break;
    case NMM_FRAME_STATE_TYPE_ID:
        errno = nmm_frame_state_write(state, (struct nmm_model const*)args, stream);
        break;
    default:
        imm_error("unknown state type id");
        errno = 1;
    }
    return errno;
}

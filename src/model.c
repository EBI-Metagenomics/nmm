#include "model.h"
#include "amino_abc.h"
#include "base_abc.h"
#include "base_lprob.h"
#include "codon_lprob.h"
#include "codon_marg.h"
#include "codon_state.h"
#include "free.h"
#include "imm/imm.h"
#include "lib/khash.h"
#include "lib/khash_ptr.h"
#include "nmm/abc_types.h"
#include "nmm/base_abc.h"
#include "nmm/base_lprob.h"
#include "nmm/codon_lprob.h"
#include "nmm/codon_marg.h"
#include "nmm/codon_state.h"
#include "nmm/frame_state.h"
#include "nmm/model.h"
#include "nmm/state_types.h"
#include <imm/model.h>

#define DEF_STRUCT(MOD)                                                                                                \
    struct MOD##_node                                                                                                  \
    {                                                                                                                  \
        uint16_t                index;                                                                                 \
        struct nmm_##MOD const* MOD;                                                                                   \
    };                                                                                                                 \
    KHASH_MAP_INIT_PTR(MOD, struct MOD##_node*)                                                                        \
    KHASH_MAP_INIT_INT64(MOD##_idx, struct MOD##_node*)

DEF_STRUCT(base_lprob)
DEF_STRUCT(codon_lprob)
DEF_STRUCT(codon_marg)

struct nmm_model
{
    struct imm_model* super;

    khash_t(base_lprob) * base_lprob_map;
    khash_t(codon_lprob) * codon_lprob_map;
    khash_t(codon_marg) * codon_marg_map;

    khash_t(base_lprob_idx) * base_lprob_idx;
    khash_t(codon_lprob_idx) * codon_lprob_idx;
    khash_t(codon_marg_idx) * codon_marg_idx;
};

static void                    create_base_lprob_map(struct nmm_model* model, struct imm_hmm_block* block);
static void                    create_codon_lprob_map(struct nmm_model* model, struct imm_hmm_block* block);
static void                    create_codon_marg_map(struct nmm_model* model, struct imm_hmm_block* block);
static void                    deep_destroy(struct nmm_model const* model);
static void                    deep_destroy_base_lprob(khash_t(base_lprob) * base_lprob_map);
static void                    deep_destroy_codon_lprob(khash_t(codon_lprob) * codon_lprob_map);
static void                    deep_destroy_codon_marg(khash_t(codon_marg) * codon_marg_map);
static void                    destroy_base_lprob_map(khash_t(base_lprob) * base_lprob_map);
static void                    destroy_codon_lprob_map(khash_t(codon_lprob) * codon_lprob_map);
static void                    destroy_codon_marg_map(khash_t(codon_marg) * codon_marg_map);
struct nmm_model*              model_new(void);
static int                     read(struct imm_model* model, FILE* stream);
static struct imm_abc const*   read_abc(FILE* stream, uint8_t type_id);
static int                     read_base_lprob(struct nmm_model* model, FILE* stream);
static int                     read_codon_lprob(struct nmm_model* model, FILE* stream);
static int                     read_codon_marg(struct nmm_model* model, FILE* stream);
static struct imm_state const* read_state(struct imm_model const* model, FILE* stream, void* args);
static int                     write(struct imm_model const* model, FILE* stream);
static int                     write_abc(struct nmm_model const* model, FILE* stream);
static int                     write_base_lprob(struct nmm_model const* model, FILE* stream);
static int                     write_codon_lprob(struct nmm_model const* model, FILE* stream);
static int                     write_codon_marg(struct nmm_model const* model, FILE* stream);
static int write_state(struct imm_model const* model, FILE* stream, struct imm_state const* state, void* args);

struct imm_abc const* nmm_model_abc(struct nmm_model const* model) { return imm_model_abc(model->super); }

struct nmm_model const* nmm_model_create(struct imm_hmm* hmm, struct imm_dp const* dp)
{
    struct nmm_model* model = malloc(sizeof(*model));

    model->super = __imm_model_create(hmm, dp, read_state, model, write_state, model);
    if (!model->super) {
        free_c(model);
        return NULL;
    }

    model->base_lprob_map = kh_init(base_lprob);
    model->codon_lprob_map = kh_init(codon_lprob);
    model->codon_marg_map = kh_init(codon_marg);

    model->base_lprob_idx = kh_init(base_lprob_idx);
    model->codon_lprob_idx = kh_init(codon_lprob_idx);
    model->codon_marg_idx = kh_init(codon_marg_idx);

    struct imm_hmm_block* block = imm_model_get_hmm_block(model->super, 0);
    create_base_lprob_map(model, block);
    create_codon_lprob_map(model, block);
    create_codon_marg_map(model, block);

    return model;
}

void nmm_model_destroy(struct nmm_model const* model)
{
    imm_model_destroy(model->super);
    destroy_base_lprob_map(model->base_lprob_map);
    destroy_codon_lprob_map(model->codon_lprob_map);
    destroy_codon_marg_map(model->codon_marg_map);

    kh_destroy(base_lprob_idx, model->base_lprob_idx);
    kh_destroy(codon_lprob_idx, model->codon_lprob_idx);
    kh_destroy(codon_marg_idx, model->codon_marg_idx);
    free_c(model);
}

uint16_t nmm_model_nbase_lprobs(struct nmm_model const* model) { return (uint16_t)kh_size(model->base_lprob_map); }

uint16_t nmm_model_ncodon_lprobs(struct nmm_model const* model) { return (uint16_t)kh_size(model->codon_lprob_map); }

uint16_t nmm_model_ncodon_margs(struct nmm_model const* model) { return (uint16_t)kh_size(model->codon_marg_map); }

#define CREATE_MODEL_INDEX_FUNC(MOD)                                                                                   \
    uint16_t model_##MOD##_index(struct nmm_model const* model, struct nmm_##MOD const* MOD)                           \
    {                                                                                                                  \
        khiter_t i = kh_get(MOD, model->MOD##_map, MOD);                                                               \
        IMM_BUG(i == kh_end(model->MOD##_map));                                                                        \
                                                                                                                       \
        struct MOD##_node* node = kh_val(model->MOD##_map, i);                                                         \
                                                                                                                       \
        return node->index;                                                                                            \
    }

CREATE_MODEL_INDEX_FUNC(base_lprob)
CREATE_MODEL_INDEX_FUNC(codon_lprob)
CREATE_MODEL_INDEX_FUNC(codon_marg)

struct nmm_model* model_new(void)
{
    struct nmm_model* model = malloc(sizeof(*model));
    model->super = __imm_model_new(read_state, model, write_state, model);

    model->base_lprob_map = kh_init(base_lprob);
    model->codon_lprob_map = kh_init(codon_lprob);
    model->codon_marg_map = kh_init(codon_marg);

    model->base_lprob_idx = kh_init(base_lprob_idx);
    model->codon_lprob_idx = kh_init(codon_lprob_idx);
    model->codon_marg_idx = kh_init(codon_marg_idx);

    return model;
}

#define CREATE_MODEL_GET_FUNC(MOD)                                                                                     \
    struct nmm_##MOD const* nmm_model_##MOD(struct nmm_model const* model, uint16_t index)                             \
    {                                                                                                                  \
        khiter_t k = kh_get(MOD##_idx, model->MOD##_idx, index);                                                       \
        if (k == kh_end(model->MOD##_idx)) {                                                                           \
            imm_error(#MOD " index not found");                                                                        \
            return NULL;                                                                                               \
        }                                                                                                              \
        return kh_val(model->MOD##_idx, k)->MOD;                                                                       \
    }

CREATE_MODEL_GET_FUNC(base_lprob)
CREATE_MODEL_GET_FUNC(codon_lprob)
CREATE_MODEL_GET_FUNC(codon_marg)

#define CREATE_MAP_FUNC(MOD, ID, NAME)                                                                                 \
    static void create_##MOD##_map(struct nmm_model* model, struct imm_hmm_block* block)                               \
    {                                                                                                                  \
        khash_t(MOD)* map = model->MOD##_map;                                                                          \
        khash_t(MOD##_idx)* idx = model->MOD##_idx;                                                                    \
                                                                                                                       \
        uint16_t j = 0;                                                                                                \
        for (uint16_t i = 0; i < imm_hmm_block_nstates(block); ++i) {                                                  \
            struct imm_state const* state = imm_hmm_block_state(block, i);                                             \
            if (imm_state_type_id(state) != NMM_##ID##_STATE_TYPE_ID)                                                  \
                continue;                                                                                              \
                                                                                                                       \
            struct nmm_##NAME##_state const* s = nmm_##NAME##_state_derived(state);                                    \
            struct nmm_##MOD const*          MOD = nmm_##NAME##_state_##MOD(s);                                        \
            khiter_t                         k = kh_get(MOD, map, MOD);                                                \
            if (k != kh_end(map))                                                                                      \
                continue;                                                                                              \
                                                                                                                       \
            struct MOD##_node* node = malloc(sizeof(*node));                                                           \
            node->index = j++;                                                                                         \
            node->MOD = MOD;                                                                                           \
            int      ret = 0;                                                                                          \
            khiter_t iter = kh_put(MOD, map, node->MOD, &ret);                                                         \
            IMM_BUG(ret == -1 || ret == 0);                                                                            \
            kh_key(map, iter) = node->MOD;                                                                             \
            kh_val(map, iter) = node;                                                                                  \
                                                                                                                       \
            iter = kh_put(MOD##_idx, model->MOD##_idx, node->index, &ret);                                             \
            IMM_BUG(ret == -1 || ret == 0);                                                                            \
            kh_key(idx, iter) = node->index;                                                                           \
            kh_val(idx, iter) = node;                                                                                  \
        }                                                                                                              \
    }

CREATE_MAP_FUNC(base_lprob, FRAME, frame)
CREATE_MAP_FUNC(codon_lprob, CODON, codon)
CREATE_MAP_FUNC(codon_marg, FRAME, frame)

#define CREATE_DESTROY_FUNC(MOD)                                                                                       \
    static void destroy_##MOD##_map(khash_t(MOD) * MOD##_map)                                                          \
    {                                                                                                                  \
        for (khint_t k = kh_begin(MOD##_map); k < kh_end(MOD##_map); ++k)                                              \
            if (kh_exist(MOD##_map, k)) {                                                                              \
                free_c(kh_val(MOD##_map, k));                                                                          \
            }                                                                                                          \
        kh_destroy(MOD, MOD##_map);                                                                                    \
    }

CREATE_DESTROY_FUNC(base_lprob)
CREATE_DESTROY_FUNC(codon_lprob)
CREATE_DESTROY_FUNC(codon_marg)

#define CREATE_DEEP_DESTROY_FUNC(MOD)                                                                                  \
    static void deep_destroy_##MOD(khash_t(MOD) * MOD##_map)                                                           \
    {                                                                                                                  \
        for (khint_t k = kh_begin(MOD##_map); k < kh_end(MOD##_map); ++k) {                                            \
            if (kh_exist(MOD##_map, k)) {                                                                              \
                struct MOD##_node const* node = kh_val(MOD##_map, k);                                                  \
                nmm_##MOD##_destroy(node->MOD);                                                                        \
                free_c(node);                                                                                          \
            }                                                                                                          \
        }                                                                                                              \
        kh_destroy(MOD, MOD##_map);                                                                                    \
    }

CREATE_DEEP_DESTROY_FUNC(base_lprob)
CREATE_DEEP_DESTROY_FUNC(codon_lprob)
CREATE_DEEP_DESTROY_FUNC(codon_marg)

static void deep_destroy(struct nmm_model const* model)
{
    __imm_model_deep_destroy(model->super);

    deep_destroy_base_lprob(model->base_lprob_map);
    deep_destroy_codon_lprob(model->codon_lprob_map);
    deep_destroy_codon_marg(model->codon_marg_map);

    free_c(model);
}

void nmm_model_append_hmm_block(struct nmm_model* model, struct imm_hmm* hmm, struct imm_dp const* dp)
{
    imm_model_append_hmm_block(model->super, hmm, dp);
}

struct imm_hmm_block* nmm_model_get_hmm_block(struct nmm_model const* model, uint8_t i)
{
    return imm_model_get_hmm_block(model->super, i);
}

uint8_t nmm_model_nhmm_blocks(struct nmm_model const* model) { return imm_model_nhmm_blocks(model->super); }

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

    if (read_base_lprob(model, stream)) {
        imm_error("could not read base_lprob");
        goto err;
    }

    if (read_codon_lprob(model, stream)) {
        imm_error("could not read codon_lprob");
        goto err;
    }

    if (read_codon_marg(model, stream)) {
        imm_error("could not read codon_marg");
        goto err;
    }

    if (__imm_model_read_hmm_blocks(model->super, stream))
        goto err;

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

#define CREATE_READ_FUNC(MOD)                                                                                          \
    static int read_##MOD(struct nmm_model* model, FILE* stream)                                                       \
    {                                                                                                                  \
        uint16_t n = 0;                                                                                                \
        if (fread(&n, sizeof(n), 1, stream) < 1) {                                                                     \
            imm_error("could not read the number of " #MOD);                                                           \
            return 1;                                                                                                  \
        }                                                                                                              \
                                                                                                                       \
        for (uint16_t i = 0; i < n; ++i) {                                                                             \
            struct nmm_base_abc const* base_abc = nmm_base_abc_derived(imm_model_abc(model->super));                   \
            if (!base_abc)                                                                                             \
                return 1;                                                                                              \
                                                                                                                       \
            struct nmm_##MOD const* MOD = MOD##_read(stream, base_abc);                                                \
            if (!MOD)                                                                                                  \
                return 1;                                                                                              \
                                                                                                                       \
            struct MOD##_node* node = malloc(sizeof(*node));                                                           \
            node->index = i;                                                                                           \
            node->MOD = MOD;                                                                                           \
                                                                                                                       \
            int      ret = 0;                                                                                          \
            khiter_t iter = kh_put(MOD, model->MOD##_map, node->MOD, &ret);                                            \
            IMM_BUG(ret == -1 || ret == 0);                                                                            \
            kh_key(model->MOD##_map, iter) = node->MOD;                                                                \
            kh_val(model->MOD##_map, iter) = node;                                                                     \
                                                                                                                       \
            iter = kh_put(MOD##_idx, model->MOD##_idx, node->index, &ret);                                             \
            IMM_BUG(ret == -1 || ret == 0);                                                                            \
            kh_key(model->MOD##_idx, iter) = node->index;                                                              \
            kh_val(model->MOD##_idx, iter) = node;                                                                     \
        }                                                                                                              \
                                                                                                                       \
        return 0;                                                                                                      \
    }

CREATE_READ_FUNC(base_lprob)
CREATE_READ_FUNC(codon_lprob)
CREATE_READ_FUNC(codon_marg)

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

    if (write_base_lprob(model, stream)) {
        imm_error("could not write_base_lprob");
        return 1;
    }

    if (write_codon_lprob(model, stream)) {
        imm_error("could not write_codon_lprob");
        return 1;
    }

    if (write_codon_marg(model, stream)) {
        imm_error("could not write_codon_marg");
        return 1;
    }

    if (__imm_model_write_hmm_blocks(model->super, stream))
        return 1;

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

#define CREATE_WRITE_FUNC(MOD)                                                                                         \
    static int write_##MOD(struct nmm_model const* model, FILE* stream)                                                \
    {                                                                                                                  \
        khash_t(MOD##_idx)* idx = model->MOD##_idx;                                                                    \
        IMM_BUG(kh_size(idx) > UINT16_MAX);                                                                            \
                                                                                                                       \
        uint16_t n = (uint16_t)kh_size(idx);                                                                           \
                                                                                                                       \
        if (fwrite(&n, sizeof(n), 1, stream) < 1) {                                                                    \
            imm_error("could not write n" #MOD);                                                                       \
            return 1;                                                                                                  \
        }                                                                                                              \
                                                                                                                       \
        for (uint16_t i = 0; i < n; ++i) {                                                                             \
            khiter_t iter = kh_get(MOD##_idx, idx, i);                                                                 \
            IMM_BUG(iter == kh_end(idx));                                                                              \
                                                                                                                       \
            struct MOD##_node const* node = kh_val(idx, iter);                                                         \
                                                                                                                       \
            if (MOD##_write(node->MOD, stream))                                                                        \
                return 1;                                                                                              \
        }                                                                                                              \
                                                                                                                       \
        return 0;                                                                                                      \
    }

CREATE_WRITE_FUNC(base_lprob)
CREATE_WRITE_FUNC(codon_lprob)
CREATE_WRITE_FUNC(codon_marg)

static int write_state(struct imm_model const* model, FILE* stream, struct imm_state const* state, void* args)
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

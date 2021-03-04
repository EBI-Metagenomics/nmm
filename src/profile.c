#include "profile.h"
#include "amino_abc.h"
#include "base_abc.h"
#include "base_lprob.h"
#include "bug.h"
#include "codon_lprob.h"
#include "codon_marg.h"
#include "free.h"
#include "lib/khash.h"
#include "lib/khash_ptr.h"
#include "nmm/nmm.h"

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

struct nmm_profile
{
    struct imm_profile* super;

    khash_t(base_lprob) * base_lprob_map;
    khash_t(codon_lprob) * codon_lprob_map;
    khash_t(codon_marg) * codon_marg_map;

    khash_t(base_lprob_idx) * base_lprob_idx;
    khash_t(codon_lprob_idx) * codon_lprob_idx;
    khash_t(codon_marg_idx) * codon_marg_idx;
};

static void                    update_base_lprob_map(struct nmm_profile* prof, struct imm_model* model);
static void                    update_codon_lprob_map(struct nmm_profile* prof, struct imm_model* model);
static void                    update_codon_marg_map(struct nmm_profile* prof, struct imm_model* model);
static void                    deep_destroy(struct nmm_profile const* prof);
static void                    deep_destroy_base_lprob(khash_t(base_lprob) * base_lprob_map);
static void                    deep_destroy_codon_lprob(khash_t(codon_lprob) * codon_lprob_map);
static void                    deep_destroy_codon_marg(khash_t(codon_marg) * codon_marg_map);
static void                    destroy_base_lprob_map(khash_t(base_lprob) * base_lprob_map);
static void                    destroy_codon_lprob_map(khash_t(codon_lprob) * codon_lprob_map);
static void                    destroy_codon_marg_map(khash_t(codon_marg) * codon_marg_map);
static int                     read(struct imm_profile* prof, FILE* stream);
static struct imm_abc const*   read_abc(FILE* stream, uint8_t type_id);
static int                     read_base_lprob(struct nmm_profile* prof, FILE* stream);
static int                     read_codon_lprob(struct nmm_profile* prof, FILE* stream);
static int                     read_codon_marg(struct nmm_profile* prof, FILE* stream);
static struct imm_state const* read_state(struct imm_profile* prof, FILE* stream);
static int                     write(struct imm_profile const* prof, FILE* stream);
static int                     write_abc(struct nmm_profile const* prof, FILE* stream);
static int                     write_base_lprob(struct nmm_profile const* prof, FILE* stream);
static int                     write_codon_lprob(struct nmm_profile const* prof, FILE* stream);
static int                     write_codon_marg(struct nmm_profile const* prof, FILE* stream);
static int                     write_state(struct imm_profile const* prof, FILE* stream, struct imm_state const* state);

struct imm_abc const* nmm_profile_abc(struct nmm_profile const* prof) { return imm_profile_abc(prof->super); }

struct nmm_profile* nmm_profile_create(struct imm_abc const* abc)
{
    struct nmm_profile* prof = malloc(sizeof(*prof));

    prof->super = __imm_profile_create(abc, (struct imm_profile_vtable){read_state, write_state}, prof);
    if (!prof->super) {
        free_c(prof);
        return NULL;
    }
    prof->base_lprob_map = kh_init(base_lprob);
    prof->codon_lprob_map = kh_init(codon_lprob);
    prof->codon_marg_map = kh_init(codon_marg);

    prof->base_lprob_idx = kh_init(base_lprob_idx);
    prof->codon_lprob_idx = kh_init(codon_lprob_idx);
    prof->codon_marg_idx = kh_init(codon_marg_idx);

    return prof;
}

void nmm_profile_destroy(struct nmm_profile const* prof, bool deep)
{
    imm_profile_destroy(prof->super, deep);

    if (deep) {
        for (uint16_t i = 0; i < nmm_profile_nbase_lprobs(prof); ++i)
            nmm_base_lprob_destroy(nmm_profile_base_lprob(prof, i));

        for (uint16_t i = 0; i < nmm_profile_ncodon_margs(prof); ++i)
            nmm_codon_marg_destroy(nmm_profile_codon_marg(prof, i));

        for (uint16_t i = 0; i < nmm_profile_ncodon_lprobs(prof); ++i)
            nmm_codon_lprob_destroy(nmm_profile_codon_lprob(prof, i));
    }

    destroy_base_lprob_map(prof->base_lprob_map);
    destroy_codon_lprob_map(prof->codon_lprob_map);
    destroy_codon_marg_map(prof->codon_marg_map);

    kh_destroy(base_lprob_idx, prof->base_lprob_idx);
    kh_destroy(codon_lprob_idx, prof->codon_lprob_idx);
    kh_destroy(codon_marg_idx, prof->codon_marg_idx);
    free_c(prof);
}

void nmm_profile_free(struct nmm_profile const* prof)
{
    imm_profile_free(prof->super);
    free_c(prof);
}

uint16_t nmm_profile_nbase_lprobs(struct nmm_profile const* prof) { return (uint16_t)kh_size(prof->base_lprob_map); }

uint16_t nmm_profile_ncodon_lprobs(struct nmm_profile const* prof) { return (uint16_t)kh_size(prof->codon_lprob_map); }

uint16_t nmm_profile_ncodon_margs(struct nmm_profile const* prof) { return (uint16_t)kh_size(prof->codon_marg_map); }

#define CREATE_PROFILE_INDEX_FUNC(MOD)                                                                                 \
    uint16_t profile_##MOD##_index(struct nmm_profile const* prof, struct nmm_##MOD const* MOD)                        \
    {                                                                                                                  \
        khiter_t i = kh_get(MOD, prof->MOD##_map, MOD);                                                                \
        BUG(i == kh_end(prof->MOD##_map));                                                                             \
                                                                                                                       \
        struct MOD##_node* node = kh_val(prof->MOD##_map, i);                                                          \
                                                                                                                       \
        return node->index;                                                                                            \
    }

CREATE_PROFILE_INDEX_FUNC(base_lprob)
CREATE_PROFILE_INDEX_FUNC(codon_lprob)
CREATE_PROFILE_INDEX_FUNC(codon_marg)

#define CREATE_MODEL_GET_FUNC(MOD)                                                                                     \
    struct nmm_##MOD const* nmm_profile_##MOD(struct nmm_profile const* prof, uint16_t index)                          \
    {                                                                                                                  \
        khiter_t k = kh_get(MOD##_idx, prof->MOD##_idx, index);                                                        \
        if (k == kh_end(prof->MOD##_idx)) {                                                                            \
            imm_error(#MOD " index not found");                                                                        \
            return NULL;                                                                                               \
        }                                                                                                              \
        return kh_val(prof->MOD##_idx, k)->MOD;                                                                        \
    }

CREATE_MODEL_GET_FUNC(base_lprob)
CREATE_MODEL_GET_FUNC(codon_lprob)
CREATE_MODEL_GET_FUNC(codon_marg)

#define CREATE_UPDATE_MAP_FUNC(MOD, ID, NAME)                                                                          \
    static void update_##MOD##_map(struct nmm_profile* prof, struct imm_model* model)                                  \
    {                                                                                                                  \
        khash_t(MOD)* map = prof->MOD##_map;                                                                           \
        khash_t(MOD##_idx)* idx = prof->MOD##_idx;                                                                     \
                                                                                                                       \
        uint16_t index = (uint16_t)kh_size(map);                                                                       \
        for (uint16_t i = 0; i < imm_model_nstates(model); ++i) {                                                      \
            struct imm_state const* state = imm_model_state(model, i);                                                 \
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
            node->index = index++;                                                                                     \
            node->MOD = MOD;                                                                                           \
            int      ret = 0;                                                                                          \
            khiter_t iter = kh_put(MOD, map, node->MOD, &ret);                                                         \
            BUG(ret == -1 || ret == 0);                                                                                \
            kh_key(map, iter) = node->MOD;                                                                             \
            kh_val(map, iter) = node;                                                                                  \
                                                                                                                       \
            iter = kh_put(MOD##_idx, prof->MOD##_idx, node->index, &ret);                                              \
            BUG(ret == -1 || ret == 0);                                                                                \
            kh_key(idx, iter) = node->index;                                                                           \
            kh_val(idx, iter) = node;                                                                                  \
        }                                                                                                              \
    }

CREATE_UPDATE_MAP_FUNC(base_lprob, FRAME, frame)
CREATE_UPDATE_MAP_FUNC(codon_lprob, CODON, codon)
CREATE_UPDATE_MAP_FUNC(codon_marg, FRAME, frame)

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

static void deep_destroy(struct nmm_profile const* prof)
{
    __imm_profile_deep_destroy(prof->super);

    deep_destroy_base_lprob(prof->base_lprob_map);
    deep_destroy_codon_lprob(prof->codon_lprob_map);
    deep_destroy_codon_marg(prof->codon_marg_map);

    free_c(prof);
}

void nmm_profile_append_model(struct nmm_profile* prof, struct imm_model* model)
{
    imm_profile_append_model(prof->super, model);
    update_base_lprob_map(prof, model);
    update_codon_lprob_map(prof, model);
    update_codon_marg_map(prof, model);
}

struct imm_model* nmm_profile_get_model(struct nmm_profile const* prof, uint8_t i)
{
    return imm_profile_get_model(prof->super, i);
}

uint8_t nmm_profile_nmodels(struct nmm_profile const* prof) { return imm_profile_nmodels(prof->super); }

struct nmm_profile* nmm_profile_read(FILE* stream)
{
    uint8_t abc_type_id = 0;
    if (fread(&abc_type_id, sizeof(abc_type_id), 1, stream) < 1) {
        imm_error("could not read abc type id");
        return NULL;
    }

    struct imm_abc const* abc = read_abc(stream, abc_type_id);
    if (!abc) {
        imm_error("could not read abc");
        return NULL;
    }

    struct nmm_profile* prof = nmm_profile_create(abc);

    __imm_profile_set_abc(prof->super, abc);

    if (read_base_lprob(prof, stream)) {
        imm_error("could not read base_lprob");
        goto err;
    }

    if (read_codon_lprob(prof, stream)) {
        imm_error("could not read codon_lprob");
        goto err;
    }

    if (read_codon_marg(prof, stream)) {
        imm_error("could not read codon_marg");
        goto err;
    }

    if (__imm_profile_read_models(prof->super, stream))
        goto err;

    return prof;
err:
    deep_destroy(prof);
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
    static int read_##MOD(struct nmm_profile* prof, FILE* stream)                                                      \
    {                                                                                                                  \
        uint16_t n = 0;                                                                                                \
        if (fread(&n, sizeof(n), 1, stream) < 1) {                                                                     \
            imm_error("could not read the number of " #MOD);                                                           \
            return 1;                                                                                                  \
        }                                                                                                              \
                                                                                                                       \
        for (uint16_t i = 0; i < n; ++i) {                                                                             \
            struct nmm_base_abc const* base_abc = nmm_base_abc_derived(imm_profile_abc(prof->super));                  \
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
            khiter_t iter = kh_put(MOD, prof->MOD##_map, node->MOD, &ret);                                             \
            BUG(ret == -1 || ret == 0);                                                                                \
            kh_key(prof->MOD##_map, iter) = node->MOD;                                                                 \
            kh_val(prof->MOD##_map, iter) = node;                                                                      \
                                                                                                                       \
            iter = kh_put(MOD##_idx, prof->MOD##_idx, node->index, &ret);                                              \
            BUG(ret == -1 || ret == 0);                                                                                \
            kh_key(prof->MOD##_idx, iter) = node->index;                                                               \
            kh_val(prof->MOD##_idx, iter) = node;                                                                      \
        }                                                                                                              \
                                                                                                                       \
        return 0;                                                                                                      \
    }

CREATE_READ_FUNC(base_lprob)
CREATE_READ_FUNC(codon_lprob)
CREATE_READ_FUNC(codon_marg)

static struct imm_state const* read_state(struct imm_profile* prof, FILE* stream)
{
    struct imm_state const* state = NULL;
    uint8_t                 type_id = 0;

    if (fread(&type_id, sizeof(type_id), 1, stream) < 1) {
        imm_error("could not read type id");
        return NULL;
    }

    switch (type_id) {
    case IMM_MUTE_STATE_TYPE_ID:
        if (!(state = imm_mute_state_read(stream, imm_profile_abc(prof))))
            imm_error("could not read mute state");
        break;
    case IMM_NORMAL_STATE_TYPE_ID:
        if (!(state = imm_normal_state_read(stream, imm_profile_abc(prof))))
            imm_error("could not read normal state");
        break;
    case IMM_TABLE_STATE_TYPE_ID:
        if (!(state = imm_table_state_read(stream, imm_profile_abc(prof))))
            imm_error("could not read table state");
        break;
    case NMM_CODON_STATE_TYPE_ID:
        if (!(state = nmm_codon_state_read(stream, __imm_profile_derived(prof))))
            imm_error("could not read codon_state");
        break;
    case NMM_FRAME_STATE_TYPE_ID:
        if (!(state = nmm_frame_state_read(stream, __imm_profile_derived(prof))))
            imm_error("could not read frame_state");
        break;
    default:
        imm_error("unknown state type id");
    }

    return state;
}

int nmm_profile_write(struct nmm_profile const* prof, FILE* stream)
{
    if (write_abc(prof, stream)) {
        imm_error("could not write abc");
        return 1;
    }

    if (write_base_lprob(prof, stream)) {
        imm_error("could not write_base_lprob");
        return 1;
    }

    if (write_codon_lprob(prof, stream)) {
        imm_error("could not write_codon_lprob");
        return 1;
    }

    if (write_codon_marg(prof, stream)) {
        imm_error("could not write_codon_marg");
        return 1;
    }

    if (__imm_profile_write_models(prof->super, stream))
        return 1;

    return 0;
}

static int write_abc(struct nmm_profile const* prof, FILE* stream)
{
    uint8_t type_id = imm_abc_type_id(imm_profile_abc(prof->super));
    if (fwrite(&type_id, sizeof(type_id), 1, stream) < 1) {
        imm_error("could not write abc type id");
        return 1;
    }

    if (imm_abc_write(imm_profile_abc(prof->super), stream)) {
        imm_error("could not write abc");
        return 1;
    }

    return 0;
}

#define CREATE_WRITE_FUNC(MOD)                                                                                         \
    static int write_##MOD(struct nmm_profile const* prof, FILE* stream)                                               \
    {                                                                                                                  \
        khash_t(MOD##_idx)* idx = prof->MOD##_idx;                                                                     \
        BUG(kh_size(idx) > UINT16_MAX);                                                                                \
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
            BUG(iter == kh_end(idx));                                                                                  \
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

static int write_state(struct imm_profile const* prof, FILE* stream, struct imm_state const* state)
{
    uint8_t type_id = imm_state_type_id(state);
    if (fwrite(&type_id, sizeof(type_id), 1, stream) < 1) {
        imm_error("could not write state type id");
        return 1;
    }

    int errno = 0;
    switch (type_id) {
    case IMM_MUTE_STATE_TYPE_ID:
        errno = imm_mute_state_write(state, prof, stream);
        break;
    case IMM_NORMAL_STATE_TYPE_ID:
        errno = imm_normal_state_write(state, prof, stream);
        break;
    case IMM_TABLE_STATE_TYPE_ID:
        errno = imm_table_state_write(state, prof, stream);
        break;
    case NMM_CODON_STATE_TYPE_ID:
        errno = nmm_codon_state_write(state, __imm_profile_derived_c(prof), stream);
        break;
    case NMM_FRAME_STATE_TYPE_ID:
        errno = nmm_frame_state_write(state, __imm_profile_derived_c(prof), stream);
        break;
    default:
        imm_error("unknown state type id");
        errno = 1;
    }
    return errno;
}

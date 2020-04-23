#include "nmm/base_abc.h"
#include "free.h"
#include "imm/imm.h"
#include "nmm/abc_types.h"
#include <stdlib.h>
#include <string.h>

static uint8_t               base_abc_type_id(struct imm_abc const* abc);
static void                  base_abc_destroy(struct imm_abc const* abc);
static struct imm_abc const* base_abc_clone(struct imm_abc const* abc);

static struct imm_abc_vtable const __vtable = {base_abc_type_id, NULL, base_abc_clone};

struct nmm_base_abc const* nmm_base_abc_create(char const* symbols, char const any_symbol)
{
    if (strlen(symbols) != NMM_BASE_ABC_SIZE) {
        imm_error("alphabet length must be %u", NMM_BASE_ABC_SIZE);
        return NULL;
    }

    struct nmm_base_abc* base_abc = malloc(sizeof(struct nmm_base_abc));
    base_abc->parent = __imm_abc_create_parent(symbols, any_symbol, __vtable, base_abc);
    return base_abc;
}

void nmm_base_abc_destroy(struct nmm_base_abc const* base_abc) { free_c(base_abc); }

struct nmm_base_abc const* nmm_base_abc_child(struct imm_abc const* abc)
{
    if (imm_abc_type_id(abc) != NMM_BASE_ABC_TYPE_ID) {
        imm_error("could not cast to base_abc");
        return NULL;
    }
    return __imm_abc_child(abc);
}

static uint8_t base_abc_type_id(struct imm_abc const* abc) { return NMM_BASE_ABC_TYPE_ID; }

static void base_abc_destroy(struct imm_abc const* abc) { free_c(__imm_abc_child(abc)); }

static struct imm_abc const* base_abc_clone(struct imm_abc const* abc)
{
    struct nmm_base_abc* base_abc = malloc(sizeof(struct nmm_base_abc));
    base_abc->parent = __imm_abc_clone_parent(abc);
    return base_abc->parent;
}

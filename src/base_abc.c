#include "base_abc.h"
#include "free.h"
#include "imm/imm.h"
#include "nmm/abc_types.h"
#include "nmm/base_abc.h"
#include <stdlib.h>
#include <string.h>

static uint8_t               base_abc_type_id(struct imm_abc const* abc);
static int                   base_write(struct imm_abc const* abc, FILE* stream);
static void                  base_abc_destroy(struct imm_abc const* abc);
static struct imm_abc const* base_abc_clone(struct imm_abc const* abc);

static struct imm_abc_vtable const __vtable = {base_abc_type_id, __imm_abc_write, base_abc_destroy,
                                               base_abc_clone};

struct nmm_base_abc const* nmm_base_abc_create(char const* symbols, char const any_symbol)
{
    if (strlen(symbols) != NMM_BASE_ABC_SIZE) {
        imm_error("alphabet length must be %u", NMM_BASE_ABC_SIZE);
        return NULL;
    }

    struct nmm_base_abc* base_abc = malloc(sizeof(*base_abc));
    struct imm_abc*      super = __imm_abc_create(symbols, any_symbol, base_abc);
    base_abc->super = super;
    super->vtable = __vtable;

    return base_abc;
}

void nmm_base_abc_destroy(struct nmm_base_abc const* base_abc)
{
    base_abc->super->vtable.destroy(base_abc->super);
}

struct nmm_base_abc const* nmm_base_abc_derived(struct imm_abc const* abc)
{
    if (imm_abc_type_id(abc) != NMM_BASE_ABC_TYPE_ID) {
        imm_error("could not cast to base_abc");
        return NULL;
    }
    return __imm_abc_derived(abc);
}

struct imm_abc const* base_abc_read(FILE* stream)
{
    struct nmm_base_abc* base_abc = malloc(sizeof(*base_abc));

    struct imm_abc* abc = __imm_abc_read(stream);
    if (!abc) {
        imm_error("could not read base_abc");
        free_c(base_abc);
        return NULL;
    }

    abc->vtable = __vtable;
    abc->derived = base_abc;

    return abc;
}

static uint8_t base_abc_type_id(struct imm_abc const* abc) { return NMM_BASE_ABC_TYPE_ID; }

static void base_abc_destroy(struct imm_abc const* abc)
{
    struct nmm_base_abc const* base_abc = __imm_abc_derived(abc);
    __imm_abc_destroy(abc);
    free_c(base_abc);
}

static struct imm_abc const* base_abc_clone(struct imm_abc const* abc)
{
    struct nmm_base_abc* base_abc = malloc(sizeof(*base_abc));
    base_abc->super = __imm_abc_clone(abc);
    return base_abc->super;
}

#include "nmm/amino_abc.h"
#include "amino_abc.h"
#include "free.h"
#include "imm/imm.h"
#include "nmm/abc_types.h"
#include <stdlib.h>
#include <string.h>

static uint8_t               amino_abc_type_id(struct imm_abc const* abc);
static int                   amino_write(struct imm_abc const* abc, FILE* stream);
static void                  amino_abc_destroy(struct imm_abc const* abc);
static struct imm_abc const* amino_abc_clone(struct imm_abc const* abc);

static struct imm_abc_vtable const __vtable = {amino_abc_type_id, __imm_abc_write,
                                               amino_abc_destroy, amino_abc_clone};

struct nmm_amino_abc const* nmm_amino_abc_create(char const* symbols, char const any_symbol)
{
    if (strlen(symbols) != NMM_AMINO_ABC_SIZE) {
        imm_error("alphabet length must be %u", NMM_AMINO_ABC_SIZE);
        return NULL;
    }

    struct nmm_amino_abc* amino_abc = malloc(sizeof(*amino_abc));
    struct imm_abc*       super = __imm_abc_create(symbols, any_symbol, amino_abc);
    amino_abc->super = super;
    super->vtable = __vtable;

    return amino_abc;
}

void nmm_amino_abc_destroy(struct nmm_amino_abc const* amino_abc)
{
    amino_abc->super->vtable.destroy(amino_abc->super);
}

struct nmm_amino_abc const* nmm_amino_abc_derived(struct imm_abc const* abc)
{
    if (imm_abc_type_id(abc) != NMM_AMINO_ABC_TYPE_ID) {
        imm_error("could not cast to amino_abc");
        return NULL;
    }
    return __imm_abc_derived(abc);
}

struct imm_abc const* amino_abc_read(FILE* stream)
{
    struct nmm_amino_abc* amino_abc = malloc(sizeof(*amino_abc));

    struct imm_abc* abc = __imm_abc_read(stream);
    if (!abc) {
        imm_error("could not read amino_abc");
        free_c(amino_abc);
        return NULL;
    }

    abc->vtable = __vtable;
    abc->derived = amino_abc;

    return abc;
}

static uint8_t amino_abc_type_id(struct imm_abc const* abc) { return NMM_AMINO_ABC_TYPE_ID; }

static void amino_abc_destroy(struct imm_abc const* abc)
{
    struct nmm_amino_abc const* amino_abc = __imm_abc_derived(abc);
    __imm_abc_destroy(abc);
    free_c(amino_abc);
}

static struct imm_abc const* amino_abc_clone(struct imm_abc const* abc)
{
    struct nmm_amino_abc* amino_abc = malloc(sizeof(*amino_abc));
    amino_abc->super = __imm_abc_clone(abc);
    return amino_abc->super;
}

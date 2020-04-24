#include "nmm/amino_abc.h"
#include "free.h"
#include "imm/imm.h"
#include "nmm/abc_types.h"
#include <stdlib.h>
#include <string.h>

static uint8_t               amino_abc_type_id(struct imm_abc const* abc);
static void                  amino_abc_destroy(struct imm_abc const* abc);
static struct imm_abc const* amino_abc_clone(struct imm_abc const* abc);

static struct imm_abc_vtable const __vtable = {amino_abc_type_id, amino_abc_destroy,
                                               amino_abc_clone};

struct nmm_amino_abc const* nmm_amino_abc_create(char const* symbols, char const any_symbol)
{
    if (strlen(symbols) != NMM_AMINO_ABC_SIZE) {
        imm_error("alphabet length must be %u", NMM_AMINO_ABC_SIZE);
        return NULL;
    }

    struct nmm_amino_abc* amino_abc = malloc(sizeof(struct nmm_amino_abc));
    amino_abc->parent = __imm_abc_create_parent(symbols, any_symbol, __vtable, amino_abc);
    return amino_abc;
}

void nmm_amino_abc_destroy(struct nmm_amino_abc const* amino_abc)
{
    struct imm_abc const* parent = amino_abc->parent;
    amino_abc_destroy(parent);
    __imm_abc_destroy_parent(parent);
}

struct nmm_amino_abc const* nmm_amino_abc_child(struct imm_abc const* abc)
{
    if (imm_abc_type_id(abc) != NMM_AMINO_ABC_TYPE_ID) {
        imm_error("could not cast to amino_abc");
        return NULL;
    }
    return __imm_abc_child(abc);
}

static uint8_t amino_abc_type_id(struct imm_abc const* abc) { return NMM_AMINO_ABC_TYPE_ID; }

static void amino_abc_destroy(struct imm_abc const* abc) { free_c(__imm_abc_child(abc)); }

static struct imm_abc const* amino_abc_clone(struct imm_abc const* abc)
{
    struct nmm_amino_abc* amino_abc = malloc(sizeof(struct nmm_amino_abc));
    amino_abc->parent = __imm_abc_clone_parent(abc);
    return amino_abc->parent;
}

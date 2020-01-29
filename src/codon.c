#include "nmm/codon.h"
#include "free.h"
#include "imm/imm.h"
#include "nmm/base.h"
#include <stdbool.h>
#include <stdlib.h>

struct nmm_codon
{
    struct nmm_base const* base;
    char                   a;
    char                   b;
    char                   c;
};

struct nmm_codon* nmm_codon_create(struct nmm_base const* base)
{
    struct nmm_codon* codon = malloc(sizeof(struct nmm_codon));
    codon->base = base;

    struct imm_abc const* abc = nmm_base_get_abc(base);
    codon->a = imm_abc_any_symbol(abc);
    codon->b = imm_abc_any_symbol(abc);
    codon->c = imm_abc_any_symbol(abc);

    return codon;
}

void nmm_codon_destroy(struct nmm_codon const* codon) { free_c(codon); }

struct nmm_base const* nmm_codon_get_base(struct nmm_codon const* codon)
{
    return codon->base;
}

int nmm_codon_set(struct nmm_codon* codon, char a, char b, char c)
{
    struct imm_abc const* abc = nmm_base_get_abc(codon->base);
    bool                  ok = imm_abc_has_symbol(abc, a) || a == imm_abc_any_symbol(abc);
    ok = ok && (imm_abc_has_symbol(abc, b) || b == imm_abc_any_symbol(abc));
    ok = ok && (imm_abc_has_symbol(abc, c) || c == imm_abc_any_symbol(abc));

    if (!ok) {
        imm_error("wrong codon symbol");
        return 1;
    }

    codon->a = a;
    codon->b = b;
    codon->c = c;

    return 0;
}

void nmm_codon_get(struct nmm_codon const* codon, char* a, char* b, char* c)
{
    *a = codon->a;
    *b = codon->b;
    *c = codon->c;
}

#ifndef NMM_BASE_STATIC_H
#define NMM_BASE_STATIC_H

struct imm_abc;

struct nmm_base
{
    struct imm_abc const* abc;
};

static inline struct imm_abc const* base_get_abc(struct nmm_base const* base)
{
    return base->abc;
}

#endif

#ifndef NMM_BASET_STATIC_H
#define NMM_BASET_STATIC_H

struct nmm_baset
{
    struct nmm_base const* base;
    double                 lprobs[4];
};

static inline struct nmm_base const* baset_get_base(struct nmm_baset const* baset)
{
    return baset->base;
}

#endif

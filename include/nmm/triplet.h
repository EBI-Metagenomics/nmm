#ifndef NMM_TRIPLET_H
#define NMM_TRIPLET_H

struct nmm_triplet
{
    char a;
    char b;
    char c;
};

static inline struct nmm_triplet NMM_TRIPLET(char a, char b, char c)
{
    return (struct nmm_triplet){a, b, c};
}

#endif

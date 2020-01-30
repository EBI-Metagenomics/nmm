#include "cartes.h"
#include <stdio.h>
#include <stdlib.h>

struct cartes
{
    char const* set;
    unsigned    set_size;
    unsigned    times;
    unsigned    iter_idx;
    char*       item;
    unsigned    nitems;
};

static unsigned ipow(unsigned base, unsigned exp);

struct cartes* cartes_create(char const* set, unsigned set_size, unsigned times)
{
    struct cartes* cartes = malloc(sizeof(struct cartes));

    cartes->set = set;
    cartes->set_size = set_size;
    cartes->times = times;
    cartes->item = malloc(sizeof(char) * (times + 1));
    cartes->item[times] = '\0';
    cartes->iter_idx = 0;
    cartes->nitems = ipow(set_size, times);

    return cartes;
}

void cartes_destroy(struct cartes const* cartes)
{
    free(cartes->item);
    free((void*)cartes);
}

char const* cartes_next(struct cartes* cartes)
{
    if (cartes->iter_idx == cartes->nitems)
        return NULL;

    char*    item = cartes->item;
    unsigned idx = cartes->iter_idx++;
    unsigned set_size = cartes->set_size;

    for (unsigned i = 0; i < cartes->times; ++i) {
        item[i] = cartes->set[(idx % ipow(set_size, i + 1)) / ipow(set_size, i)];
    }

    return item;
}

static unsigned ipow(unsigned base, unsigned exp)
{
    unsigned result = 1;
    for (;;) {
        if (exp & 1)
            result *= base;
        exp >>= 1;
        if (!exp)
            break;
        base *= base;
    }

    return result;
}

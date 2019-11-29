#include "cartes.h"
#include <stdio.h>
#include <stdlib.h>

struct cartes
{
    char const* set;
    int         set_size;
    int         times;
    int         iter_idx;
    char*       item;
    int         nitems;
};

static int ipow(int base, int exp);

struct cartes* cartes_create(char const* set, int set_size, int times)
{
    struct cartes* cartes = malloc(sizeof(struct cartes));

    cartes->set = set;
    cartes->set_size = set_size;
    cartes->times = times;
    cartes->item = malloc(sizeof(char) * (set_size + 1));
    cartes->item[set_size] = '\0';
    cartes->iter_idx = 0;
    cartes->nitems = ipow(set_size, times);

    return cartes;
}

void cartes_destroy(struct cartes* cartes)
{
    if (!cartes) {
        fprintf(stderr, "cartes should not be NULL");
        exit(1);
        return;
    }

    free(cartes->item);
    free(cartes);
}

char const* cartes_next(struct cartes* cartes)
{
    if (cartes->iter_idx == cartes->nitems)
        return NULL;

    char* item = cartes->item;
    int   idx = cartes->iter_idx++;
    int   set_size = cartes->set_size;

    for (int i = 0; i < cartes->times; ++i) {
        item[i] = cartes->set[(idx % ipow(set_size, i + 1)) / ipow(set_size, i)];
    }

    return item;
}

static int ipow(int base, int exp)
{
    int result = 1;
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

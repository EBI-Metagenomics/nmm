#ifndef CARTES_H
#define CARTES_H

struct cartes
{
    char const* set;
    int         set_size;
    int         times;
    int         iter_idx;
    char*       item;
    int         nitems;
};

struct cartes* cartes_create(char const* set, int set_size, int times);
void           cartes_destroy(struct cartes* cartes);
char const*    cartes_next(struct cartes* cartes);

/* int  cartes_nitems(int seq_len, int times); */
/* void cartes_item(char const* seq, int seq_len, int times, int idx, char* item); */

#endif

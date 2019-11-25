#ifndef CARTES_H
#define CARTES_H

struct cartes;

struct cartes* cartes_create(char const* set, int set_size, int times);
void           cartes_destroy(struct cartes* cartes);
char const*    cartes_next(struct cartes* cartes);

#endif

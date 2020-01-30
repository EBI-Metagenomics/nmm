#ifndef CARTES_H
#define CARTES_H

struct cartes;

struct cartes* cartes_create(char const* set, unsigned set_size, unsigned times);
void           cartes_destroy(struct cartes const* cartes);
char const*    cartes_next(struct cartes* cartes);

#endif

#ifndef IO_H
#define IO_H

#include <inttypes.h>

struct nmm_base_table;
struct nmm_codon_lprob;
struct nmm_codon_table;
struct nmm_io;

uint32_t io_baset_index(struct nmm_io const* io, struct nmm_base_table const* baset);
struct nmm_base_table const* io_get_baset(struct nmm_io const* io, uint32_t index);
uint32_t io_codonp_index(struct nmm_io const* io, struct nmm_codon_lprob const* codonp);
struct nmm_codon_lprob const* io_get_codonp(struct nmm_io const* io, uint32_t index);
uint32_t io_codont_index(struct nmm_io const* io, struct nmm_codon_table const* codont);
struct nmm_codon_table const* io_get_codont(struct nmm_io const* io, uint32_t index);

#endif

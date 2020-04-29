#ifndef IO_H
#define IO_H

#include <inttypes.h>

struct nmm_base_table;
struct nmm_codon_lprob;
struct nmm_codon_table;
struct nmm_io;

uint32_t io_baset_index(struct nmm_io const* io, struct nmm_base_table const* baset);
uint32_t io_codonp_index(struct nmm_io const* io, struct nmm_codon_lprob const* codonp);
uint32_t io_codont_index(struct nmm_io const* io, struct nmm_codon_table const* codont);

#endif

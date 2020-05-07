#ifndef NMM_IO_H
#define NMM_IO_H

#include "nmm/export.h"
#include <inttypes.h>
#include <stdio.h>

struct imm_dp;
struct imm_hmm;
struct imm_io;
struct nmm_io;
struct nmm_base_table;
struct nmm_codon_table;

NMM_API struct nmm_io const* nmm_io_create(struct imm_hmm* hmm, struct imm_dp const* dp);
NMM_API struct nmm_io const* nmm_io_create_from_file(FILE* stream);
NMM_API struct nmm_io const* nmm_io_derived(struct imm_io const* io);
NMM_API void                 nmm_io_destroy(struct nmm_io const* io);
NMM_API struct imm_io const* nmm_io_super(struct nmm_io const* io);
NMM_API int                  nmm_io_write(struct nmm_io const* io, FILE* stream);

NMM_API uint32_t                      nmm_io_nbase_tables(struct nmm_io const* io);
NMM_API struct nmm_base_table const*  nmm_io_base_table(struct nmm_io const* io, uint32_t index);
NMM_API uint32_t                      nmm_io_ncodon_tables(struct nmm_io const* io);
NMM_API struct nmm_codon_table const* nmm_io_codon_table(struct nmm_io const* io,
                                                            uint32_t             index);
NMM_API uint32_t                      nmm_io_ncodon_lprobs(struct nmm_io const* io);
NMM_API struct nmm_codon_lprob const* nmm_io_codon_lprob(struct nmm_io const* io,
                                                            uint32_t             index);

#endif

#ifndef NMM_IO_H
#define NMM_IO_H

#include "nmm/export.h"
#include <stdio.h>

struct imm_dp;
struct imm_hmm;
struct imm_io;
struct nmm_io;

NMM_EXPORT struct nmm_io const* nmm_io_create(struct imm_hmm* hmm, struct imm_dp const* dp);
NMM_EXPORT struct nmm_io const* nmm_io_create_from_file(FILE* stream);
NMM_EXPORT void                 nmm_io_destroy(struct nmm_io const* io);
NMM_EXPORT int                  nmm_io_write(struct nmm_io const* io, FILE* stream);
NMM_EXPORT struct nmm_io const* nmm_io_derived(struct imm_io const* io);

#endif

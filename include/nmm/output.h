#ifndef NMM_OUTPUT_H
#define NMM_OUTPUT_H

#include "nmm/export.h"
#include <stdio.h>

struct nmm_model;
struct nmm_output;

NMM_API int                nmm_output_close(struct nmm_output* output);
NMM_API struct nmm_output* nmm_output_create(char const* filepath);
NMM_API int                nmm_output_destroy(struct nmm_output* output);
NMM_API int nmm_output_write(struct nmm_output* output, struct nmm_model const* model);

#endif
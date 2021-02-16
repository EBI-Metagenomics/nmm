#ifndef NMM_INPUT_H
#define NMM_INPUT_H

#include "nmm/export.h"
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>

struct nmm_model;

struct nmm_input;

NMM_API int                     nmm_input_close(struct nmm_input* input);
NMM_API struct nmm_input*       nmm_input_create(char const* filepath);
NMM_API int                     nmm_input_destroy(struct nmm_input* input);
NMM_API bool                    nmm_input_eof(struct nmm_input const* input);
NMM_API int                     nmm_input_fseek(struct nmm_input* input, int64_t offset);
NMM_API int64_t                 nmm_input_ftell(struct nmm_input* input);
NMM_API struct nmm_model const* nmm_input_read(struct nmm_input* input);

#endif

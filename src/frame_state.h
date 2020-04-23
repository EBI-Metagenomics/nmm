#ifndef FRAME_STATE_H
#define FRAME_STATE_H

#include <stdio.h>

struct imm_state;

struct imm_state const* frame_state_read(FILE* stream);

#endif

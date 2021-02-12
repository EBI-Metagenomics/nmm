#include "nmm/input.h"
#include "free.h"
#include "imm/imm.h"
#include "model.h"
#include "nmm/io.h"
#include "nmm/model.h"
#include <stdlib.h>
#include <string.h>

struct nmm_input
{
    FILE*       stream;
    char const* filepath;

    bool eof;
    bool closed;
};

static struct nmm_model const* read_block(struct nmm_input* input, uint8_t block_type);

int nmm_input_close(struct nmm_input* input)
{
    if (input->closed)
        return 0;

    int errno = 0;

    if (fclose(input->stream)) {
        imm_error("failed to close file %s", input->filepath);
        errno = 1;
    }

    input->closed = true;
    return errno;
}

struct nmm_input* nmm_input_create(char const* filepath)
{
    FILE* stream = fopen(filepath, "r");
    if (!stream) {
        imm_error("could not open file %s for reading", filepath);
        return NULL;
    }

    struct nmm_input* input = malloc(sizeof(*input));
    input->stream = stream;
    input->filepath = strdup(filepath);
    input->eof = false;
    input->closed = false;

    return input;
}

int nmm_input_destroy(struct nmm_input* input)
{
    int errno = nmm_input_close(input);
    free_c(input->filepath);
    free_c(input);
    return errno;
}

bool nmm_input_eof(struct nmm_input const* input) { return input->eof; }

int nmm_input_fseek(struct nmm_input* input, long offset)
{
    return fseek(input->stream, offset, SEEK_SET);
}

long nmm_input_ftell(struct nmm_input* input) { return ftell(input->stream); }

struct nmm_model const* nmm_input_read(struct nmm_input* input)
{
    uint8_t block_type = 0x00;

    if (fread(&block_type, sizeof(block_type), 1, input->stream) < 1) {
        imm_error("could not read block type");
        return NULL;
    }

    return read_block(input, block_type);
}

static struct nmm_model const* read_block(struct nmm_input* input, uint8_t block_type)
{
    if (block_type == IMM_IO_BLOCK_EOF) {
        input->eof = true;
        return NULL;
    }

    if (block_type != NMM_IO_BLOCK_MODEL) {
        imm_error("unknown block type");
        return NULL;
    }

    struct nmm_model const* model = NULL;
    if (!(model = nmm_model_read(input->stream))) {
        imm_error("failed to read file %s", input->filepath);
        return NULL;
    }
    return model;
}

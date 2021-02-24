#include "free.h"
#include "nmm/nmm.h"
#include "profile.h"
#include <stdlib.h>
#include <string.h>

struct nmm_input
{
    FILE*       stream;
    bool        own_stream;
    char const* filepath;
    bool        eof;
};

static struct nmm_input*         input_screate(char const* filepath, FILE* restrict stream, bool own_stream);
static struct nmm_profile const* read_block(struct nmm_input* input, uint8_t block_type);

int nmm_input_close(struct nmm_input* input)
{
    if (!input->stream)
        return 0;

    if (fclose(input->stream)) {
        imm_error("failed to close file %s", input->filepath);
        input->stream = NULL;
        return 1;
    }

    input->stream = NULL;
    return 0;
}

struct nmm_input* nmm_input_create(char const* filepath)
{
    FILE* stream = fopen(filepath, "r");
    if (!stream) {
        imm_error("could not open file %s for reading", filepath);
        return NULL;
    }
    return input_screate(filepath, stream, true);
}

int nmm_input_destroy(struct nmm_input* input)
{
    int errno = 0;
    if (input->own_stream)
        errno = nmm_input_close(input);
    free_c(input->filepath);
    free_c(input);
    return errno;
}

bool nmm_input_eof(struct nmm_input const* input) { return input->eof; }

int nmm_input_fseek(struct nmm_input* input, int64_t offset) { return imm_file_seek(input->stream, offset, SEEK_SET); }

int64_t nmm_input_ftell(struct nmm_input* input) { return imm_file_tell(input->stream); }

struct nmm_profile const* nmm_input_read(struct nmm_input* input)
{
    uint8_t block_type = 0x00;

    if (fread(&block_type, sizeof(block_type), 1, input->stream) < 1) {
        imm_error("could not read block type");
        return NULL;
    }

    return read_block(input, block_type);
}

static struct nmm_profile const* read_block(struct nmm_input* input, uint8_t block_type)
{
    if (block_type == IMM_IO_BLOCK_EOF) {
        input->eof = true;
        return NULL;
    }

    if (block_type != NMM_IO_BLOCK_MODEL) {
        imm_error("unknown block type");
        return NULL;
    }

    struct nmm_profile const* prof = NULL;
    if (!(prof = nmm_profile_read(input->stream))) {
        imm_error("failed to read file %s", input->filepath);
        return NULL;
    }
    return prof;
}

struct nmm_input* nmm_input_screate(char const* filepath, FILE* restrict stream)
{
    return input_screate(filepath, stream, false);
}

static struct nmm_input* input_screate(char const* filepath, FILE* restrict stream, bool own_stream)
{
    struct nmm_input* input = malloc(sizeof(*input));
    input->stream = stream;
    input->own_stream = own_stream;
    input->filepath = strdup(filepath);
    input->eof = false;
    return input;
}

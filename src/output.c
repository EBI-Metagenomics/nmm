#include "free.h"
#include "nmm/nmm.h"
#include "profile.h"
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

struct nmm_output
{
    FILE*       stream;
    FILE*       stream_idx;
    char const* filepath;
    char const* filepath_idx;
    bool        closed;
};

int nmm_output_close(struct nmm_output* output)
{
    if (output->closed)
        return 0;

    uint8_t block_type = IMM_IO_BLOCK_EOF;
    int     errno = 0;

    if (fwrite(&block_type, sizeof(block_type), 1, output->stream) < 1) {
        imm_error("could not write block type");
        errno = 1;
    }

    if (fclose(output->stream)) {
        imm_error("failed to close file %s", output->filepath);
        errno = 1;
    }

    if (fclose(output->stream_idx)) {
        imm_error("failed to close file %s", output->filepath_idx);
        errno = 1;
    }

    output->closed = true;
    return errno;
}

struct nmm_output* nmm_output_create(char const* filepath)
{
    FILE* stream = fopen(filepath, "w");
    if (!stream) {
        imm_error("could not open file %s for writing", filepath);
        return NULL;
    }

    char* filepath_idx = malloc(sizeof(char) * (strlen(filepath) + 5));
    strcpy(filepath_idx, filepath);
    strcat(filepath_idx, ".idx");
    FILE* stream_idx = fopen(filepath_idx, "w");
    if (!stream_idx) {
        imm_error("could not open file %s for writing", filepath_idx);
        free_c(filepath_idx);
        return NULL;
    }

    struct nmm_output* output = malloc(sizeof(*output));
    output->stream = stream;
    output->stream_idx = stream_idx;
    output->filepath = strdup(filepath);
    output->filepath_idx = filepath_idx;
    output->closed = false;

    return output;
}

int nmm_output_destroy(struct nmm_output* output)
{
    int errno = nmm_output_close(output);
    free_c(output->filepath);
    free_c(output->filepath_idx);
    free_c(output);
    return errno;
}

int64_t nmm_output_ftell(struct nmm_output* output) { return imm_file_tell(output->stream); }

int nmm_output_write(struct nmm_output* output, struct nmm_profile const* prof)
{
    uint8_t block_type = NMM_IO_BLOCK_MODEL;

    long offset = imm_file_tell(output->stream);
    if (offset < -1) {
        imm_error("could not ftell");
        return 1;
    }

    if (fprintf(output->stream_idx, "%ld\n", offset) < 0) {
        imm_error("could not write offset");
        return 1;
    }

    if (fwrite(&block_type, sizeof(block_type), 1, output->stream) < 1) {
        imm_error("could not write block type");
        return 1;
    }

    if (nmm_profile_write(prof, output->stream)) {
        imm_error("could not write nmm profile");
        return 1;
    }
    return 0;
}

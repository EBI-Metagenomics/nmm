#include "nmm/output.h"
#include "free.h"
#include "imm/imm.h"
#include "model.h"
#include "nmm/io.h"
#include "nmm/model.h"
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

struct nmm_output
{
    FILE*       stream;
    char const* filepath;
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

    struct nmm_output* output = malloc(sizeof(*output));
    output->stream = stream;
    output->filepath = strdup(filepath);
    output->closed = false;

    return output;
}

int nmm_output_destroy(struct nmm_output* output)
{
    int errno = nmm_output_close(output);
    free_c(output->filepath);
    free_c(output);
    return errno;
}

int nmm_output_write(struct nmm_output* output, struct nmm_model const* model)
{
    uint8_t block_type = NMM_IO_BLOCK_MODEL;

    if (fwrite(&block_type, sizeof(block_type), 1, output->stream) < 1) {
        imm_error("could not write block type");
        return 1;
    }

    if (nmm_model_write(model, output->stream)) {
        imm_error("could not write nmm model");
        return 1;
    }
    return 0;
}

#include "nmm/output.h"
#include "free.h"
#include "imm/imm.h"
#include "model.h"
#include "nmm/io.h"
#include "nmm/model.h"
#include <stdlib.h>
#include <string.h>

struct nmm_output
{
    struct imm_output* super;
};

struct nmm_output* nmm_output_create(char const* filepath)
{
    struct nmm_output* output = malloc(sizeof(*output));

    if (!(output->super = imm_output_create(filepath))) {
        free_c(output);
        return NULL;
    }
    return output;
}

int nmm_output_destroy(struct nmm_output const* output)
{
    int errno = imm_output_destroy(output->super);
    free_c(output);
    return errno;
}

int nmm_output_write(struct nmm_output* output, struct nmm_model const* model)
{
    uint8_t block_type = NMM_IO_BLOCK_MODEL;

    if (fwrite(&block_type, sizeof(block_type), 1, output->super->stream) < 1) {
        imm_error("could not write block type");
        return 1;
    }

    if (nmm_model_write(model, output->super->stream)) {
        imm_error("could not write nmm model");
        return 1;
    }
    return 0;
}

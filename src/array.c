#include "nmm/array.h"
#include <inttypes.h>
#include <limits.h>

int nmm_array3d_write(struct nmm_array3d const* arr, FILE* stream)
{
    if (arr->strides[0] > UINT16_MAX || arr->strides[1] > UINT16_MAX ||
        arr->strides[2] > UINT16_MAX) {
        imm_error("strides not within the uint16_t range");
        return 1;
    }

    uint16_t const strides[3] = {(uint16_t)arr->strides[0], (uint16_t)arr->strides[1],
                                 (uint16_t)arr->strides[2]};

    if (fwrite(strides, sizeof(*strides), 3, stream) < 3) {
        imm_error("could not write array strides");
        return 1;
    }

    unsigned len = nmm_array3d_length(arr);
    if (fwrite(arr->values, sizeof(*arr->values), len, stream) < len) {
        imm_error("could not write array values");
        return 1;
    }

    return 0;
}

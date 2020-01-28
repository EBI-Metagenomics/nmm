#ifndef NMM_BASE_H
#define NMM_BASE_H

#include "nmm/api.h"

struct imm_abc;
struct nmm_base;

NMM_API struct nmm_base const* nmm_base_create(struct imm_abc const* abc);
NMM_API void                   nmm_base_destroy(struct nmm_base const* base);
NMM_API struct imm_abc const*  nmm_base_get_abc(struct nmm_base const* base);

#endif

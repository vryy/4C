#ifndef POST_FLUID3_IS_H
#define POST_FLUID3_IS_H


#ifdef D_FLUID3_IS

#include "../post_common/post_common.h"
#include "../output/gid.h"

void fluid3_is_write_gauss(GIDSET* gid);

void fluid3_is_write_mesh(FIELD_DATA *field, GIDSET* gid, INT* first_mesh);

#endif
#endif

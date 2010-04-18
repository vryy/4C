#ifndef POST_FLUID2_IS_H
#define POST_FLUID2_IS_H


#ifdef D_FLUID2_IS

#include "../post_common/post_common.h"
#include "../output/gid.h"

void fluid2_is_write_gauss(GIDSET* gid);

void fluid2_is_write_mesh(FIELD_DATA *field, GIDSET* gid, INT* first_mesh);

#endif
#endif

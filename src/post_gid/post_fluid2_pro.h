#ifndef POST_FLUID2_PRO_H
#define POST_FLUID2_PRO_H


#ifdef D_FLUID2_PRO

#include "../post_common/post_common.h"
#include "../output/gid.h"

void fluid2_pro_write_gauss(GIDSET* gid);

void fluid2_pro_write_mesh(FIELD_DATA *field, GIDSET* gid, INT* first_mesh);

#endif
#endif

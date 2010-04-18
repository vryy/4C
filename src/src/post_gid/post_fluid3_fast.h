
#ifndef POST_FLUID3_FAST_H
#define POST_FLUID3_FAST_H

#ifdef D_FLUID3_F

#include "../post_common/post_common.h"
#include "../output/gid.h"

void fluid3_fast_write_domain(FIELD_DATA *field, GIDSET* gid, CHUNK_DATA* chunk);

void fluid3_fast_write_stress(FIELD_DATA *field, GIDSET* gid, CHUNK_DATA* chunk, DOUBLE time, INT step);

void fluid3_fast_write_gauss(GIDSET* gid);

void fluid3_fast_write_mesh(FIELD_DATA *field, GIDSET* gid, INT* first_mesh);


#endif
#endif

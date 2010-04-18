
#ifndef POST_SOLID3_H
#define POST_SOLID3_H

#ifdef D_SOLID3

#include "../post_common/post_common.h"
#include "../output/gid.h"


void solid3_write_domain(FIELD_DATA *field, GIDSET* gid, CHUNK_DATA* chunk);


void solid3_write_stress(FIELD_DATA *field, GIDSET* gid, CHUNK_DATA* chunk, DOUBLE time, INT step);


void solid3_write_gauss(GIDSET* gid);


void solid3_write_mesh(FIELD_DATA *field, GIDSET* gid, INT* first_mesh);


#endif
#endif

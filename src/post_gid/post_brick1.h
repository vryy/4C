
#ifndef POST_BRICK1_H
#define POST_BRICK1_H

#ifdef D_BRICK1

#include "../post_common/post_common.h"
#include "../output/gid.h"


void brick1_write_domain(FIELD_DATA *field, GIDSET* gid, CHUNK_DATA* chunk);


void brick1_write_stress(FIELD_DATA *field, GIDSET* gid, CHUNK_DATA* chunk, DOUBLE time, INT step);


void brick1_write_gauss(GIDSET* gid);


void brick1_write_mesh(FIELD_DATA *field, GIDSET* gid, INT* first_mesh);


#endif
#endif

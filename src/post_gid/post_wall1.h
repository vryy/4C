#ifndef POST_WALL1_H
#define POST_WALL1_H


#ifdef D_WALL1

#include "../post_common/post_common.h"
#include "../output/gid.h"


void wall1_write_domain(FIELD_DATA *field, GIDSET* gid, CHUNK_DATA* chunk);


void wall1_write_stress(FIELD_DATA *field, GIDSET* gid, CHUNK_DATA* chunk, DOUBLE time, INT step);


void wall1_write_gauss(GIDSET* gid);


void wall1_write_mesh(FIELD_DATA *field, GIDSET* gid, INT* first_mesh);


#endif
#endif

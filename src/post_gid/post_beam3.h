#ifndef POST_BEAM3_H
#define POST_BEAM3_H


#ifdef D_BEAM3

#include "../post_common/post_common.h"
#include "../output/gid.h"

void beam3_write_domain(FIELD_DATA *field, GIDSET* gid, CHUNK_DATA* chunk);


void beam3_write_stress(FIELD_DATA *field, GIDSET* gid, CHUNK_DATA* chunk, DOUBLE time, INT step);


void beam3_write_gauss(GIDSET* gid);


void beam3_write_mesh(FIELD_DATA *field, GIDSET* gid, INT* first_mesh);

#endif
#endif

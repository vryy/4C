#ifndef POST_INTERF_H
#define POST_INTERF_H

#ifdef D_INTERF

#include "../post_common/post_common.h"
#include "../output/gid.h"


void interf_write_stress(FIELD_DATA *field, GIDSET* gid, CHUNK_DATA* chunk, DOUBLE time, INT step);

void interf_write_gauss(GIDSET* gid);

void interf_write_mesh(FIELD_DATA *field, GIDSET* gid, INT* first_mesh);


#endif
#endif

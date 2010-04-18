#ifndef POST_WALLGE_H
#define POST_WALLGE_H


#ifdef D_WALLGE

#include "../post_common/post_common.h"
#include "../output/gid.h"

void wallge_write_stress(FIELD_DATA *field, GIDSET* gid, CHUNK_DATA* chunk, DOUBLE time, INT step);

void wallge_write_gauss(GIDSET* gid);

void wallge_write_mesh(FIELD_DATA *field, GIDSET* gid, INT* first_mesh);

#endif
#endif

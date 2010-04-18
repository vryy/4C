
#ifndef POST_SHELL9_H
#define POST_SHELL9_H

#ifdef D_SHELL9

#include "../post_common/post_common.h"
#include "../output/gid.h"


void shell9_write_domain(FIELD_DATA *field, GIDSET* gid, CHUNK_DATA* chunk);


void shell9_write_stress(FIELD_DATA *field, GIDSET* gid, CHUNK_DATA* chunk, DOUBLE time, INT step);


void shell9_write_gauss(GIDSET* gid);


void shell9_write_mesh(FIELD_DATA *field, GIDSET* gid, INT* first_mesh);


void shell9_write_displacement(FIELD_DATA *field, RESULT_DATA* result);


#endif
#endif

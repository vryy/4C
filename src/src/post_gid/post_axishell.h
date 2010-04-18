#ifndef POST_AXISHELL_H
#define POST_AXISHELL_H

#ifdef D_AXISHELL

#include "../post_common/post_common.h"
#include "../output/gid.h"

void axishell_write_stress(FIELD_DATA *field, GIDSET* gid, CHUNK_DATA* chunk, DOUBLE time, INT step);

void axishell_write_gauss(GIDSET* gid);

void axishell_write_mesh(FIELD_DATA *field, GIDSET* gid, INT* first_mesh);


#endif
#endif

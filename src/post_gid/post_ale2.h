#ifndef POST_ALE2_H
#define POST_ALE2_H

#ifdef D_ALE

#include "../post_common/post_common.h"
#include "../output/gid.h"

void ale2_write_domain(FIELD_DATA *field, GIDSET* gid, CHUNK_DATA* chunk);

void ale2_write_gauss(GIDSET* gid);

void ale2_write_mesh(FIELD_DATA *field, GIDSET* gid, INT* first_mesh);

#endif
#endif


#ifndef POST_SHELL8_H
#define POST_SHELL8_H

#ifdef D_SHELL8

#include "../post_common/post_common.h"
#include "../output/gid.h"


void shell8_write_domain(FIELD_DATA *field, GIDSET* gid, CHUNK_DATA* chunk);


void shell8_write_stress(FIELD_DATA *field, GIDSET* gid, CHUNK_DATA* chunk, DOUBLE time, INT step);


void shell8_write_gauss(GIDSET* gid);


void shell8_write_mesh(FIELD_DATA *field, GIDSET* gid, INT* first_mesh);



/*----------------------------------------------------------------------*/
/*!
  \brief Write displacements for shell8 elements as bricks.

  The normal functions are sufficient to visualize the middle
  layer. However, there is a third dimension.

  \author u.kue
  \date 12/04
*/
/*----------------------------------------------------------------------*/
void shell8_write_displacement(FIELD_DATA *field, RESULT_DATA* result);


#endif
#endif

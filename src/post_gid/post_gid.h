/*!
\file
\brief Postprocessing utility that takes ccarat output and produces a
GiD input file.

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>

Filters like this one are special inhabitants of the ccarat
world. They are always single processor applications yet they share
some code with ccarat and are closely linked to ccarat internals.

The general idea is that we cannot load the whole result data into
memory at once.

\author u.kue
\date 09/04

*/

#ifndef POST_GID_H
#define POST_GID_H

#include "../post_common/post_common.h"
#include "../output/gid.h"


/*----------------------------------------------------------------------*/
/*!
  \brief Convert an array of field local node ids to global node ids
  in GiD (fortran) style.

  This is an utility function called when meshes are written. There we
  read the mesh chunk that contains the field local node ids. But we
  have to give global node ids to GiD. Furthermore the internal ids
  are counted from zero but GiD requires ids to be counted from one.

  \param field      (i) the field we write
  \param local_ids  (i) array of known field local ids
  \param global_ids (o) array of gid ids to be found
  \param numnp      (i) number of node ids in these arrays

  \author u.kue
  \date 12/04
*/
/*----------------------------------------------------------------------*/
void get_gid_node_ids(FIELD_DATA* field, INT* local_ids, INT* global_ids, INT numnp);


/*----------------------------------------------------------------------*/
/*!
  \brief Write coordinates.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
void write_coords(FIELD_DATA* field, GIDSET* gid);


#endif

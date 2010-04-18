/*!
  \brief Simple octtree.

  The octree depends on the node structure.
 */

#ifndef POST_OCTTREE_H
#define POST_OCTTREE_H

#include "post_common.h"

typedef struct _NODELIST
{
  struct _NODE       *node;        /* pointer to the node */
  struct _NODELIST   *next;        /* pointer to the next element in the NODELIST */

} NODELIST;

typedef struct _OCTREE
{
  INT                 layeroct;                   /* number of layer that I am on (numbering starts with 0)  */
  INT                 numnpoct;                   /* number of nodes belonging to me */
  DOUBLE              xoct[6];                    /* coordinte vector of my geometry(x_min, x_max, y_min, y_max, z_min, z_max) */

  struct _NODELIST    *nodelist;                  /* pointer to the list of nodes that belong to me */
  struct _OCTREE      *next[2];                   /* pointers to next octree elements */
} OCTREE;

void post_create_octree(POST_DISCRETIZATION *discret,
                        OCTREE* octree);
void post_destroy_octree(OCTREE* octree);
NODE* post_octree_find_node(OCTREE *octree, DOUBLE x, DOUBLE y, DOUBLE z);
OCTREE* post_octree_insert_node(OCTREE* octree, NODE* n);

/* maximum number of nodes allowed in one partition of the octree */
#define MAXNODEPART     100

void post_octree_partition(OCTREE *octreetemp);

#ifdef D_FSI

void post_fsi_initcoupling(POST_DISCRETIZATION* struct_discret,
                           POST_DISCRETIZATION* fluid_discret,
                           POST_DISCRETIZATION* ale_discret,
                           INT *fluid_ale_connect);


void post_fsi_divide_leaf (OCTREE   *octree,
                          OCTREE   *octree0,
                          OCTREE   *octree1,
                          DOUBLE   xplushalf,
                          int      geo);

void init_post_node_check(OCTREE* octree, POST_DISCRETIZATION* discret);

#endif

#endif

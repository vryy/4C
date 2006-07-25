#ifndef POST_OCTTREE_H
#define POST_OCTTREE_H

#include "post_common.h"

#ifdef D_FSI

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


/* maximum number of nodes allowed in one partition of the octree */
#define MAXNODEPART     100


void post_fsi_create_octree(OCTREE *octreetemp);

void post_fsi_initcoupling(POST_DISCRETIZATION* struct_discret,
                           POST_DISCRETIZATION* fluid_discret,
                           POST_DISCRETIZATION* ale_discret,
                           INT *fluid_ale_connect);

NODE* post_fsi_sort_ale_node(OCTREE *octree,NODE *actanode);

void post_fsi_sort_nodes (OCTREE   *octree,
                          OCTREE   *octree0,
                          OCTREE   *octree1,
                          DOUBLE   xplushalf,
                          int      geo);

void init_post_node_check(OCTREE* octree, POST_DISCRETIZATION* discret);

#endif

#endif

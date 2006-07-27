
#include "post_octtree.h"

#ifdef D_FSI

/*-----------------------------------------------------------------------*/
/*!
  \brief Split the octree leaf

  \param octree    (i) the current leaf
  \param lhs       (o) new left leaf
  \param rhs       (o) new right leaf
  \param xplushalf (i) dividing line
  \param direction (i) direction in which split occures

*/
/*-----------------------------------------------------------------------*/
void post_fsi_divide_leaf(OCTREE   *octree,
                          OCTREE   *lhs,
                          OCTREE   *rhs,
                          DOUBLE   xplushalf,
                          int      direction)
{
  NODELIST  *nodelist;
  NODELIST  *nodelist0;
  NODELIST  *nodelist1;
  INT       i;

#ifdef DEBUG
  dstrc_enter("post_fsi_divide_leaf");
#endif

  /*initialisation*/
  nodelist=octree->nodelist;
  nodelist0=lhs->nodelist;
  nodelist1=rhs->nodelist;

  for (i=0;i<(octree->numnpoct);i++)
  {
    /* plausibility check */
    /*dsassert(octree->nodelist==NULL,"No nodes in this octree partition\n");*/

    /*if node lies in right octree partition*/
    if ((nodelist->node->x[direction])>xplushalf)
    {
      if (rhs->numnpoct==0)    /*if adding first node*/
      {
        rhs->nodelist=nodelist;
        nodelist1=nodelist;
        rhs->numnpoct++;
      }
      else                       /*if adding any other node but the first one*/
      {
        nodelist1->next=nodelist;
        nodelist1=nodelist1->next;
        rhs->numnpoct++;
      }
    }

    /*if node lies in left octree partition*/
    else
    {
      if (lhs->numnpoct==0)    /*if adding first node*/
      {
        lhs->nodelist=nodelist;
        nodelist0=nodelist;
        lhs->numnpoct++;
      }
      else                        /*if adding any other node but the first one*/
      {
        nodelist0->next=nodelist;
        nodelist0=nodelist0->next;
        lhs->numnpoct++;
      }
    }
    nodelist=nodelist->next;  /*next node of the octree*/
  }
  if (nodelist0!=NULL) nodelist0->next=NULL;
  if (nodelist1!=NULL) nodelist1->next=NULL;
  octree->nodelist=NULL;
  octree->numnpoct=0;

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*-----------------------------------------------------------------------*/
/*!
  \brief create the octree structure

  recursive function to subdivide the partitions of the octree
  structure if there are still too many nodes in this partition
*/
/*-----------------------------------------------------------------------*/
void post_octree_partition(OCTREE *octreetemp)
{
  OCTREE                *octreemem;               /*temporary pointer for memory allocation*/
  DOUBLE                x1, x2, x3, xplushalf;    /*length of the edges of the octree partition, half the length of the edge*/
  INT                   geo;                      /* flag (0,1,2)*/

#ifdef DEBUG
  dstrc_enter("post_octree_partition");
#endif

  octreemem=NULL;
  if (octreetemp->numnpoct>MAXNODEPART)
  {
    /*Initialisation of the new octree elements*/
    octreemem=(OCTREE*)CCACALLOC(2,sizeof(OCTREE));

    octreemem[0].layeroct=(octreetemp->layeroct)+1;
    octreemem[1].layeroct=(octreetemp->layeroct)+1;
    octreemem[0].numnpoct=0;
    octreemem[1].numnpoct=0;

    /*Calculation of the coordinates of the octree partitions */
    x1=(octreetemp->xoct[1])-(octreetemp->xoct[0]);
    x2=(octreetemp->xoct[3])-(octreetemp->xoct[2]);
    x3=(octreetemp->xoct[5])-(octreetemp->xoct[4]);

    if (x1>=x2 && x1>=x3)    /*x-edge is the longest*/
    {
      geo=0;
      xplushalf=(octreetemp->xoct[0])+x1*.5;
      /*coordinates the edge that is  divided*/
      octreemem[0].xoct[0]=octreetemp->xoct[0];
      octreemem[0].xoct[1]=xplushalf;
      octreemem[1].xoct[0]=xplushalf;
      octreemem[1].xoct[1]=octreetemp->xoct[1];
      /*coordinates of edges that are not divided*/
      octreemem[0].xoct[2]=octreetemp->xoct[2];
      octreemem[0].xoct[3]=octreetemp->xoct[3];
      octreemem[0].xoct[4]=octreetemp->xoct[4];
      octreemem[0].xoct[5]=octreetemp->xoct[5];
      octreemem[1].xoct[2]=octreetemp->xoct[2];
      octreemem[1].xoct[3]=octreetemp->xoct[3];
      octreemem[1].xoct[4]=octreetemp->xoct[4];
      octreemem[1].xoct[5]=octreetemp->xoct[5];

      post_fsi_divide_leaf(octreetemp,&octreemem[0],&octreemem[1], xplushalf, geo);
    }
    else if (x2>x1 && x2>=x3)     /*y-edge is the longest*/
    {
      geo=1;
      xplushalf=(octreetemp->xoct[2])+x2*.5;
      /*coordinates the edge that is  divided*/
      octreemem[0].xoct[2]=octreetemp->xoct[2];
      octreemem[0].xoct[3]=xplushalf;
      octreemem[1].xoct[2]=xplushalf;
      octreemem[1].xoct[3]=octreetemp->xoct[3];
      /*coordinates of edges that are not divided*/
      octreemem[0].xoct[0]=octreetemp->xoct[0];
      octreemem[0].xoct[1]=octreetemp->xoct[1];
      octreemem[0].xoct[4]=octreetemp->xoct[4];
      octreemem[0].xoct[5]=octreetemp->xoct[5];
      octreemem[1].xoct[0]=octreetemp->xoct[0];
      octreemem[1].xoct[1]=octreetemp->xoct[1];
      octreemem[1].xoct[4]=octreetemp->xoct[4];
      octreemem[1].xoct[5]=octreetemp->xoct[5];

      post_fsi_divide_leaf(octreetemp, &octreemem[0], &octreemem[1], xplushalf, geo);
    }
    else if (x3>x1 && x3>x2)      /*z-edge is the longest*/
    {
      geo=2;
      xplushalf=(octreetemp->xoct[4])+x3*.5;
      /*coordinates the edge that is  divided*/
      octreemem[0].xoct[4]=octreetemp->xoct[4];
      octreemem[0].xoct[5]=xplushalf;
      octreemem[1].xoct[4]=xplushalf;
      octreemem[1].xoct[5]=octreetemp->xoct[5];
      /*coordinates of edges that are not divided*/
      octreemem[0].xoct[0]=octreetemp->xoct[0];
      octreemem[0].xoct[1]=octreetemp->xoct[1];
      octreemem[0].xoct[2]=octreetemp->xoct[2];
      octreemem[0].xoct[3]=octreetemp->xoct[3];
      octreemem[1].xoct[0]=octreetemp->xoct[0];
      octreemem[1].xoct[1]=octreetemp->xoct[1];
      octreemem[1].xoct[2]=octreetemp->xoct[2];
      octreemem[1].xoct[3]=octreetemp->xoct[3];

      post_fsi_divide_leaf(octreetemp, &octreemem[0], &octreemem[1], xplushalf, geo);
    }
    else
    {
      dserror("Impossible error!");
    }
    /*adding the new partitions to the octree structure*/
    octreetemp->next[0]=&octreemem[0];
    octreetemp->next[1]=&octreemem[1];

    /*recursive call of the fsi_create_octree function*/
    if (octreemem!=NULL)
    {
      post_octree_partition((octreetemp->next[0]));
      post_octree_partition((octreetemp->next[1]));
    }
  }
  /*else: no further partition required*/

#ifdef DEBUG
  dstrc_exit();
#endif
}


static void post_octree_leaf_add_node(OCTREE *octree, NODE *actnode)
{
  NODELIST *nodelistmem;
#ifdef DEBUG
  dstrc_enter("post_octree_leaf_add_node");
#endif

  octree->xoct[0] = MIN(actnode->x[0], octree->xoct[0]);
  octree->xoct[1] = MAX(actnode->x[0], octree->xoct[1]);
  octree->xoct[2] = MIN(actnode->x[1], octree->xoct[2]);
  octree->xoct[3] = MAX(actnode->x[1], octree->xoct[3]);
  octree->xoct[4] = MIN(actnode->x[2], octree->xoct[4]);
  octree->xoct[5] = MAX(actnode->x[2], octree->xoct[5]);

  /*Create the initial nodelist containing all fluid nodes*/
  nodelistmem=(NODELIST*)CCACALLOC(1,sizeof(NODELIST));
  nodelistmem->node=actnode;
  nodelistmem->next=NULL;

  if (octree->numnpoct==0)
  {
    /* if adding first node */
    octree->nodelist = nodelistmem;
  }
  else
  {
    /* if adding any other node but the first one */
    nodelistmem->next = octree->nodelist;
    octree->nodelist = nodelistmem;
  }
  octree->numnpoct++;

#ifdef DEBUG
  dstrc_exit();
#endif
}


void post_create_octree(POST_DISCRETIZATION *discret,
                        OCTREE* octree)
{
  INT i;
  NODE *actnode;
  NODELIST *nodelistmem;

#ifdef DEBUG
  dstrc_enter("post_create_octree");
#endif

  /*Initialize the octree structure*/
  octree->layeroct = 0;
  octree->numnpoct = 0;
  /* the initial values of the octree coordinates are arbitrarily set
   * to those of the first fluid node*/
  actnode=&(discret->node[0]);
  octree->xoct[0]=actnode->x[0];
  octree->xoct[1]=actnode->x[0];
  octree->xoct[2]=actnode->x[1];
  octree->xoct[3]=actnode->x[1];
  octree->xoct[4]=actnode->x[2];
  octree->xoct[5]=actnode->x[2];
  octree->nodelist=NULL;
  octree->next[0]=NULL;
  octree->next[1]=NULL;

  /* loop all fluid nodes to look for size of the FSI interface */
  for (i=0;i<discret->field->numnp;i++)
  {
    actnode = &(discret->node[i]);
    post_octree_leaf_add_node(octree, actnode);
  }
  /*create the octree*/
  post_octree_partition(octree);

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
void post_destroy_octree(OCTREE* octree)
{
  NODELIST* nodelist;

#ifdef DEBUG
  dstrc_enter("post_destroy_octree");
#endif

  nodelist = octree->nodelist;
  while (nodelist!=NULL)
  {
    NODELIST* n;
    n = nodelist;
    nodelist = nodelist->next;
    CCAFREE(n);
  }

  if (octree->next[0]!=NULL)
  {
    post_destroy_octree(octree->next[0]);
    post_destroy_octree(octree->next[1]);

    CCAFREE(octree->next[0]);
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*-----------------------------------------------------------------------*/
/*!
  \brief find the node at the specified position (if any)
*/
/*-----------------------------------------------------------------------*/
NODE* post_octree_find_node(OCTREE *octree,DOUBLE x, DOUBLE y, DOUBLE z)
{
  NODE   *return_node;

#ifdef DEBUG
  dstrc_enter("post_octree_find_node");
#endif

  return_node=NULL;

  /* see if the node could belong to me (with some tolerance) */
  if ((octree->xoct[0]-1e-6 <= x) && (x <= octree->xoct[1]+1e-6) &&
      (octree->xoct[2]-1e-6 <= y) && (y <= octree->xoct[3]+1e-6) &&
      (octree->xoct[4]-1e-6 <= z) && (z <= octree->xoct[5]+1e-6))
  {
    if (octree->nodelist==NULL)
    {
      OCTREE *lhs;
      OCTREE *rhs;

      lhs = octree->next[0];
      rhs = octree->next[1];

      /*if fluid node is in the left octree partition*/
      if (lhs!=NULL)
        return_node = post_octree_find_node(lhs, x, y, z);
      if (return_node==NULL)
        if (rhs!=NULL)
          return_node = post_octree_find_node(rhs, x, y, z);
    }
    else
    {
      INT i;
      NODE* n;
      NODELIST* nl;
      nl = octree->nodelist;
      for (i=0; i<octree->numnpoct; ++i)
      {
        n = nl->node;
        if (fabs(n->x[0]-x) < 1e-6 &&
            fabs(n->x[1]-y) < 1e-6 &&
            fabs(n->x[2]-z) < 1e-6)
        {
          return_node = n;
          break;
        }

        nl = nl->next;
      }
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return return_node;
}


/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
OCTREE* post_octree_insert_node(OCTREE* octree, NODE* n)
{
  OCTREE *leaf;

#ifdef DEBUG
  dstrc_enter("post_octree_insert_node");
#endif

  leaf = NULL;

  /* see if the node could belong to me (with some tolerance) */
  if ((octree->xoct[0] <= n->x[0]) && (n->x[0] <= octree->xoct[1]) &&
      (octree->xoct[2] <= n->x[1]) && (n->x[1] <= octree->xoct[3]) &&
      (octree->xoct[4] <= n->x[2]) && (n->x[2] <= octree->xoct[5]))
  {
    if (octree->nodelist==NULL)
    {
      OCTREE *lhs;
      OCTREE *rhs;

      lhs = octree->next[0];
      rhs = octree->next[1];

      /*if fluid node is in the left octree partition*/
      if (lhs!=NULL)
        leaf = post_octree_insert_node(lhs, n);
      if (leaf==NULL)
        if (rhs!=NULL)
          leaf = post_octree_insert_node(rhs, n);
    }
    else
    {
      leaf = octree;
      post_octree_leaf_add_node(octree, n);
      post_octree_partition(octree);
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return leaf;
}


/*-----------------------------------------------------------------------*/
/*!
  \brief initialise fsi coupling conditions

  \param *structfield   FIELD         (i)      structure field
  \param *fluidfield    FIELD         (i)      fluid field
  \param *alefield      FIELD         (i)      ale field

 */
/*-----------------------------------------------------------------------*/
void post_fsi_initcoupling(POST_DISCRETIZATION* struct_discret,
                           POST_DISCRETIZATION* fluid_discret,
                           POST_DISCRETIZATION* ale_discret,
                           INT *fluid_ale_connect)
{
  INT     numfnp, numsnp, numanp;              /* number of nodes         */
  INT     i;                                   /* simply some counters    */

  NODE   *actfnode=NULL,  *actanode=NULL;  /* actual nodes */
  OCTREE octree;                            /* head of the octree structure */
  NODELIST *nodelisttemp,*nodelistmem;        /* temporary pointers */

#ifdef DEBUG
  dstrc_enter("post_fsi_initcoupling");
#endif

  /* find number of nodes in different fields */
  numsnp  = struct_discret->field->numnp;
  numfnp  = fluid_discret->field->numnp;
  numanp  = ale_discret->field->numnp;

  /*initialization of the fluid_ale_connect array*/
  for (i=0;i<numfnp;i++)
  {
    fluid_ale_connect[i]=-1;
  }

  post_create_octree(fluid_discret, &octree);

  /*find the corresponding fluid node to every ale node using the octree*/
  for (i=0;i<numanp;i++)
  {
    actanode=&(ale_discret->node[i]);
    actfnode=post_octree_find_node(&octree, actanode->x[0], actanode->x[1], actanode->x[2]);
    if (actfnode!=NULL)
    {
      fluid_ale_connect[actfnode->Id_loc]=i;
    }
  }

  post_destroy_octree(&octree);

#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}

#endif

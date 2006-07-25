
#include "post_octtree.h"

#ifdef D_FSI


/*-----------------------------------------------------------------------*/
/*!
  \brief sort the FSI Interface structure nodes into the different
         partitions of the octree structure

  -sorting the structure nodes into the octree partitions
  -no tolerance implemented yet
*/
/*-----------------------------------------------------------------------*/
void post_fsi_sort_nodes (OCTREE   *octree,
                          OCTREE   *octree0,
                          OCTREE   *octree1,
                          DOUBLE   xplushalf,
                          int      geo)
{
  NODELIST  *nodelist;
  NODELIST  *nodelist0;
  NODELIST  *nodelist1;
  INT       i;

  /*initialisation*/
  nodelist=octree->nodelist;
  nodelist0=octree0->nodelist;
  nodelist1=octree1->nodelist;

  for (i=0;i<(octree->numnpoct);i++)
  {
    /* plausibility check */
    /*dsassert(octree->nodelist==NULL,"No nodes in this octree partition\n");*/

    /*if node lies in right octree partition*/
    if ((nodelist->node->x[geo])>xplushalf)
    {
      if (octree1->numnpoct==0)    /*if adding first node*/
      {
        octree1->nodelist=nodelist;
        nodelist1=nodelist;
        octree1->numnpoct++;
      }
      else                       /*if adding any other node but the first one*/
      {
        nodelist1->next=nodelist;
        nodelist1=nodelist1->next;
        octree1->numnpoct++;
      }
    }

    /*if node lies in left octree partition*/
    else
    {
      if (octree0->numnpoct==0)    /*if adding first node*/
      {
        octree0->nodelist=nodelist;
        nodelist0=nodelist;
        octree0->numnpoct++;
      }
      else                        /*if adding any other node but the first one*/
      {
        nodelist0->next=nodelist;
        nodelist0=nodelist0->next;
        octree0->numnpoct++;
      }
    }
    nodelist=nodelist->next;  /*next node of the octree*/
  }
  if (nodelist0!=NULL) nodelist0->next=NULL;
  if (nodelist1!=NULL) nodelist1->next=NULL;
  octree->nodelist=NULL;
  octree->numnpoct=0;
}


/*-----------------------------------------------------------------------*/
/*!
  \brief sort the FSI Interface fluid nodes into the different
         partitions of the octree structure

  sorting the fluid nodes into the octree partitions
*/
/*-----------------------------------------------------------------------*/
NODE* post_fsi_sort_ale_node(OCTREE *octree,NODE *actanode)
{
  NODE   *return_node;

#ifdef DEBUG
  dstrc_enter("post_fsi_sort_ale_node");
#endif

  return_node=NULL;

  /* see if the node could belong to me (with some tolerance) */
  if ((octree->xoct[0]-1e-6 <= actanode->x[0]) && (actanode->x[0] <= octree->xoct[1]+1e-6) &&
      (octree->xoct[2]-1e-6 <= actanode->x[1]) && (actanode->x[1] <= octree->xoct[3]+1e-6) &&
      (octree->xoct[4]-1e-6 <= actanode->x[2]) && (actanode->x[2] <= octree->xoct[5]+1e-6))
  {
    if (octree->nodelist==NULL)
    {
      OCTREE *lhs;
      OCTREE *rhs;

      lhs = octree->next[0];
      rhs = octree->next[1];

      /*if fluid node is in the left octree partition*/
      if (lhs!=NULL)
        return_node = post_fsi_sort_ale_node(lhs, actanode);
      if (return_node==NULL)
        if (rhs!=NULL)
          return_node = post_fsi_sort_ale_node(rhs, actanode);
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
        if (fabs(n->x[0]-actanode->x[0]) < 1e-6 &&
            fabs(n->x[1]-actanode->x[1]) < 1e-6 &&
            fabs(n->x[2]-actanode->x[2]) < 1e-6)
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
/*!
  \brief create the octree structure for the FSI interface nodes

  recursive function to subdivide the partitions of the octree structure if
  there are still too many nodes in this partition
*/
/*-----------------------------------------------------------------------*/
void post_fsi_create_octree (OCTREE *octreetemp)
{
  OCTREE                *octreemem;               /*temporary pointer for memory allocation*/
  DOUBLE                x1, x2, x3, xplushalf;    /*length of the edges of the octree partition, half the length of the edge*/
  INT                   geo;                      /* flag (0,1,2)*/


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

      post_fsi_sort_nodes(octreetemp,&octreemem[0],&octreemem[1], xplushalf, geo);
    }
    if (x2>x1 && x2>=x3)     /*y-edge is the longest*/
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

      post_fsi_sort_nodes(octreetemp, &octreemem[0], &octreemem[1], xplushalf, geo);
    }
    if (x3>x1 && x3>x2)      /*z-edge is the longest*/
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

      post_fsi_sort_nodes(octreetemp, &octreemem[0], &octreemem[1], xplushalf, geo);
    }
    /*adding the new partitions to the octree structure*/
    octreetemp->next[0]=&octreemem[0];
    octreetemp->next[1]=&octreemem[1];


    /*recursive call of the fsi_create_octree function*/
    if (octreemem!=NULL)
    {
      post_fsi_create_octree((octreetemp->next[0]));

      post_fsi_create_octree((octreetemp->next[1]));

    }
  }
  /*else: no further partition required*/
}

/*-----------------------------------------------------------------------*/

/*-----------------------------------------------------------------------*/
/*!
  \brief initialise fsi coupling conditions

  create mulitfield solution history and set pointers to the corresponding
  nodes of the other fields

  \param *structfield   FIELD         (i)      structure field
  \param *fluidfield    FIELD         (i)      fluid field
  \param *alefield      FIELD         (i)      ale field

  \return void

  \author genk
  \date   09/02

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



/**********************new************************************************************/
  OCTREE                octree;                            /* head of the octree structure */
  NODELIST              *nodelisttemp,*nodelistmem;        /* temporary pointers */

/*************new***********************************************************************/

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

  /*Initialize the octree structure*/
  octree.layeroct = 0;
  octree.numnpoct = 0;
  /*the initial values of the octree coordinates are arbitrarily set to those of the first fluid node*/
  actfnode=&(fluid_discret->node[0]);
  octree.xoct[0]=actfnode->x[0];
  octree.xoct[1]=actfnode->x[0];
  octree.xoct[2]=actfnode->x[1];
  octree.xoct[3]=actfnode->x[1];
  octree.xoct[4]=actfnode->x[2];
  octree.xoct[5]=actfnode->x[2];
  octree.nodelist=NULL;
  octree.next[0]=NULL;
  octree.next[1]=NULL;

  /*loop all fluid nodes to look for size of the FSI interface*/
  for(i=0;i<numfnp;i++)
  {
    actfnode = &(fluid_discret->node[i]);
    if (actfnode->x[0]<octree.xoct[0]) octree.xoct[0]=actfnode->x[0];
    if (actfnode->x[0]>octree.xoct[1]) octree.xoct[1]=actfnode->x[0];
    if (actfnode->x[1]<octree.xoct[2]) octree.xoct[2]=actfnode->x[1];
    if (actfnode->x[1]>octree.xoct[3]) octree.xoct[3]=actfnode->x[1];
    if (actfnode->x[2]<octree.xoct[4]) octree.xoct[4]=actfnode->x[2];
    if (actfnode->x[2]>octree.xoct[5]) octree.xoct[5]=actfnode->x[2];

    /*Create the initial nodelist containing all fluid nodes*/
    nodelistmem=(NODELIST*)CCACALLOC(1,sizeof(NODELIST));
    nodelistmem->node=actfnode;
    nodelistmem->next=NULL;

    if (octree.numnpoct==0) /*if adding first node*/
    {
      octree.nodelist=nodelistmem;
      nodelisttemp=octree.nodelist;
    }
    else          /*if adding any other node but the first one*/
    {
      nodelisttemp->next=nodelistmem;
      nodelisttemp=nodelisttemp->next;
    }
    octree.numnpoct++;
  }
  /*create the octree*/
  post_fsi_create_octree(&octree);

  /*find the corresponding fluid node to every ale node using the octree*/
  for (i=0;i<numanp;i++)
  {
    actanode=&(ale_discret->node[i]);
    actfnode=post_fsi_sort_ale_node(&octree, actanode);
    if (actfnode!=NULL)
    {
      fluid_ale_connect[actfnode->Id_loc]=i;
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}

#endif  /* ifdef D_FSI */

/*-----------------------------------------------------------------------*/
/*!
\file
\brief coupling conditions for fsi problems


<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

 */
/*-----------------------------------------------------------------------*/

#ifndef CCADISCRET

/*!
\addtogroup FSI
*//*! @{ (documentation module open)*/


#ifdef D_FSI
#ifdef FSI_NONMATCH

#include "../headers/standardtypes.h"
#include "../solver/solver.h"
#include "../fluid2/fluid2.h"
#include "../fluid3/fluid3.h"
#include "../fluid3/fluid3_prototypes.h"
#include "fsi_prototypes.h"


/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;


/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | ALLDYNA               *alldyn;                                       |
 *----------------------------------------------------------------------*/
extern ALLDYNA      *alldyn;


/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate design if needed                                 |
 | defined in global_control.c                                          |
 *----------------------------------------------------------------------*/
extern struct _DESIGN *design;

/**************************new******************************************/
/*-------------------------------------------------------------------------------------*
 | struct NODELIST contains all nodes belonging to a specific octree partition (leaf)  |
 *-------------------------------------------------------------------------------------*/
  typedef struct _NODELIST
  {
    struct _NODE       *node;        /* pointer to the node */
    struct _NODELIST   *next;        /* pointer to the next element in the NODELIST */

  } NODELIST;

/*----------------------------------------------------------------------------------------*
 | struct OCTREE contains all information about a specific partition (leaf) of the octree |
 *----------------------------------------------------------------------------------------*/
  typedef struct _OCTREE
  {
    INT                 layeroct;    /* number of layer that I am on (numbering starts with 0)  */
    INT                 numnpoct;    /* number of nodes belonging to me */
    DOUBLE              xoct[6];     /* coordinate vector of my geometry(x_min, x_max, y_min, y_max, z_min, z_max) */

    struct _NODELIST    *nodelist;   /* pointer to the list of nodes belonging to me */
    struct _OCTREE      *next[2];    /* pointers to next 2 octree elements */

  } OCTREE;

#define MAXNODEPART   16      /* maximum number of nodes allowed in one partition of the octree (leaf)*/
#define DISTLIMIT     3       /* maximum distance of a node from an element allowed for pairing the two */
                              /* if the distance is greater -> remesh*/
#define EPS           1e-7    /* tolerance for host element search*/
#define ERROR         1e-6    /* error for host element serach (solving for local coordinates)*/
/****************************new*****************************************/

/*-----------------------------------------------------------------------*/
/*!
  \brief create fsi coupling conditions

  In this routine the fsi-coupling conditions from the input file are
  evaluated and transformed into Neumann conditions for structure
  and into Dirichlet conditions for fluid and ale.

  \return void

  \author genk
  \date   09/02

 */
/*-----------------------------------------------------------------------*/
void fsi_createfsicoup()
{

  INT i,j;                                  /* simply some counters       */
  INT hasdirich,hascouple,hasfsi,hasneum;   /* different flags            */

  FIELDTYP  fieldtyp;

  DNODE    *actdnode;
  DLINE    *actdline;
  DSURF    *actdsurf;


#ifdef DEBUG
  dstrc_enter("fsi_createfsicoup");
#endif


  /* loop dsurfs */
  /* ----------- */
  for (i=0; i<design->ndsurf; i++)
  {
    hasdirich=0;
    hascouple=0;
    hasfsi   =0;
    hasneum  =0;
    actdsurf = &(design->dsurf[i]);


    /* check for conditions */
    if (actdsurf->dirich!=NULL) hasdirich++;
    if (actdsurf->couple!=NULL) hascouple++;
    if (actdsurf->fsicouple!=NULL) hasfsi++;
    if (actdsurf->neum!=NULL) hasneum++;



    /* surface is a fsi coupling surface */
    if (hasfsi!=0)
    {

      fieldtyp=actdsurf->fsicouple->fieldtyp;
      switch (fieldtyp)
      {

        case structure:
          dsassert(hasdirich==0,
              "dirich- and fsi-coupling condition defined on same DSURF\n");
          dsassert(hascouple==0,
              "coupling- and fsi-coupling condition defined on same DSURF\n");

          if (hasneum==0)
          {
            actdsurf->neum = (NEUM_CONDITION*)CCACALLOC(1,sizeof(NEUM_CONDITION));
            actdsurf->neum->neum_type = neum_FSI;
          }
          else if (hasneum>0 && actdsurf->neum->neum_type==neum_live)
          {
            /* combination of live load and FSI load */
            actdsurf->neum->neum_type = neum_live_FSI;
          }
          else if (hasneum>0 && actdsurf->neum->neum_type==neum_orthopressure)
          {
            /* combination of live load and FSI load */
            actdsurf->neum->neum_type = neum_opres_FSI;
          }
          else
            dserror("neumann- and fsi-couping condition defined on same DSRUF\n");
          break;


        case fluid:
          dsassert(hasdirich==0,
              "dirich- and fsi-coupling condition defined on same DSURF\n");
          dsassert(hascouple==0,
              "coupling- and fsi-coupling condition defined on same DSURF\n");
          dsassert(hasneum==0,
              "neumann- and fsi-coupling condition defined on same DSURF\n");


          /* allocate space for a dirichlet condition in this dsurf */
          actdsurf->dirich = (DIRICH_CONDITION*)CCACALLOC(1,sizeof(DIRICH_CONDITION));
          amdef("onoff",&(actdsurf->dirich->dirich_onoff),MAXDOFPERNODE,1,"IV");
          amdef("val",&(actdsurf->dirich->dirich_val),MAXDOFPERNODE,1,"DV");
          amdef("curve",&(actdsurf->dirich->curve),MAXDOFPERNODE,1,"IV");
          amdef("function",&(actdsurf->dirich->funct),MAXDOFPERNODE,1,"IV");
          amzero(&(actdsurf->dirich->dirich_onoff));
          amzero(&(actdsurf->dirich->dirich_val));
          amzero(&(actdsurf->dirich->curve));
          amzero(&(actdsurf->dirich->funct));

          /* initialise for fsi-coupling */
          for (j=0;j<genprob.ndim;j++)
            actdsurf->dirich->dirich_onoff.a.iv[j] = 1;
          actdsurf->dirich->dirich_type=dirich_FSI;


          /* also create the dirichlet condition for the ale */
          if (genprob.create_ale ==1)
          {
            /* allocate space for a dirichlet condition in this dsurf */
            actdsurf->ale_dirich =
              (DIRICH_CONDITION*)CCACALLOC(1,sizeof(DIRICH_CONDITION));
            amdef("onoff",    &(actdsurf->ale_dirich->dirich_onoff),MAXDOFPERNODE, 1,"IV");
            amdef("val",      &(actdsurf->ale_dirich->dirich_val),  MAXDOFPERNODE, 1,"DV");
            amdef("curve",    &(actdsurf->ale_dirich->curve),       MAXDOFPERNODE, 1,"IV");
            amdef("function", &(actdsurf->ale_dirich->funct),       MAXDOFPERNODE, 1,"IV");
            amzero(&(actdsurf->ale_dirich->dirich_onoff));
            amzero(&(actdsurf->ale_dirich->dirich_val));
            amzero(&(actdsurf->ale_dirich->funct));
            amzero(&(actdsurf->ale_dirich->curve));

            /* initialise for fsi-coupling */
            for (j=0;j<genprob.ndim;j++)
              actdsurf->ale_dirich->dirich_onoff.a.iv[j] = 1;
            actdsurf->ale_dirich->dirich_type=dirich_FSI;
          }
          break;


        case ale:
          dsassert(hasdirich==0,
              "dirich- and fsi-coupling condition defined on same DSURF\n");
          dsassert(hascouple==0,
              "coupling- and fsi-coupling condition defined on same DSURF\n");
          dsassert(hasneum==0,
              "neumann- and fsi-coupling condition defined on same DSURF\n");

          /* allocate space for a dirichlet condition in this dsurf */
          actdsurf->dirich = (DIRICH_CONDITION*)CCACALLOC(1,sizeof(DIRICH_CONDITION));
          amdef("onoff",&(actdsurf->dirich->dirich_onoff),MAXDOFPERNODE,1,"IV");
          amdef("val",&(actdsurf->dirich->dirich_val),MAXDOFPERNODE,1,"DV");
          amdef("curve",&(actdsurf->dirich->curve),MAXDOFPERNODE,1,"IV");
          amdef("function",&(actdsurf->dirich->funct),MAXDOFPERNODE,1,"IV");
          amzero(&(actdsurf->dirich->dirich_onoff));
          amzero(&(actdsurf->dirich->dirich_val));
          amzero(&(actdsurf->dirich->funct));
          amzero(&(actdsurf->dirich->curve));

          /* initialise for fsi-coupling */
          for (j=0;j<genprob.ndim;j++)
            actdsurf->dirich->dirich_onoff.a.iv[j] = 1;
          actdsurf->dirich->dirich_type=dirich_FSI;
          break;


        default:
          dserror("fieldtyp unknown!\n");

      } /* switch(fieldtyp) */

    }  /* if (hasfsi!=0) */


    /* surface is not a fsi coupling surface */
    else
    {

#if 1
      /* surface is an exterior surface */
      if (actdsurf->ndvol == 1 )
      {

        /* create dirichlet condition for ale */
        actdsurf->ale_dirich =
          (DIRICH_CONDITION*)CCACALLOC(1,sizeof(DIRICH_CONDITION));
        amdef("onoff",    &(actdsurf->ale_dirich->dirich_onoff),MAXDOFPERNODE, 1,"IV");
        amdef("val",      &(actdsurf->ale_dirich->dirich_val),  MAXDOFPERNODE, 1,"DV");
        amdef("curve",    &(actdsurf->ale_dirich->curve),       MAXDOFPERNODE, 1,"IV");
        amdef("function", &(actdsurf->ale_dirich->funct),       MAXDOFPERNODE, 1,"IV");
        amzero(&(actdsurf->ale_dirich->dirich_onoff));
        amzero(&(actdsurf->ale_dirich->dirich_val));
        amzero(&(actdsurf->ale_dirich->funct));
        amzero(&(actdsurf->ale_dirich->curve));

        /* initialise for ale boundary */
        for (j=0;j<genprob.ndim;j++)
          actdsurf->ale_dirich->dirich_onoff.a.iv[j] = 1;
        actdsurf->ale_dirich->dirich_type = dirich_none;
      }
#endif


    }  /* else */

  }  /* for (i=0; i<design->ndsurf; i++) */



  /* loop dlines */
  /* ----------- */
  for (i=0; i<design->ndline; i++)
  {
    hasdirich=0;
    hascouple=0;
    hasfsi   =0;
    hasneum  =0;
    actdline = &(design->dline[i]);

    /* check for conditions */
    if (actdline->dirich!=NULL) hasdirich++;
    if (actdline->couple!=NULL) hascouple++;
    if (actdline->fsicouple!=NULL) hasfsi++;
    if (actdline->neum!=NULL) hasneum++;


    /* surface is a fsi coupling surface */
    if (hasfsi!=0)
    {

      fieldtyp=actdline->fsicouple->fieldtyp;
      switch (fieldtyp)
      {

        case structure:
          if (hasdirich!=0) dswarning(1,7);
          dsassert(hasdirich==0,
              "dirich- and fsi-coupling condition defined on same DLINE\n");
          dsassert(hascouple==0,
              "coupling- and fsi-coupling condition defined on same DLINE\n");

          actdline->neum = (NEUM_CONDITION*)CCACALLOC(1,sizeof(NEUM_CONDITION));
          if (!actdline->neum) dserror("Allocation of memory failed");
          actdline->neum->neum_type = neum_FSI;
          break;


        case fluid:
          dsassert(hasdirich==0,
              "dirich- and fsi-coupling condition defined on same DLINE\n");
          dsassert(hascouple==0,
              "coupling- and fsi-coupling condition defined on same DLINE\n");
          dsassert(hasneum==0,
              "neumann- and fsi-coupling condition defined on same DLINE\n");

          /* allocate space for a dirichlet condition in this dline */
          actdline->dirich = (DIRICH_CONDITION*)CCACALLOC(1,sizeof(DIRICH_CONDITION));
          amdef("onoff",&(actdline->dirich->dirich_onoff),MAXDOFPERNODE,1,"IV");
          amdef("val",&(actdline->dirich->dirich_val),MAXDOFPERNODE,1,"DV");
          amdef("curve",&(actdline->dirich->curve),MAXDOFPERNODE,1,"IV");
          amdef("function",&(actdline->dirich->funct),MAXDOFPERNODE,1,"IV");
          amzero(&(actdline->dirich->dirich_onoff));
          amzero(&(actdline->dirich->dirich_val));
          amzero(&(actdline->dirich->curve));
          amzero(&(actdline->dirich->funct));

          /* initialise for fsi-coupling */
          for (j=0;j<genprob.ndim;j++)
            actdline->dirich->dirich_onoff.a.iv[j] = 1;
          actdline->dirich->dirich_type=dirich_FSI;


          /* also create the dirichlet condition for the ale */
          if (genprob.create_ale ==1)
          {
            /* allocate space for a dirichlet condition in this dsurf */
            actdline->ale_dirich =
              (DIRICH_CONDITION*)CCACALLOC(1,sizeof(DIRICH_CONDITION));
            amdef("onoff",    &(actdline->ale_dirich->dirich_onoff), MAXDOFPERNODE,1,"IV");
            amdef("val",      &(actdline->ale_dirich->dirich_val),   MAXDOFPERNODE,1,"DV");
            amdef("curve",    &(actdline->ale_dirich->curve),        MAXDOFPERNODE,1,"IV");
            amdef("function", &(actdline->ale_dirich->funct),        MAXDOFPERNODE,1,"IV");
            amzero(&(actdline->ale_dirich->dirich_onoff));
            amzero(&(actdline->ale_dirich->dirich_val));
            amzero(&(actdline->ale_dirich->funct));
            amzero(&(actdline->ale_dirich->curve));

            /* initialise for fsi-coupling */
            for (j=0;j<genprob.ndim;j++)
              actdline->ale_dirich->dirich_onoff.a.iv[j] = 1;
            actdline->ale_dirich->dirich_type=dirich_FSI;
          }
          break;


        case ale:
          dsassert(hasdirich==0,
              "dirich- and fsi-coupling condition defined on same DLINE\n");
          dsassert(hascouple==0,
              "coupling- and fsi-coupling condition defined on same DLINE\n");
          dsassert(hasneum==0,
              "neumann- and fsi-coupling condition defined on same DLINE\n");

          /* allocate space for a dirichlet condition in this dline */
          actdline->dirich = (DIRICH_CONDITION*)CCACALLOC(1,sizeof(DIRICH_CONDITION));
          amdef("onoff",&(actdline->dirich->dirich_onoff),MAXDOFPERNODE,1,"IV");
          amdef("val",&(actdline->dirich->dirich_val),MAXDOFPERNODE,1,"DV");
          amdef("curve",&(actdline->dirich->curve),MAXDOFPERNODE,1,"IV");
          amdef("function",&(actdline->dirich->funct),MAXDOFPERNODE,1,"IV");
          amzero(&(actdline->dirich->dirich_onoff));
          amzero(&(actdline->dirich->dirich_val));
          amzero(&(actdline->dirich->curve));
          amzero(&(actdline->dirich->funct));

          /* initialise for fsi-coupling */
          for (j=0;j<genprob.ndim;j++)
            actdline->dirich->dirich_onoff.a.iv[j] = 1;
          actdline->dirich->dirich_type=dirich_FSI;
          break;


        default:
          dserror("fieldtyp unknown!\n");

      } /* switch(fieldtyp) */

    }  /* if (hasfsi!=0) */


    /* line is not a fsi coupling line */
    else
    {

      /* line is an exterior line */
      if (actdline->ndsurf == 1 )
      {

        /* create dirichlet condition for ale */
        actdline->ale_dirich =
          (DIRICH_CONDITION*)CCACALLOC(1,sizeof(DIRICH_CONDITION));
        amdef("onoff",    &(actdline->ale_dirich->dirich_onoff),MAXDOFPERNODE, 1,"IV");
        amdef("val",      &(actdline->ale_dirich->dirich_val),  MAXDOFPERNODE, 1,"DV");
        amdef("curve",    &(actdline->ale_dirich->curve),       MAXDOFPERNODE, 1,"IV");
        amdef("function", &(actdline->ale_dirich->funct),       MAXDOFPERNODE, 1,"IV");
        amzero(&(actdline->ale_dirich->dirich_onoff));
        amzero(&(actdline->ale_dirich->dirich_val));
        amzero(&(actdline->ale_dirich->funct));
        amzero(&(actdline->ale_dirich->curve));

        /* initialise for ale boundary */
        for (j=0;j<genprob.ndim;j++)
          actdline->ale_dirich->dirich_onoff.a.iv[j] = 1;
        actdline->ale_dirich->dirich_type = dirich_none;
      }


    }  /* else */

  }  /* for (i=0; i<design->ndline; i++) */



  /* loop dnodes */
  /* ----------- */
  for (i=0; i<design->ndnode; i++)
  {
    hasdirich=0;
    hascouple=0;
    hasfsi   =0;
    hasneum  =0;
    actdnode = &(design->dnode[i]);

    /* check for conditions */
    if (actdnode->dirich!=NULL) hasdirich++;
    if (actdnode->couple!=NULL) hascouple++;
    if (actdnode->fsicouple!=NULL) hasfsi++;
    if (actdnode->neum!=NULL) hasneum++;
    if (hasfsi==0) continue;


    fieldtyp=actdnode->fsicouple->fieldtyp;
    switch (fieldtyp)
    {

      case structure:
        dsassert(hasdirich==0,
            "dirich- and fsi-coupling condition defined on same DNODE\n");
        dsassert(hascouple==0,
            "coupling- and fsi-coupling condition defined on same DNODE\n");

        actdnode->neum = (NEUM_CONDITION*)CCACALLOC(1,sizeof(NEUM_CONDITION));
        if (!actdnode->neum) dserror("Allocation of memory failed");
        actdnode->neum->neum_type = neum_FSI;
        break;


      case fluid:
        dsassert(hasdirich==0,
            "dirich- and fsi-coupling condition defined on same DNODE\n");
        dsassert(hascouple==0,
            "coupling- and fsi-coupling condition defined on same DNODE\n");
        dsassert(hasneum==0,
            "neumann- and fsi-coupling condition defined on same DNODE\n");

        /* allocate space for a dirichlet condition in this dnode */
        actdnode->dirich = (DIRICH_CONDITION*)CCACALLOC(1,sizeof(DIRICH_CONDITION));
        amdef("onoff",&(actdnode->dirich->dirich_onoff),MAXDOFPERNODE,1,"IV");
        amdef("val",&(actdnode->dirich->dirich_val),MAXDOFPERNODE,1,"DV");
        amdef("curve",&(actdnode->dirich->curve),MAXDOFPERNODE,1,"IV");
        amdef("function",&(actdnode->dirich->funct),MAXDOFPERNODE,1,"IV");
        amzero(&(actdnode->dirich->dirich_onoff));
        amzero(&(actdnode->dirich->dirich_val));
        amzero(&(actdnode->dirich->curve));
        amzero(&(actdnode->dirich->funct));

        /* initialise for fsi-coupling */
        for (j=0;j<genprob.ndim;j++)
          actdnode->dirich->dirich_onoff.a.iv[j] = 1;
        actdnode->dirich->dirich_type=dirich_FSI;


        /* also create the dirichlet condition for the ale */
        if (genprob.create_ale ==1)
        {
          /* allocate space for a dirichlet condition in this dsurf */
          actdnode->ale_dirich =
            (DIRICH_CONDITION*)CCACALLOC(1,sizeof(DIRICH_CONDITION));
          amdef("onoff",   &(actdnode->ale_dirich->dirich_onoff),MAXDOFPERNODE,1,"IV");
          amdef("val",     &(actdnode->ale_dirich->dirich_val),  MAXDOFPERNODE,1,"DV");
          amdef("curve",   &(actdnode->ale_dirich->curve),       MAXDOFPERNODE,1,"IV");
          amdef("function",&(actdnode->ale_dirich->funct),       MAXDOFPERNODE,1,"IV");
          amzero(&(actdnode->ale_dirich->dirich_onoff));
          amzero(&(actdnode->ale_dirich->dirich_val));
          amzero(&(actdnode->ale_dirich->funct));
          amzero(&(actdnode->ale_dirich->curve));

          /* initialise for fsi-coupling */
          for (j=0;j<genprob.ndim;j++)
            actdnode->ale_dirich->dirich_onoff.a.iv[j] = 1;
          actdnode->ale_dirich->dirich_type=dirich_FSI;
        }
        break;


      case ale:
        dsassert(hasdirich==0,
            "dirich- and fsi-coupling condition defined on same DNODE\n");
        dsassert(hascouple==0,
            "coupling- and fsi-coupling condition defined on same DNODE\n");
        dsassert(hasneum==0,
            "neumann- and fsi-coupling condition defined on same DNODE\n");

        /* allocate space for a dirichlet condition in this dnode */
        actdnode->dirich = (DIRICH_CONDITION*)CCACALLOC(1,sizeof(DIRICH_CONDITION));
        amdef("onoff",&(actdnode->dirich->dirich_onoff),MAXDOFPERNODE,1,"IV");
        amdef("val",&(actdnode->dirich->dirich_val),MAXDOFPERNODE,1,"DV");
        amdef("curve",&(actdnode->dirich->curve),MAXDOFPERNODE,1,"IV");
        amdef("function",&(actdnode->dirich->funct),MAXDOFPERNODE,1,"IV");
        amzero(&(actdnode->dirich->dirich_onoff));
        amzero(&(actdnode->dirich->dirich_val));
        amzero(&(actdnode->dirich->curve));
        amzero(&(actdnode->dirich->funct));

        /* initialise for fsi-coupling */
        for (j=0;j<genprob.ndim;j++)
          actdnode->dirich->dirich_onoff.a.iv[j] = 1;
        actdnode->dirich->dirich_type=dirich_FSI;
        break;


      default:
        dserror("fieldtyp unknown!\n");


    } /* switch(fieldtyp) */

  }  /* for (i=0; i<design->ndnode; i++) */


#ifdef DEBUG
  dstrc_exit();
#endif



  return;
} /* end of fsi_creatcoup*/

/*-----------------------------------------------------------------------*/
/*!
  \brief sort the FSI Interface structure nodes into the different
         partitions of the octree

  The mother partition has to be divided and the two empty children already
  exist. -> Each structure node in the mother has to be sorted into the children
  octree partitions. Two new children nodelists are created by resetting the
  pointers to the nodes and attaching them to the nodelists of the children.

  for example:

  mother: node0->node1->node2->node3->node4
  child1: node0->node2->node3
  child2: node1->node4

  \return void

  \author henke
  \date   05/06
*/
/*-----------------------------------------------------------------------*/
void fsi_sort_struct_nodes (OCTREE   *octree,
                            OCTREE   *octree0,
                            OCTREE   *octree1,
                            DOUBLE   xdivide,
                            int      geo)
{
  NODELIST  *nodelist;
  NODELIST  *nodelist0;
  NODELIST  *nodelist1;
  INT       i;

  /*initialisation*/
  nodelist=octree->nodelist;          /*mother*/
  nodelist0=octree0->nodelist;        /*child1*/
  nodelist1=octree1->nodelist;        /*child2*/

  for (i=0;i<(octree->numnpoct);i++)
  {
    /*if node lies in right octree partition*/
    if ((nodelist->node->x[geo])>xdivide)

      /*if adding first node*/
       if (octree1->numnpoct==0)
       {
         octree1->nodelist=nodelist;
         nodelist1=nodelist;
         octree1->numnpoct++;
       }
       /*if adding any other node but the first one*/
       else
       {
         nodelist1->next=nodelist;
         nodelist1=nodelist1->next;
         octree1->numnpoct++;
       }
    /*if node lies in left octree partition*/
    else

      /*if adding first node*/
       if (octree0->numnpoct==0)
       {
         octree0->nodelist=nodelist;
         nodelist0=nodelist;
         octree0->numnpoct++;
       }
       /*if adding any other node but the first one*/
       else
       {
         nodelist0->next=nodelist;
         nodelist0=nodelist0->next;
         octree0->numnpoct++;
       }
    /*next node of the mother octree partition*/
    nodelist=nodelist->next;
  }
  /*separate the nodelists from each other*/
  if (nodelist0!=NULL) nodelist0->next=NULL;
  if (nodelist1!=NULL) nodelist1->next=NULL;

  /*delete nodelist of mother octree partition*/
  octree->nodelist=NULL;
  octree->numnpoct=0;
}

/*------------------------------------------------------------------------------*/
/*!
  \brief coupling of fluid and structure field
         find a structure element (host element) for every FSI fluid node

   - sort fluid node into the octree
   - find structure host element for fluid node


   The fluid node is sorted into the octree.

   example: MAXNODEPART=2, 5 FSI struct nodes

                                     layer0
                                    /      \
                              layer1        layer1
                             (2 nodes)     /      \
                              (leaf)    layer2   layer2
                                       (1 node) (2 nodes)
                                        (leaf)   (leaf)

   Elements attached to the FSI structure nodes in this octree partition are
   checked for containing this fluid node. If no structure element is found
   that contains the fluid node, e.g. due to curved surfaces or numerical
   inaccuracy, the closest struct element is chosen to be the host element.
   In the end the fluid node will know its host element (see: struct GNODE)
   and the struct element will know the local coordinates (r,s,t) of the
   fluid nodes located inside it (see: struct ELEMENT).

   A Newton-Raphson-Algorithm is used to solve the nonlinear equation system for
   the local coordinates (r,s,t).

   example: Hex8 element

   isoparametric approach for geometric mapping:
      _        _          _
      x = Sum( N(r,s,t) * x )

   -> x1 = c11*r + c12*s + c13*t + c14*r*s + c15*r*t + c16*s*t + c17*r*s*t
      x2 = c21*r + c22*s + c23*t + c24*r*s + c25*r*t + c26*s*t + c27*r*s*t
      x3 = c31*r + c32*s + c33*t + c34*r*s + c35*r*t + c36*s*t + c37*r*s*t
                                                                               _
   The coefficients c are functions of the global nodal coordinates (coeff[][](x)) and
   are explicitly implemented for certain element types. The global coordinates x1,x2,x3
   are given as the nodal coordinates of the fluid node.

   -> f1(r,s,t) = x1 - c11*r + c12*s + c13*t + c14*r*s + c15*r*t + c16*s*t + c17*r*s*t = 0
      f2(r,s,t) = x2 - c21*r + c22*s + c23*t + c24*r*s + c25*r*t + c26*s*t + c27*r*s*t = 0
      f3(r,s,t) = x3 - c31*r + c32*s + c33*t + c34*r*s + c35*r*t + c36*s*t + c37*r*s*t = 0

   This nonlinear equation system is linearized introducing the Jacobian matrix:

   -> f1(r,s,t) = f1(r0,s0,t0) + df1/dr*deltar + df1/ds*deltas + df1/dt*deltat + higher order terms
      f2(r,s,t) = f2(r0,s0,t0) + df2/dr*deltar + df2/ds*deltas + df2/dt*deltat + higher order terms
      f3(r,s,t) = f3(r0,s0,t0) + df3/dr*deltar + df3/ds*deltas + df3/dt*deltat + higher order terms

   _   |df1/dr  df1/ds  df1/dt|
   J = |df2/dr  df2/ds  df2/dt|
       |df3/dr  df3/ds  df3/dt|

   The elements of the Jacobian matrix (jacobi[][]) are explicitely implemented for certain element types.
      _          _             _   _
   -> residual = f(r0,s0,t0) + J * delta

   This linear equation system is solved explicitely for the delta vector. The values for r0,s0,t0 are
   advanced by deltar,deltas,deltat. This iterative procedure is repeated until the residual is smaller
   than a specified error. The values of r0,s0,t0 have now converged to the wanted values r,s,t. The
   position of the fluid node (r,s,t) in relation to the checked struct element is used to decide whether
   the element is the host element or not.

  \return void

  \author henke
  \date   05/06
*/
/*-----------------------------------------------------------------------*/
void fsi_find_hostelement(OCTREE *octree,NODE *actfnode)
{
  OCTREE     *octreetemp;         /*pointer to the current octree partition*/
  NODELIST   *nodelist;           /**/
  NODE       *actsnode;           /*pointer to the current struct node*/
  ELEMENT    *hostele;            /*pointer to the structure host element paired with the fluid node (host element)*/
  DOUBLE     r,s,t;               /*local coordinates of the fluid node inside a struct element*/
  GNODE      *actfgnode;          /*pointer to the GNODE of actfnode*/

  ARRAY      *x;                  /*array for the coordinates of the nodes of the element: x[dimension][number of nodes of the element]*/
  ARRAY      *coeff;              /*array for the coefficients needed to calculate the local coordinates via the shape functions*/
  ARRAY      *jacobi;             /*array for the linearized residuals*/
  INT        i,j,k,m;             /*counters*/

  INT        max=5;                              /*number of fluid nodes in one host element for which space is allocated*/
  DOUBLE     error;                              /*error*/
  DOUBLE     res[3];                             /*residual*/
  DOUBLE     deltar,deltas,deltat;               /*components of delta vevtor*/
  DOUBLE     actdistance=DISTLIMIT, distance=0;  /*distance of a fluid node from a structure element*/
  ELEMENT    *acthostele=NULL;                   /*pointer to the actually closest structure element*/
  DOUBLE     actr=0,acts=0,actt=0;               /*temporary storing of local coordinates of the fluid node inside a struct element*/
  INT        found=0;                            /*flag indicating whether pairing was successfull*/

  actfgnode=actfnode->gnode;
  octreetemp=octree;

  if (actfgnode->fsicouple==NULL)
    dserror("fluid node %d is not an interface node\n", actfnode->Id);

/*printf("\nEinsortieren von Fluid Knoten %d\n", actfnode->Id);*/
/*printf("FLUID node %d: %f %f %f\n",actfnode->Id,actfnode->x[0],actfnode->x[1],actfnode->x[2]);*/
/*printf("Octreelayer: %d\n", octreetemp->layeroct);*/
/*printf("Octreegrenzen:\n");*/
/*printf("xoct[0]: %f\n", octreetemp->xoct[0]);*/
/*printf("xoct[1]: %f\n", octreetemp->xoct[1]);*/
/*printf("xoct[2]: %f\n", octreetemp->xoct[2]);*/
/*printf("xoct[3]: %f\n", octreetemp->xoct[3]);*/
/*printf("xoct[4]: %f\n", octreetemp->xoct[4]);*/
/*printf("xoct[5]: %f\n", octreetemp->xoct[5]);*/

  /*plausibility check: the fluid node has to lie within the boundaries of octree (FSI interface)*/
  if (!((octreetemp->xoct[0]<=actfnode->x[0])&&(actfnode->x[0]<=octreetemp->xoct[1])&&
        (octreetemp->xoct[2]<=actfnode->x[1])&&(actfnode->x[1]<=octreetemp->xoct[3])&&
        (octreetemp->xoct[4]<=actfnode->x[2])&&(actfnode->x[2]<=octreetemp->xoct[5])))
  {
    dserror("fluid node %d does not lie within the boundaries of the octree (FSI interface)\n",actfnode->Id);
  }
  /*********************************/
  /*sort fluid node into the octree*/
  /*********************************/

  /*loop until the lowest level of the octree is reached*/
  while (octreetemp->nodelist==NULL)
  {
    /*if fluid node is in the left octree partition*/
    if ((octreetemp->next[0]->xoct[0]<=actfnode->x[0])&&(actfnode->x[0]<=octreetemp->next[0]->xoct[1])&&
        (octreetemp->next[0]->xoct[2]<=actfnode->x[1])&&(actfnode->x[1]<=octreetemp->next[0]->xoct[3])&&
        (octreetemp->next[0]->xoct[4]<=actfnode->x[2])&&(actfnode->x[2]<=octreetemp->next[0]->xoct[5]))
    {
       /*this way is not heading towards a leaf*/
       if((octreetemp->next[0]->next[0]!=NULL) && (octreetemp->next[0]->next[1]!=NULL))
       {
         octreetemp=octreetemp->next[0];
       }
       /*this way is heading towards a leaf*/
       else
       {
         /*this is a leaf -> done!*/
	 if(octreetemp->next[0]->nodelist!=NULL)
         {
           octreetemp=octreetemp->next[0];
  	 }
         /*this case is theoretically impossible, but - how knows?*/
	 else
	 {
	   dserror("Bad octree!\n");
	   break;
	 }
       }
    }
    /*if fluid node is in the right octree partition*/
    else
    {
      /*this way is not heading towards a leaf*/
      if((octreetemp->next[1]->next[0]!=NULL) && (octreetemp->next[1]->next[1]!=NULL))
      {
	octreetemp=octreetemp->next[1];
      }
      /*this way is heading towards a leaf*/
      else
      {
        /*this is a leaf -> done!*/
	if(octreetemp->next[1]->nodelist!=NULL)
        {
          octreetemp=octreetemp->next[1];
        }
        /*this case is theoretically impossible, but - how knows?*/
        else
	{
	  dserror("Bad octree!\n");
	  break;
	}
      }
    }
  }

/*if (actfnode->Id == 656)*/
/*{*/
/*nodelist0=octreetemp->nodelist;*/
/*printf("%d struct nodes in Octree!\n", octreetemp->numnpoct);*/
/*printf("Grenzen des Octrees\n");*/
/*printf("xunten %f xoben %f\n",octreetemp->xoct[0], octreetemp->xoct[1]);*/
/*printf("yunten %f yoben %f\n",octreetemp->xoct[2], octreetemp->xoct[3]);*/
/*printf("zunten %f zoben %f\n",octreetemp->xoct[4], octreetemp->xoct[5]);*/
/*for (i=0;i<octreetemp->numnpoct;i++)*/
/*{*/
/*  printf("Knoten %d hat Koordinaten %f %f %f\n",nodelist0->node->Id,nodelist0->node->x[0],nodelist0->node->x[1],nodelist0->node->x[2]);*/
/*  nodelist0=nodelist0->next;*/
/*}*/
/*printf("fluid node %d has coordinates %f %f %f\n", actfnode->Id,actfnode->x[0],actfnode->x[1],actfnode->x[2]);*/
/*}*/

  /********************************************/
  /*find structure host element for fluid node*/
  /********************************************/

  nodelist=octreetemp->nodelist;
  actsnode=nodelist->node;

  /*loop all structure nodes in octree partition*/
  for (i=0;i<octreetemp->numnpoct;i++)
  {
    /*plausibility check: the struct node has to lie within the boundaries of octree (FSI interface)*/
    if (!((octreetemp->xoct[0]-EPS<=actsnode->x[0])&&(actsnode->x[0]<=octreetemp->xoct[1]+EPS)&&
          (octreetemp->xoct[2]-EPS<=actsnode->x[1])&&(actsnode->x[1]<=octreetemp->xoct[3]+EPS)&&
          (octreetemp->xoct[4]-EPS<=actsnode->x[2])&&(actsnode->x[2]<=octreetemp->xoct[5]+EPS)))
    {
      dserror("structure node %d does not lie within the boundaries of the octree (FSI interface)\n",actsnode->Id);
    }

    /*Note: It is necessary that all FSI elements can access their *coupleptr later on.
     * Therefore all elements have to be looped prior to the host element search.*/

    /*loop all elements of this node*/
    for (j=0;j<actsnode->numele;j++)
    {
      hostele=actsnode->element[j];

      /*if space for the coupling information was not allocated yet (elements belong to several nodes!)*/
      if (hostele->coupleptr==NULL)
      {
        /*allocate space for the coupling information*/
        hostele->coupleptr=(FSI_COUPLE_NODE*)CCACALLOC(1,sizeof(FSI_COUPLE_NODE));
        hostele->coupleptr->numnp=0;
        hostele->coupleptr->couplenode=(NODE**)CCACALLOC(max,sizeof(NODE*));
        for (k=0;k<max;k++) hostele->coupleptr->couplenode[k]=NULL;
      }
    }

    /*loop all elements of this node*/
    for (j=0;j<actsnode->numele;j++)
    {
      hostele=actsnode->element[j];
      found=0;

      /*checking the struct element type*/
      if (!(hostele->distyp==hex8 ||
            hostele->distyp==hex27 ||
            hostele->distyp==tet4 ||
            hostele->distyp==tet10 ||
            hostele->distyp==tri3))
         dserror("Invalid element type: only hex8,hex27,tet4,tet10 allowed\n");

      /*if there is no more space for fluid nodes in the *couplenode array of the struct element -> reallocate*/
      if (hostele->coupleptr->numnp==max)
      {
        max+=max;
        hostele->coupleptr->couplenode=(NODE**)CCAREALLOC(hostele->coupleptr->couplenode,max*sizeof(NODE*));
      }

      /*allocate space for the coordinates of the nodes belonging to the element*/
      x=(ARRAY*)CCACALLOC(1,sizeof(ARRAY));
      amdef("x",x,3,hostele->numnp,"DA");
      amzero(x);

      /*write coordinates of nodes belonging to element into array*/
      for (k=0;k<hostele->numnp;k++)
      {
        for (m=0;m<3;m++)   /*loop all dimensions*/
        {
          x->a.da[m][k]=hostele->node[k]->x[m];
        }
      }

      /*allocate space for the coefficients*/
      coeff=(ARRAY*)CCACALLOC(1,sizeof(ARRAY));
      amdef("coeff",coeff,3,hostele->numnp,"DA");
      amzero(coeff);

      /*allocate space for the Jacobian matrix*/
      jacobi=(ARRAY*)CCACALLOC(1,sizeof(ARRAY));
      amdef("jacobi",jacobi,3,3,"DA");
      amzero(jacobi);

      /*set initial values for the local coordinates and the error*/
      /*Arbitrarily the center of the element is chosen.*/
      r=0.0;
      s=0.0;
      t=0.0;
      error=1;

      switch (hostele->distyp)
      {
        /*HEX8 element*/
      case hex8:
      {
        for (m=0;m<3;m++)     /*loop all dimensions*/
        {
          coeff->a.da[m][0]= .125*( x->a.da[m][0]+x->a.da[m][1]+x->a.da[m][2]+x->a.da[m][3]+x->a.da[m][4]+x->a.da[m][5]+x->a.da[m][6]+x->a.da[m][7]);
	  coeff->a.da[m][1]= .125*( x->a.da[m][0]+x->a.da[m][1]-x->a.da[m][2]-x->a.da[m][3]+x->a.da[m][4]+x->a.da[m][5]-x->a.da[m][6]-x->a.da[m][7]);
	  coeff->a.da[m][2]= .125*(-x->a.da[m][0]+x->a.da[m][1]+x->a.da[m][2]-x->a.da[m][3]-x->a.da[m][4]+x->a.da[m][5]+x->a.da[m][6]-x->a.da[m][7]);
	  coeff->a.da[m][3]= .125*(-x->a.da[m][0]-x->a.da[m][1]-x->a.da[m][2]-x->a.da[m][3]+x->a.da[m][4]+x->a.da[m][5]+x->a.da[m][6]+x->a.da[m][7]);
	  coeff->a.da[m][4]= .125*(-x->a.da[m][0]+x->a.da[m][1]-x->a.da[m][2]+x->a.da[m][3]-x->a.da[m][4]+x->a.da[m][5]-x->a.da[m][6]+x->a.da[m][7]);
	  coeff->a.da[m][5]= .125*(-x->a.da[m][0]-x->a.da[m][1]+x->a.da[m][2]+x->a.da[m][3]+x->a.da[m][4]+x->a.da[m][5]-x->a.da[m][6]-x->a.da[m][7]);
	  coeff->a.da[m][6]= .125*( x->a.da[m][0]-x->a.da[m][1]-x->a.da[m][2]+x->a.da[m][3]-x->a.da[m][4]+x->a.da[m][5]+x->a.da[m][6]-x->a.da[m][7]);
	  coeff->a.da[m][7]= .125*( x->a.da[m][0]-x->a.da[m][1]+x->a.da[m][2]-x->a.da[m][3]-x->a.da[m][4]+x->a.da[m][5]-x->a.da[m][6]+x->a.da[m][7]);
        }

        /*solve the nonlinear equation system for the natural coordinats of the fluid node*/
        /*loop until the error is smaller than the specified value*/
        while (error>ERROR)
        {
          /*initialize the Jacobian matrix (linearization)*/
          jacobi->a.da[0][0]=coeff->a.da[0][1]+coeff->a.da[0][4]*s+coeff->a.da[0][5]*t+coeff->a.da[0][7]*s*t;
          jacobi->a.da[0][1]=coeff->a.da[0][2]+coeff->a.da[0][4]*r+coeff->a.da[0][6]*t+coeff->a.da[0][7]*r*t;
	  jacobi->a.da[0][2]=coeff->a.da[0][3]+coeff->a.da[0][5]*r+coeff->a.da[0][6]*s+coeff->a.da[0][7]*r*s;
          jacobi->a.da[1][0]=coeff->a.da[1][1]+coeff->a.da[1][4]*s+coeff->a.da[1][5]*t+coeff->a.da[1][7]*s*t;
          jacobi->a.da[1][1]=coeff->a.da[1][2]+coeff->a.da[1][4]*r+coeff->a.da[1][6]*t+coeff->a.da[1][7]*r*t;
          jacobi->a.da[1][2]=coeff->a.da[1][3]+coeff->a.da[1][5]*r+coeff->a.da[1][6]*s+coeff->a.da[1][7]*r*s;
          jacobi->a.da[2][0]=coeff->a.da[2][1]+coeff->a.da[2][4]*s+coeff->a.da[2][5]*t+coeff->a.da[2][7]*s*t;
          jacobi->a.da[2][1]=coeff->a.da[2][2]+coeff->a.da[2][4]*r+coeff->a.da[2][6]*t+coeff->a.da[2][7]*r*t;
          jacobi->a.da[2][2]=coeff->a.da[2][3]+coeff->a.da[2][5]*r+coeff->a.da[2][6]*s+coeff->a.da[2][7]*r*s;

          /*calculate the residual*/
          for (k=0;k<3;k++)   /*loop all dimensions*/
          {
            res[k]=actfnode->x[k]-(coeff->a.da[k][0]+coeff->a.da[k][1]*r+coeff->a.da[k][2]*s+coeff->a.da[k][3]*t+coeff->a.da[k][4]*r*s+coeff->a.da[k][5]*r*t+coeff->a.da[k][6]*s*t+coeff->a.da[k][7]*r*s*t);
          }

          /*solving the linear equation system explicitly for the delta vector*/
          deltar=((jacobi->a.da[0][1]*jacobi->a.da[1][2]*res[2]-jacobi->a.da[0][1]*res[1]*jacobi->a.da[2][2]+jacobi->a.da[0][2]*jacobi->a.da[2][1]*res[1]-jacobi->a.da[0][2]*jacobi->a.da[1][1]*res[2]+res[0]*jacobi->a.da[1][1]*jacobi->a.da[2][2]-res[0]*jacobi->a.da[2][1]*jacobi->a.da[1][2])/
                  (jacobi->a.da[0][0]*jacobi->a.da[1][1]*jacobi->a.da[2][2]-jacobi->a.da[0][0]*jacobi->a.da[2][1]*jacobi->a.da[1][2]-jacobi->a.da[1][0]*jacobi->a.da[0][1]*jacobi->a.da[2][2]+jacobi->a.da[2][1]*jacobi->a.da[1][0]*jacobi->a.da[0][2]-jacobi->a.da[1][1]*jacobi->a.da[2][0]*jacobi->a.da[0][2]+jacobi->a.da[2][0]*jacobi->a.da[0][1]*jacobi->a.da[1][2]));
	  deltas=-((jacobi->a.da[0][0]*jacobi->a.da[1][2]*res[2]-jacobi->a.da[0][0]*res[1]*jacobi->a.da[2][2]-jacobi->a.da[1][0]*jacobi->a.da[0][2]*res[2]-jacobi->a.da[1][2]*jacobi->a.da[2][0]*res[0]+res[1]*jacobi->a.da[2][0]*jacobi->a.da[0][2]+jacobi->a.da[1][0]*res[0]*jacobi->a.da[2][2])/
                   (jacobi->a.da[0][0]*jacobi->a.da[1][1]*jacobi->a.da[2][2]-jacobi->a.da[0][0]*jacobi->a.da[2][1]*jacobi->a.da[1][2]-jacobi->a.da[1][0]*jacobi->a.da[0][1]*jacobi->a.da[2][2]+jacobi->a.da[2][1]*jacobi->a.da[1][0]*jacobi->a.da[0][2]-jacobi->a.da[1][1]*jacobi->a.da[2][0]*jacobi->a.da[0][2]+jacobi->a.da[2][0]*jacobi->a.da[0][1]*jacobi->a.da[1][2]));
	  deltat=((jacobi->a.da[2][1]*jacobi->a.da[1][0]*res[0]-jacobi->a.da[0][0]*jacobi->a.da[2][1]*res[1]+jacobi->a.da[0][0]*jacobi->a.da[1][1]*res[2]-jacobi->a.da[1][1]*jacobi->a.da[2][0]*res[0]-jacobi->a.da[1][0]*jacobi->a.da[0][1]*res[2]+jacobi->a.da[2][0]*jacobi->a.da[0][1]*res[1])/
                  (jacobi->a.da[0][0]*jacobi->a.da[1][1]*jacobi->a.da[2][2]-jacobi->a.da[0][0]*jacobi->a.da[2][1]*jacobi->a.da[1][2]-jacobi->a.da[1][0]*jacobi->a.da[0][1]*jacobi->a.da[2][2]+jacobi->a.da[2][1]*jacobi->a.da[1][0]*jacobi->a.da[0][2]-jacobi->a.da[1][1]*jacobi->a.da[2][0]*jacobi->a.da[0][2]+jacobi->a.da[2][0]*jacobi->a.da[0][1]*jacobi->a.da[1][2]));

          /*advance the values for r,s,t approaching the solution*/
          r=r+deltar;
	  s=s+deltas;
	  t=t+deltat;
          /*calculate the error*/
          error=sqrt(res[0]*res[0]+res[1]*res[1]+res[2]*res[2]);

/*printf("für das element %d\n", hostele->Id);*/
/*for (k=0;k<hostele->numnp;k++)*/
/*{*/
/*  printf("node %d x[0,%d] = %f x[1,%d] = %f x[2,%d] = %f\n",hostele->node[k]->Id,k,x->a.da[0][k],k,x->a.da[1][k],k,x->a.da[2][k]);*/
/*}*/
/*printf("r %f,s %f,t %f\n", r, s, t);*/
/*printf("deltar %f,deltas %f,deltat %f\n", deltar, deltas, deltat);*/
/*printf("|res|=%f\n", error);*/
        }

        /*The local coordinates r,s,t of the fluid node related to the struct element are now known.
         *Use them to find the position in relation to the struct element*/

        /*criterion for hex8: -1< r,s,t <1 -> fluid node inside struct element*/
        if(r>=-1-EPS && s>=-1-EPS && t>=-1-EPS &&r<=1+EPS && s<=1+EPS &&t<=1+EPS)
	{
          /*host element found!*/
          found=1;
/*printf("fluid node %d coupled with struct element %d\n", actfnode->Id, hostele->Id);*/
/*printf("r %f,s %f,t %f\n\n", r, s, t);*/
/*printf("coordinates of fluid node %d: %f %f %f\n",actfnode->Id,actfnode->x[0],actfnode->x[1],actfnode->x[2]);*/
        }
        else
	{
          /*struct element is not the host element*/
          found=0;
/*printf("no host element found\n");*/
        }
        break;
      }/*end of if (HEX8)*/

      /*HEX27 element*/
      case hex27:
      {
        for (m=0;m<3;m++)     /*loop all dimensions*/
        {
          coeff->a.da[m][0]=x->a.da[m][26];
          coeff->a.da[m][1]=.125*(4*x->a.da[m][21]-4*x->a.da[m][23]);
          coeff->a.da[m][2]=.125*(-4*x->a.da[m][22]+4*x->a.da[m][24]);
          coeff->a.da[m][3]=.125*(4*x->a.da[m][20]-4*x->a.da[m][25]);
          coeff->a.da[m][4]=.125*(-2*x->a.da[m][15]+2*x->a.da[m][14]+2*x->a.da[m][12]-2*x->a.da[m][13]);
          coeff->a.da[m][5]=.125*(2*x->a.da[m][8]-2*x->a.da[m][16]-2*x->a.da[m][10]+2*x->a.da[m][18]);
          coeff->a.da[m][6]=.125*(2*x->a.da[m][17]-2*x->a.da[m][19]+2*x->a.da[m][11]-2*x->a.da[m][9]);
          coeff->a.da[m][7]=.125*(-x->a.da[m][3]+x->a.da[m][7]+x->a.da[m][2]-x->a.da[m][1]-x->a.da[m][4]+x->a.da[m][0]-x->a.da[m][6]+x->a.da[m][5]);
          coeff->a.da[m][8]=.125*(-8*x->a.da[m][26]+4*x->a.da[m][21]+4*x->a.da[m][23]);
          coeff->a.da[m][9]=.125*(4*x->a.da[m][24]+4*x->a.da[m][22]-8*x->a.da[m][26]);
          coeff->a.da[m][10]=.125*(4*x->a.da[m][25]+4*x->a.da[m][20]-8*x->a.da[m][26]);
          coeff->a.da[m][11]=.125*(2*x->a.da[m][15]+4*x->a.da[m][22]-2*x->a.da[m][13]-2*x->a.da[m][14]-4*x->a.da[m][24]+2*x->a.da[m][12]);
          coeff->a.da[m][12]=.125*(-2*x->a.da[m][18]+2*x->a.da[m][8]+2*x->a.da[m][10]-4*x->a.da[m][20]-2*x->a.da[m][16]+4*x->a.da[m][25]);
          coeff->a.da[m][13]=.125*(2*x->a.da[m][13]-2*x->a.da[m][14]+4*x->a.da[m][23]-4*x->a.da[m][21]+2*x->a.da[m][12]-2*x->a.da[m][15]);
          coeff->a.da[m][14]=.125*(-4*x->a.da[m][20]-2*x->a.da[m][19]+2*x->a.da[m][11]+2*x->a.da[m][9]+4*x->a.da[m][25]-2*x->a.da[m][17]);
          coeff->a.da[m][15]=.125*(2*x->a.da[m][8]+4*x->a.da[m][23]-2*x->a.da[m][10]+2*x->a.da[m][16]-4*x->a.da[m][21]-2*x->a.da[m][18]);
          coeff->a.da[m][16]=.125*(-2*x->a.da[m][17]+2*x->a.da[m][19]+4*x->a.da[m][22]-2*x->a.da[m][9]+2*x->a.da[m][11]-4*x->a.da[m][24]);
          coeff->a.da[m][17]=.125*(2*x->a.da[m][13]+2*x->a.da[m][12]+2*x->a.da[m][15]-4*x->a.da[m][22]-4*x->a.da[m][23]-4*x->a.da[m][24]+2*x->a.da[m][14]+8*x->a.da[m][26]-4*x->a.da[m][21]);
          coeff->a.da[m][18]=.125*(8*x->a.da[m][26]-4*x->a.da[m][21]-4*x->a.da[m][23]+2*x->a.da[m][16]-4*x->a.da[m][25]-4*x->a.da[m][20]+2*x->a.da[m][18]+2*x->a.da[m][10]+2*x->a.da[m][8]);
          coeff->a.da[m][19]=.125*(-4*x->a.da[m][24]+2*x->a.da[m][11]+8*x->a.da[m][26]-4*x->a.da[m][25]+2*x->a.da[m][17]+2*x->a.da[m][9]-4*x->a.da[m][22]-4*x->a.da[m][20]+2*x->a.da[m][19]);
          coeff->a.da[m][20]=.125*(-x->a.da[m][1]+x->a.da[m][6]+x->a.da[m][5]+x->a.da[m][0]+2*x->a.da[m][9]-x->a.da[m][7]-2*x->a.da[m][11]+x->a.da[m][3]-2*x->a.da[m][17]-x->a.da[m][2]+2*x->a.da[m][19]-x->a.da[m][4]);
          coeff->a.da[m][21]=.125*(-2*x->a.da[m][18]+2*x->a.da[m][10]-x->a.da[m][5]-x->a.da[m][3]+x->a.da[m][1]+2*x->a.da[m][16]+x->a.da[m][0]-x->a.da[m][2]+x->a.da[m][7]-x->a.da[m][4]+x->a.da[m][6]-2*x->a.da[m][8]);
          coeff->a.da[m][22]=.125*(x->a.da[m][6]+2*x->a.da[m][13]+2*x->a.da[m][15]-x->a.da[m][3]-x->a.da[m][7]-2*x->a.da[m][14]-x->a.da[m][5]+x->a.da[m][2]-x->a.da[m][1]+x->a.da[m][4]-2*x->a.da[m][12]+x->a.da[m][0]);
          coeff->a.da[m][23]=.125*(x->a.da[m][5]-x->a.da[m][7]-2*x->a.da[m][13]+2*x->a.da[m][10]+x->a.da[m][0]+2*x->a.da[m][14]+x->a.da[m][1]-2*x->a.da[m][16]-2*x->a.da[m][12]-x->a.da[m][2]+2*x->a.da[m][18]-4*x->a.da[m][23]+4*x->a.da[m][21]+x->a.da[m][4]-x->a.da[m][3]-x->a.da[m][6]+2*x->a.da[m][15]-2*x->a.da[m][8]);
          coeff->a.da[m][24]=.125*(-2*x->a.da[m][15]-x->a.da[m][6]-x->a.da[m][1]+2*x->a.da[m][9]-x->a.da[m][5]+4*x->a.da[m][24]+2*x->a.da[m][13]-2*x->a.da[m][12]-2*x->a.da[m][11]+x->a.da[m][3]-4*x->a.da[m][22]+2*x->a.da[m][17]+x->a.da[m][7]-x->a.da[m][2]+x->a.da[m][0]-2*x->a.da[m][19]+2*x->a.da[m][14]+x->a.da[m][4]);
          coeff->a.da[m][25]=.125*(-x->a.da[m][5]+x->a.da[m][3]-x->a.da[m][7]-x->a.da[m][4]-2*x->a.da[m][10]+4*x->a.da[m][20]-2*x->a.da[m][9]-x->a.da[m][6]+2*x->a.da[m][16]+x->a.da[m][0]-2*x->a.da[m][11]+x->a.da[m][1]+2*x->a.da[m][18]-4*x->a.da[m][25]+2*x->a.da[m][17]+x->a.da[m][2]+2*x->a.da[m][19]-2*x->a.da[m][8]);
          coeff->a.da[m][26]=.125*(x->a.da[m][0]+x->a.da[m][1]+x->a.da[m][2]+x->a.da[m][3]+x->a.da[m][4]+x->a.da[m][5]+x->a.da[m][6]+x->a.da[m][7]-2*x->a.da[m][8]-2*x->a.da[m][9]-2*x->a.da[m][10]-2*x->a.da[m][11]-2*x->a.da[m][12]-2*x->a.da[m][13]-2*x->a.da[m][14]-2*x->a.da[m][15]-2*x->a.da[m][16]-2*x->a.da[m][17]-2*x->a.da[m][18]-2*x->a.da[m][19]+4*x->a.da[m][20]+4*x->a.da[m][21]+4*x->a.da[m][22]+4*x->a.da[m][23]+4*x->a.da[m][24]+4*x->a.da[m][25]-8*x->a.da[m][26]);
        }

        /*solve the nonlinear equation system for the natural coordinats of the fluid node*/
        /*loop until the error is smaller than the specified value*/
        while (error>ERROR)   /*How large is the error?*/
        {
          /*initialize the Jacobian matrix*/
          jacobi->a.da[0][0]=coeff->a.da[0][1]+coeff->a.da[0][4]*s+coeff->a.da[0][5]*t+coeff->a.da[0][7]*s*t+coeff->a.da[0][8]*2*r+coeff->a.da[0][11]*2*r*s+coeff->a.da[0][12]*2*r*t+coeff->a.da[0][13]*s*s+coeff->a.da[0][15]*t*t+coeff->a.da[0][17]*2*r*s*s+coeff->a.da[0][18]*2*r*t*t+coeff->a.da[0][20]*2*r*s*t+coeff->a.da[0][21]*s*s*t+coeff->a.da[0][22]*s*t*t+coeff->a.da[0][23]*s*s*t*t+coeff->a.da[0][24]*2*r*s*t*t+coeff->a.da[0][25]*2*r*s*s*t+coeff->a.da[0][26]*2*r*s*s*t*t;
          jacobi->a.da[0][1]=coeff->a.da[0][2]+coeff->a.da[0][4]*r+coeff->a.da[0][6]*t+coeff->a.da[0][7]*r*t+coeff->a.da[0][9]*2*s+coeff->a.da[0][11]*r*r+coeff->a.da[0][13]*r*2*s+coeff->a.da[0][14]*2*s*t+coeff->a.da[0][16]*2*t+coeff->a.da[0][17]*r*r*2*s+coeff->a.da[0][19]*2*s+t*t+coeff->a.da[0][20]*r*r*t+coeff->a.da[0][21]*r*2*s*t+coeff->a.da[0][22]*r*t*t+coeff->a.da[0][23]*r*2*s*t*t+coeff->a.da[0][24]*r*r*t*t+coeff->a.da[0][25]*r*r*2*s*t+coeff->a.da[0][26]*r*r*2*s*t*t;
          jacobi->a.da[0][2]=coeff->a.da[0][3]+coeff->a.da[0][5]*r+coeff->a.da[0][6]*s+coeff->a.da[0][7]*r*t+coeff->a.da[0][10]*2*t+coeff->a.da[0][12]*r*r+coeff->a.da[0][14]*s*s+coeff->a.da[0][15]*r*2*t+coeff->a.da[0][16]*s*2*t+coeff->a.da[0][18]*r*r*2*t+coeff->a.da[0][19]*s*s*2*t+coeff->a.da[0][20]*r*r*s+coeff->a.da[0][21]*r*s*s+coeff->a.da[0][22]*r*s*2*t+coeff->a.da[0][23]*r*s*s*2*t+coeff->a.da[0][24]*r*r*s*2*t+coeff->a.da[0][25]*r*r*s*s+coeff->a.da[0][26]*r*r*s*s*2*t;
          jacobi->a.da[1][0]=coeff->a.da[1][1]+coeff->a.da[1][4]*s+coeff->a.da[1][5]*t+coeff->a.da[1][7]*s*t+coeff->a.da[1][8]*2*r+coeff->a.da[1][11]*2*r*s+coeff->a.da[1][12]*2*r*t+coeff->a.da[1][13]*s*s+coeff->a.da[1][15]*t*t+coeff->a.da[1][17]*2*r*s*s+coeff->a.da[1][18]*2*r*t*t+coeff->a.da[1][20]*2*r*s*t+coeff->a.da[1][21]*s*s*t+coeff->a.da[1][22]*s*t*t+coeff->a.da[1][23]*s*s*t*t+coeff->a.da[1][24]*2*r*s*t*t+coeff->a.da[1][25]*2*r*s*s*t+coeff->a.da[1][26]*2*r*s*s*t*t;
          jacobi->a.da[1][1]=coeff->a.da[1][2]+coeff->a.da[1][4]*r+coeff->a.da[1][6]*t+coeff->a.da[1][7]*r*t+coeff->a.da[1][9]*2*s+coeff->a.da[1][11]*r*r+coeff->a.da[1][13]*r*2*s+coeff->a.da[1][14]*2*s*t+coeff->a.da[1][16]*2*t+coeff->a.da[1][17]*r*r*2*s+coeff->a.da[1][19]*2*s+t*t+coeff->a.da[1][20]*r*r*t+coeff->a.da[1][21]*r*2*s*t+coeff->a.da[1][22]*r*t*t+coeff->a.da[1][23]*r*2*s*t*t+coeff->a.da[1][24]*r*r*t*t+coeff->a.da[1][25]*r*r*2*s*t+coeff->a.da[1][26]*r*r*2*s*t*t;
          jacobi->a.da[1][2]=coeff->a.da[1][3]+coeff->a.da[1][5]*r+coeff->a.da[1][6]*s+coeff->a.da[1][7]*r*t+coeff->a.da[1][10]*2*t+coeff->a.da[1][12]*r*r+coeff->a.da[1][14]*s*s+coeff->a.da[1][15]*r*2*t+coeff->a.da[1][16]*s*2*t+coeff->a.da[1][18]*r*r*2*t+coeff->a.da[1][19]*s*s*2*t+coeff->a.da[1][20]*r*r*s+coeff->a.da[1][21]*r*s*s+coeff->a.da[1][22]*r*s*2*t+coeff->a.da[1][23]*r*s*s*2*t+coeff->a.da[1][24]*r*r*s*2*t+coeff->a.da[1][25]*r*r*s*s+coeff->a.da[1][26]*r*r*s*s*2*t;
          jacobi->a.da[2][0]=coeff->a.da[2][1]+coeff->a.da[2][4]*s+coeff->a.da[2][5]*t+coeff->a.da[2][7]*s*t+coeff->a.da[2][8]*2*r+coeff->a.da[2][11]*2*r*s+coeff->a.da[2][12]*2*r*t+coeff->a.da[2][13]*s*s+coeff->a.da[2][15]*t*t+coeff->a.da[2][17]*2*r*s*s+coeff->a.da[2][18]*2*r*t*t+coeff->a.da[2][20]*2*r*s*t+coeff->a.da[2][21]*s*s*t+coeff->a.da[2][22]*s*t*t+coeff->a.da[2][23]*s*s*t*t+coeff->a.da[2][24]*2*r*s*t*t+coeff->a.da[2][25]*2*r*s*s*t+coeff->a.da[2][26]*2*r*s*s*t*t;
          jacobi->a.da[2][1]=coeff->a.da[2][2]+coeff->a.da[2][4]*r+coeff->a.da[2][6]*t+coeff->a.da[2][7]*r*t+coeff->a.da[2][9]*2*s+coeff->a.da[2][11]*r*r+coeff->a.da[2][13]*r*2*s+coeff->a.da[2][14]*2*s*t+coeff->a.da[2][16]*2*t+coeff->a.da[2][17]*r*r*2*s+coeff->a.da[2][19]*2*s+t*t+coeff->a.da[2][20]*r*r*t+coeff->a.da[2][21]*r*2*s*t+coeff->a.da[2][22]*r*t*t+coeff->a.da[2][23]*r*2*s*t*t+coeff->a.da[2][24]*r*r*t*t+coeff->a.da[2][25]*r*r*2*s*t+coeff->a.da[2][26]*r*r*2*s*t*t;
          jacobi->a.da[2][2]=coeff->a.da[2][3]+coeff->a.da[2][5]*r+coeff->a.da[2][6]*s+coeff->a.da[2][7]*r*t+coeff->a.da[2][10]*2*t+coeff->a.da[2][12]*r*r+coeff->a.da[2][14]*s*s+coeff->a.da[2][15]*r*2*t+coeff->a.da[2][16]*s*2*t+coeff->a.da[2][18]*r*r*2*t+coeff->a.da[2][19]*s*s*2*t+coeff->a.da[2][20]*r*r*s+coeff->a.da[2][21]*r*s*s+coeff->a.da[2][22]*r*s*2*t+coeff->a.da[2][23]*r*s*s*2*t+coeff->a.da[2][24]*r*r*s*2*t+coeff->a.da[2][25]*r*r*s*s+coeff->a.da[2][26]*r*r*s*s*2*t;

          /*calculate the residual*/
          for (k=0;k<3;k++)    /*loop all dimensions*/
          {
            res[k]=actfnode->x[k]-(coeff->a.da[k][0] + coeff->a.da[k][1]*r + coeff->a.da[k][2]*s + coeff->a.da[k][3]*t + coeff->a.da[k][4]*r*s + coeff->a.da[k][5]*r*t + coeff->a.da[k][6]*s*t + coeff->a.da[k][7]*r*s*t + coeff->a.da[k][8]*r*r + coeff->a.da[k][9]*s*s + coeff->a.da[k][10]*t*t + coeff->a.da[k][11]*r*r*s + coeff->a.da[k][12]*r*r*t + coeff->a.da[k][13]*r*s*s + coeff->a.da[k][14]*s*s*t + coeff->a.da[k][15]*r*t*t + coeff->a.da[k][16]*s*t*t + coeff->a.da[k][17]*r*r*s*s + coeff->a.da[k][18]*r*r*t*t + coeff->a.da[k][19]*s*s*t*t + coeff->a.da[k][20]*r*r*s*t + coeff->a.da[k][21]*r*s*s*t + coeff->a.da[k][22]*r*s*t*t + coeff->a.da[k][23]*r*s*s*t*t + coeff->a.da[k][24]*r*r*s*t*t + coeff->a.da[k][25]*r*r*s*s*t + coeff->a.da[k][26]*r*r*s*s*t*t);
          }

          /*solving the linear equation system explicitly for the delta vector*/
          deltar=(jacobi->a.da[0][1]*jacobi->a.da[1][2]*res[2]-jacobi->a.da[0][1]*res[1]*jacobi->a.da[2][2]+jacobi->a.da[0][2]*jacobi->a.da[2][1]*res[1]-jacobi->a.da[0][2]*jacobi->a.da[1][1]*res[2]+res[0]*jacobi->a.da[1][1]*jacobi->a.da[2][2]-res[0]*jacobi->a.da[2][1]*jacobi->a.da[1][2])/
            (jacobi->a.da[0][0]*jacobi->a.da[1][1]*jacobi->a.da[2][2]-jacobi->a.da[0][0]*jacobi->a.da[2][1]*jacobi->a.da[1][2]-jacobi->a.da[1][0]*jacobi->a.da[0][1]*jacobi->a.da[2][2]+jacobi->a.da[2][1]*jacobi->a.da[1][0]*jacobi->a.da[0][2]-jacobi->a.da[1][1]*jacobi->a.da[2][0]*jacobi->a.da[0][2]+jacobi->a.da[2][0]*jacobi->a.da[0][1]*jacobi->a.da[1][2]);
	  deltas=-(jacobi->a.da[0][0]*jacobi->a.da[1][2]*res[2]-jacobi->a.da[0][0]*res[1]*jacobi->a.da[2][2]-jacobi->a.da[1][0]*jacobi->a.da[0][2]*res[2]-jacobi->a.da[1][2]*jacobi->a.da[2][0]*res[0]+res[1]*jacobi->a.da[2][0]*jacobi->a.da[0][2]+jacobi->a.da[1][0]*res[0]*jacobi->a.da[2][2])/
            (jacobi->a.da[0][0]*jacobi->a.da[1][1]*jacobi->a.da[2][2]-jacobi->a.da[0][0]*jacobi->a.da[2][1]*jacobi->a.da[1][2]-jacobi->a.da[1][0]*jacobi->a.da[0][1]*jacobi->a.da[2][2]+jacobi->a.da[2][1]*jacobi->a.da[1][0]*jacobi->a.da[0][2]-jacobi->a.da[1][1]*jacobi->a.da[2][0]*jacobi->a.da[0][2]+jacobi->a.da[2][0]*jacobi->a.da[0][1]*jacobi->a.da[1][2]);
	  deltat=(jacobi->a.da[2][1]*jacobi->a.da[1][0]*res[0]-jacobi->a.da[0][0]*jacobi->a.da[2][1]*res[1]+jacobi->a.da[0][0]*jacobi->a.da[1][1]*res[2]-jacobi->a.da[1][1]*jacobi->a.da[2][0]*res[0]-jacobi->a.da[1][0]*jacobi->a.da[0][1]*res[2]+jacobi->a.da[2][0]*jacobi->a.da[0][1]*res[1])/
            (jacobi->a.da[0][0]*jacobi->a.da[1][1]*jacobi->a.da[2][2]-jacobi->a.da[0][0]*jacobi->a.da[2][1]*jacobi->a.da[1][2]-jacobi->a.da[1][0]*jacobi->a.da[0][1]*jacobi->a.da[2][2]+jacobi->a.da[2][1]*jacobi->a.da[1][0]*jacobi->a.da[0][2]-jacobi->a.da[1][1]*jacobi->a.da[2][0]*jacobi->a.da[0][2]+jacobi->a.da[2][0]*jacobi->a.da[0][1]*jacobi->a.da[1][2]);

          /*advance the values for r,s,t approaching the solution*/
          r=r+deltar;
	  s=s+deltas;
	  t=t+deltat;
          /*calculate the error*/
	  error=sqrt(res[0]*res[0]+res[1]*res[1]+res[2]*res[2]);
        }

        /*The local coordinates r,s,t of the fluid node related to the struct element are now known.
         *Use them to find the position in relation to the struct element*/

        /*criterion for hex27: -1< r,s,t <1 -> fluid node inside struct element*/
        if(r>=-1 && s>=-1 && t>=-1 &&r<=1 && s<=1 &&t<=1)
	{
          /*host element found!*/
	  found=1;
       	}
	else
	{
          /*struct element is not the host element*/
          found=0;
	}
        break;
      }/*end of if(HEX27)*/

      /*TET4 element*/
      case tet4:
      {
        for (m=0;m<3;m++)     /*loop all dimensions*/
        {
          coeff->a.da[m][0]= x->a.da[m][0];
          coeff->a.da[m][0]=-x->a.da[m][0]+x->a.da[m][1];
          coeff->a.da[m][0]=-x->a.da[m][0]+x->a.da[m][2];
          coeff->a.da[m][0]=-x->a.da[m][0]+x->a.da[m][3];
        }

        /*solve the nonlinear equation system for the natural coordinats of the fluid node*/
        /*loop until the error is smaller than the specified value*/
        while (error>ERROR)
        {
          /*initialize the Jacobian matrix*/
          jacobi->a.da[0][0]=coeff->a.da[0][1]+coeff->a.da[0][4]*s+coeff->a.da[0][5]*t+coeff->a.da[0][7]*s*t;
          jacobi->a.da[0][1]=coeff->a.da[0][2]+coeff->a.da[0][4]*r+coeff->a.da[0][6]*t+coeff->a.da[0][7]*r*t;
          jacobi->a.da[0][2]=coeff->a.da[0][3]+coeff->a.da[0][5]*r+coeff->a.da[0][6]*s+coeff->a.da[0][7]*r*s;
          jacobi->a.da[1][0]=coeff->a.da[1][1]+coeff->a.da[1][4]*s+coeff->a.da[1][5]*t+coeff->a.da[1][7]*s*t;
          jacobi->a.da[1][1]=coeff->a.da[1][2]+coeff->a.da[1][4]*r+coeff->a.da[1][6]*t+coeff->a.da[1][7]*r*t;
          jacobi->a.da[1][2]=coeff->a.da[1][3]+coeff->a.da[1][5]*r+coeff->a.da[1][6]*s+coeff->a.da[1][7]*r*s;
          jacobi->a.da[2][0]=coeff->a.da[2][1]+coeff->a.da[2][4]*s+coeff->a.da[2][5]*t+coeff->a.da[2][7]*s*t;
          jacobi->a.da[2][1]=coeff->a.da[2][2]+coeff->a.da[2][4]*r+coeff->a.da[2][6]*t+coeff->a.da[2][7]*r*t;
          jacobi->a.da[2][2]=coeff->a.da[2][3]+coeff->a.da[2][5]*r+coeff->a.da[2][6]*s+coeff->a.da[2][7]*r*s;

          /*calculate the residual*/
          for (k=0;k<3;k++)    /*loop all dimensions*/
          {
            res[k]=actfnode->x[k]-(coeff->a.da[m][0]+coeff->a.da[m][1]*r+coeff->a.da[m][2]*s+coeff->a.da[m][3]*t+coeff->a.da[m][4]*r*s+coeff->a.da[m][5]*r*t+coeff->a.da[m][6]*s*t+coeff->a.da[m][7]*r*s*t);
          }

          /*solving the linear equation system explicitly for the delta vector*/
          deltar=(jacobi->a.da[0][1]*jacobi->a.da[1][2]*res[2]-jacobi->a.da[0][1]*res[1]*jacobi->a.da[2][2]+jacobi->a.da[0][2]*jacobi->a.da[2][1]*res[1]-jacobi->a.da[0][2]*jacobi->a.da[1][1]*res[2]+res[0]*jacobi->a.da[1][1]*jacobi->a.da[2][2]-res[0]*jacobi->a.da[2][1]*jacobi->a.da[1][2])/
            (jacobi->a.da[0][0]*jacobi->a.da[1][1]*jacobi->a.da[2][2]-jacobi->a.da[0][0]*jacobi->a.da[2][1]*jacobi->a.da[1][2]-jacobi->a.da[1][0]*jacobi->a.da[0][1]*jacobi->a.da[2][2]+jacobi->a.da[2][1]*jacobi->a.da[1][0]*jacobi->a.da[0][2]-jacobi->a.da[1][1]*jacobi->a.da[2][0]*jacobi->a.da[0][2]+jacobi->a.da[2][0]*jacobi->a.da[0][1]*jacobi->a.da[1][2]);
	  deltas=-(jacobi->a.da[0][0]*jacobi->a.da[1][2]*res[2]-jacobi->a.da[0][0]*res[1]*jacobi->a.da[2][2]-jacobi->a.da[1][0]*jacobi->a.da[0][2]*res[2]-jacobi->a.da[1][2]*jacobi->a.da[2][0]*res[0]+res[1]*jacobi->a.da[2][0]*jacobi->a.da[0][2]+jacobi->a.da[1][0]*res[0]*jacobi->a.da[2][2])/
            (jacobi->a.da[0][0]*jacobi->a.da[1][1]*jacobi->a.da[2][2]-jacobi->a.da[0][0]*jacobi->a.da[2][1]*jacobi->a.da[1][2]-jacobi->a.da[1][0]*jacobi->a.da[0][1]*jacobi->a.da[2][2]+jacobi->a.da[2][1]*jacobi->a.da[1][0]*jacobi->a.da[0][2]-jacobi->a.da[1][1]*jacobi->a.da[2][0]*jacobi->a.da[0][2]+jacobi->a.da[2][0]*jacobi->a.da[0][1]*jacobi->a.da[1][2]);
	  deltat=(jacobi->a.da[2][1]*jacobi->a.da[1][0]*res[0]-jacobi->a.da[0][0]*jacobi->a.da[2][1]*res[1]+jacobi->a.da[0][0]*jacobi->a.da[1][1]*res[2]-jacobi->a.da[1][1]*jacobi->a.da[2][0]*res[0]-jacobi->a.da[1][0]*jacobi->a.da[0][1]*res[2]+jacobi->a.da[2][0]*jacobi->a.da[0][1]*res[1])/
            (jacobi->a.da[0][0]*jacobi->a.da[1][1]*jacobi->a.da[2][2]-jacobi->a.da[0][0]*jacobi->a.da[2][1]*jacobi->a.da[1][2]-jacobi->a.da[1][0]*jacobi->a.da[0][1]*jacobi->a.da[2][2]+jacobi->a.da[2][1]*jacobi->a.da[1][0]*jacobi->a.da[0][2]-jacobi->a.da[1][1]*jacobi->a.da[2][0]*jacobi->a.da[0][2]+jacobi->a.da[2][0]*jacobi->a.da[0][1]*jacobi->a.da[1][2]);

          /*advance the values for r,s,t approaching the solution*/
	  r=r+deltar;
	  s=s+deltas;
	  t=t+deltat;
          /*calculate the error*/
	  error=sqrt(res[0]*res[0]+res[1]*res[1]+res[2]*res[2]);
        }

        /*The local coordinates r,s,t of the fluid node related to the struct element are now known.
         *Use them to find the position in relation to the struct element*/

        /*criterion for tet4 r,s,t>0 && r+s+t<=1 -> fluid node inside struct element*/
        if (r>=-EPS && s >=-EPS && t>=-EPS && (r+s+t)<=1+EPS)
        {
          /*host element found!*/
          found=1;
        }
        else
        {
          /*struct element is not the host element*/
          found=0;
        }
        break;
      }/*end of if(TET4)*/

      /*TET10 element*/
      case tet10:
      {
        for (m=0;m<3;m++)     /*loop all dimensions*/
        {
          coeff->a.da[m][0]=x->a.da[m][0];
          coeff->a.da[m][1]=-x->a.da[m][1]+4*x->a.da[m][4]-3*x->a.da[m][0];
          coeff->a.da[m][2]=-x->a.da[m][2]+4*x->a.da[m][6]-3*x->a.da[m][0];
          coeff->a.da[m][3]=-3*x->a.da[m][0]-x->a.da[m][3]+4*x->a.da[m][9];
          coeff->a.da[m][4]=-4*x->a.da[m][6]+4*x->a.da[m][5]-4*x->a.da[m][4]+4*x->a.da[m][0];
          coeff->a.da[m][5]=4*x->a.da[m][7]-4*x->a.da[m][4]+4*x->a.da[m][0]-4*x->a.da[m][9];
          coeff->a.da[m][6]=4*x->a.da[m][0]+4*x->a.da[m][8]-4*x->a.da[m][9]-4*x->a.da[m][6];
          coeff->a.da[m][7]=2*x->a.da[m][1]-4*x->a.da[m][4]+2*x->a.da[m][0];
          coeff->a.da[m][8]=2*x->a.da[m][2]-4*x->a.da[m][6]+2*x->a.da[m][0];
          coeff->a.da[m][9]=2*x->a.da[m][0]+2*x->a.da[m][3]-4*x->a.da[m][9];
        }

        /*solve the nonlinear equation system for the natural coordinats of the fluid node*/
        /*loop until the error is smaller than the specified value*/
        while (error>ERROR)
        {
          /*initialize the Jacobian matrix*/
          jacobi->a.da[0][0]=coeff->a.da[0][1]+coeff->a.da[0][4]*s+coeff->a.da[0][5]*t+coeff->a.da[0][7]*2*r;
          jacobi->a.da[0][1]=coeff->a.da[0][2]+coeff->a.da[0][4]*r+coeff->a.da[0][6]*t+coeff->a.da[0][8]*2*s;
          jacobi->a.da[0][2]=coeff->a.da[0][3]+coeff->a.da[0][5]*r+coeff->a.da[0][6]*s+coeff->a.da[0][9]*2*t;
          jacobi->a.da[1][0]=coeff->a.da[1][1]+coeff->a.da[1][4]*s+coeff->a.da[1][5]*t+coeff->a.da[1][7]*2*r;
          jacobi->a.da[1][1]=coeff->a.da[1][2]+coeff->a.da[1][4]*r+coeff->a.da[1][6]*t+coeff->a.da[1][8]*2*s;
          jacobi->a.da[1][2]=coeff->a.da[1][3]+coeff->a.da[1][5]*r+coeff->a.da[1][6]*s+coeff->a.da[1][9]*2*t;
          jacobi->a.da[2][0]=coeff->a.da[2][1]+coeff->a.da[2][4]*s+coeff->a.da[2][5]*t+coeff->a.da[2][7]*2*r;
          jacobi->a.da[2][1]=coeff->a.da[2][2]+coeff->a.da[2][4]*r+coeff->a.da[2][6]*t+coeff->a.da[2][8]*2*s;
          jacobi->a.da[2][2]=coeff->a.da[2][3]+coeff->a.da[2][5]*r+coeff->a.da[2][6]*s+coeff->a.da[2][9]*2*t;

          /*calculate the residual*/
          for (k=0;k<3;k++)    /*loop all dimensions*/
          {
            res[k]=actfnode->x[k]-(coeff->a.da[k][0]+coeff->a.da[k][1]*r+coeff->a.da[k][2]*s+coeff->a.da[k][3]*t+coeff->a.da[k][4]*r*s+coeff->a.da[k][5]*r*t+coeff->a.da[k][6]*s*t+coeff->a.da[k][7]*r*r+coeff->a.da[k][8]*s*s+coeff->a.da[k][9]*t*t);
          }
          /*solving the linear equation system explicitly for the delta vector*/
          deltar =((-jacobi->a.da[0][1]*jacobi->a.da[1][2]*res[2]+jacobi->a.da[0][1]*jacobi->a.da[2][2]*res[1]-jacobi->a.da[0][2]*jacobi->a.da[2][1]*res[1]+jacobi->a.da[0][2]*res[2]*jacobi->a.da[1][1]-res[0]*jacobi->a.da[2][2]*jacobi->a.da[1][1]+res[0]*jacobi->a.da[1][2]*jacobi->a.da[2][1])/
                   (-jacobi->a.da[0][0]*jacobi->a.da[2][2]*jacobi->a.da[1][1]+jacobi->a.da[2][0]*jacobi->a.da[0][2]*jacobi->a.da[1][1]-jacobi->a.da[1][2]*jacobi->a.da[2][0]*jacobi->a.da[0][1]-jacobi->a.da[1][0]*jacobi->a.da[0][2]*jacobi->a.da[2][1]+jacobi->a.da[2][2]*jacobi->a.da[1][0]*jacobi->a.da[0][1]+jacobi->a.da[0][0]*jacobi->a.da[1][2]*jacobi->a.da[2][1]));
          deltas =-((-jacobi->a.da[0][0]*jacobi->a.da[1][2]*res[2]+jacobi->a.da[0][0]*jacobi->a.da[2][2]*res[1]-jacobi->a.da[2][2]*jacobi->a.da[1][0]*res[0]-jacobi->a.da[2][0]*jacobi->a.da[0][2]*res[1]+jacobi->a.da[1][2]*jacobi->a.da[2][0]*res[0]+jacobi->a.da[1][0]*jacobi->a.da[0][2]*res[2])/
                    (-jacobi->a.da[0][0]*jacobi->a.da[2][2]*jacobi->a.da[1][1]+jacobi->a.da[2][0]*jacobi->a.da[0][2]*jacobi->a.da[1][1]-jacobi->a.da[1][2]*jacobi->a.da[2][0]*jacobi->a.da[0][1]-jacobi->a.da[1][0]*jacobi->a.da[0][2]*jacobi->a.da[2][1]+jacobi->a.da[2][2]*jacobi->a.da[1][0]*jacobi->a.da[0][1]+jacobi->a.da[0][0]*jacobi->a.da[1][2]*jacobi->a.da[2][1]));
          deltat =-((-jacobi->a.da[0][0]*jacobi->a.da[2][1]*res[1]+jacobi->a.da[0][0]*res[2]*jacobi->a.da[1][1]+jacobi->a.da[2][1]*jacobi->a.da[1][0]*res[0]+jacobi->a.da[2][0]*jacobi->a.da[0][1]*res[1]-jacobi->a.da[2][0]*res[0]*jacobi->a.da[1][1]-res[2]*jacobi->a.da[1][0]*jacobi->a.da[0][1])/
                    (-jacobi->a.da[0][0]*jacobi->a.da[2][2]*jacobi->a.da[1][1]+jacobi->a.da[2][0]*jacobi->a.da[0][2]*jacobi->a.da[1][1]-jacobi->a.da[1][2]*jacobi->a.da[2][0]*jacobi->a.da[0][1]-jacobi->a.da[1][0]*jacobi->a.da[0][2]*jacobi->a.da[2][1]+jacobi->a.da[2][2]*jacobi->a.da[1][0]*jacobi->a.da[0][1]+jacobi->a.da[0][0]*jacobi->a.da[1][2]*jacobi->a.da[2][1]));

          /*advance the values for r,s,t approaching the solution*/
	  r=r+deltar;
	  s=s+deltas;
	  t=t+deltat;
          /*calculate the error*/
	  error=sqrt(res[0]*res[0]+res[1]*res[1]+res[2]*res[2]);
        }

        /*The local coordinates r,s,t of the fluid node related to the struct element are now known.
         *Use them to find the position in relation to the struct element*/

        /*criterion for tet10: r,s,t>0 && r+s+t<=1 -> fluid node inside struct element*/
        if (r>=0 && s >=0 && t>=0 && (r+s+t)<=1)
        {
          /*host element found!*/
          found=1;
        }
        else
        {
          /*struct element is not the host element*/
          found=0;
        }
        break;
      }/*end of if(TET10)*/

      case tri3:
      {
        INT i, j;
        DOUBLE A[3][3];
        DOUBLE rhs[3];
        DOUBLE det;
        DOUBLE b00,b01,b02,b10,b11,b12,b20,b21,b22;
        DOUBLE detinv;

        /* use the first two dimensions to build a linear system */

        A[0][0] = x->a.da[0][1] - x->a.da[0][0];
        A[0][1] = x->a.da[0][2] - x->a.da[0][0];
        A[0][2] = 1;
        A[1][0] = x->a.da[1][1] - x->a.da[1][0];
        A[1][1] = x->a.da[1][2] - x->a.da[1][0];
        A[1][2] = 1;
        A[2][0] = x->a.da[2][1] - x->a.da[2][0];
        A[2][1] = x->a.da[2][2] - x->a.da[2][0];
        A[2][2] = 1;

        rhs[0] = actfnode->x[0] - x->a.da[0][0];
        rhs[1] = actfnode->x[1] - x->a.da[1][0];
        rhs[2] = actfnode->x[2] - x->a.da[2][0];

        b00 = A[0][0];
        b01 = A[0][1];
        b02 = A[0][2];
        b10 = A[1][0];
        b11 = A[1][1];
        b12 = A[1][2];
        b20 = A[2][0];
        b21 = A[2][1];
        b22 = A[2][2];

        A[0][0] =   b11*b22 - b21*b12;
        A[1][0] = - b10*b22 + b20*b12;
        A[2][0] =   b10*b21 - b20*b11;
        A[0][1] = - b01*b22 + b21*b02;
        A[1][1] =   b00*b22 - b20*b02;
        A[2][1] = - b00*b21 + b20*b01;
        A[0][2] =   b01*b12 - b11*b02;
        A[1][2] = - b00*b12 + b10*b02;
        A[2][2] =   b00*b11 - b10*b01;

        det = b00*A[0][0]+b01*A[1][0]+b02*A[2][0];
        detinv = 1.0/det;
        for (i=0; i<3; i++)
          for (j=0; j<3; j++)
            A[i][j] *= detinv;

        /* solve linear system */

        r = rhs[0]*A[0][0] + rhs[1]*A[0][1] + rhs[2]*A[0][2];
        s = rhs[0]*A[1][0] + rhs[1]*A[1][1] + rhs[2]*A[1][2];
        t = rhs[0]*A[2][0] + rhs[1]*A[2][1] + rhs[2]*A[2][2];

        /* check third equation and ranges */
        if ((r >= -EPS) &&
            (s >= -EPS) &&
            (r+s <= 1+2*EPS) &&
            (t >= -EPS) &&
            (t <= EPS))
        {
          found=1;
        }
        else
        {
          found=0;
        }
        break;
      }

      default:
        dserror("unsupported element type '%d'", hostele->distyp);
      }

      /*free space for arrays*/
      amdel(x);
      amdel(coeff);
      amdel(jacobi);

      /*if no host element was found, determine the distance to the current one*/
      if (found==0)
      {
        distance=sqrt(r*r+s*s+t*t);

        /*if this element is closer than the one before*/
        if (distance<actdistance)
        {
          actdistance=distance;
          /*store data in temporary variables*/
          acthostele=hostele;
          actr=r;
          acts=s;
          actt=t;
        }
      }

      if (found==1) break;
    }/*end of lopp over elements belonging to struct node*/

    /*take the next node in the nodelist of the octree partition*/
    nodelist=nodelist->next;
    /*if not the last node in the octree partition*/
    if (nodelist!=NULL)
       actsnode=nodelist->node;

    if (found==1) break;
  }/*end of loop over structure nodes in octree*/

  /*if there was no host element found, take the closest one*/
  if (found==0 && actdistance<DISTLIMIT)
  {
    /*copy data stored in temporary variables*/
    hostele=acthostele;
    r=actr;
    s=acts;
    t=actt;

    /*host element found!*/
    found=1;
/*printf("fluid node %d coupled with closest struct element %d\n", actfnode->Id, hostele->Id);*/
/*printf("r %f,s %f,t %f\n", r, s, t);*/
/*dserror("fluid node %d coupled with closest struct element %d\n", actfnode->Id, hostele->Id);*/
  }

  /*if a host element was found for this fluid node*/
  if (found==1)
  {
    /*write the coupling information for this fluid node*/
    actfnode->gnode->coupleptr->hostele=hostele;
    actfnode->gnode->coupleptr->x[0]=r;
    actfnode->gnode->coupleptr->x[1]=s;
    actfnode->gnode->coupleptr->x[2]=t;
    /*write the coupling information for the structure host element*/
    hostele->coupleptr->couplenode[hostele->coupleptr->numnp]=actfnode;
    hostele->coupleptr->numnp++;
  }

  /*if no host element could be found for this fluid node*/
  if (found==0 && actdistance==DISTLIMIT)
    dserror("Fluid node %d could not be paired with an element - distance: %f too large -> remesh!\n", actfnode->Id, actdistance);
}

/*-----------------------------------------------------------------------*/
/*!
  \brief create the octree for the FSI structure nodes

  recursive function to subdivide the partitions of the octree until the maximum
  number of nodes per partition (MAXNODEPART) is arrived


  The current octree partition is subdivided into two partitions by dividing the
  longest edge of the mother partition. The position on the edge where the partition
  is divided is calculated as the average of the coordinates of all contained
  structure nodes in this direction. So there will never be an empty child partition
  created. A function is called to sort the struct nodes into the children partitions.
  This function will call itself until the number of struct nodes in all leaves is
  smaller than the given limit MAXNODEPART.

  \return void

  \author henke
  \date   05/06
*/
/*-----------------------------------------------------------------------*/
void fsi_create_octree (OCTREE *octreetemp)
{
  OCTREE        *octreemem;        /*temporary pointer for memory allocation*/
  DOUBLE        x1,x2,x3;          /*length of the edges of the octree partition*/
  DOUBLE        xdivide;           /*position on the longest edge where the partition is divided*/
  INT           i;                 /*counter*/
  INT           geo;               /* flag (0,1,2)*/
  NODELIST      *nodelisttemp;     /*temporary pointer*/

  /*if there are too many nodes in this octree partition -> subdivide!*/
  if (octreetemp->numnpoct>MAXNODEPART)
  {
    /*allocate space for two children partitions*/
    octreemem=(OCTREE*)CCACALLOC(2,sizeof(OCTREE));

    /*initialisation of the children octree partitions*/
    octreemem[0].layeroct=(octreetemp->layeroct)+1;
    octreemem[1].layeroct=(octreetemp->layeroct)+1;
    octreemem[0].numnpoct=0;
    octreemem[1].numnpoct=0;

    /*calculate the lengths of the edges of the mother partition*/
    x1=(octreetemp->xoct[1])-(octreetemp->xoct[0]);
    x2=(octreetemp->xoct[3])-(octreetemp->xoct[2]);
    x3=(octreetemp->xoct[5])-(octreetemp->xoct[4]);

    xdivide=0;
    nodelisttemp=octreetemp->nodelist;

    /*x-edge is the longest*/
    if (x1>=x2 && x1>=x3)
    {
      geo=0;
      /*calculate the average coordinate of all nodes in octree partition in direction of the longest edge*/
      for (i=0;i<octreetemp->numnpoct;i++)
      {
        xdivide += nodelisttemp->node->x[geo];
        nodelisttemp=nodelisttemp->next;
      }
      xdivide=xdivide/(octreetemp->numnpoct);

      /*pausibility check: xdivide has to be between left and right border of partition*/
      if (!((octreetemp->xoct[0]<=xdivide)&&(xdivide<=octreetemp->xoct[1])))
      {
        dserror("invalid division of octree partition\n");
      }

      /*coordinates of edge that is divided + tolerance*/
      octreemem[0].xoct[0]=octreetemp->xoct[0]-EPS;
      octreemem[0].xoct[1]=xdivide+EPS;
      octreemem[1].xoct[0]=xdivide-EPS;
      octreemem[1].xoct[1]=octreetemp->xoct[1]+EPS;
      /*coordinates of edges that are not divided + tolerace*/
      octreemem[0].xoct[2]=octreetemp->xoct[2]-EPS;
      octreemem[0].xoct[3]=octreetemp->xoct[3]+EPS;
      octreemem[0].xoct[4]=octreetemp->xoct[4]-EPS;
      octreemem[0].xoct[5]=octreetemp->xoct[5]+EPS;
      octreemem[1].xoct[2]=octreetemp->xoct[2]-EPS;
      octreemem[1].xoct[3]=octreetemp->xoct[3]+EPS;
      octreemem[1].xoct[4]=octreetemp->xoct[4]-EPS;
      octreemem[1].xoct[5]=octreetemp->xoct[5]+EPS;

      /*sort struct nodes into children partitions*/
      fsi_sort_struct_nodes(octreetemp,&octreemem[0],&octreemem[1],xdivide,geo);
    }
    /*y-edge is the longest*/
    if (x2>x1 && x2>=x3)
    {
      geo=1;
      /*calculate the average coordinate of all nodes in octree partition indirection of the longest edge*/
      for (i=0;i<octreetemp->numnpoct;i++)
      {
        xdivide += nodelisttemp->node->x[geo];
        nodelisttemp=nodelisttemp->next;
      }
      xdivide=xdivide/(octreetemp->numnpoct);

      /*pausibility check: xdivide has to be between left and right border of partition*/
      if (!((octreetemp->xoct[2]<=xdivide)&&(xdivide<=octreetemp->xoct[3])))
        dserror("invalid division of octree partition\n");

      /*coordinates of edge that is divided + tolerance*/
      octreemem[0].xoct[2]=octreetemp->xoct[2]-EPS;
      octreemem[0].xoct[3]=xdivide+EPS;
      octreemem[1].xoct[2]=xdivide-EPS;
      octreemem[1].xoct[3]=octreetemp->xoct[3]+EPS;
      /*coordinates of edges that are not divided + tolerance*/
      octreemem[0].xoct[0]=octreetemp->xoct[0]-EPS;
      octreemem[0].xoct[1]=octreetemp->xoct[1]+EPS;
      octreemem[0].xoct[4]=octreetemp->xoct[4]-EPS;
      octreemem[0].xoct[5]=octreetemp->xoct[5]+EPS;
      octreemem[1].xoct[0]=octreetemp->xoct[0]-EPS;
      octreemem[1].xoct[1]=octreetemp->xoct[1]+EPS;
      octreemem[1].xoct[4]=octreetemp->xoct[4]-EPS;
      octreemem[1].xoct[5]=octreetemp->xoct[5]+EPS;

      /*sort struct nodes into children partitions*/
      fsi_sort_struct_nodes(octreetemp, &octreemem[0], &octreemem[1],xdivide,geo);
    }
    /*z-edge is the longest*/
    if (x3>x1 && x3>x2)
    {
      geo=2;
      /*calculate the average coordinate of all nodes in octree partition in direction of the longest edge*/
      for (i=0;i<octreetemp->numnpoct;i++)
      {
        xdivide += nodelisttemp->node->x[geo];
        nodelisttemp=nodelisttemp->next;
      }
      xdivide=xdivide/(octreetemp->numnpoct);

      /*pausibility check: xdivide has to be between left and right border of partition*/
      if (!((octreetemp->xoct[4]<=xdivide)&&(xdivide<=octreetemp->xoct[5])))
        dserror("invalid division of octree partition\n");

      /*coordinates of edge that is divided + tolerance*/
      octreemem[0].xoct[4]=octreetemp->xoct[4]-EPS;
      octreemem[0].xoct[5]=xdivide+EPS;
      octreemem[1].xoct[4]=xdivide-EPS;
      octreemem[1].xoct[5]=octreetemp->xoct[5]+EPS;
      /*coordinates of edges that are not divided + tolerace*/
      octreemem[0].xoct[0]=octreetemp->xoct[0]-EPS;
      octreemem[0].xoct[1]=octreetemp->xoct[1]+EPS;
      octreemem[0].xoct[2]=octreetemp->xoct[2]-EPS;
      octreemem[0].xoct[3]=octreetemp->xoct[3]+EPS;
      octreemem[1].xoct[0]=octreetemp->xoct[0]-EPS;
      octreemem[1].xoct[1]=octreetemp->xoct[1]+EPS;
      octreemem[1].xoct[2]=octreetemp->xoct[2]-EPS;
      octreemem[1].xoct[3]=octreetemp->xoct[3]+EPS;

      /*sort struct nodes into children partitions*/
      fsi_sort_struct_nodes(octreetemp, &octreemem[0], &octreemem[1],xdivide,geo);
    }
    /*adding the children partitions to the octree structure*/
    octreetemp->next[0]=&octreemem[0];
    octreetemp->next[1]=&octreemem[1];

    /*recursive call of the fsi_create_octree function*/
    fsi_create_octree((octreetemp->next[0]));
    fsi_create_octree((octreetemp->next[1]));
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
void fsi_initcoupling(
    FIELD              *structfield,
    INT                 disnum_s,
    FIELD              *fluidfield,
    INT                 disnum_f,
    FIELD              *alefield,
    INT                 disnum_a
    )

{

  INT     numfnp, numsnp, numanp;              /* number of nodes         */
  INT     numdf;                               /* number of dofs          */
  INT     numc=0;                              /* number of columns in mf */
  INT     numaf,numsf,numff;
  INT     i,j;                                 /* simply some counters    */
  INT     ierr;                                /* flag                    */
  INT     dim;                                 /* dimension of problem    */
  INT     sfound,afound;                       /* flag                    */
  INT     is_ale;


  DOUBLE  tol=EPS4;                            /* tolerance for node dist */

  NODE   *actfnode=NULL, *actsnode=NULL, *actanode=NULL;  /* actual nodes */
  GNODE  *actfgnode=NULL,*actsgnode=NULL,*actagnode=NULL; /* actual gnodes*/

  ELEMENT *actele;

  ARRAY   aindex_a;
  INT    *aindex;

  FSI_DYNAMIC    *fsidyn;

/*********************start*new************************************************************/
  OCTREE                octree;                            /* head of the octree structure */
  NODELIST              *nodelisttemp,*nodelistmem;        /* temporary pointers */
  INT                   match=0;                           /* flag for plausibility check */
  DOUBLE                difference=0;                      /* variable for plausibility check */
/************end*new***********************************************************************/

#ifdef DEBUG
  dstrc_enter("fsi_initcoupling");
#endif


  fsidyn= alldyn[3].fsidyn;


  /* find number of nodes in different fields */
  numsnp  = structfield->dis[disnum_s].numnp;
  numfnp  = fluidfield->dis[disnum_f].numnp;
  numanp  = alefield->dis[disnum_a].numnp;
  numaf   = genprob.numaf;
  numff   = genprob.numff;
  numsf   = genprob.numsf;
  dim     = genprob.ndim;

#ifdef PERF
perf_begin(53);
#endif



  /* allocate space for mulitfield solution history:
   * -----------------------------------------------
   *
   * the following data are necassary:
   * ale needs displacements of the structural fsi coupling nodes
   * fluid needs velocity of all ale nodes
   * structure needs stresses of the fluid fsi coupling nodes
   */


  /* multifield solution history of fluid nodes:
   * -------------------------------------------
   *   actfnode->sol_mf.a.da[0][i]: velocity transfered to ale
   *   actfnode->sol_mf.a.da[1][i]: stresses transfered to structure
   *     2D: SIGMA11,SIGMA22,SIGMA12
   *     3D: SIGMA11,SIGMA22,SIGMA33,SIGMA12,SIGMA13,SIGMA23
   */


  /* loop fluid nodes */
  for (i=0;i<numfnp;i++)
  {
    actfnode  = &(fluidfield->dis[disnum_f].node[i]);
    actfgnode = actfnode->gnode;
    numdf = actfnode->numdf;

    if (numdf==3 || numdf==5) numc=numdf;
    else if (numdf==4 || numdf==7) numc=IMAX(6,numdf);
    else  dserror("number of fluid dofs not possible!\n");

    amdef("sol_mf",&actfnode->sol_mf,2,numc,"DA");
    amzero(&(actfnode->sol_mf));

    actfgnode->mfcpnode=(NODE**)CCACALLOC(3,sizeof(NODE*));

    if (actfgnode->mfcpnode==NULL)
      dserror("Allocation of coupling node pointers failed");

    for (j=0;j<3;j++) actfgnode->mfcpnode[j]=NULL;

  } /* end of loop over fluid nodes */



  /* multifield solution history of structural nodes:
   * ------------------------------------------------
   * actsnode->sol_mf.a.da[0][i]: actual displacements transfered to ale
   * actsnode->sol_mf.a.da[1][i]: displacements of old iteration step
   * actsnode->sol_mf.a.da[2][i]: displacements of old time step
   * actsnode->sol_mf.a.da[3][i]: dispi
   * actsnode->sol_mf.a.da[4][i]: couplingforces at the end of time step
   * actsnode->sol_mf.a.da[5][i]: coupling forces at the beginning of time step
   */

  /* loop structure nodes */
  for (i=0;i<numsnp;i++)
  {
    actsnode  = &(structfield->dis[disnum_s].node[i]);
    actsgnode = actsnode->gnode;
    numdf = actsnode->numdf;

    amdef("sol_mf",&actsnode->sol_mf,6,numdf,"DA");
    amzero(&(actsnode->sol_mf));

    actsgnode->mfcpnode=(NODE**)CCACALLOC(3,sizeof(NODE*));

    if (actsgnode->mfcpnode==NULL)
      dserror("Allocation of coupling node pointers failed");

    for (j=0;j<3;j++) actsgnode->mfcpnode[j]=NULL;

  } /* end of loop over struct nodes */



  /* multifield solution history of ale nodes:
   * -----------------------------------------
   * actanode->sol_mf.a.da[0][i]: displacements at (n)
   * actanode->sol_mf.a.da[1][i]: displacements at (n+1)
   */

  /* loop ale nodes */
  for (i=0;i<numanp;i++)
  {
    actanode  = &(alefield->dis[disnum_a].node[i]);
    actagnode = actanode->gnode;
    numdf = actanode->numdf;

    amdef("sol_mf",&actanode->sol_mf,3,numdf,"DA");
    amzero(&(actanode->sol_mf));

    actagnode->mfcpnode=(NODE**)CCACALLOC(3,sizeof(NODE*));

    if (actagnode->mfcpnode==NULL)
      dserror("Allocation of coupling node pointers failed");

    for (j=0;j<3;j++) actagnode->mfcpnode[j]=NULL;

  } /* end of loop over ale nodes */


#ifdef PERF
perf_end(53);
#endif

#ifdef PERF
perf_begin(54);
#endif

/********************start*new*******************************/

    /*initialize the octree structure*/
    octree.layeroct = 0;
    octree.numnpoct = 0;
    /*the initial values of the octree coordinates are arbitrarily
     * set to those of the first FSI structure node*/
    for (i=0;i<numsnp;i++)
    {
      actsnode=&(structfield->dis[disnum_s].node[i]);
      actsgnode=actsnode->gnode;
      if (actsgnode->fsicouple!=NULL) break;
    }
    octree.xoct[0]=actsnode->x[0];
    octree.xoct[1]=actsnode->x[0];
    octree.xoct[2]=actsnode->x[1];
    octree.xoct[3]=actsnode->x[1];
    octree.xoct[4]=actsnode->x[2];
    octree.xoct[5]=actsnode->x[2];
    octree.nodelist=NULL;
    octree.next[0]=NULL;
    octree.next[1]=NULL;

    /*loop all fluid nodes to look for size of the FSI interface*/
    for(i=0;i<numfnp;i++)
    {
      actfnode = &(fluidfield->dis[disnum_f].node[i]);
      actfgnode = actfnode->gnode;

      /*no FSI node*/
      if (actfgnode->fsicouple==NULL) continue;

      /*Find the min. and max. coordinates of the FSI interface*/
      if (actfnode->x[0]<octree.xoct[0]) octree.xoct[0]=actfnode->x[0];
      if (actfnode->x[0]>octree.xoct[1]) octree.xoct[1]=actfnode->x[0];
      if (actfnode->x[1]<octree.xoct[2]) octree.xoct[2]=actfnode->x[1];
      if (actfnode->x[1]>octree.xoct[3]) octree.xoct[3]=actfnode->x[1];
      if (actfnode->x[2]<octree.xoct[4]) octree.xoct[4]=actfnode->x[2];
      if (actfnode->x[2]>octree.xoct[5]) octree.xoct[5]=actfnode->x[2];
    }

    /*loop all structure nodes to look for size of the FSI interface*/
    for(i=0;i<numsnp;i++)
    {
      actsnode = &(structfield->dis[disnum_s].node[i]);
      actsgnode = actsnode->gnode;

      /*no FSI node*/
      if (actsgnode->fsicouple==NULL) continue;

      /*Find the min. and max. coordinates of the FSI interface*/
      if (actsnode->x[0]<octree.xoct[0]) octree.xoct[0]=actsnode->x[0];
      if (actsnode->x[0]>octree.xoct[1]) octree.xoct[1]=actsnode->x[0];
      if (actsnode->x[1]<octree.xoct[2]) octree.xoct[2]=actsnode->x[1];
      if (actsnode->x[1]>octree.xoct[3]) octree.xoct[3]=actsnode->x[1];
      if (actsnode->x[2]<octree.xoct[4]) octree.xoct[4]=actsnode->x[2];
      if (actsnode->x[2]>octree.xoct[5]) octree.xoct[5]=actsnode->x[2];

      /*Create the initial nodelist containing all structure nodes*/
      nodelistmem=(NODELIST*)CCACALLOC(1,sizeof(NODELIST));
      nodelistmem->node=actsnode;
      nodelistmem->next=NULL;

      /*if adding first node*/
      if (octree.numnpoct==0)
      {
        octree.nodelist=nodelistmem;
        nodelisttemp=octree.nodelist;
      }
      /*if adding any other node but the first one*/
      else
      {
        nodelisttemp->next=nodelistmem;
        nodelisttemp=nodelisttemp->next;
      }

      octree.numnpoct++;
    }
    /*add a tolerance to the FSI interface*/
    octree.xoct[0]=octree.xoct[0]-EPS;
    octree.xoct[1]=octree.xoct[1]+EPS;
    octree.xoct[2]=octree.xoct[2]-EPS;
    octree.xoct[3]=octree.xoct[3]+EPS;
    octree.xoct[4]=octree.xoct[4]-EPS;
    octree.xoct[5]=octree.xoct[5]+EPS;

/* printf("Gebiet\n");
   printf("x-Koord unten: %f\n",octree.xoct[0]);
   printf("x-Koord oben: %f\n",octree.xoct[1]);
   printf("y-Koord unten: %f\n",octree.xoct[2]);
   printf("y-Koord oben: %f\n",octree.xoct[3]);
   printf("z-Koord unten: %f\n",octree.xoct[4]);
   printf("z-Koord oben: %f\n",octree.xoct[5]);*/

#if 0
    /*plausibility check*/
    dsassert(octree.numnpoct==0,"No nodes in the first octree partition\n");
    dsassert(octree.nodelist==NULL,"There was no nodelist created\n");
#endif

    /*build the octree for all FSI structure nodes*/
    fsi_create_octree(&octree);

/*printf("Octree successfully built!\n");*/

/*********************end*new*********************************/

  /* create and initialise index arrays */
  aindex = amdef("aindex",&aindex_a,numanp,1,"IV");
  for (i=0;i<numanp;i++) aindex[i]=1;
/*   sindex = amdef("sindex",&sindex_a,numsnp,1,"IV"); */
/*   for (i=0;i<numsnp;i++) sindex[i]=1; */

  /* loop fluid nodes  */
  for (i=0;i<numfnp;i++)
  {
    actfnode  = &(fluidfield->dis[disnum_f].node[i]);
    actfgnode = actfnode->gnode;
    is_ale=0;

    for (j=0;j<actfnode->numele;j++)
    {

      actele= actfnode->element[j];

      switch (actele->eltyp)
      {

#ifdef D_FLUID2
        case el_fluid2:
          if (actele->e.f2->is_ale>0) is_ale++;
          break;
#endif


#ifdef D_FLUID3
        case el_fluid3:
          if (actele->e.f3->is_ale>0) is_ale++;
          break;

#endif


#ifdef D_FLUID3_F
        case el_fluid3_fast:
          if (actele->e.f3->is_ale>0) is_ale++;
          break;
#endif

        default:
          dserror("eltyp unknow\n");

      }  /* switch (actele->eltyp) */

      if (is_ale>0) break;

    }  /* for (j=0;j<actfnode->numele;j++) */
       /*falls am aktuellen fluid knoten ein element haengt, das teil des ALEs ist, dann steht jetzt is_ale auf 1*/
       /*das heisst, dass ich jetzt weiss, ob ich einen fluid knoten habe, der teil des Ales ist*/
       /*printf("fluid node %d: is_ale: %d\n", actfnode->Id, is_ale);*/
    sfound=0;
    afound=0;


#ifdef PERF
perf_begin(56);
#endif

/************************************start*new*******************************/
    /*couple FSI fluid node with a structure element*/
    if (actfgnode->fsicouple != NULL)
    {
      /*allocate space for the coupling information*/
      actfgnode->coupleptr=(FSI_COUPLE_ELEMENT*)CCACALLOC(1,sizeof(FSI_COUPLE_ELEMENT));
      actfgnode->coupleptr->hostele=NULL;

      /*initialization with arbitrary values for local coordintes r,s,t*/
      for (j=0;j<3;j++) actfgnode->coupleptr->x[j]=2;

      /*find the structure host element for this fluid node*/
      fsi_find_hostelement(&octree,actfnode);
    }
/**********************************end*new************************************/

#ifdef PERF
    perf_end(56);
#endif


#ifdef PERF
    perf_begin(57);
#endif

    if (genprob.create_ale != 1 || fluidfield->subdivide > 0)  /*??????*/
    {
      /* loop ale nodes and find corresponding node */
      if (is_ale>0)
        for (j=0;j<numanp;j++)
        {

          if (aindex[j]==0) continue;  /*falls der index == 0 ist, dann hab ich schoen einen coupling node gefunden*/

          actanode  = &(alefield->dis[disnum_a].node[j]);

          cheque_distance(&(actfnode->x[0]),&(actanode->x[0]),tol,&ierr);
          if(ierr==0) continue;   /*falls die distance zu groß ist -> nächster ale node*/

          afound++;
          actagnode = actanode->gnode;
          aindex[j]=0;

          break;

        } /* end of loop over ale nodes */
    }

    else
    {
      if (is_ale>0)
      {
        for (j=0; j<actfnode->element[0]->numnp; j++)
        {
          if (actfnode->Id == actfnode->element[0]->node[j]->Id )
          {
            switch (actfnode->element[0]->eltyp)
            {
#ifdef D_FLUID3
              case el_fluid3:
              case el_fluid3_fast:
                actanode = actfnode->element[0]->e.f3->ale_ele->node[j];
                break;
#endif

#ifdef D_FLUID2
              case el_fluid2:
                actanode = actfnode->element[0]->e.f2->ale_ele->node[j];
                break;
#endif

              default:
                dserror("");
                break;
            }  /* switch (actfnode->element[0]->eltyp) */

            break;  /* break the loop over the nodes of the element */

          }  /* if (actfnode->Id == actfnode->element[0]->node[j]->Id ) */
        }


#ifdef DEBUG
        cheque_distance(&(actfnode->x[0]),&(actanode->x[0]),tol,&ierr);
        if(ierr==0) dserror("Wrong node matching for fluid <-> ale!!!");
#endif

        afound++;
        actagnode = actanode->gnode;

      }
    }

#ifdef PERF
perf_end(57);
#endif


/*********************start*change****************************/
    /* set pointers to corresponding nodes */
    /* so far all mfcouple fields are initialized with NULL,
     * therefore they are simply overwritten where necessary */
                   actfgnode->mfcpnode[numff]=actfnode;
    if(afound>0)   actfgnode->mfcpnode[numaf]=actanode;
    if(afound>0)   actagnode->mfcpnode[numff]=actfnode;
    if(afound>0)   actagnode->mfcpnode[numaf]=actanode;
/*********************end*change******************************/

/*****************start*new************************************/
    /*copy the coupling information from the FSI fluid node to the matching FSI ale node*/
    if (actfgnode->fsicouple != NULL)
    {
      if (afound>0)
      {
        actagnode->coupleptr=(FSI_COUPLE_ELEMENT*)CCACALLOC(1,sizeof(FSI_COUPLE_ELEMENT));
        actagnode->coupleptr->hostele=NULL;
        actagnode->coupleptr->hostele=actfgnode->coupleptr->hostele;
        for (j=0;j<3;j++)
        {
          actagnode->coupleptr->x[j]=actfgnode->coupleptr->x[j];
        }
      }
      else
      {
        dserror("no FSI coupling information for ale node, because there is no ale node for fluid node %d\n", actfnode->Id);
      }
    }
/*****************end*new**************************************/

  }  /* for (i=0;i<numfnp;i++) */

/* loop structure nodes */
/* (i=0;i<numsnp;i++)*/
/*{*/
/*  actsnode  = &(structfield->dis[disnum_s].node[i]);*/
/*  actsgnode = actsnode->gnode;*/
/*  numdf = actsnode->numdf;*/
/*  if (actsgnode->fsicouple!=NULL)*/
/*  {*/
/*    if (actsnode->Id==2359)*/
/*    {*/
/*      for (j=0;j<actsnode->numele;j++)*/
/*      {*/
/*        printf("element %d: \n", actsnode->element[j]->Id);*/
/*        for (m=0;m<actsnode->element[j]->coupleptr->numnp;m++)*/
/*        {*/
/*          printf("fluid node in this element: %d\n",actsnode->element[j]->coupleptr->couplenode[m]->Id);*/
/*        }*/
/*      }*/
/*    }*/
/*  }*/
/*}*/

  /*plausibility check: are ale and fluid mesh identical (for automtic generation CA 1)?
   *                    are ale and fluid mesh identical within a tolerance (for separate generation)?*/

  /* loop ale nodes */
  for (i=0;i<numanp;i++)
  {
    actanode = &(alefield->dis[disnum_a].node[i]);
    actagnode = actanode->gnode;
    numdf=actanode->numdf;
    match=0;

    /*loop fluid nodes*/
    for (j=0;j<numfnp;j++)
    {
      actfnode = &(fluidfield->dis[disnum_f].node[j]);
      actfgnode = actfnode->gnode;

      if ((-EPS<=(actanode->x[0]-actfnode->x[0])&&(actanode->x[0]-actfnode->x[0])<=EPS)&&
          (-EPS<=(actanode->x[1]-actfnode->x[1])&&(actanode->x[1]-actfnode->x[1])<=EPS)&&
          (-EPS<=(actanode->x[2]-actfnode->x[2])&&(actanode->x[2]-actfnode->x[2])<=EPS))
      {
/*#ifdef: automatic generation (CA 1)
 *
 *       difference=actanode->x[0]-actfnode->x[0]+actanode->x[1]-actfnode->x[1]+actanode->x[2]-actanode->x[2];
 *       if (difference!=0) dserror("ale and fluid mesh are not identical!\n");
 *       difference=0;
 *
 *#endif: automatic generation (CA 1)*/

        match=1;
      }
    }
    if (match==0) dserror("ale and fluid mesh are not identical!\n");
  }

#ifdef PERF
perf_end(54);
#endif

  amdel(&aindex_a);

  if (genprob.visual>0) goto end;
  if (fsidyn->coupmethod == 0) goto end;    /* INT    coupmethod; !< flag, 0=mortar , 1=conforming */


#ifdef PERF
perf_begin(55);
#endif



  /* plausibility checks */
  /* ------------------- */
  for (i=0;i<numfnp;i++)
  {
    actfnode  = &(fluidfield->dis[disnum_f].node[i]);
    actfgnode = actfnode->gnode;

    if (actfgnode->fsicouple==NULL) continue;

    /*actsnode = actfgnode->mfcpnode[numsf];*/
    actanode = actfgnode->mfcpnode[numaf];
    /*actsgnode = actsnode->gnode;*/
    actagnode = actanode->gnode;


/* if (actsgnode == NULL)   */
/*  dserror("No structure node for fluid node: #%d", actfnode->Id);  */


   if (actagnode == NULL)
      dserror("No ale node for fluid node: #%d", actfnode->Id);


    /* check locsys */                         /*  Was ist eine local system ID? */
    if (actfnode->locsysId!=0)
      dserror("No locsys at FSI coupling node: #%d", actfnode->Id);
/*    if (actsnode->locsysId!=0)*/
/*      dserror("No locsys at FSI coupling node: #%d", actsnode->Id);*/
    if (actanode->locsysId!=0)
      dserror("No locsys at FSI coupling node: #%d", actanode->Id);


    /* check coupling Ids */
/*    if (actfgnode->fsicouple->fsi_coupleId!=actsgnode->fsicouple->fsi_coupleId)*/
/*      dserror("FSI Coupling Condition Fluid-Struct: wrong coupleId, fluid #%d, struct #%d",*/
/*          actfnode->Id, actsnode->Id);*/

    if (actfgnode->fsicouple->fsi_coupleId!=actagnode->fsicouple->fsi_coupleId)
      dserror("FSI Coupling Condition Fluid-Ale: wrong coupleId, fluid #%d, ale #%d",
          actfnode->Id, actanode->Id);




    /* check coupling conditions */ /*Was genau bringt dieser Check? Conforming - Nonconforming?Oder Die Art der Felder paßt nicht zusammen?*/
/*    if (actsgnode->fsicouple==NULL)    */
/*    {   */
/*      dserror("FSI Coupling Condition Fluid-Struct not the same: struct node #%d",   */
/*          actsnode->Id);         */
/*    }      */

    if (actagnode->fsicouple==NULL)
      dserror("FSI Coupling Condition Fluid-Ale not the same: ale node #%d",
          actanode->Id);





    /* check mesh */
/*    if (actfgnode->fsicouple->fsi_mesh!=actsgnode->fsicouple->fsi_mesh)*/
/*      dserror("FSI Coupling Condition Fluid-Struct: wrong mesh, fluid #%d, struct #%d",*/
/*          actfnode->Id, actsnode->Id);*/

    if (actfgnode->fsicouple->fsi_mesh!=actagnode->fsicouple->fsi_mesh)
      dserror("FSI Coupling Condition Fluid-Ale: wrong mesh, fluid #%d, ale #%d",
          actfnode->Id, actanode->Id);

    /***********************start*new*********************/
    if (actfgnode->coupleptr->hostele==NULL)
      dserror("fluid node %d was not coupled with any struct element\n", actfnode->Id);
    /**********************end*new***********************/

    /* check dirich conds of fluid node */
    if (actfgnode->dirich==NULL)
      dserror("No dirich condition for fsi-coupling fluid node #%d",actfnode->Id);
    if(actfgnode->dirich->dirich_type!=dirich_FSI) continue;
    if (actfnode->numdf!=dim+1)
      dserror("numdf not possible for fluid node #%d at FSI interface!", actfnode->Id);
    for (j=0;j<dim;j++)
      if (actfgnode->dirich->dirich_onoff.a.iv[j]!=1)
        dserror("onoff(%d)=%d at fluid node #%d",
            j,actfgnode->dirich->dirich_onoff.a.iv[j],actfnode->Id);
    if (actfgnode->dirich->dirich_onoff.a.iv[dim]!=0)
      dserror("onoff(%d)=%d at fluid node #%d",
          dim,actfgnode->dirich->dirich_onoff.a.iv[j],actfnode->Id);

  }  /* for (i=0;i<numfnp;i++) */


  for (i=0;i<numanp;i++)
  {
    actanode  = &(alefield->dis[disnum_a].node[i]);
    actagnode = actanode->gnode;

    if (actagnode->fsicouple==NULL) continue;
    {
      /***********************start*new*********************/
      if (actagnode->coupleptr->hostele==NULL)
        dserror("ale node %d was not coupled with any struct  element\n", actanode->Id);
      /**********************end*new***********************/
    }

    if (actagnode->dirich==NULL)
      dserror("No dirich condition for fsi-coupling ale node #%d",actanode->Id);

    if(actagnode->dirich->dirich_type!=dirich_FSI) continue;

    for(j=0;j<actanode->numdf;j++)
      if(actagnode->dirich->dirich_onoff.a.iv[j]!=1)
        dserror("wrong onoff() at ale node #%d",actanode->Id);

    for (j=actanode->numdf;j<MAXDOFPERNODE;j++)
      actagnode->dirich->dirich_onoff.a.iv[j]=0;

  }  /* for (i=0;i<numanp;i++) */



#ifdef PERF
perf_end(55);
#endif

  /* print out coupling */
  if (genprob.visual==0)
    out_fsi(fluidfield);

end:


#ifdef DEBUG
  dstrc_exit();
#endif


  return;
} /* end of fsi_initcoup*/







/*-----------------------------------------------------------------------*/
/*!
  \brief determine structural fsi interface dofs

  create array sid (structural interface dofs)

  \param *structfield   FIELD         (i)      structure field
  \param *fsidyn        FSI_DYNAMIC   (i)

  \return void

  \author genk
  \date   01/03

 */
/*-----------------------------------------------------------------------*/
void fsi_struct_intdofs(
    FIELD              *structfield,
    INT                 disnum
    )

{

  INT    i,j;                    /* some counters                         */
  INT    numnp_total;            /* total number of structure dofs        */
  INT    dof;                    /* actual dof                            */
  INT    counter=0;
  INT    numaf;
  INT   *sid;                    /* structural interface dofs             */

  NODE  *actsnode /*,*actanode*/;    /* actual nodes                          */
  GNODE *actsgnode/*,*actagnode*/;  /* actual gnodes                         */

  FSI_DYNAMIC *fsidyn;


#ifdef DEBUG
  dstrc_enter("fsi_struct_intdofs");
#endif


  fsidyn = alldyn[3].fsidyn;

  numnp_total = structfield->dis[disnum].numnp;
  numaf       = genprob.numaf;

  sid = amdef("sid",&fsidyn->sid,structfield->dis[disnum].numdf,1,"IV");
  amzero(&fsidyn->sid);


  /* loop structure nodes */
  for (i=0;i<numnp_total;i++)
  {
    actsnode  = &(structfield->dis[disnum].node[i]);
    actsgnode = actsnode->gnode;

    /*****start*new****************************/
    if (actsgnode->fsicouple==NULL) continue;
    for (j=0;j<actsnode->numdf;j++)
    {
      dof = actsnode->dof[j];
      sid[dof]=1;
      counter++;
    }
    /*****end*new******************************/

    /*****start*replace****************************/
    /*actanode  = actsgnode->mfcpnode[numaf];*/
    /*if (actanode == NULL) continue;*/     /*falls kein Eintrag im mfcpnode*/
    /*actagnode = actanode->gnode;*/

    /* check for coupling nodes */
    /*if(actagnode->dirich == NULL)*/         /*falls keine Dirichlet Bedingung*/
    /*  dserror("no dirich condition for coupled ALE node #%d",actanode->Id);*/

    /*if (actagnode->dirich->dirich_type != dirich_FSI) continue;*/    /*falls keine FSI >Dirichlet Bedingung*/

    /*for (j=0;j<actanode->numdf;j++)*/          /*jetzt muss es ein FSI Node sein*/
    /*{*/
    /*  dof = actsnode->dof[j];*/
    /*  sid[dof]=1;*/
    /*  counter++;*/
    /*}*/
    /*****end*replace****************************/

  } /* end of loop over nodes */

  fsidyn->numsid=counter;


#ifdef DEBUG
  dstrc_exit();
#endif


  return;

} /* end of fsi_struct_intdofs*/



#endif
#endif  /* ifdef D_FSI */

/*! @} (documentation module close)*/

#endif

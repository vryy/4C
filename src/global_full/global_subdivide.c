/*-----------------------------------------------------------------------*/
/*!
\file
\brief Brief description.

  Very detailed description.

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

 */
/*-----------------------------------------------------------------------*/

/*!
\addtogroup Global
*//*! @{ (documentation module open)*/

#ifdef SUBDIV



#include "../headers/standardtypes.h"

#if defined(D_FLUID3) || defined(D_FLUID3_F)
#include "../fluid3/fluid3.h"
#endif

#if defined(D_FLUID2)
#include "../fluid2/fluid2.h"
#endif

#ifdef D_ALE
#include "../ale2/ale2.h"
#include "../ale3/ale3.h"
#endif

#ifdef D_SHELL8
#include "../shell8/shell8.h"
#endif

#ifdef D_WALL1
#include "../wall1/wall1.h"
#endif

#ifdef D_BRICK1
#include "../brick1/brick1.h"
#endif



/*----------------------------------------------------------------------*
  |                                                       m.gee 06/01    |
  | general problem data                                                 |
  | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;


/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate design if needed                                 |
 | defined in global_control.c                                          |
 *----------------------------------------------------------------------*/
extern struct _DESIGN *design;


/*!----------------------------------------------------------------------
\brief ranks and communicators

<pre>                                                         m.gee 8/00
This structure struct _PAR par; is defined in main_ccarat.c
and the type is in partition.h
</pre>

*----------------------------------------------------------------------*/
 extern struct _PAR   par;



/*
 * prototypes for this file
 */
void global_subdivide_hex(
    FIELD        *actfield
    );

void global_subdivide_quad(
    FIELD        *actfield
    );

void global_subdivide_tet(
    FIELD        *actfield
    );

void global_subdivide_tri(
    FIELD        *actfield
    );


/*-----------------------------------------------------------------------*/
/*!
  \brief Brief description

  Very detailed description what the function is doing why.

  \param b      DOUBLE  (o) explanation

  \return void

  \author mn
  \date   09/05

 */
/*-----------------------------------------------------------------------*/
void global_subdivide(
    FIELD              *actfield
    )

{

  DISCRET       *io_dis;


#ifdef DEBUG
  dstrc_enter("global_subdivide");
#endif

  if (actfield->subdivide > MAX_DIVIDE)
    dserror_args(__FILE__, __LINE__, "MAX_DIVIDE is too small: MAX_DIVIDE >= %2i",
        actfield->subdivide);


  if (par.myrank==0)
  {
    switch (actfield->fieldtyp)
    {
      case fluid:
        printf("\nFluid-field: ");
        break;

      case ale:
        printf("\nAle-field: ");
        break;

      case structure:
        printf("\nStructure-field: ");
        break;


      case none:
      default:
        dserror("Unknown field type!!");
        break;
    }
    printf("Subdivide = %d \n",actfield->subdivide);

  }  /* if (par.myrank==0) */


  io_dis    = &(actfield->dis[0]);


  /* determine element type in this dis */
  /* must be a uniform discretization */
  switch (io_dis->element[0].distyp)
  {
    case quad4:
      global_subdivide_quad(actfield);
      break;

    case quad8:
    case quad9:
      dserror("Subdivision only possible for linear elements: Use quad4 for io!!");
      break;

    case tri3:
      global_subdivide_tri(actfield);
      break;

    case tri6:
      dserror("Subdivision only possible for linear elements: Use tri3 for io!!");
      break;

    case hex8:
      global_subdivide_hex(actfield);
      break;

    case hex20:
    case hex27:
      dserror("Subdivision only possible for linear elements: Use hex8 for io!!");
      break;

    case tet4:
      global_subdivide_tet(actfield);
      break;

    case tet10:
      dserror("Subdivision only possible for linear elements: Use tet4 for io!!");
      break;

    case line2:
    case line3:
      /* TODO */
      dserror("Subdivision not yet possible for line elements!!");
      break;

    case dis_none:
    default:
      dserror("Unknown distype for this element!!");
      break;
  }


#ifdef DEBUG
  dstrc_exit();
#endif


  return;

}  /* END of global_subdivide */





/*-----------------------------------------------------------------------*/
/*!
  \brief Brief description

  Very detailed description what the function is doing why.

  \param b      DOUBLE  (o) explanation

  \return void

  \author mn
  \date   09/05

 */
/*-----------------------------------------------------------------------*/
void global_subdivide_hex(
    FIELD        *actfield
    )

{

  DISCRET       *io_dis;
  DISCRET       *cal_dis;

  INT            subdivide;
  INT            nodeid;
  INT            eleid;

  INT            numnd;
  INT            numele;
  INT            i,j;
  INT            node_counter;
  INT            ele_counter;

  GNODE         *actgnode = NULL;
  GLINE         *actgline = NULL;
  GSURF         *actgsurf = NULL;
  GVOL          *actgvol  = NULL;

  ELEMENT       *actele, *oldele;

  INT            id;

  INT            tmp[24];
  INT            counter;
  INT            counter2;
  INT            k,l,m,g;

  INT            node_matrix[MAX_DIVIDE+1][MAX_DIVIDE+1][MAX_DIVIDE+1];

  INT            pro_total=0, pro_fac=0, p=0, p_counter=0;

#ifdef DEBUG
  dstrc_enter("global_subdivide_hex");
#endif


  io_dis    = &(actfield->dis[0]);
  cal_dis   = &(actfield->dis[1]);

  subdivide = actfield->subdivide;
  nodeid    = genprob.maxnode;
  eleid     = genprob.nele;


  /* calculate number of nodes on second dis */
  numnd = io_dis->ngnode +
          io_dis->ngline * (subdivide - 1) +
          io_dis->ngsurf * (subdivide - 1) * (subdivide - 1) +
          io_dis->ngvol  * (subdivide - 1) * (subdivide - 1) * (subdivide - 1);



  /* allocate the nodes to the dis */
  cal_dis->node  = (NODE*)CCACALLOC(numnd,sizeof(NODE));
  cal_dis->numnp = numnd;
  node_counter = 0;



  /* calculate number of elements in second dis */
  numele = io_dis->ngvol * subdivide * subdivide * subdivide;


  /* allocate elements in second dis */
  cal_dis->element = (ELEMENT*)CCACALLOC(numele,sizeof(ELEMENT));
  cal_dis->numele  = numele;
  ele_counter = 0;






  /* loop gnodes and create nodes on second dis */
  for (i=0; i<io_dis->ngnode; i++)
  {
    actgnode = &(io_dis->gnode[i]);

    /* set Id of the new node */
    cal_dis->node[node_counter].Id = nodeid;

    /* set coords of the new node */
    for (j=0; j<3; j++)
      cal_dis->node[node_counter].x[j] = actgnode->node->x[j];

    /* store slave_node Id in master gnode */
    actgnode->slave_node = nodeid;

    /* create pointer from slave to master node */
    cal_dis->node[node_counter].slave_node  = NULL;
    cal_dis->node[node_counter].master_node = actgnode->node;

    /* create pointer from master to slave node */
    actgnode->node->slave_node  = &(cal_dis->node[node_counter]);
    actgnode->node->master_node = NULL;


    node_counter++;
    nodeid++;
  }



  /* loop glines and create nodes */
  for (i=0; i<io_dis->ngline; i++)
  {

    DOUBLE x1[3],dx[3];

    actgline = &(io_dis->gline[i]);

    /* coordinates of the beginning of this gline */
    x1[0] = actgline->gnode[0]->node->x[0];
    x1[1] = actgline->gnode[0]->node->x[1];
    x1[2] = actgline->gnode[0]->node->x[2];

    /* lenght of this gline in the three directions divided by subdivide */
    dx[0] = (actgline->gnode[1]->node->x[0]-actgline->gnode[0]->node->x[0])/subdivide;
    dx[1] = (actgline->gnode[1]->node->x[1]-actgline->gnode[0]->node->x[1])/subdivide;
    dx[2] = (actgline->gnode[1]->node->x[2]-actgline->gnode[0]->node->x[2])/subdivide;

    for (j=0; j<subdivide-1; j++)
    {
      cal_dis->node[node_counter].Id = nodeid;

      /* store slave_node Id in gline */
      actgline->slave_node[j]   = nodeid;

      /* set coords of the new node */
      for (k=0; k<3; k++)
        cal_dis->node[node_counter].x[k] =  x1[k] + (j+1)*dx[k];


      node_counter++;
      nodeid++;
    }

  }  /* for (i=0; i<io_dis->ngline; i++) */











  /* loop gsurfs and create nodes */
  for (i=0; i<io_dis->ngsurf; i++)
  {
    GNODE  *gnode0=NULL, *gnode1=NULL, *gnode2=NULL, *gnode3=NULL;

    DOUBLE x0[3],x1[3],x2[3],x3[3];

    actgsurf = &(io_dis->gsurf[i]);


    /* specify node 0 */
    gnode0 = actgsurf->gnode[0];
    actgsurf->node_ind[0] = 0;




    /* find the line 0 */
    for (j=0; j<actgsurf->ngline; j++)
    {
      /* this line is in same direction */
      if (gnode0 == actgsurf->gline[j]->gnode[0])
      {
        actgline = actgsurf->gline[j];
        actgsurf->line_ind[0][0] = j;
        actgsurf->line_ind[1][0] = 0;

        gnode1 = actgline->gnode[1];

        /* find the node 1 */
        for (k=0; k<actgsurf->ngnode; k++)
        {
          if (gnode1 == actgsurf->gnode[k])
          {
            actgsurf->node_ind[1] = k;
            break;
          }
        }

        break;
      }

      /* this line is in opposite direction */
      if (gnode0 == actgsurf->gline[j]->gnode[1])
      {
        actgline = actgsurf->gline[j];
        actgsurf->line_ind[0][0] = j;
        actgsurf->line_ind[1][0] = 1;

        gnode1 = actgline->gnode[0];

        /* find the node 1 */
        for (k=0; k<actgsurf->ngnode; k++)
        {
          if (gnode1 == actgsurf->gnode[k])
          {
            actgsurf->node_ind[1] = k;
            break;
          }
        }

        break;
      }
    }  /* for (j=1; j<gsurf->ngline; j++) */



    /* find the line 1 */
    for (j=0; j<actgsurf->ngline; j++)
    {

      if (actgsurf->gline[j] == actgline)
        continue;

      /* this line is in same direction */
      if (gnode1 == actgsurf->gline[j]->gnode[0])
      {
        actgline = actgsurf->gline[j];
        actgsurf->line_ind[0][1] = j;
        actgsurf->line_ind[1][1] = 0;

        gnode2 = actgline->gnode[1];

        /* find the node 2 */
        for (k=0; k<actgsurf->ngnode; k++)
        {
          if (gnode2 == actgsurf->gnode[k])
          {
            actgsurf->node_ind[2] = k;
            break;
          }
        }

        break;
      }

      /* this line is in opposite direction */
      if (gnode1 == actgsurf->gline[j]->gnode[1])
      {
        actgline = actgsurf->gline[j];
        actgsurf->line_ind[0][1] = j;
        actgsurf->line_ind[1][1] = 1;

        gnode2 = actgline->gnode[0];

        /* find the node 2 */
        for (k=0; k<actgsurf->ngnode; k++)
        {
          if (gnode2 == actgsurf->gnode[k])
          {
            actgsurf->node_ind[2] = k;
            break;
          }
        }

        break;
      }
    }  /* for (j=1; j<gsurf->ngline; j++) */


    /* find the line 2 */
    for (j=0; j<actgsurf->ngline; j++)
    {

      if (actgsurf->gline[j] == actgline)
        continue;

      /* this line is in same direction */
      if (gnode2 == actgsurf->gline[j]->gnode[1])
      {
        actgline = actgsurf->gline[j];
        actgsurf->line_ind[0][2] = j;
        actgsurf->line_ind[1][2] = 0;

        gnode3 = actgline->gnode[0];

        /* find the node 3 */
        for (k=0; k<actgsurf->ngnode; k++)
        {
          if (gnode3 == actgsurf->gnode[k])
          {
            actgsurf->node_ind[3] = k;
            break;
          }
        }

        break;
      }

      /* this line is in opposite direction */
      if (gnode2 == actgsurf->gline[j]->gnode[0])
      {
        actgline = actgsurf->gline[j];
        actgsurf->line_ind[0][2] = j;
        actgsurf->line_ind[1][2] = 1;

        gnode3 = actgline->gnode[1];

        /* find the node 3 */
        for (k=0; k<actgsurf->ngnode; k++)
        {
          if (gnode3 == actgsurf->gnode[k])
          {
            actgsurf->node_ind[3] = k;
            break;
          }
        }

        break;
      }
    }  /* for (j=1; j<gsurf->ngline; j++) */





    /* find the line 3 */
    for (j=0; j<actgsurf->ngline; j++)
    {

      if (actgsurf->gline[j] == actgline)
        continue;

      /* this line is in same direction */
      if (gnode3 == actgsurf->gline[j]->gnode[1])
      {
        actgline = actgsurf->gline[j];
        actgsurf->line_ind[0][3] = j;
        actgsurf->line_ind[1][3] = 0;

        actgnode = actgline->gnode[0];

        if (actgnode != gnode0)
          dserror("NODE MIXUP!!\n");

        break;
      }

      /* this line is in opposite direction */
      if (gnode3 == actgsurf->gline[j]->gnode[0])
      {
        actgline = actgsurf->gline[j];
        actgsurf->line_ind[0][3] = j;
        actgsurf->line_ind[1][3] = 1;

        actgnode = actgline->gnode[1];

        if (actgnode != gnode0)
          dserror("NODE MIXUP!!\n");

        break;
      }
    }  /* for (j=1; j<gsurf->ngline; j++) */







    /*
       printf("GSURF %2i:  line_ind: %1i(%1i) %1i(%1i) %1i(%1i) %1i(%1i)\n",
       actgsurf->Id,
       actgsurf->line_ind[0][0],actgsurf->line_ind[1][0],
       actgsurf->line_ind[0][1],actgsurf->line_ind[1][1],
       actgsurf->line_ind[0][2],actgsurf->line_ind[1][2],
       actgsurf->line_ind[0][3],actgsurf->line_ind[1][3]);

       printf("GSURF %2i:  actgsurf->node_ind: %1i %1i %1i %1i\n\n",
       actgsurf->Id,
       actgsurf->node_ind[0],actgsurf->node_ind[1],
       actgsurf->node_ind[2],actgsurf->node_ind[3]);
       */


    /* set the nodal coordinates */
    x0[0] = (actgsurf->gnode[actgsurf->node_ind[0]]->node->x[0])/(subdivide*subdivide);
    x0[1] = (actgsurf->gnode[actgsurf->node_ind[0]]->node->x[1])/(subdivide*subdivide);
    x0[2] = (actgsurf->gnode[actgsurf->node_ind[0]]->node->x[2])/(subdivide*subdivide);

    x1[0] = (actgsurf->gnode[actgsurf->node_ind[1]]->node->x[0])/(subdivide*subdivide);
    x1[1] = (actgsurf->gnode[actgsurf->node_ind[1]]->node->x[1])/(subdivide*subdivide);
    x1[2] = (actgsurf->gnode[actgsurf->node_ind[1]]->node->x[2])/(subdivide*subdivide);

    x2[0] = (actgsurf->gnode[actgsurf->node_ind[2]]->node->x[0])/(subdivide*subdivide);
    x2[1] = (actgsurf->gnode[actgsurf->node_ind[2]]->node->x[1])/(subdivide*subdivide);
    x2[2] = (actgsurf->gnode[actgsurf->node_ind[2]]->node->x[2])/(subdivide*subdivide);

    x3[0] = (actgsurf->gnode[actgsurf->node_ind[3]]->node->x[0])/(subdivide*subdivide);
    x3[1] = (actgsurf->gnode[actgsurf->node_ind[3]]->node->x[1])/(subdivide*subdivide);
    x3[2] = (actgsurf->gnode[actgsurf->node_ind[3]]->node->x[2])/(subdivide*subdivide);



    /* generate the nodes */
    for (j=0; j<subdivide-1; j++)
    {
      for (k=0; k<subdivide-1; k++)
      {
        cal_dis->node[node_counter].Id  = nodeid;

        /* store slave_node Id in gsurf */
        actgsurf->slave_node[j][k] = nodeid;

        /* set coords of the new node */
        for (l=0; l<3; l++)
          cal_dis->node[node_counter].x[l] =
            (x1[l]*(j+1) + x0[l]*(subdivide-j-1)) * (subdivide-k-1) +
            (x2[l]*(j+1) + x3[l]*(subdivide-j-1)) * (k+1);

        node_counter++;
        nodeid++;
      }
    }

  }  /* for (i=0; i<io_dis->ngsurf; i++) */




  if (par.myrank==0)
  {
    printf("  creating subelements on %3i GVOLS\n",io_dis->ngvol);

    pro_fac   = io_dis->ngvol/50 + 1;

    pro_total = io_dis->ngvol/pro_fac;

    printf("  0");
    for (p=0; p<pro_total; p++)
      printf(" ");
    printf("%3i\n",io_dis->ngvol);

    printf("  |");
    for (p=0; p<pro_total; p++)
      printf("-");
    printf("|\n");


    printf("  |");
    p_counter = 0;
  }

  /* loop gvols and create nodes */
  for (g=0; g<io_dis->ngvol; g++)
  {

    GNODE  *gnode0=NULL, *gnode1=NULL, *gnode2=NULL, *gnode3=NULL;
    GNODE  *gnode4=NULL, *gnode5=NULL, *gnode6=NULL, *gnode7=NULL;

    DOUBLE x0[3],x1[3],x2[3],x3[3];
    DOUBLE x4[3],x5[3],x6[3],x7[3];

    actgvol = &(io_dis->gvol[g]);


    /* LINE 0: 0-1 */
    /*=============*/

    actgline = actgvol->gline[0];
    actgvol->line_ind[0][0] = 0;
    actgvol->line_ind[1][0] = 0;



    /* coordinates of node 0 = beginning of line 1 */
    gnode0 = actgline->gnode[0];

    x0[0] = gnode0->node->x[0]/(subdivide*subdivide*subdivide);
    x0[1] = gnode0->node->x[1]/(subdivide*subdivide*subdivide);
    x0[2] = gnode0->node->x[2]/(subdivide*subdivide*subdivide);


    /* coordinates of node 1 = end of line 1 */
    gnode1 = actgline->gnode[1];

    x1[0] = gnode1->node->x[0]/(subdivide*subdivide*subdivide);
    x1[1] = gnode1->node->x[1]/(subdivide*subdivide*subdivide);
    x1[2] = gnode1->node->x[2]/(subdivide*subdivide*subdivide);


    /* LINE 1: 1-2 */
    /*=============*/

    for (j=1; j<actgvol->ngline; j++)
    {
      /* this line is in same direction */
      if (gnode1 == actgvol->gline[j]->gnode[0])
      {
        actgline = actgvol->gline[j];
        actgvol->line_ind[0][1] = j;
        actgvol->line_ind[1][1] = 0;

        gnode2 = actgline->gnode[1];

        /* coordinates of the third node */
        x2[0] = gnode2->node->x[0]/(subdivide*subdivide*subdivide);
        x2[1] = gnode2->node->x[1]/(subdivide*subdivide*subdivide);
        x2[2] = gnode2->node->x[2]/(subdivide*subdivide*subdivide);

        break;
      }

      /* this line is in opposite direction */
      if (gnode1 == actgvol->gline[j]->gnode[1])
      {
        actgline = actgvol->gline[j];
        actgvol->line_ind[0][1] = j;
        actgvol->line_ind[1][1] = 1;

        gnode2 = actgline->gnode[0];

        /* coordinates of the third node */
        x2[0] = gnode2->node->x[0]/(subdivide*subdivide*subdivide);
        x2[1] = gnode2->node->x[1]/(subdivide*subdivide*subdivide);
        x2[2] = gnode2->node->x[2]/(subdivide*subdivide*subdivide);

        break;
      }

    }  /* for (j=1; j<actgvol->ngline; j++) */


    actgsurf = NULL;



    /* GSURF 0 */
    /*=========*/

    for (j=0; j<actgvol->ngsurf; j++)
    {
      INT found = 0;


      for (k=0; k<actgvol->gsurf[j]->ngline; k++)
        if (actgvol->gsurf[j]->gline[k] == actgvol->gline[actgvol->line_ind[0][0]])
          found++;

      for (k=0; k<actgvol->gsurf[j]->ngline; k++)
        if (actgvol->gsurf[j]->gline[k] == actgvol->gline[actgvol->line_ind[0][1]])
          found++;

      if (found == 2)
      {
        actgsurf = actgvol->gsurf[j];
        actgvol->surf_ind[0] = j;
        break;
      }

    }  /* for (j=0; j<actgvol->ngsurf; j++) */

    if (actgsurf == NULL)
    {
      dserror("Surface 0 not found!!\n");
    }



    /* LINE 6: 2-6 */
    /*=============*/

    for (j=0; j<actgvol->ngline; j++)
    {

      INT found = 0;

      /* skip all lines of the first surface */
      for (k=0; k<actgsurf->ngline; k++)
        if (actgvol->gline[j] == actgsurf->gline[k])
        {
          found = 1;
          break;
        }

      if (found == 1)
        continue;


      /* this line is in same direction */
      if (gnode2 == actgvol->gline[j]->gnode[0])
      {
        actgline = actgvol->gline[j];
        actgvol->line_ind[0][6] = j;
        actgvol->line_ind[1][6] = 0;

        gnode6 = actgline->gnode[1];

        /* coordinates of the seventh node */
        x6[0] = gnode6->node->x[0]/(subdivide*subdivide*subdivide);
        x6[1] = gnode6->node->x[1]/(subdivide*subdivide*subdivide);
        x6[2] = gnode6->node->x[2]/(subdivide*subdivide*subdivide);

        break;
      }

      /* this line is in opposite direction */
      if (gnode2 == actgvol->gline[j]->gnode[1])
      {
        actgline = actgvol->gline[j];
        actgvol->line_ind[0][6] = j;
        actgvol->line_ind[1][6] = 1;

        gnode6 = actgline->gnode[0];

        /* coordinates of the seventh node */
        x6[0] = gnode6->node->x[0]/(subdivide*subdivide*subdivide);
        x6[1] = gnode6->node->x[1]/(subdivide*subdivide*subdivide);
        x6[2] = gnode6->node->x[2]/(subdivide*subdivide*subdivide);

        break;
      }
    }  /* for (j=1; j<gsurf->ngline; j++) */




    /* LINE 2: 3-2 */
    /*=============*/

    for (j=0; j<actgvol->ngline; j++)
    {
      INT   found = 0;

      if (actgvol->gline[j] == actgvol->gline[actgvol->line_ind[0][0]])
        continue;

      if (actgvol->gline[j] == actgvol->gline[actgvol->line_ind[0][1]])
        continue;

      for (k=0; k<actgsurf->ngline; k++)
        if (actgvol->gline[j] == actgsurf->gline[k])
          found = 1;

      if (found == 0)
        continue;


      /* this line is in same direction */
      if (gnode2 == actgvol->gline[j]->gnode[1])
      {
        actgline = actgvol->gline[j];
        actgvol->line_ind[0][2] = j;
        actgvol->line_ind[1][2] = 0;

        gnode3 = actgline->gnode[0];

        /* coordinates of the fourth node */
        x3[0] = gnode3->node->x[0]/(subdivide*subdivide*subdivide);
        x3[1] = gnode3->node->x[1]/(subdivide*subdivide*subdivide);
        x3[2] = gnode3->node->x[2]/(subdivide*subdivide*subdivide);

        break;
      }

      /* this line is in opposite direction */
      if (gnode2 == actgvol->gline[j]->gnode[0])
      {
        actgline = actgvol->gline[j];
        actgvol->line_ind[0][2] = j;
        actgvol->line_ind[1][2] = 1;

        gnode3 = actgline->gnode[1];

        /* coordinates of the fourth node */
        x3[0] = gnode3->node->x[0]/(subdivide*subdivide*subdivide);
        x3[1] = gnode3->node->x[1]/(subdivide*subdivide*subdivide);
        x3[2] = gnode3->node->x[2]/(subdivide*subdivide*subdivide);

        break;
      }
    }  /* for (j=1; j<gsurf->ngline; j++) */



    /* LINE 3: 0-3 */
    /*=============*/

    for (j=0; j<actgvol->ngline; j++)
    {
      INT   found = 0;

      if (actgvol->gline[j] == actgvol->gline[actgvol->line_ind[0][0]])
        continue;

      if (actgvol->gline[j] == actgvol->gline[actgvol->line_ind[0][1]])
        continue;

      if (actgvol->gline[j] == actgvol->gline[actgvol->line_ind[0][2]])
        continue;

      for (k=0; k<actgsurf->ngline; k++)
        if (actgvol->gline[j] == actgsurf->gline[k])
          found = 1;

      if (found == 0)
        continue;


      /* this line is in same direction */
      if (gnode3 == actgvol->gline[j]->gnode[1])
      {
        actgline = actgvol->gline[j];

        actgvol->line_ind[0][3] = j;
        actgvol->line_ind[1][3] = 0;
        break;
      }

      /* this line is in opposite direction */
      if (gnode3 == actgvol->gline[j]->gnode[0])
      {
        actgline = actgvol->gline[j];

        actgvol->line_ind[0][3] = j;
        actgvol->line_ind[1][3] = 1;
      }
    }  /* for (j=1; j<gsurf->ngline; j++) */



    /* LINE 7: 3-7 */
    /*=============*/

    for (j=0; j<actgvol->ngline; j++)
    {

      INT found = 0;

      /* skip all lines of the first surface */
      for (k=0; k<actgsurf->ngline; k++)
        if (actgvol->gline[j] == actgsurf->gline[k])
        {
          found = 1;
          break;
        }

      if (found == 1)
        continue;


      /* this line is in same direction */
      if (gnode3 == actgvol->gline[j]->gnode[0])
      {
        actgline = actgvol->gline[j];
        actgvol->line_ind[0][7] = j;
        actgvol->line_ind[1][7] = 0;

        gnode7 = actgline->gnode[1];

        /* coordinates of the eigth node */
        x7[0] = gnode7->node->x[0]/(subdivide*subdivide*subdivide);
        x7[1] = gnode7->node->x[1]/(subdivide*subdivide*subdivide);
        x7[2] = gnode7->node->x[2]/(subdivide*subdivide*subdivide);

        break;
      }

      /* this line is in opposite direction */
      if (gnode3 == actgvol->gline[j]->gnode[1])
      {
        actgline = actgvol->gline[j];
        actgvol->line_ind[0][7] = j;
        actgvol->line_ind[1][7] = 1;

        gnode7 = actgline->gnode[0];

        /* coordinates of the eigth node */
        x7[0] = gnode7->node->x[0]/(subdivide*subdivide*subdivide);
        x7[1] = gnode7->node->x[1]/(subdivide*subdivide*subdivide);
        x7[2] = gnode7->node->x[2]/(subdivide*subdivide*subdivide);

        break;
      }
    }  /* for (j=1; j<gsurf->ngline; j++) */




    /* LINE 5: 1-5 */
    /*=============*/

    for (j=0; j<actgvol->ngline; j++)
    {

      INT found = 0;

      /* skip all lines of the first surface */
      for (k=0; k<actgsurf->ngline; k++)
        if (actgvol->gline[j] == actgsurf->gline[k])
        {
          found = 1;
          break;
        }

      if (found == 1)
        continue;


      /* this line is in same direction */
      if (gnode1 == actgvol->gline[j]->gnode[0])
      {
        actgline = actgvol->gline[j];
        actgvol->line_ind[0][5] = j;
        actgvol->line_ind[1][5] = 0;

        gnode5 = actgline->gnode[1];

        /* coordinates of the sixth node */
        x5[0] = gnode5->node->x[0]/(subdivide*subdivide*subdivide);
        x5[1] = gnode5->node->x[1]/(subdivide*subdivide*subdivide);
        x5[2] = gnode5->node->x[2]/(subdivide*subdivide*subdivide);

        break;
      }

      /* this line is in opposite direction */
      if (gnode1 == actgvol->gline[j]->gnode[1])
      {
        actgline = actgvol->gline[j];
        actgvol->line_ind[0][5] = j;
        actgvol->line_ind[1][5] = 1;

        gnode5 = actgline->gnode[0];


        /* coordinates of the sixth node */
        x5[0] = gnode5->node->x[0]/(subdivide*subdivide*subdivide);
        x5[1] = gnode5->node->x[1]/(subdivide*subdivide*subdivide);
        x5[2] = gnode5->node->x[2]/(subdivide*subdivide*subdivide);

        break;
      }
    }  /* for (j=1; j<gsurf->ngline; j++) */




    /* LINE 4: 0-4 */
    /*=============*/

    for (j=0; j<actgvol->ngline; j++)
    {

      INT found = 0;

      /* skip all lines of the first surface */
      for (k=0; k<actgsurf->ngline; k++)
        if (actgvol->gline[j] == actgsurf->gline[k])
        {
          found = 1;
          break;
        }

      if (found == 1)
        continue;


      /* this line is in same direction */
      if (gnode0 == actgvol->gline[j]->gnode[0])
      {
        actgline = actgvol->gline[j];
        actgvol->line_ind[0][4] = j;
        actgvol->line_ind[1][4] = 0;

        gnode4 = actgline->gnode[1];

        /* coordinates of the fifth node */
        x4[0] = gnode4->node->x[0]/(subdivide*subdivide*subdivide);
        x4[1] = gnode4->node->x[1]/(subdivide*subdivide*subdivide);
        x4[2] = gnode4->node->x[2]/(subdivide*subdivide*subdivide);

        break;
      }

      /* this line is in opposite direction */
      if (gnode0 == actgvol->gline[j]->gnode[1])
      {
        actgline = actgvol->gline[j];
        actgvol->line_ind[0][4] = j;
        actgvol->line_ind[1][4] = 1;

        gnode4 = actgline->gnode[0];

        /* coordinates of the fifth node */
        x4[0] = gnode4->node->x[0]/(subdivide*subdivide*subdivide);
        x4[1] = gnode4->node->x[1]/(subdivide*subdivide*subdivide);
        x4[2] = gnode4->node->x[2]/(subdivide*subdivide*subdivide);

        break;
      }
    }  /* for (j=1; j<gsurf->ngline; j++) */






    /* LINE 8: 4-5 */
    /*=============*/

    for (j=0; j<actgvol->ngline; j++)
    {
      /* this line is in same direction */
      if (gnode4 == actgvol->gline[j]->gnode[0] &&
          gnode5 == actgvol->gline[j]->gnode[1])
      {
        actgline = actgvol->gline[j];
        actgvol->line_ind[0][8] = j;
        actgvol->line_ind[1][8] = 0;
        break;
      }

      /* this line is in opposite direction */
      if (gnode4 == actgvol->gline[j]->gnode[1] &&
          gnode5 == actgvol->gline[j]->gnode[0])
      {
        actgline = actgvol->gline[j];
        actgvol->line_ind[0][8] = j;
        actgvol->line_ind[1][8] = 1;
        break;
      }
    }  /* for (j=1; j<gsurf->ngline; j++) */




    /* LINE 9: 5-6 */
    /*=============*/

    for (j=0; j<actgvol->ngline; j++)
    {
      /* this line is in same direction */
      if (gnode5 == actgvol->gline[j]->gnode[0] &&
          gnode6 == actgvol->gline[j]->gnode[1])
      {
        actgline = actgvol->gline[j];
        actgvol->line_ind[0][9] = j;
        actgvol->line_ind[1][9] = 0;
        break;
      }

      /* this line is in opposite direction */
      if (gnode5 == actgvol->gline[j]->gnode[1] &&
          gnode6 == actgvol->gline[j]->gnode[0])
      {
        actgline = actgvol->gline[j];
        actgvol->line_ind[0][9] = j;
        actgvol->line_ind[1][9] = 1;
        break;
      }
    }  /* for (j=1; j<gsurf->ngline; j++) */





    /* LINE 10: 7-6 */
    /*==============*/

    for (j=0; j<actgvol->ngline; j++)
    {
      /* this line is in same direction */
      if (gnode6 == actgvol->gline[j]->gnode[1] &&
          gnode7 == actgvol->gline[j]->gnode[0])
      {
        actgline = actgvol->gline[j];
        actgvol->line_ind[0][10] = j;
        actgvol->line_ind[1][10] = 0;
        break;
      }

      /* this line is in opposite direction */
      if (gnode6 == actgvol->gline[j]->gnode[0] &&
          gnode7 == actgvol->gline[j]->gnode[1])
      {
        actgline = actgvol->gline[j];
        actgvol->line_ind[0][10] = j;
        actgvol->line_ind[1][10] = 1;
        break;
      }
    }  /* for (j=1; j<gsurf->ngline; j++) */




    /* LINE 11: 4-7 */
    /*==============*/

    for (j=0; j<actgvol->ngline; j++)
    {
      /* this line is in same direction */
      if (gnode7 == actgvol->gline[j]->gnode[1] &&
          gnode4 == actgvol->gline[j]->gnode[0])
      {
        actgline = actgvol->gline[j];
        actgvol->line_ind[0][11] = j;
        actgvol->line_ind[1][11] = 0;
        break;
      }

      /* this line is in opposite direction */
      if (gnode7 == actgvol->gline[j]->gnode[0] &&
          gnode4 == actgvol->gline[j]->gnode[1])
      {
        actgline = actgvol->gline[j];
        actgvol->line_ind[0][11] = j;
        actgvol->line_ind[1][11] = 1;
        break;
      }
    }  /* for (j=0; j<gsurf->ngline; j++) */



    /* GSURF 1 */
    /*=========*/

    for (j=0; j<actgvol->ngsurf; j++)
    {
      INT found = 0;


      for (k=0; k<actgvol->gsurf[j]->ngline; k++)
        if (actgvol->gsurf[j]->gline[k] == actgvol->gline[actgvol->line_ind[0][0]])
          found++;

      for (k=0; k<actgvol->gsurf[j]->ngline; k++)
        if (actgvol->gsurf[j]->gline[k] == actgvol->gline[actgvol->line_ind[0][5]])
          found++;

      if (found == 2)
      {
        actgsurf = actgvol->gsurf[j];
        actgvol->surf_ind[1] = j;
        break;
      }

    }  /* for (j=0; j<actgvol->ngsurf; j++) */

    if (actgsurf == NULL)
    {
      dserror("Surface 1 not found!!\n");
    }


    /* GSURF 2 */
    /*=========*/

    for (j=0; j<actgvol->ngsurf; j++)
    {
      INT found = 0;


      for (k=0; k<actgvol->gsurf[j]->ngline; k++)
        if (actgvol->gsurf[j]->gline[k] == actgvol->gline[actgvol->line_ind[0][5]])
          found++;

      for (k=0; k<actgvol->gsurf[j]->ngline; k++)
        if (actgvol->gsurf[j]->gline[k] == actgvol->gline[actgvol->line_ind[0][1]])
          found++;

      if (found == 2)
      {
        actgsurf = actgvol->gsurf[j];
        actgvol->surf_ind[2] = j;
        break;
      }

    }  /* for (j=0; j<actgvol->ngsurf; j++) */

    if (actgsurf == NULL)
    {
      dserror("Surface 2 not found!!\n");
    }


    /* GSURF 3 */
    /*=========*/

    for (j=0; j<actgvol->ngsurf; j++)
    {
      INT found = 0;


      for (k=0; k<actgvol->gsurf[j]->ngline; k++)
        if (actgvol->gsurf[j]->gline[k] == actgvol->gline[actgvol->line_ind[0][2]])
          found++;

      for (k=0; k<actgvol->gsurf[j]->ngline; k++)
        if (actgvol->gsurf[j]->gline[k] == actgvol->gline[actgvol->line_ind[0][6]])
          found++;

      if (found == 2)
      {
        actgsurf = actgvol->gsurf[j];
        actgvol->surf_ind[3] = j;
        break;
      }

    }  /* for (j=0; j<actgvol->ngsurf; j++) */

    if (actgsurf == NULL)
    {
      dserror("Surface 3 not found!!\n");
    }


    /* GSURF 4 */
    /*=========*/

    for (j=0; j<actgvol->ngsurf; j++)
    {
      INT found = 0;


      for (k=0; k<actgvol->gsurf[j]->ngline; k++)
        if (actgvol->gsurf[j]->gline[k] == actgvol->gline[actgvol->line_ind[0][3]])
          found++;

      for (k=0; k<actgvol->gsurf[j]->ngline; k++)
        if (actgvol->gsurf[j]->gline[k] == actgvol->gline[actgvol->line_ind[0][4]])
          found++;

      if (found == 2)
      {
        actgsurf = actgvol->gsurf[j];
        actgvol->surf_ind[4] = j;
        break;
      }

    }  /* for (j=0; j<actgvol->ngsurf; j++) */

    if (actgsurf == NULL)
    {
      dserror("Surface 4 not found!!\n");
    }


    /* GSURF 5 */
    /*=========*/

    for (j=0; j<actgvol->ngsurf; j++)
    {
      INT found = 0;


      for (k=0; k<actgvol->gsurf[j]->ngline; k++)
        if (actgvol->gsurf[j]->gline[k] == actgvol->gline[actgvol->line_ind[0][8]])
          found++;

      for (k=0; k<actgvol->gsurf[j]->ngline; k++)
        if (actgvol->gsurf[j]->gline[k] == actgvol->gline[actgvol->line_ind[0][9]])
          found++;

      if (found == 2)
      {
        actgsurf = actgvol->gsurf[j];
        actgvol->surf_ind[5] = j;
        break;
      }

    }  /* for (j=0; j<actgvol->ngsurf; j++) */

    if (actgsurf == NULL)
    {
      dserror("Surface 5 not found!!\n");
    }




    /* check gnode - gline */
    for (j=0; j<12; j++)
    {
      INT links = 0, rechts = 1;
      GNODE *gnode_l=NULL, *gnode_r=NULL;

      INT line_node_ind[12][2];

      line_node_ind[0][0] = 0;
      line_node_ind[0][1] = 1;
      line_node_ind[1][0] = 1;
      line_node_ind[1][1] = 2;
      line_node_ind[2][0] = 3;
      line_node_ind[2][1] = 2;
      line_node_ind[3][0] = 0;
      line_node_ind[3][1] = 3;
      line_node_ind[4][0] = 0;
      line_node_ind[4][1] = 4;
      line_node_ind[5][0] = 1;
      line_node_ind[5][1] = 5;
      line_node_ind[6][0] = 2;
      line_node_ind[6][1] = 6;
      line_node_ind[7][0] = 3;
      line_node_ind[7][1] = 7;
      line_node_ind[8][0] = 4;
      line_node_ind[8][1] = 5;
      line_node_ind[9][0] = 5;
      line_node_ind[9][1] = 6;
      line_node_ind[10][0] = 7;
      line_node_ind[10][1] = 6;
      line_node_ind[11][0] = 4;
      line_node_ind[11][1] = 7;


      if(actgvol->line_ind[1][j] == 0)
      {
        links = 0;
        rechts = 1;
      }
      else
      {
        links = 1;
        rechts = 0;
      }

      switch (line_node_ind[j][0])
      {
        case 0: gnode_l = gnode0; break;
        case 1: gnode_l = gnode1; break;
        case 2: gnode_l = gnode2; break;
        case 3: gnode_l = gnode3; break;
        case 4: gnode_l = gnode4; break;
        case 5: gnode_l = gnode5; break;
        case 6: gnode_l = gnode6; break;
        case 7: gnode_l = gnode7; break;
      }

      switch (line_node_ind[j][1])
      {
        case 0: gnode_r = gnode0; break;
        case 1: gnode_r = gnode1; break;
        case 2: gnode_r = gnode2; break;
        case 3: gnode_r = gnode3; break;
        case 4: gnode_r = gnode4; break;
        case 5: gnode_r = gnode5; break;
        case 6: gnode_r = gnode6; break;
        case 7: gnode_r = gnode7; break;
      }

      if (actgvol->gline[actgvol->line_ind[0][j]]->gnode[links] != gnode_l)
        dserror_args(__FILE__, __LINE__, "Gline %1i, links != gnode %1i\n",
            j,line_node_ind[j][0]);

      if (actgvol->gline[actgvol->line_ind[0][j]]->gnode[rechts] != gnode_r)
        dserror_args(__FILE__, __LINE__, "Gline %1i, rechts != gnode %1i\n",
            j,line_node_ind[j][1]);
    }





    /* check gnodes */
    for (j=0; j<8; j++)
    {
      for (k=0;k <8; k++)
      {
        GNODE *gnode_l=NULL, *gnode_r=NULL;

        switch (j)
        {
          case 0: gnode_l = gnode0; break;
          case 1: gnode_l = gnode1; break;
          case 2: gnode_l = gnode2; break;
          case 3: gnode_l = gnode3; break;
          case 4: gnode_l = gnode4; break;
          case 5: gnode_l = gnode5; break;
          case 6: gnode_l = gnode6; break;
          case 7: gnode_l = gnode7; break;
        }
        switch (k)
        {
          case 0: gnode_r = gnode0; break;
          case 1: gnode_r = gnode1; break;
          case 2: gnode_r = gnode2; break;
          case 3: gnode_r = gnode3; break;
          case 4: gnode_r = gnode4; break;
          case 5: gnode_r = gnode5; break;
          case 6: gnode_r = gnode6; break;
          case 7: gnode_r = gnode7; break;
        }

        if (gnode_r == gnode_l && j!=k)
        {
          dserror_args(__FILE__, __LINE__, "GNODE %1i == GNODE %1i !!\n",j,k);
        }
      }
    }


    /*
       printf("GVOL %2i:  line_ind: ", actgvol->Id);
       for (j=0; j<actgvol->ngline; j++)
       printf("%1i(%1i) ",actgvol->line_ind[0][j],actgvol->line_ind[1][j]);
       printf("\n\n");
       */


    /* generate the nodes */
    for (j=0; j<subdivide-1; j++)
    {
      for (k=0; k<subdivide-1; k++)
      {
        for (l=0; l<subdivide-1; l++)
        {
          cal_dis->node[node_counter].Id    = nodeid;

          /* store slave_node Id in gsurf */
          actgvol->slave_node[j][k][l] = nodeid;

          /* write nodeid into node_matrix */
          node_matrix[j+1][k+1][l+1]   = nodeid;

          /* set coords of the new node */
          for (m=0; m<3; m++)
            cal_dis->node[node_counter].x[m] =
              x0[m]*(subdivide-j-1) * (subdivide-k-1) * (subdivide-l-1) +
              x1[m]*(j+1)           * (subdivide-k-1) * (subdivide-l-1) +
              x2[m]*(j+1)           * (k+1)           * (subdivide-l-1) +
              x3[m]*(subdivide-j-1) * (k+1)           * (subdivide-l-1) +

              x4[m]*(subdivide-j-1) * (subdivide-k-1) * (l+1) +
              x5[m]*(j+1)           * (subdivide-k-1) * (l+1) +
              x6[m]*(j+1)           * (k+1)           * (l+1) +
              x7[m]*(subdivide-j-1) * (k+1)           * (l+1);

          node_counter++;
          nodeid++;
        }
      }
    }




    /* write corner nodes into node_matrix */
    node_matrix[        0][        0][        0] = gnode0->slave_node;
    node_matrix[subdivide][        0][        0] = gnode1->slave_node;
    node_matrix[subdivide][subdivide][        0] = gnode2->slave_node;
    node_matrix[        0][subdivide][        0] = gnode3->slave_node;

    node_matrix[        0][        0][subdivide] = gnode4->slave_node;
    node_matrix[subdivide][        0][subdivide] = gnode5->slave_node;
    node_matrix[subdivide][subdivide][subdivide] = gnode6->slave_node;
    node_matrix[        0][subdivide][subdivide] = gnode7->slave_node;




    /* write line nodeid into node_matrix */
    for (j=0;j<subdivide-1; j++)
    {
      /* line 0 */
      if (actgvol->line_ind[1][0] == 0)
        node_matrix[j+1][0][0] =
          actgvol->gline[actgvol->line_ind[0][0]]->slave_node[j];
      else
        node_matrix[subdivide-1-j][0][0] =
          actgvol->gline[actgvol->line_ind[0][0]]->slave_node[j];

      /* line 1 */
      if (actgvol->line_ind[1][1] == 0)
        node_matrix[subdivide][j+1][0] =
          actgvol->gline[actgvol->line_ind[0][1]]->slave_node[j];
      else
        node_matrix[subdivide][subdivide-1-j][0] =
          actgvol->gline[actgvol->line_ind[0][1]]->slave_node[j];

      /* line 2 */
      if (actgvol->line_ind[1][2] == 0)
        node_matrix[j+1][subdivide][0] =
          actgvol->gline[actgvol->line_ind[0][2]]->slave_node[j];
      else
        node_matrix[subdivide-1-j][subdivide][0] =
          actgvol->gline[actgvol->line_ind[0][2]]->slave_node[j];

      /* line 3 */
      if (actgvol->line_ind[1][3] == 0)
        node_matrix[0][j+1][0] =
          actgvol->gline[actgvol->line_ind[0][3]]->slave_node[j];
      else
        node_matrix[0][subdivide-1-j][0] =
          actgvol->gline[actgvol->line_ind[0][3]]->slave_node[j];

      /* line 4 */
      if (actgvol->line_ind[1][4] == 0)
        node_matrix[0][0][j+1] =
          actgvol->gline[actgvol->line_ind[0][4]]->slave_node[j];
      else
        node_matrix[0][0][subdivide-1-j] =
          actgvol->gline[actgvol->line_ind[0][4]]->slave_node[j];

      /* line 5 */
      if (actgvol->line_ind[1][5] == 0)
        node_matrix[subdivide][0][j+1] =
          actgvol->gline[actgvol->line_ind[0][5]]->slave_node[j];
      else
        node_matrix[subdivide][0][subdivide-1-j] =
          actgvol->gline[actgvol->line_ind[0][5]]->slave_node[j];

      /* line 6 */
      if (actgvol->line_ind[1][6] == 0)
        node_matrix[subdivide][subdivide][j+1] =
          actgvol->gline[actgvol->line_ind[0][6]]->slave_node[j];
      else
        node_matrix[subdivide][subdivide][subdivide-1-j] =
          actgvol->gline[actgvol->line_ind[0][6]]->slave_node[j];

      /* line 7 */
      if (actgvol->line_ind[1][7] == 0)
        node_matrix[0][subdivide][j+1] =
          actgvol->gline[actgvol->line_ind[0][7]]->slave_node[j];
      else
        node_matrix[0][subdivide][subdivide-1-j] =
          actgvol->gline[actgvol->line_ind[0][7]]->slave_node[j];

      /* line 8 */
      if (actgvol->line_ind[1][8] == 0)
        node_matrix[j+1][0][subdivide] =
          actgvol->gline[actgvol->line_ind[0][8]]->slave_node[j];
      else
        node_matrix[subdivide-1-j][0][subdivide] =
          actgvol->gline[actgvol->line_ind[0][8]]->slave_node[j];

      /* line 9 */
      if (actgvol->line_ind[1][9] == 0)
        node_matrix[subdivide][j+1][subdivide] =
          actgvol->gline[actgvol->line_ind[0][9]]->slave_node[j];
      else
        node_matrix[subdivide][subdivide-1-j][subdivide] =
          actgvol->gline[actgvol->line_ind[0][9]]->slave_node[j];

      /* line 10 */
      if (actgvol->line_ind[1][10] == 0)
        node_matrix[j+1][subdivide][subdivide] =
          actgvol->gline[actgvol->line_ind[0][10]]->slave_node[j];
      else
        node_matrix[subdivide-1-j][subdivide][subdivide] =
          actgvol->gline[actgvol->line_ind[0][10]]->slave_node[j];

      /* line 11 */
      if (actgvol->line_ind[1][11] == 0)
        node_matrix[0][j+1][subdivide] =
          actgvol->gline[actgvol->line_ind[0][11]]->slave_node[j];
      else
        node_matrix[0][subdivide-1-j][subdivide] =
          actgvol->gline[actgvol->line_ind[0][11]]->slave_node[j];
    }



    /* write the surface nodes in node_matrix */
    for (i=0; i<subdivide-1;i++)
    {
      for (j=0; j<subdivide-1;j++)
      {

        GSURF *gsurf0 =  actgvol->gsurf[actgvol->surf_ind[0]];
        GSURF *gsurf1 =  actgvol->gsurf[actgvol->surf_ind[1]];
        GSURF *gsurf2 =  actgvol->gsurf[actgvol->surf_ind[2]];
        GSURF *gsurf3 =  actgvol->gsurf[actgvol->surf_ind[3]];
        GSURF *gsurf4 =  actgvol->gsurf[actgvol->surf_ind[4]];
        GSURF *gsurf5 =  actgvol->gsurf[actgvol->surf_ind[5]];


        /* SURFACE 0 */
        if (gsurf0->gline[0]->gnode[0] == gnode0 && gsurf0->gline[0]->gnode[1] == gnode1)
          node_matrix[i+1][j+1][0] = gsurf0->slave_node[i][j];

        if (gsurf0->gline[0]->gnode[0] == gnode1 && gsurf0->gline[0]->gnode[1] == gnode2)
          node_matrix[subdivide-1-j][i+1][0] = gsurf0->slave_node[i][j];

        if (gsurf0->gline[0]->gnode[0] == gnode2 && gsurf0->gline[0]->gnode[1] == gnode3)
          node_matrix[subdivide-1-i][subdivide-1-j][0] = gsurf0->slave_node[i][j];

        if (gsurf0->gline[0]->gnode[0] == gnode3 && gsurf0->gline[0]->gnode[1] == gnode0)
          node_matrix[j][subdivide-1-i][0] = gsurf0->slave_node[i][j];

        /* up side down */
        if (gsurf0->gline[0]->gnode[0] == gnode1 && gsurf0->gline[0]->gnode[1] == gnode0)
          node_matrix[subdivide-1-i][j+1][0] = gsurf0->slave_node[i][j];

        if (gsurf0->gline[0]->gnode[0] == gnode2 && gsurf0->gline[0]->gnode[1] == gnode1)
          node_matrix[subdivide-1-j][subdivide-1-i][0] = gsurf0->slave_node[i][j];

        if (gsurf0->gline[0]->gnode[0] == gnode3 && gsurf0->gline[0]->gnode[1] == gnode2)
          node_matrix[i+1][subdivide-1-j][0] = gsurf0->slave_node[i][j];

        if (gsurf0->gline[0]->gnode[0] == gnode0 && gsurf0->gline[0]->gnode[1] == gnode3)
          node_matrix[j+1][i+1][0] = gsurf0->slave_node[i][j];



        /* SURFACE 1 */
        if (gsurf1->gline[0]->gnode[0] == gnode0 && gsurf1->gline[0]->gnode[1] == gnode1)
          node_matrix[i+1][0][j+1] = gsurf1->slave_node[i][j];

        if (gsurf1->gline[0]->gnode[0] == gnode1 && gsurf1->gline[0]->gnode[1] == gnode5)
          node_matrix[subdivide-1-j][0][i+1] = gsurf1->slave_node[i][j];

        if (gsurf1->gline[0]->gnode[0] == gnode5 && gsurf1->gline[0]->gnode[1] == gnode4)
          node_matrix[subdivide-1-i][0][subdivide-1-j] = gsurf1->slave_node[i][j];

        if (gsurf1->gline[0]->gnode[0] == gnode4 && gsurf1->gline[0]->gnode[1] == gnode0)
          node_matrix[j+1][0][subdivide-1-i] = gsurf1->slave_node[i][j];

        /* up side down */
        if (gsurf1->gline[0]->gnode[0] == gnode1 && gsurf1->gline[0]->gnode[1] == gnode0)
          node_matrix[subdivide-1-i][0][j+1] = gsurf1->slave_node[i][j];

        if (gsurf1->gline[0]->gnode[0] == gnode5 && gsurf1->gline[0]->gnode[1] == gnode1)
          node_matrix[subdivide-1-j][0][subdivide-1-i] = gsurf1->slave_node[i][j];

        if (gsurf1->gline[0]->gnode[0] == gnode4 && gsurf1->gline[0]->gnode[1] == gnode5)
          node_matrix[i+1][0][subdivide-1-j] = gsurf1->slave_node[i][j];

        if (gsurf1->gline[0]->gnode[0] == gnode0 && gsurf1->gline[0]->gnode[1] == gnode4)
          node_matrix[j+1][0][i+1] = gsurf1->slave_node[i][j];


        /* SURFACE 2 */
        if (gsurf2->gline[0]->gnode[0] == gnode1 && gsurf2->gline[0]->gnode[1] == gnode2)
          node_matrix[subdivide][i+1][j+1] = gsurf2->slave_node[i][j];

        if (gsurf2->gline[0]->gnode[0] == gnode2 && gsurf2->gline[0]->gnode[1] == gnode6)
          node_matrix[subdivide][subdivide-1-j][i+1] = gsurf2->slave_node[i][j];

        if (gsurf2->gline[0]->gnode[0] == gnode6 && gsurf2->gline[0]->gnode[1] == gnode5)
          node_matrix[subdivide][subdivide-1-i][subdivide-1-j] = gsurf2->slave_node[i][j];

        if (gsurf2->gline[0]->gnode[0] == gnode5 && gsurf2->gline[0]->gnode[1] == gnode1)
          node_matrix[subdivide][j+1][subdivide-1-i] = gsurf2->slave_node[i][j];

        /* up side down */
        if (gsurf2->gline[0]->gnode[0] == gnode2 && gsurf2->gline[0]->gnode[1] == gnode1)
          node_matrix[subdivide][subdivide-1-i][j+1] = gsurf2->slave_node[i][j];

        if (gsurf2->gline[0]->gnode[0] == gnode6 && gsurf2->gline[0]->gnode[1] == gnode2)
          node_matrix[subdivide][subdivide-1-j][subdivide-1-i] = gsurf2->slave_node[i][j];

        if (gsurf2->gline[0]->gnode[0] == gnode5 && gsurf2->gline[0]->gnode[1] == gnode6)
          node_matrix[subdivide][i+1][subdivide-1-j] = gsurf2->slave_node[i][j];

        if (gsurf2->gline[0]->gnode[0] == gnode1 && gsurf2->gline[0]->gnode[1] == gnode5)
          node_matrix[subdivide][j+1][i+1] = gsurf2->slave_node[i][j];

        /* SURFACE 3 */
        if (gsurf3->gline[0]->gnode[0] == gnode2 && gsurf3->gline[0]->gnode[1] == gnode3)
          node_matrix[subdivide-1-i][subdivide][j+1] = gsurf3->slave_node[i][j];

        if (gsurf3->gline[0]->gnode[0] == gnode3 && gsurf3->gline[0]->gnode[1] == gnode7)
          node_matrix[j+1][subdivide][i+1] = gsurf3->slave_node[i][j];

        if (gsurf3->gline[0]->gnode[0] == gnode7 && gsurf3->gline[0]->gnode[1] == gnode6)
          node_matrix[i+1][subdivide][subdivide-1-j] = gsurf3->slave_node[i][j];

        if (gsurf3->gline[0]->gnode[0] == gnode6 && gsurf3->gline[0]->gnode[1] == gnode2)
          node_matrix[subdivide-1-j][subdivide][subdivide-1-i] = gsurf3->slave_node[i][j];

        /* up side down */
        if (gsurf3->gline[0]->gnode[0] == gnode3 && gsurf3->gline[0]->gnode[1] == gnode2)
          node_matrix[i+1][subdivide][j+1] = gsurf3->slave_node[i][j];

        if (gsurf3->gline[0]->gnode[0] == gnode7 && gsurf3->gline[0]->gnode[1] == gnode3)
          node_matrix[j+1][subdivide][subdivide-1-i] = gsurf3->slave_node[i][j];

        if (gsurf3->gline[0]->gnode[0] == gnode6 && gsurf3->gline[0]->gnode[1] == gnode7)
          node_matrix[subdivide-1-i][subdivide][subdivide-1-j] = gsurf3->slave_node[i][j];

        if (gsurf3->gline[0]->gnode[0] == gnode2 && gsurf3->gline[0]->gnode[1] == gnode6)
          node_matrix[subdivide-1-j][subdivide][i+1] = gsurf3->slave_node[i][j];


        /* SURFACE 4 */
        if (gsurf4->gline[0]->gnode[0] == gnode3 && gsurf4->gline[0]->gnode[1] == gnode0)
          node_matrix[0][subdivide-1-i][j+1] = gsurf4->slave_node[i][j];

        if (gsurf4->gline[0]->gnode[0] == gnode0 && gsurf4->gline[0]->gnode[1] == gnode4)
          node_matrix[0][j+1][i+1] = gsurf4->slave_node[i][j];

        if (gsurf4->gline[0]->gnode[0] == gnode4 && gsurf4->gline[0]->gnode[1] == gnode7)
          node_matrix[0][i+1][subdivide-1-j] = gsurf4->slave_node[i][j];

        if (gsurf4->gline[0]->gnode[0] == gnode7 && gsurf4->gline[0]->gnode[1] == gnode3)
          node_matrix[0][j+1][subdivide-1-i] = gsurf4->slave_node[i][j];

        /* up side down */
        if (gsurf4->gline[0]->gnode[0] == gnode0 && gsurf4->gline[0]->gnode[1] == gnode3)
          node_matrix[0][i+1][j+1] = gsurf4->slave_node[i][j];

        if (gsurf4->gline[0]->gnode[0] == gnode4 && gsurf4->gline[0]->gnode[1] == gnode0)
          node_matrix[0][j+1][subdivide-1-i] = gsurf4->slave_node[i][j];

        if (gsurf4->gline[0]->gnode[0] == gnode7 && gsurf4->gline[0]->gnode[1] == gnode4)
          node_matrix[0][subdivide-1-i][subdivide-1-j] = gsurf4->slave_node[i][j];

        if (gsurf4->gline[0]->gnode[0] == gnode3 && gsurf4->gline[0]->gnode[1] == gnode7)
          node_matrix[0][subdivide-1-j][i+1] = gsurf4->slave_node[i][j];


        /* SURFACE 5 */
        if (gsurf5->gline[0]->gnode[0] == gnode4 && gsurf5->gline[0]->gnode[1] == gnode5)
          node_matrix[i+1][j+1][subdivide] = gsurf5->slave_node[i][j];

        if (gsurf5->gline[0]->gnode[0] == gnode5 && gsurf5->gline[0]->gnode[1] == gnode6)
          node_matrix[subdivide-1-j][i+1][subdivide] = gsurf5->slave_node[i][j];

        if (gsurf5->gline[0]->gnode[0] == gnode6 && gsurf5->gline[0]->gnode[1] == gnode7)
          node_matrix[subdivide-1-i][subdivide-1-j][subdivide] = gsurf5->slave_node[i][j];

        if (gsurf5->gline[0]->gnode[0] == gnode7 && gsurf5->gline[0]->gnode[1] == gnode4)
          node_matrix[j+1][subdivide-1-i][subdivide] = gsurf5->slave_node[i][j];

        /* up side down */
        if (gsurf5->gline[0]->gnode[0] == gnode5 && gsurf5->gline[0]->gnode[1] == gnode4)
          node_matrix[subdivide-1-i][j+1][subdivide] = gsurf5->slave_node[i][j];

        if (gsurf5->gline[0]->gnode[0] == gnode6 && gsurf5->gline[0]->gnode[1] == gnode5)
          node_matrix[subdivide-1-j][subdivide-1-i][subdivide] = gsurf5->slave_node[i][j];

        if (gsurf5->gline[0]->gnode[0] == gnode7 && gsurf5->gline[0]->gnode[1] == gnode6)
          node_matrix[i+1][subdivide-1-j][subdivide] = gsurf5->slave_node[i][j];

        if (gsurf5->gline[0]->gnode[0] == gnode4 && gsurf5->gline[0]->gnode[1] == gnode7)
          node_matrix[j+1][i+1][subdivide] = gsurf5->slave_node[i][j];

      }  /* for (j=0; j<subdivide-1;j++) */

    }  /* for (i=0; i<subdivide-1;i++) */






#if 0
    /* print the node_matrix */
    /* only for debugging */
    for (k=subdivide; k>=0;k--)
    {
      for (l=subdivide; l>=0;l--)
      {
        for (j=0; j<=subdivide;j++)
        {
          printf(" %5i ",node_matrix[j][k][l]);
        }
        printf("\n");
      }
      printf("\n");
    }
#endif




    /* ==========================================*
     *                                           *
     * generate the new elements                 *
     *                                           *
     * ==========================================*/


    oldele = actgvol->element;
    oldele->master_ele = NULL;
    oldele->slave_ele  = NULL;


    if (oldele->distyp != hex8)
      dserror("Subdivision only possible for uniform discretizations!!");


    for (l=0; l<subdivide; l++)
    {

      for (k=0; k<subdivide; k++)
      {

        for (j=0; j<subdivide; j++)
        {

          actele = &(cal_dis->element[ele_counter]);


          /* set general properties */
          actele->Id      = eleid;
          actele->numnp   = 8;
          actele->distyp  = hex8;

          actele->master_ele = oldele;
          actele->slave_ele  = NULL;



          /* allocate lm */
          actele->lm = (INT*)CCACALLOC(actele->numnp,sizeof(INT));
          if (actele->lm==NULL) dserror("Allocation of lm in ELEMENT failed\n");


          /* fill lm */
          if (subdivide == 1)
            for (m=0; m<actele->numnp; m++)
            {
              actele->lm[m] =  oldele->node[m]->gnode->slave_node;
            }
          else
          {
            actele->lm[0] = node_matrix[j  ][k  ][l  ];
            actele->lm[1] = node_matrix[j+1][k  ][l  ];
            actele->lm[2] = node_matrix[j+1][k+1][l  ];
            actele->lm[3] = node_matrix[j  ][k+1][l  ];
            actele->lm[4] = node_matrix[j  ][k  ][l+1];
            actele->lm[5] = node_matrix[j+1][k  ][l+1];
            actele->lm[6] = node_matrix[j+1][k+1][l+1];
            actele->lm[7] = node_matrix[j  ][k+1][l+1];
          }



          switch (oldele->eltyp)
          {

#ifdef D_ALE
            case el_ale3:

              actele->eltyp   = el_ale3;

              /* allocate the element */
              actele->e.ale3 = (ALE3*)CCACALLOC(1,sizeof(ALE3));
              if (actele->e.ale3==NULL) dserror("Allocation of element ALE3 failed");

              /* copy element properties */
              actele->mat                = oldele->mat;
              actele->e.ale3->nGP[0]     = oldele->e.ale3->nGP[0];
              actele->e.ale3->nGP[1]     = oldele->e.ale3->nGP[1];
              actele->e.ale3->nGP[2]     = oldele->e.ale3->nGP[2];
              actele->e.ale3->jacobi     = oldele->e.ale3->jacobi;
              break;
#endif

#ifdef D_FLUID3
            case el_fluid3:

              actele->eltyp   = el_fluid3;

              /* allocate the element */
              actele->e.f3 = (FLUID3*)CCACALLOC(1,sizeof(FLUID3));
              if (actele->e.f3==NULL) dserror("Allocation of element FLUID3 failed");

              /* copy element properties */
              actele->mat              = oldele->mat;
              actele->e.f3->nGP[0]     = oldele->e.f3->nGP[0];
              actele->e.f3->nGP[1]     = oldele->e.f3->nGP[1];
              actele->e.f3->nGP[2]     = oldele->e.f3->nGP[2];
              actele->e.f3->is_ale     = oldele->e.f3->is_ale;
              actele->e.f3->create_ale = oldele->e.f3->create_ale;
              actele->e.f3->fs_on      = 0;
              break;
#endif


#ifdef D_FLUID3_F
            case el_fluid3_fast:

              actele->eltyp   = el_fluid3_fast;

              /* allocate the element */
              actele->e.f3 = (FLUID3*)CCACALLOC(1,sizeof(FLUID3));
              if (actele->e.f3==NULL) dserror("Allocation of element FLUID3 failed");

              /* copy element properties */
              actele->mat              = oldele->mat;
              actele->e.f3->nGP[0]     = oldele->e.f3->nGP[0];
              actele->e.f3->nGP[1]     = oldele->e.f3->nGP[1];
              actele->e.f3->nGP[2]     = oldele->e.f3->nGP[2];
              actele->e.f3->is_ale     = oldele->e.f3->is_ale;
              actele->e.f3->create_ale = oldele->e.f3->create_ale;
              actele->e.f3->fs_on      = 0;
              break;
#endif


#ifdef D_BRICK1
            case el_brick1:

              actele->eltyp   = el_brick1;

              /* allocate the element */
              actele->e.c1 = (BRICK1*)CCACALLOC(1,sizeof(BRICK1));
              if (actele->e.c1==NULL) dserror("Allocation of element BRICK1 failed");

              /* copy element properties */
              actele->mat              = oldele->mat;
              actele->e.c1->nGP[0]     = oldele->e.c1->nGP[0];
              actele->e.c1->nGP[1]     = oldele->e.c1->nGP[1];
              actele->e.c1->nGP[2]     = oldele->e.c1->nGP[2];
              actele->e.c1->nhyb       = oldele->e.c1->nhyb;
              actele->e.c1->form       = oldele->e.c1->form;
              actele->e.c1->stresstyp  = oldele->e.c1->stresstyp;
              break;
#endif

            case el_none:
            default:
              dserror("Unknown element type!!");
              break;
          }



#if 0
          /* print the location vector */
          /* only for debugging */
          for (m=0; m<actele->numnp;m++)
          {
            printf(" %5i ",actele->lm[m]);
          }
          printf("\n");
#endif


          ele_counter++;
          eleid++;
        }
      }
    }


    if (par.myrank==0)
    {
      p_counter++;
      if (p_counter == pro_fac)
      {
        printf("*");
        fflush(stdout);
        p_counter = 0;
      }
    }
  }  /* for (g=0; g<io_dis->ngvol; g++) */

  if (par.myrank==0)
    printf("|\n");



  /* check the number of created nodes */
  if (node_counter != numnd )
    dserror_args(__FILE__, __LINE__,
        "Wrong number of nodes created!! numnd = %d;  node_counter = %d",
        numnd,node_counter);


  genprob.maxnode = nodeid + 1;
  genprob.nnode += numnd;




  /* check the number of created elelments */
  if (ele_counter != numele )
    dserror_args(__FILE__, __LINE__,
        "Wrong number of elements created!! numnd = %d;  ele_counter = %d",
        numnd,ele_counter);


  genprob.nele += numele;



  /* build pointers between nodes and elements */
  inp_topology(cal_dis);



  /* create all gnodes... and build the pointers between them all */
  inp_detailed_topology(cal_dis);



#if 0
  design->ndnode_fenode = (INT*)CCACALLOC(design->ndnode,sizeof(INT));
  design->dnode_fenode2 = (INT**)CCACALLOC(design->ndnode,sizeof(INT*));

  design->ndline_fenode = (INT*)CCACALLOC(design->ndline,sizeof(INT));;
  design->dline_fenode2 = (INT**)CCACALLOC(design->ndline,sizeof(INT*));;

  design->ndsurf_fenode = (INT*)CCACALLOC(design->ndsurf,sizeof(INT));;
  design->dsurf_fenode2 = (INT**)CCACALLOC(design->ndsurf,sizeof(INT*));;

  design->ndvol_fenode  = (INT*)CCACALLOC(design->ndvol,sizeof(INT));;
  design->dvol_fenode2  = (INT**)CCACALLOC(design->ndvol,sizeof(INT*));;
#endif



  for (i=0; i<design->ndnode; i++)
    design->ndnode_fenode[i] = 0;

  for (i=0; i<design->ndline; i++)
    design->ndline_fenode[i] = 0;

  for (i=0; i<design->ndsurf; i++)
    design->ndsurf_fenode[i] = 0;

  for (i=0; i<design->ndvol; i++)
    design->ndvol_fenode[i] = 0;


  /* dnodes */
  /*========*/


  /* count number of gnodes on each design node */
  for (i=0; i<io_dis->ngnode; i++)
  {
    actgnode = &(io_dis->gnode[i]);
    if (actgnode->ondesigntyp == ondnode)
    {
      id = actgnode->d.dnode->Id;
      design->ndnode_fenode[id]+=1;
    }
  }


  /* allocate dnode_fenode[i] */
  for (i=0; i<design->ndnode; i++)
  {
    if (design->ndnode_fenode[i] > 0)
      design->dnode_fenode2[i] = (INT*)CCAREALLOC(design->dnode_fenode2[i],
          design->ndnode_fenode[i]*sizeof(INT));

    if (design->ndnode_fenode[i] > 1 )
      dserror("Only one fe_node possible on each dnode !!");
  }


  /* fill dnode_fenode */
  for (i=0; i<io_dis->ngnode; i++)
  {
    actgnode = &(io_dis->gnode[i]);
    if (actgnode->ondesigntyp == ondnode)
    {
      id = actgnode->d.dnode->Id;
      design->dnode_fenode2[id][0] = actgnode->slave_node;
    }
  }



  /* dlines */
  /*========*/


  /* count number of gnodes on each design line */
  for (i=0; i<io_dis->ngline; i++)
  {
    actgline = &(io_dis->gline[i]);
    if (actgline->dline != NULL)
    {
      id = actgline->dline->Id;
      design->ndline_fenode[id] += (subdivide + 1);
    }
  }



  /* allocate dline_fenode[i] */
  for (i=0; i<design->ndline; i++)
  {
    if (design->ndline_fenode[i] > 0)
      design->dline_fenode2[i] = (INT*)CCAREALLOC(design->dline_fenode2[i],
          design->ndline_fenode[i]*sizeof(INT));
    counter = 0;


    /* fill dline_fenode */
    for (j=0; j<io_dis->ngline; j++)
    {
      actgline = &(io_dis->gline[j]);
      if (actgline->dline != NULL)
      {
        if (actgline->dline->Id != i)
          continue;
        design->dline_fenode2[i][counter++] = actgline->gnode[0]->slave_node;
        design->dline_fenode2[i][counter++] = actgline->gnode[1]->slave_node;
        /* Mittelknoten */
        for (k=0; k<subdivide-1; k++)
          design->dline_fenode2[i][counter++] = actgline->slave_node[k];
      }
    }

    if (counter != design->ndline_fenode[i])
      dserror_args(__FILE__, __LINE__,
          "Irgendetwas stimmt hier nicht !! counter = %d;  ndline_fenode[i] =%d",
          counter,design->ndline_fenode[i] );

  }




  /* dsurf */
  /*========*/


  /* count number of gnodes on each design surf */
  for (i=0; i<io_dis->ngsurf; i++)
  {
    actgsurf = &(io_dis->gsurf[i]);
    if (actgsurf->dsurf != NULL)
    {
      id = actgsurf->dsurf->Id;
      design->ndsurf_fenode[id] += (subdivide + 1) * (subdivide + 1);
    }
  }



  /* allocate dsurf_fenode[i] */
  for (i=0; i<design->ndsurf; i++)
  {
    if (design->ndsurf_fenode[i] > 0)
      design->dsurf_fenode2[i] = (INT*)CCAREALLOC(design->dsurf_fenode2[i],
          design->ndsurf_fenode[i]*sizeof(INT));
    counter = 0;


    /* fill dsurf_fenode */
    for (j=0; j<io_dis->ngsurf; j++)
    {
      actgsurf = &(io_dis->gsurf[j]);
      if (actgsurf->dsurf != NULL)
      {
        if (actgsurf->dsurf->Id != i)
          continue;
        design->dsurf_fenode2[i][counter++] = actgsurf->gnode[0]->slave_node;
        design->dsurf_fenode2[i][counter++] = actgsurf->gnode[1]->slave_node;
        design->dsurf_fenode2[i][counter++] = actgsurf->gnode[2]->slave_node;
        design->dsurf_fenode2[i][counter++] = actgsurf->gnode[3]->slave_node;
        /* Mittelknoten der Flaeche */
        for (k=0; k<subdivide-1; k++)
          for (l=0; l<subdivide-1; l++)
            design->dsurf_fenode2[i][counter++] = actgsurf->slave_node[k][l];
        /* Mittelknoten der Linien */
        for (k=0; k<actgsurf->ngline; k++)
          for (l=0; l<subdivide-1; l++)
            design->dsurf_fenode2[i][counter++] = actgsurf->gline[k]->slave_node[l];

      }
    }

    if (counter != design->ndsurf_fenode[i])
      dserror_args(__FILE__, __LINE__,
          "Irgendetwas stimmt hier nicht !! counter = %d;  ndsurf_fenode[i] =%d",
          counter,design->ndsurf_fenode[i] );

  }




  /* dvol */
  /*========*/


  /* counter number of gnodes on each design vol */
  for (i=0; i<io_dis->ngvol; i++)
  {
    actgvol = &(io_dis->gvol[i]);
    if (actgvol->dvol != NULL)
    {
      id = actgvol->dvol->Id;
      design->ndvol_fenode[id] += (subdivide + 1) * (subdivide + 1) * (subdivide + 1);
    }
  }



  /* allocate dvol_fenode[i] */
  for (i=0; i<design->ndvol; i++)
  {
    if (design->ndvol_fenode[i] > 0)
      design->dvol_fenode2[i] = (INT*)CCAREALLOC(design->dvol_fenode2[i],
          design->ndvol_fenode[i]*sizeof(INT));
    counter = 0;


    /* fill dvol_fenode */
    for (j=0; j<io_dis->ngvol; j++)
    {
      actgvol = &(io_dis->gvol[j]);
      if (actgvol->dvol != NULL)
      {
        if (actgvol->dvol->Id != i)
          continue;

        counter2 = 0;
        for (k=0; k<6; k++)
          for (l=0; l<4; l++)
            tmp[counter2++] = actgvol->gsurf[k]->gnode[l]->slave_node;

        for (k=0; k<24; k++)
          for (l=0; l<24; l++)
            if (tmp[k] == tmp[l] && k!=l)
              tmp[l] = -1;


        for (k=0; k<24; k++)
          if (tmp[k] != -1)
            design->dvol_fenode2[i][counter++] = tmp[k];


        /* Mittelknoten des Vols */
        for (k=0; k<subdivide-1; k++)
          for (l=0; l<subdivide-1; l++)
            for (m=0; m<subdivide-1; m++)
              design->dvol_fenode2[i][counter++] = actgvol->slave_node[k][l][m];

        /* Mittelknoten der Flaechen */
        for (k=0; k<actgvol->ngsurf; k++)
          for (l=0; l<subdivide-1; l++)
            for (m=0; m<subdivide-1; m++)
              design->dvol_fenode2[i][counter++] = actgvol->gsurf[k]->slave_node[l][m];
        /* Mittelknoten der Linien */
        for (k=0; k<actgvol->ngline; k++)
          for (l=0; l<subdivide-1; l++)
            design->dvol_fenode2[i][counter++] = actgvol->gline[k]->slave_node[l];
      }
    }

    if (counter != design->ndvol_fenode[i])
      dserror_args(__FILE__, __LINE__,
          "Irgendetwas stimmt hier nicht !! counter = %d;  ndvol_fenode[i] =%d",
          counter,design->ndvol_fenode[i] );

  }



  /* create the pointers between gnodes... and the design */
  inpdesign_topology_fe(cal_dis);


#ifdef DEBUG
  dstrc_exit();
#endif


  return;

}  /* END of global_subdivide_hex */















/*-----------------------------------------------------------------------*/
/*!
  \brief Brief description

  Very detailed description what the function is doing why.

  \param b      DOUBLE  (o) explanation

  \return void

  \author mn
  \date   09/05

 */
/*-----------------------------------------------------------------------*/
void global_subdivide_tet(
    FIELD        *actfield
    )

{

  DISCRET       *io_dis;
  DISCRET       *cal_dis;

  INT            subdivide;
  INT            nodeid;
  INT            eleid;

  INT            numnd;
  INT            numele;
  INT            i,j;
  INT            node_counter;
  INT            ele_counter;

  GNODE         *actgnode = NULL;
  GLINE         *actgline = NULL;
  GSURF         *actgsurf = NULL;
  GVOL          *actgvol  = NULL;

  ELEMENT       *actele, *oldele;

  INT            id;

  INT            tmp[24];
  INT            counter;
  INT            counter2;
  INT            k,l,m,g;

  INT            node_matrix[MAX_DIVIDE+1][MAX_DIVIDE+1][MAX_DIVIDE+1];

  INT            pro_total=0, pro_fac=0, p=0, p_counter=0;

#ifdef DEBUG
  dstrc_enter("global_subdivide_tet");
#endif


  io_dis    = &(actfield->dis[0]);
  cal_dis   = &(actfield->dis[1]);

  subdivide = actfield->subdivide;
  nodeid    = genprob.maxnode;
  eleid     = genprob.nele;


  if (subdivide != 2)
  {
    dserror("Subdivision of tets only possible for subdivide = 2!!!");
  }


  /* calculate number of nodes on second dis */
  numnd = io_dis->ngnode + io_dis->ngline + io_dis->ngsurf + io_dis->ngvol;



  /* allocate the nodes to the dis */
  cal_dis->node  = (NODE*)CCACALLOC(numnd,sizeof(NODE));
  cal_dis->numnp = numnd;
  node_counter = 0;



  /* calculate number of elements in second dis */
  numele = io_dis->ngvol * 4;


  /* allocate elements in second dis */
  cal_dis->element = (ELEMENT*)CCACALLOC(numele,sizeof(ELEMENT));
  cal_dis->numele  = numele;
  ele_counter = 0;






  /* loop gnodes and create nodes on second dis */
  for (i=0; i<io_dis->ngnode; i++)
  {
    actgnode = &(io_dis->gnode[i]);

    /* set Id of the new node */
    cal_dis->node[node_counter].Id = nodeid;

    /* set coords of the new node */
    for (j=0; j<3; j++)
      cal_dis->node[node_counter].x[j] = actgnode->node->x[j];

    /* store slave_node Id in master gnode */
    actgnode->slave_node = nodeid;

    /* create pointer from slave to master node */
    cal_dis->node[node_counter].slave_node  = NULL;
    cal_dis->node[node_counter].master_node = actgnode->node;

    /* create pointer from master to slave node */
    actgnode->node->slave_node  = &(cal_dis->node[node_counter]);
    actgnode->node->master_node = NULL;


    node_counter++;
    nodeid++;
  }



  /* loop glines and create nodes */
  for (i=0; i<io_dis->ngline; i++)
  {

    DOUBLE x1[3],dx[3];

    actgline = &(io_dis->gline[i]);

    /* coordinates of the beginning of this gline */
    x1[0] = actgline->gnode[0]->node->x[0];
    x1[1] = actgline->gnode[0]->node->x[1];
    x1[2] = actgline->gnode[0]->node->x[2];

    /* lenght of this gline in the three directions divided by subdivide */
    dx[0] = (actgline->gnode[1]->node->x[0] - actgline->gnode[0]->node->x[0])/2;
    dx[1] = (actgline->gnode[1]->node->x[1] - actgline->gnode[0]->node->x[1])/2;
    dx[2] = (actgline->gnode[1]->node->x[2] - actgline->gnode[0]->node->x[2])/2;

    cal_dis->node[node_counter].Id = nodeid;

    /* store slave_node Id in gline */
    actgline->slave_node[0]   = nodeid;

    /* set coords of the new node */
    for (k=0; k<3; k++)
      cal_dis->node[node_counter].x[k] =  x1[k] + dx[k];

    node_counter++;
    nodeid++;

  }  /* for (i=0; i<io_dis->ngline; i++) */











  /* loop gsurfs and create nodes */
  for (i=0; i<io_dis->ngsurf; i++)
  {
    GNODE  *gnode0, *gnode1, *gnode2;

    actgsurf = &(io_dis->gsurf[i]);


    /* specify the nodes */
    gnode0 = actgsurf->gnode[0];
    gnode1 = actgsurf->gnode[1];
    gnode2 = actgsurf->gnode[2];


    /* generate the node */
    cal_dis->node[node_counter].Id  = nodeid;

    /* store slave_node Id in gsurf */
    actgsurf->slave_node[0][0] = nodeid;

    /* set coords of the new node */
    for (l=0; l<3; l++)
      cal_dis->node[node_counter].x[l] =
        (gnode0->node->x[l] + gnode1->node->x[l] + gnode2->node->x[l]) / 3;


    node_counter++;
    nodeid++;

  }  /* for (i=0; i<io_dis->ngsurf; i++) */






  if (par.myrank==0)
  {
    printf("  creating subelements on %3i GVOLS\n",io_dis->ngvol);

    pro_fac   = io_dis->ngvol/50 + 1;

    pro_total = io_dis->ngvol/pro_fac;

    printf("  0");
    for (p=0; p<pro_total; p++)
      printf(" ");
    printf("%3i\n",io_dis->ngvol);

    printf("  |");
    for (p=0; p<pro_total; p++)
      printf("-");
    printf("|\n");


    printf("  |");
    p_counter = 0;
  }

  /* loop gvols and create nodes */
  for (g=0; g<io_dis->ngvol; g++)
  {

    NODE    *node0,  *node1,  *node2,  *node3;
    GNODE  *gnode0, *gnode1, *gnode2, *gnode3;


    actgvol = &(io_dis->gvol[g]);



#if 0
    /* print the node_matrix */
    /* only for debugging */
    for (j=2; j>=0;j--)
    {
      for (k=2; k>=0;k--)
      {
        for (i=0; i<=2;i++)
        {
         node_matrix[i][j][k] = -17;
        }
      }
    }
#endif



    node0  = actgvol->element->node[0];
    gnode0 = node0->gnode;
    node1  = actgvol->element->node[1];
    gnode1 = node1->gnode;
    node2  = actgvol->element->node[2];
    gnode2 = node2->gnode;
    node3  = actgvol->element->node[3];
    gnode3 = node3->gnode;


    /* generate the node */
    cal_dis->node[node_counter].Id    = nodeid;

    /* store slave_node Id in gsurf */
    actgvol->slave_node[0][0][0] = nodeid;

    /* write nodeid into node_matrix */
    node_matrix[1][1][1]         = nodeid;

    /* set coords of the new node */
    for (m=0; m<3; m++)
      cal_dis->node[node_counter].x[m] =
        (node0->x[m] + node1->x[m] + node2->x[m] + node3->x[m]) / 4;

    node_counter++;
    nodeid++;



    /* identify the lines */


    /* LINE 0: 0-1 */
    /*==============*/
    for (j=0; j<actgvol->ngline; j++)
    {
      /* this line is in same direction */
      if (gnode0 == actgvol->gline[j]->gnode[0] &&
          gnode1 == actgvol->gline[j]->gnode[1])
      {
        actgline = actgvol->gline[j];
        actgvol->line_ind[0][0] = j;
        actgvol->line_ind[1][0] = 0;
        break;
      }

      /* this line is in opposite direction */
      if (gnode1 == actgvol->gline[j]->gnode[0] &&
          gnode0 == actgvol->gline[j]->gnode[1])
      {
        actgline = actgvol->gline[j];
        actgvol->line_ind[0][0] = j;
        actgvol->line_ind[1][0] = 1;
        break;
      }
    }  /* for (j=0; j<gsurf->ngline; j++) */




    /* LINE 1: 1-2 */
    /*==============*/
    for (j=0; j<actgvol->ngline; j++)
    {
      /* this line is in same direction */
      if (gnode1 == actgvol->gline[j]->gnode[0] &&
          gnode2 == actgvol->gline[j]->gnode[1])
      {
        actgline = actgvol->gline[j];
        actgvol->line_ind[0][1] = j;
        actgvol->line_ind[1][1] = 0;
        break;
      }

      /* this line is in opposite direction */
      if (gnode2 == actgvol->gline[j]->gnode[0] &&
          gnode1 == actgvol->gline[j]->gnode[1])
      {
        actgline = actgvol->gline[j];
        actgvol->line_ind[0][1] = j;
        actgvol->line_ind[1][1] = 1;
        break;
      }
    }  /* for (j=0; j<gsurf->ngline; j++) */



    /* LINE 2: 2-0 */
    /*==============*/
    for (j=0; j<actgvol->ngline; j++)
    {
      /* this line is in same direction */
      if (gnode2 == actgvol->gline[j]->gnode[0] &&
          gnode0 == actgvol->gline[j]->gnode[1])
      {
        actgline = actgvol->gline[j];
        actgvol->line_ind[0][2] = j;
        actgvol->line_ind[1][2] = 0;
        break;
      }

      /* this line is in opposite direction */
      if (gnode0 == actgvol->gline[j]->gnode[0] &&
          gnode2 == actgvol->gline[j]->gnode[1])
      {
        actgline = actgvol->gline[j];
        actgvol->line_ind[0][2] = j;
        actgvol->line_ind[1][2] = 1;
        break;
      }
    }  /* for (j=0; j<gsurf->ngline; j++) */



    /* LINE 3: 0-3 */
    /*==============*/
    for (j=0; j<actgvol->ngline; j++)
    {
      /* this line is in same direction */
      if (gnode0 == actgvol->gline[j]->gnode[0] &&
          gnode3 == actgvol->gline[j]->gnode[1])
      {
        actgline = actgvol->gline[j];
        actgvol->line_ind[0][3] = j;
        actgvol->line_ind[1][3] = 0;
        break;
      }

      /* this line is in opposite direction */
      if (gnode3 == actgvol->gline[j]->gnode[0] &&
          gnode0 == actgvol->gline[j]->gnode[1])
      {
        actgline = actgvol->gline[j];
        actgvol->line_ind[0][3] = j;
        actgvol->line_ind[1][3] = 1;
        break;
      }
    }  /* for (j=0; j<gsurf->ngline; j++) */



    /* LINE 4: 1-3 */
    /*==============*/
    for (j=0; j<actgvol->ngline; j++)
    {
      /* this line is in same direction */
      if (gnode1 == actgvol->gline[j]->gnode[0] &&
          gnode3 == actgvol->gline[j]->gnode[1])
      {
        actgline = actgvol->gline[j];
        actgvol->line_ind[0][4] = j;
        actgvol->line_ind[1][4] = 0;
        break;
      }

      /* this line is in opposite direction */
      if (gnode3 == actgvol->gline[j]->gnode[0] &&
          gnode1 == actgvol->gline[j]->gnode[1])
      {
        actgline = actgvol->gline[j];
        actgvol->line_ind[0][4] = j;
        actgvol->line_ind[1][4] = 1;
        break;
      }
    }  /* for (j=0; j<gsurf->ngline; j++) */



    /* LINE 5: 2-3 */
    /*==============*/
    for (j=0; j<actgvol->ngline; j++)
    {
      /* this line is in same direction */
      if (gnode2 == actgvol->gline[j]->gnode[0] &&
          gnode3 == actgvol->gline[j]->gnode[1])
      {
        actgline = actgvol->gline[j];
        actgvol->line_ind[0][5] = j;
        actgvol->line_ind[1][5] = 0;
        break;
      }

      /* this line is in opposite direction */
      if (gnode3 == actgvol->gline[j]->gnode[0] &&
          gnode2 == actgvol->gline[j]->gnode[1])
      {
        actgline = actgvol->gline[j];
        actgvol->line_ind[0][5] = j;
        actgvol->line_ind[1][5] = 1;
        break;
      }
    }  /* for (j=0; j<gsurf->ngline; j++) */




    /* identify the surfaces */

    /* GSURF 0 */
    /*=========*/
    for (j=0; j<actgvol->ngsurf; j++)
    {
      INT found = 0;

      for (k=0; k<actgvol->gsurf[j]->ngline; k++)
        if (actgvol->gsurf[j]->gline[k] == actgvol->gline[actgvol->line_ind[0][0]])
          found++;

      for (k=0; k<actgvol->gsurf[j]->ngline; k++)
        if (actgvol->gsurf[j]->gline[k] == actgvol->gline[actgvol->line_ind[0][1]])
          found++;

      if (found == 2)
      {
        actgsurf = actgvol->gsurf[j];
        actgvol->surf_ind[0] = j;
        break;
      }
    }  /* for (j=0; j<actgvol->ngsurf; j++) */
    if (actgsurf == NULL)
      dserror("Surface 0 not found!!\n");



    /* GSURF 1 */
    /*=========*/
    for (j=0; j<actgvol->ngsurf; j++)
    {
      INT found = 0;

      for (k=0; k<actgvol->gsurf[j]->ngline; k++)
        if (actgvol->gsurf[j]->gline[k] == actgvol->gline[actgvol->line_ind[0][0]])
          found++;

      for (k=0; k<actgvol->gsurf[j]->ngline; k++)
        if (actgvol->gsurf[j]->gline[k] == actgvol->gline[actgvol->line_ind[0][4]])
          found++;

      if (found == 2)
      {
        actgsurf = actgvol->gsurf[j];
        actgvol->surf_ind[1] = j;
        break;
      }
    }  /* for (j=0; j<actgvol->ngsurf; j++) */
    if (actgsurf == NULL)
      dserror("Surface 1 not found!!\n");



    /* GSURF 2 */
    /*=========*/
    for (j=0; j<actgvol->ngsurf; j++)
    {
      INT found = 0;

      for (k=0; k<actgvol->gsurf[j]->ngline; k++)
        if (actgvol->gsurf[j]->gline[k] == actgvol->gline[actgvol->line_ind[0][4]])
          found++;

      for (k=0; k<actgvol->gsurf[j]->ngline; k++)
        if (actgvol->gsurf[j]->gline[k] == actgvol->gline[actgvol->line_ind[0][1]])
          found++;

      if (found == 2)
      {
        actgsurf = actgvol->gsurf[j];
        actgvol->surf_ind[2] = j;
        break;
      }
    }  /* for (j=0; j<actgvol->ngsurf; j++) */
    if (actgsurf == NULL)
      dserror("Surface 2 not found!!\n");



    /* GSURF 3 */
    /*=========*/
    for (j=0; j<actgvol->ngsurf; j++)
    {
      INT found = 0;

      for (k=0; k<actgvol->gsurf[j]->ngline; k++)
        if (actgvol->gsurf[j]->gline[k] == actgvol->gline[actgvol->line_ind[0][5]])
          found++;

      for (k=0; k<actgvol->gsurf[j]->ngline; k++)
        if (actgvol->gsurf[j]->gline[k] == actgvol->gline[actgvol->line_ind[0][2]])
          found++;

      if (found == 2)
      {
        actgsurf = actgvol->gsurf[j];
        actgvol->surf_ind[3] = j;
        break;
      }
    }  /* for (j=0; j<actgvol->ngsurf; j++) */
    if (actgsurf == NULL)
      dserror("Surface 3 not found!!\n");




    /* write corner nodes into node_matrix */
    node_matrix[0][0][0] = gnode0->slave_node;
    node_matrix[2][0][0] = gnode1->slave_node;
    node_matrix[0][2][0] = gnode2->slave_node;
    node_matrix[0][0][2] = gnode3->slave_node;





    /* write line nodeid into node_matrix */

    /* line 0 */
    node_matrix[1][0][0] =
      actgvol->gline[actgvol->line_ind[0][0]]->slave_node[0];

    /* line 1 */
    node_matrix[2][1][0] =
      actgvol->gline[actgvol->line_ind[0][1]]->slave_node[0];

    /* line 2 */
    node_matrix[0][1][0] =
      actgvol->gline[actgvol->line_ind[0][2]]->slave_node[0];

    /* line 3 */
    node_matrix[0][0][1] =
      actgvol->gline[actgvol->line_ind[0][3]]->slave_node[0];

    /* line 4 */
    node_matrix[2][0][1] =
      actgvol->gline[actgvol->line_ind[0][4]]->slave_node[0];

    /* line 5 */
    node_matrix[0][2][1] =
      actgvol->gline[actgvol->line_ind[0][5]]->slave_node[0];




    /* write the surface nodes in node_matrix */

    /* SURFACE 0 */
    node_matrix[1][1][0] = actgvol->gsurf[actgvol->surf_ind[0]]->slave_node[0][0];

    /* SURFACE 1 */
    node_matrix[1][0][1] = actgvol->gsurf[actgvol->surf_ind[1]]->slave_node[0][0];

    /* SURFACE 2 */
    node_matrix[2][1][1] = actgvol->gsurf[actgvol->surf_ind[2]]->slave_node[0][0];

    /* SURFACE 3 */
    node_matrix[0][1][1] = actgvol->gsurf[actgvol->surf_ind[3]]->slave_node[0][0];




#if 0
    /* print the node_matrix */
    /* only for debugging */
    for (j=2; j>=0;j--)
    {
      for (k=2; k>=0;k--)
      {
        for (i=0; i<=2;i++)
        {
          printf(" %5i ",node_matrix[i][j][k]);
        }
        printf("\n");
      }
      printf("\n");
    }
#endif



    /* ==========================================*
     *                                           *
     * generate the new elements                 *
     *                                           *
     * ==========================================*/


    oldele = actgvol->element;

    oldele->master_ele = NULL;
    oldele->slave_ele  = NULL;


    if (oldele->distyp != tet4)
      dserror("Subdivision only possible for uniform discretizations!!");


    for (l=0; l<4; l++)
    {

      actele = &(cal_dis->element[ele_counter]);


      /* set general properties */
      actele->Id      = eleid;
      actele->numnp   = 8;
      actele->distyp  = hex8;

      actele->master_ele = oldele;
      actele->slave_ele  = NULL;


      /* allocate lm */
      actele->lm = (INT*)CCACALLOC(actele->numnp,sizeof(INT));
      if (actele->lm==NULL) dserror("Allocation of lm in ELEMENT failed\n");


      /* fill lm */
      switch (l)
      {
        case 0:
          actele->lm[0] = node_matrix[0][0][0];
          actele->lm[1] = node_matrix[1][0][0];
          actele->lm[2] = node_matrix[1][1][0];
          actele->lm[3] = node_matrix[0][1][0];
          actele->lm[4] = node_matrix[0][0][1];
          actele->lm[5] = node_matrix[1][0][1];
          actele->lm[6] = node_matrix[1][1][1];
          actele->lm[7] = node_matrix[0][1][1];
          break;

        case 1:
          actele->lm[0] = node_matrix[2][0][0];
          actele->lm[1] = node_matrix[2][1][0];
          actele->lm[2] = node_matrix[1][1][0];
          actele->lm[3] = node_matrix[1][0][0];
          actele->lm[4] = node_matrix[2][0][1];
          actele->lm[5] = node_matrix[2][1][1];
          actele->lm[6] = node_matrix[1][1][1];
          actele->lm[7] = node_matrix[1][0][1];
          break;

        case 2:
          actele->lm[0] = node_matrix[0][2][0];
          actele->lm[1] = node_matrix[0][1][0];
          actele->lm[2] = node_matrix[1][1][0];
          actele->lm[3] = node_matrix[2][1][0];
          actele->lm[4] = node_matrix[0][2][1];
          actele->lm[5] = node_matrix[0][1][1];
          actele->lm[6] = node_matrix[1][1][1];
          actele->lm[7] = node_matrix[2][1][1];
          break;

        case 3:
          actele->lm[0] = node_matrix[0][0][2];
          actele->lm[1] = node_matrix[0][0][1];
          actele->lm[2] = node_matrix[0][1][1];
          actele->lm[3] = node_matrix[0][2][1];
          actele->lm[4] = node_matrix[2][0][1];
          actele->lm[5] = node_matrix[1][0][1];
          actele->lm[6] = node_matrix[1][1][1];
          actele->lm[7] = node_matrix[2][1][1];
          break;

        default:
          dserror("This is NOT possible!!!");
          break;
      }





      switch (oldele->eltyp)
      {

#ifdef D_ALE
        case el_ale3:

          actele->eltyp   = el_ale3;

          /* allocate the element */
          actele->e.ale3 = (ALE3*)CCACALLOC(1,sizeof(ALE3));
          if (actele->e.ale3==NULL) dserror("Allocation of element ALE3 failed");

          /* copy element properties */
          actele->mat                = oldele->mat;
          actele->e.ale3->jacobi     = oldele->e.ale3->jacobi;

          switch (oldele->e.ale3->nGP[0])
          {
            case 1:
              actele->e.ale3->nGP[0]     = 1;
              actele->e.ale3->nGP[1]     = 1;
              actele->e.ale3->nGP[2]     = 1;
              break;

            case 4:
              actele->e.ale3->nGP[0]     = 2;
              actele->e.ale3->nGP[1]     = 2;
              actele->e.ale3->nGP[2]     = 2;
              break;

            default:
              dserror("Unknown number of Gauss Points for tet element!!");
              break;
          }
          break;
#endif

        case el_brick1:
          dserror("Subdivison not yet implemented for brick1");
          break;

#ifdef D_FLUID3
        case el_fluid3:

          actele->eltyp   = el_fluid3;

          /* allocate the element */
          actele->e.f3 = (FLUID3*)CCACALLOC(1,sizeof(FLUID3));
          if (actele->e.f3==NULL) dserror("Allocation of element FLUID3 failed");

          /* copy element properties */
          actele->mat              = oldele->mat;
          actele->e.f3->is_ale     = oldele->e.f3->is_ale;
          actele->e.f3->create_ale = oldele->e.f3->create_ale;
          actele->e.f3->fs_on      = 0;

          switch (oldele->e.f3->nGP[0])
          {
            case 1:
              actele->e.f3->nGP[0]     = 1;
              actele->e.f3->nGP[1]     = 1;
              actele->e.f3->nGP[2]     = 1;
              break;

            case 4:
              actele->e.f3->nGP[0]     = 2;
              actele->e.f3->nGP[1]     = 2;
              actele->e.f3->nGP[2]     = 2;
              break;

            default:
              dserror("Unknown number of Gauss Points for tet element!!");
              break;
          }
          break;
#endif


#ifdef D_FLUID3_F
        case el_fluid3_fast:

          actele->eltyp   = el_fluid3_fast;

          /* allocate the element */
          actele->e.f3 = (FLUID3*)CCACALLOC(1,sizeof(FLUID3));
          if (actele->e.f3==NULL) dserror("Allocation of element FLUID3 failed");

          /* copy element properties */
          actele->mat              = oldele->mat;
          actele->e.f3->is_ale     = oldele->e.f3->is_ale;
          actele->e.f3->create_ale = oldele->e.f3->create_ale;
          actele->e.f3->fs_on      = 0;

          switch (oldele->e.f3->nGP[0])
          {
            case 1:
              actele->e.f3->nGP[0]     = 1;
              actele->e.f3->nGP[1]     = 1;
              actele->e.f3->nGP[2]     = 1;
              break;

            case 4:
              actele->e.f3->nGP[0]     = 2;
              actele->e.f3->nGP[1]     = 2;
              actele->e.f3->nGP[2]     = 2;
              break;

            default:
              dserror("Unknown number of Gauss Points for tet element!!");
              break;
          }
          break;
#endif

        case el_none:
        default:
          dserror("Unknown element type!!");
          break;
      }



#if 0
      /* print the location vector */
      /* only for debugging */
      for (m=0; m<actele->numnp;m++)
      {
        printf(" %5i ",actele->lm[m]);
      }
      printf("\n");
#endif


      ele_counter++;
      eleid++;
    }

    if (par.myrank==0)
    {
      p_counter++;
      if (p_counter == pro_fac)
      {
        printf("*");
        p_counter = 0;
      }
    }
  }  /* for (g=0; g<io_dis->ngvol; g++) */

  if (par.myrank==0)
    printf("|\n");




  /* check the number of created nodes */
  if (node_counter != numnd )
    dserror_args(__FILE__, __LINE__,
        "Wrong number of nodes created!! numnd = %d;  node_counter = %d",
        numnd,node_counter);


  genprob.maxnode = nodeid + 1;
  genprob.nnode += numnd;




  /* check the number of created elelments */
  if (ele_counter != numele )
    dserror_args(__FILE__, __LINE__,
        "Wrong number of elements created!! numnd = %d;  ele_counter = %d",
        numnd,ele_counter);


  genprob.nele += numele;



  /* build pointers between nodes and elements */
  inp_topology(cal_dis);



  /* create all gnodes... and build the pointers between them all */
  inp_detailed_topology(cal_dis);



  design->ndnode_fenode = (INT*)CCACALLOC(design->ndnode,sizeof(INT));
  design->dnode_fenode2 = (INT**)CCACALLOC(design->ndnode,sizeof(INT*));

  design->ndline_fenode = (INT*)CCACALLOC(design->ndline,sizeof(INT));;
  design->dline_fenode2 = (INT**)CCACALLOC(design->ndline,sizeof(INT*));;

  design->ndsurf_fenode = (INT*)CCACALLOC(design->ndsurf,sizeof(INT));;
  design->dsurf_fenode2 = (INT**)CCACALLOC(design->ndsurf,sizeof(INT*));;

  design->ndvol_fenode  = (INT*)CCACALLOC(design->ndvol,sizeof(INT));;
  design->dvol_fenode2  = (INT**)CCACALLOC(design->ndvol,sizeof(INT*));;



  for (i=0; i<design->ndnode; i++)
    design->ndnode_fenode[i] = 0;

  for (i=0; i<design->ndline; i++)
    design->ndline_fenode[i] = 0;

  for (i=0; i<design->ndsurf; i++)
    design->ndsurf_fenode[i] = 0;

  for (i=0; i<design->ndvol; i++)
    design->ndvol_fenode[i] = 0;


  /* dnodes */
  /*========*/


  /* count number of gnodes on each design node */
  for (i=0; i<io_dis->ngnode; i++)
  {
    actgnode = &(io_dis->gnode[i]);
    if (actgnode->ondesigntyp == ondnode)
    {
      id = actgnode->d.dnode->Id;
      design->ndnode_fenode[id]+=1;
    }
  }


  /* allocate dnode_fenode[i] */
  for (i=0; i<design->ndnode; i++)
  {
    if (design->ndnode_fenode[i] > 0)
      design->dnode_fenode2[i] = (INT*)CCAREALLOC(design->dnode_fenode2[i],
          design->ndnode_fenode[i]*sizeof(INT));

    if (design->ndnode_fenode[i] > 1 )
      dserror("Only one fe_node possible on each dnode !!");

  }


  /* fill dnode_fenode */
  for (i=0; i<io_dis->ngnode; i++)
  {
    actgnode = &(io_dis->gnode[i]);
    if (actgnode->ondesigntyp == ondnode)
    {
      id = actgnode->d.dnode->Id;
      design->dnode_fenode2[id][0] = actgnode->slave_node;
    }
  }



  /* dlines */
  /*========*/


  /* count number of gnodes on each design line */
  for (i=0; i<io_dis->ngline; i++)
  {
    actgline = &(io_dis->gline[i]);
    if (actgline->dline != NULL)
    {
      id = actgline->dline->Id;
      design->ndline_fenode[id] += 3;
    }
  }



  /* allocate dline_fenode[i] */
  for (i=0; i<design->ndline; i++)
  {
    if (design->ndline_fenode[i] > 0)
      design->dline_fenode2[i] = (INT*)CCAREALLOC(design->dline_fenode2[i],
          design->ndline_fenode[i]*sizeof(INT));
    counter = 0;


    /* fill dline_fenode */
    for (j=0; j<io_dis->ngline; j++)
    {
      actgline = &(io_dis->gline[j]);
      if (actgline->dline != NULL)
      {
        if (actgline->dline->Id != i)
          continue;
        design->dline_fenode2[i][counter++] = actgline->gnode[0]->slave_node;
        design->dline_fenode2[i][counter++] = actgline->gnode[1]->slave_node;
        /* Mittelknoten */
        design->dline_fenode2[i][counter++] = actgline->slave_node[0];
      }
    }

    if (counter != design->ndline_fenode[i])
      dserror_args(__FILE__, __LINE__,
          "Irgendetwas stimmt hier nicht !! counter = %d;  ndline_fenode[i] =%d",
          counter,design->ndline_fenode[i] );

  }




  /* dsurf */
  /*========*/


  /* count number of gnodes on each design surf */
  for (i=0; i<io_dis->ngsurf; i++)
  {
    actgsurf = &(io_dis->gsurf[i]);
    if (actgsurf->dsurf != NULL)
    {
      id = actgsurf->dsurf->Id;
      design->ndsurf_fenode[id] += 7;
    }
  }



  /* allocate dsurf_fenode[i] */
  for (i=0; i<design->ndsurf; i++)
  {
    if (design->ndsurf_fenode[i] > 0)
      design->dsurf_fenode2[i] = (INT*)CCAREALLOC(design->dsurf_fenode2[i],
          design->ndsurf_fenode[i]*sizeof(INT));
    counter = 0;


    /* fill dsurf_fenode */
    for (j=0; j<io_dis->ngsurf; j++)
    {
      actgsurf = &(io_dis->gsurf[j]);
      if (actgsurf->dsurf != NULL)
      {
        if (actgsurf->dsurf->Id != i)
          continue;
        design->dsurf_fenode2[i][counter++] = actgsurf->gnode[0]->slave_node;
        design->dsurf_fenode2[i][counter++] = actgsurf->gnode[1]->slave_node;
        design->dsurf_fenode2[i][counter++] = actgsurf->gnode[2]->slave_node;
        /* Mittelknoten der Flaeche */
        design->dsurf_fenode2[i][counter++] = actgsurf->slave_node[0][0];
        /* Mittelknoten der Linien */
        for (k=0; k<actgsurf->ngline; k++)
          design->dsurf_fenode2[i][counter++] = actgsurf->gline[k]->slave_node[0];

      }
    }

    if (counter != design->ndsurf_fenode[i])
      dserror_args(__FILE__, __LINE__,
          "Irgendetwas stimmt hier nicht !! counter = %d;  ndsurf_fenode[i] =%d",
          counter,design->ndsurf_fenode[i] );

  }




  /* dvol */
  /*========*/


  /* counter number of gnodes on each design vol */
  for (i=0; i<io_dis->ngvol; i++)
  {
    actgvol = &(io_dis->gvol[i]);
    if (actgvol->dvol != NULL)
    {
      id = actgvol->dvol->Id;
      design->ndvol_fenode[id] += 15;
    }
  }



  /* allocate dvol_fenode[i] */
  for (i=0; i<design->ndvol; i++)
  {
    if (design->ndvol_fenode[i] > 0)
      design->dvol_fenode2[i] = (INT*)CCAREALLOC(design->dvol_fenode2[i],
          design->ndvol_fenode[i]*sizeof(INT));
    counter = 0;


    /* fill dvol_fenode */
    for (j=0; j<io_dis->ngvol; j++)
    {
      actgvol = &(io_dis->gvol[j]);
      if (actgvol->dvol != NULL)
      {
        if (actgvol->dvol->Id != i)
          continue;

        counter2 = 0;
        for (k=0; k<4; k++)
          for (l=0; l<3; l++)
            tmp[counter2++] = actgvol->gsurf[k]->gnode[l]->slave_node;

        for (k=0; k<12; k++)
          for (l=0; l<12; l++)
            if (tmp[k] == tmp[l] && k!=l)
              tmp[l] = -1;


        for (k=0; k<12; k++)
          if (tmp[k] != -1)
            design->dvol_fenode2[i][counter++] = tmp[k];


        /* Mittelknoten des Vols */
        design->dvol_fenode2[i][counter++] = actgvol->slave_node[0][0][0];

        /* Mittelknoten der Flaechen */
        for (k=0; k<actgvol->ngsurf; k++)
          design->dvol_fenode2[i][counter++] = actgvol->gsurf[k]->slave_node[0][0];
        /* Mittelknoten der Linien */
        for (k=0; k<actgvol->ngline; k++)
          design->dvol_fenode2[i][counter++] = actgvol->gline[k]->slave_node[0];
      }
    }

    if (counter != design->ndvol_fenode[i])
      dserror_args(__FILE__, __LINE__,
          "Irgendetwas stimmt hier nicht !! counter = %d;  ndvol_fenode[i] =%d",
          counter,design->ndvol_fenode[i] );

  }



  /* create the pointers between gnodes... and the design */
  inpdesign_topology_fe(cal_dis);


#ifdef DEBUG
  dstrc_exit();
#endif


  return;

}  /* END of global_subdivide_tet */














/*-----------------------------------------------------------------------*/
/*!
  \brief Brief description

  Very detailed description what the function is doing why.

  \param b      DOUBLE  (o) explanation

  \return void

  \author mn
  \date   09/05

 */
/*-----------------------------------------------------------------------*/
void global_subdivide_quad(
    FIELD        *actfield
    )

{

  DISCRET       *io_dis;
  DISCRET       *cal_dis;

  INT            subdivide;
  INT            nodeid;
  INT            eleid;

  INT            numnd;
  INT            numele;
  INT            i,j;
  INT            node_counter;
  INT            ele_counter;

  GNODE         *actgnode = NULL;
  GLINE         *actgline = NULL;
  GSURF         *actgsurf = NULL;

  ELEMENT       *actele, *oldele;

  INT            id;

  INT            counter;
  INT            k,l,m,g;

  INT            node_matrix[MAX_DIVIDE+1][MAX_DIVIDE+1];

  INT            pro_total=0, pro_fac=0, p=0, p_counter=0;

#ifdef DEBUG
  dstrc_enter("global_subdivide_quad");
#endif


  io_dis    = &(actfield->dis[0]);
  cal_dis   = &(actfield->dis[1]);

  subdivide = actfield->subdivide;
  nodeid    = genprob.maxnode;
  eleid     = genprob.nele;


  /* calculate number of nodes on second dis */
  numnd = io_dis->ngnode +
          io_dis->ngline * (subdivide - 1) +
          io_dis->ngsurf * (subdivide - 1) * (subdivide - 1);



  /* allocate the nodes to the dis */
  cal_dis->node  = (NODE*)CCACALLOC(numnd,sizeof(NODE));
  cal_dis->numnp = numnd;
  node_counter = 0;



  /* calculate number of elements in second dis */
  numele = io_dis->ngsurf * subdivide * subdivide;


  /* allocate elements in second dis */
  cal_dis->element = (ELEMENT*)CCACALLOC(numele,sizeof(ELEMENT));
  cal_dis->numele  = numele;
  ele_counter = 0;






  /* loop gnodes and create nodes on second dis */
  for (i=0; i<io_dis->ngnode; i++)
  {
    actgnode = &(io_dis->gnode[i]);

    /* set Id of the new node */
    cal_dis->node[node_counter].Id = nodeid;

    /* set coords of the new node */
    for (j=0; j<3; j++)
      cal_dis->node[node_counter].x[j] = actgnode->node->x[j];

    /* store slave_node Id in master gnode */
    actgnode->slave_node = nodeid;

    /* create pointer from slave to master node */
    cal_dis->node[node_counter].slave_node  = NULL;
    cal_dis->node[node_counter].master_node = actgnode->node;

    /* create pointer from master to slave node */
    actgnode->node->slave_node  = &(cal_dis->node[node_counter]);
    actgnode->node->master_node = NULL;


    node_counter++;
    nodeid++;
  }



  /* loop glines and create nodes */
  for (i=0; i<io_dis->ngline; i++)
  {

    DOUBLE x1[3],dx[3];

    actgline = &(io_dis->gline[i]);

    /* coordinates of the beginning of this gline */
    x1[0] = actgline->gnode[0]->node->x[0];
    x1[1] = actgline->gnode[0]->node->x[1];
    x1[2] = actgline->gnode[0]->node->x[2];

    /* lenght of this gline in the three directions divided by subdivide */
    dx[0] = (actgline->gnode[1]->node->x[0] - actgline->gnode[0]->node->x[0])/subdivide;
    dx[1] = (actgline->gnode[1]->node->x[1] - actgline->gnode[0]->node->x[1])/subdivide;
    dx[2] = (actgline->gnode[1]->node->x[2] - actgline->gnode[0]->node->x[2])/subdivide;

    for (j=0; j<subdivide-1; j++)
    {
      cal_dis->node[node_counter].Id = nodeid;

      /* store slave_node Id in gline */
      actgline->slave_node[j]   = nodeid;

      /* set coords of the new node */
      for (k=0; k<3; k++)
        cal_dis->node[node_counter].x[k] =  x1[k] + (j+1)*dx[k];


      node_counter++;
      nodeid++;
    }

  }  /* for (i=0; i<io_dis->ngline; i++) */










  if (par.myrank==0)
  {
    printf("  creating subelements on %3i GSURFS\n",io_dis->ngsurf);

    pro_fac   = io_dis->ngsurf/50 + 1;

    pro_total = io_dis->ngsurf/pro_fac;

    printf("  0");
    for (p=0; p<pro_total; p++)
      printf(" ");
    printf("%3i\n",io_dis->ngsurf);

    printf("  |");
    for (p=0; p<pro_total; p++)
      printf("-");
    printf("|\n");


    printf("  |");
    p_counter = 0;
  }

  /* loop gsurfs and create nodes */
  for (g=0; g<io_dis->ngsurf; g++)
  {
    GNODE  *gnode0=NULL, *gnode1=NULL, *gnode2=NULL, *gnode3=NULL;

    DOUBLE x0[3],x1[3],x2[3],x3[3];

    actgsurf = &(io_dis->gsurf[g]);


    /* specify node 0 */
    gnode0 = actgsurf->gnode[0];
    actgsurf->node_ind[0] = 0;



    /* find the line 0 */
    for (j=0; j<actgsurf->ngline; j++)
    {
      /* this line is in same direction */
      if (gnode0 == actgsurf->gline[j]->gnode[0])
      {
        actgline = actgsurf->gline[j];
        actgsurf->line_ind[0][0] = j;
        actgsurf->line_ind[1][0] = 0;

        gnode1 = actgline->gnode[1];

        /* find the node 1 */
        for (k=0; k<actgsurf->ngnode; k++)
        {
          if (gnode1 == actgsurf->gnode[k])
          {
            actgsurf->node_ind[1] = k;
            break;
          }
        }

        break;
      }

      /* this line is in opposite direction */
      if (gnode0 == actgsurf->gline[j]->gnode[1])
      {
        actgline = actgsurf->gline[j];
        actgsurf->line_ind[0][0] = j;
        actgsurf->line_ind[1][0] = 1;

        gnode1 = actgline->gnode[0];

        /* find the node 1 */
        for (k=0; k<actgsurf->ngnode; k++)
        {
          if (gnode1 == actgsurf->gnode[k])
          {
            actgsurf->node_ind[1] = k;
            break;
          }
        }

        break;
      }
    }  /* for (j=1; j<gsurf->ngline; j++) */




    /* find the line 1 */
    for (j=0; j<actgsurf->ngline; j++)
    {

      if (actgsurf->gline[j] == actgline)
        continue;

      /* this line is in same direction */
      if (gnode1 == actgsurf->gline[j]->gnode[0])
      {
        actgline = actgsurf->gline[j];
        actgsurf->line_ind[0][1] = j;
        actgsurf->line_ind[1][1] = 0;

        gnode2 = actgline->gnode[1];

        /* find the node 2 */
        for (k=0; k<actgsurf->ngnode; k++)
        {
          if (gnode2 == actgsurf->gnode[k])
          {
            actgsurf->node_ind[2] = k;
            break;
          }
        }

        break;
      }

      /* this line is in opposite direction */
      if (gnode1 == actgsurf->gline[j]->gnode[1])
      {
        actgline = actgsurf->gline[j];
        actgsurf->line_ind[0][1] = j;
        actgsurf->line_ind[1][1] = 1;

        gnode2 = actgline->gnode[0];

        /* find the node 2 */
        for (k=0; k<actgsurf->ngnode; k++)
        {
          if (gnode2 == actgsurf->gnode[k])
          {
            actgsurf->node_ind[2] = k;
            break;
          }
        }

        break;
      }
    }  /* for (j=1; j<gsurf->ngline; j++) */


    /* find the line 2 */
    for (j=0; j<actgsurf->ngline; j++)
    {

      if (actgsurf->gline[j] == actgline)
        continue;

      /* this line is in same direction */
      if (gnode2 == actgsurf->gline[j]->gnode[1])
      {
        actgline = actgsurf->gline[j];
        actgsurf->line_ind[0][2] = j;
        actgsurf->line_ind[1][2] = 0;

        gnode3 = actgline->gnode[0];

        /* find the node 3 */
        for (k=0; k<actgsurf->ngnode; k++)
        {
          if (gnode3 == actgsurf->gnode[k])
          {
            actgsurf->node_ind[3] = k;
            break;
          }
        }

        break;
      }

      /* this line is in opposite direction */
      if (gnode2 == actgsurf->gline[j]->gnode[0])
      {
        actgline = actgsurf->gline[j];
        actgsurf->line_ind[0][2] = j;
        actgsurf->line_ind[1][2] = 1;

        gnode3 = actgline->gnode[1];

        /* find the node 3 */
        for (k=0; k<actgsurf->ngnode; k++)
        {
          if (gnode3 == actgsurf->gnode[k])
          {
            actgsurf->node_ind[3] = k;
            break;
          }
        }

        break;
      }
    }  /* for (j=1; j<gsurf->ngline; j++) */





    /* find the line 3 */
    for (j=0; j<actgsurf->ngline; j++)
    {

      if (actgsurf->gline[j] == actgline)
        continue;

      /* this line is in same direction */
      if (gnode3 == actgsurf->gline[j]->gnode[1])
      {
        actgline = actgsurf->gline[j];
        actgsurf->line_ind[0][3] = j;
        actgsurf->line_ind[1][3] = 0;

        actgnode = actgline->gnode[0];

        if (actgnode != gnode0)
          dserror("NODE MIXUP!!\n");

        break;
      }

      /* this line is in opposite direction */
      if (gnode3 == actgsurf->gline[j]->gnode[0])
      {
        actgline = actgsurf->gline[j];
        actgsurf->line_ind[0][3] = j;
        actgsurf->line_ind[1][3] = 1;

        actgnode = actgline->gnode[1];

        if (actgnode != gnode0)
          dserror("NODE MIXUP!!\n");

        break;
      }
    }  /* for (j=1; j<gsurf->ngline; j++) */


    /*
       printf("GSURF %2i:  line_ind: %1i(%1i) %1i(%1i) %1i(%1i) %1i(%1i)\n",
       actgsurf->Id,
       actgsurf->line_ind[0][0],actgsurf->line_ind[1][0],
       actgsurf->line_ind[0][1],actgsurf->line_ind[1][1],
       actgsurf->line_ind[0][2],actgsurf->line_ind[1][2],
       actgsurf->line_ind[0][3],actgsurf->line_ind[1][3]);

       printf("GSURF %2i:  node_ind: %1i %1i %1i %1i\n\n",
       actgsurf->Id,
       actgsurf->node_ind[0],actgsurf->node_ind[1],
       actgsurf->node_ind[2],actgsurf->node_ind[3]);
       */


    /* set the nodal coordinates */
    x0[0] = (actgsurf->gnode[actgsurf->node_ind[0]]->node->x[0])/(subdivide*subdivide);
    x0[1] = (actgsurf->gnode[actgsurf->node_ind[0]]->node->x[1])/(subdivide*subdivide);
    x0[2] = (actgsurf->gnode[actgsurf->node_ind[0]]->node->x[2])/(subdivide*subdivide);

    x1[0] = (actgsurf->gnode[actgsurf->node_ind[1]]->node->x[0])/(subdivide*subdivide);
    x1[1] = (actgsurf->gnode[actgsurf->node_ind[1]]->node->x[1])/(subdivide*subdivide);
    x1[2] = (actgsurf->gnode[actgsurf->node_ind[1]]->node->x[2])/(subdivide*subdivide);

    x2[0] = (actgsurf->gnode[actgsurf->node_ind[2]]->node->x[0])/(subdivide*subdivide);
    x2[1] = (actgsurf->gnode[actgsurf->node_ind[2]]->node->x[1])/(subdivide*subdivide);
    x2[2] = (actgsurf->gnode[actgsurf->node_ind[2]]->node->x[2])/(subdivide*subdivide);

    x3[0] = (actgsurf->gnode[actgsurf->node_ind[3]]->node->x[0])/(subdivide*subdivide);
    x3[1] = (actgsurf->gnode[actgsurf->node_ind[3]]->node->x[1])/(subdivide*subdivide);
    x3[2] = (actgsurf->gnode[actgsurf->node_ind[3]]->node->x[2])/(subdivide*subdivide);



    /* generate the nodes */
    for (j=0; j<subdivide-1; j++)
    {
      for (k=0; k<subdivide-1; k++)
      {
        cal_dis->node[node_counter].Id  = nodeid;

        /* store slave_node Id in gsurf */
        actgsurf->slave_node[j][k] = nodeid;

        /* set coords of the new node */
        for (l=0; l<3; l++)
          cal_dis->node[node_counter].x[l] =
            (x1[l]*(j+1) + x0[l]*(subdivide-j-1)) * (subdivide-k-1) +
            (x2[l]*(j+1) + x3[l]*(subdivide-j-1)) * (k+1);

        node_counter++;
        nodeid++;
      }
    }










    /* check gnode - gline */
    for (j=0; j<4; j++)
    {
      INT links = 0, rechts = 1;
      GNODE *gnode_l=NULL, *gnode_r=NULL;

      INT line_node_ind[4][2];

      line_node_ind[0][0] = 0;
      line_node_ind[0][1] = 1;
      line_node_ind[1][0] = 1;
      line_node_ind[1][1] = 2;
      line_node_ind[2][0] = 3;
      line_node_ind[2][1] = 2;
      line_node_ind[3][0] = 0;
      line_node_ind[3][1] = 3;


      if(actgsurf->line_ind[1][j] == 0)
      {
        links = 0;
        rechts = 1;
      }
      else
      {
        links = 1;
        rechts = 0;
      }

      switch (line_node_ind[j][0])
      {
        case 0: gnode_l = gnode0; break;
        case 1: gnode_l = gnode1; break;
        case 2: gnode_l = gnode2; break;
        case 3: gnode_l = gnode3; break;
      }

      switch (line_node_ind[j][1])
      {
        case 0: gnode_r = gnode0; break;
        case 1: gnode_r = gnode1; break;
        case 2: gnode_r = gnode2; break;
        case 3: gnode_r = gnode3; break;
      }

      if (actgsurf->gline[actgsurf->line_ind[0][j]]->gnode[links] != gnode_l)
        dserror_args(__FILE__, __LINE__, "Gline %1i, links != gnode %1i\n",
            j,line_node_ind[j][0]);

      if (actgsurf->gline[actgsurf->line_ind[0][j]]->gnode[rechts] != gnode_r)
        dserror_args(__FILE__, __LINE__, "Gline %1i, rechts != gnode %1i\n",
            j,line_node_ind[j][1]);
    }





    /* check gnodes */
    for (j=0; j<4; j++)
    {
      for (k=0;k <4; k++)
      {
        GNODE *gnode_l=NULL, *gnode_r=NULL;

        switch (j)
        {
          case 0: gnode_l = gnode0; break;
          case 1: gnode_l = gnode1; break;
          case 2: gnode_l = gnode2; break;
          case 3: gnode_l = gnode3; break;
        }
        switch (k)
        {
          case 0: gnode_r = gnode0; break;
          case 1: gnode_r = gnode1; break;
          case 2: gnode_r = gnode2; break;
          case 3: gnode_r = gnode3; break;
        }

        if (gnode_r == gnode_l && j!=k)
        {
          dserror_args(__FILE__, __LINE__, "GNODE %1i == GNODE %1i !!\n",j,k);
        }
      }
    }


    /*
       printf("GVOL %2i:  line_ind: ", actgsurf->Id);
       for (j=0; j<actgsurf->ngline; j++)
       printf("%1i(%1i) ",
           actgsurf->actgsurf->line_ind[0][j],actgsurf->actgsurf->line_ind[1][j]);
       printf("\n\n");
       */




    /* write corner nodes into node_matrix */
    node_matrix[        0][        0] = gnode0->slave_node;
    node_matrix[subdivide][        0] = gnode1->slave_node;
    node_matrix[subdivide][subdivide] = gnode2->slave_node;
    node_matrix[        0][subdivide] = gnode3->slave_node;




    /* write line nodeid into node_matrix */
    for (j=0;j<subdivide-1; j++)
    {
      /* line 0 */
      if (actgsurf->line_ind[1][0] == 0)
        node_matrix[j+1][0] =
          actgsurf->gline[actgsurf->line_ind[0][0]]->slave_node[j];
      else
        node_matrix[subdivide-1-j][0] =
          actgsurf->gline[actgsurf->line_ind[0][0]]->slave_node[j];

      /* line 1 */
      if (actgsurf->line_ind[1][1] == 0)
        node_matrix[subdivide][j+1] =
          actgsurf->gline[actgsurf->line_ind[0][1]]->slave_node[j];
      else
        node_matrix[subdivide][subdivide-1-j] =
          actgsurf->gline[actgsurf->line_ind[0][1]]->slave_node[j];

      /* line 2 */
      if (actgsurf->line_ind[1][2] == 0)
        node_matrix[j+1][subdivide] =
          actgsurf->gline[actgsurf->line_ind[0][2]]->slave_node[j];
      else
        node_matrix[subdivide-1-j][subdivide] =
          actgsurf->gline[actgsurf->line_ind[0][2]]->slave_node[j];

      /* line 3 */
      if (actgsurf->line_ind[1][3] == 0)
        node_matrix[0][j+1] =
          actgsurf->gline[actgsurf->line_ind[0][3]]->slave_node[j];
      else
        node_matrix[0][subdivide-1-j] =
          actgsurf->gline[actgsurf->line_ind[0][3]]->slave_node[j];

    }



    /* write the surface nodes in node_matrix */
    for (i=0; i<subdivide-1;i++)
    {
      for (j=0; j<subdivide-1;j++)
      {

        GSURF *gsurf0 =  actgsurf;

        /* SURFACE 0 */
        if (gsurf0->gnode[1] == gnode1)
          node_matrix[i+1][j+1] = gsurf0->slave_node[i][j];

        else if (gsurf0->gnode[1] == gnode3)
          node_matrix[j+1][i+1] = gsurf0->slave_node[i][j];

        else
          dserror("Surface orientation error!!!");

      }  /* for (j=0; j<subdivide-1;j++) */

    }  /* for (i=0; i<subdivide-1;i++) */






#if 0
    /* print the node_matrix */
    /* only for debugging */
    for (k=subdivide; k>=0;k--)
    {
      for (l=subdivide; l>=0;l--)
      {
        for (j=0; j<=subdivide;j++)
        {
          printf(" %5i ",node_matrix[j][k][l]);
        }
        printf("\n");
      }
      printf("\n");
    }
#endif




    /* ==========================================*
     *                                           *
     * generate the new elements                 *
     *                                           *
     * ==========================================*/


    oldele = actgsurf->element;
    oldele->master_ele = NULL;
    oldele->slave_ele  = NULL;


    if (oldele->distyp != quad4)
      dserror("Subdivision only possible for uniform discretizations!!");


      for (k=0; k<subdivide; k++)
      {

        for (j=0; j<subdivide; j++)
        {

          INT nhyb = 0;

          actele = &(cal_dis->element[ele_counter]);


          /* set general properties */
          actele->Id      = eleid;
          actele->numnp   = 4;
          actele->distyp  = quad4;

          actele->master_ele = oldele;
          actele->slave_ele  = NULL;


          /* allocate lm */
          actele->lm = (INT*)CCACALLOC(actele->numnp,sizeof(INT));
          if (actele->lm==NULL) dserror("Allocation of lm in ELEMENT failed\n");


          /* fill lm */
          if (subdivide == 1)
            for (m=0; m<actele->numnp; m++)
            {
              actele->lm[m] =  oldele->node[m]->gnode->slave_node;
            }
          else
          {
            actele->lm[0] = node_matrix[j  ][k  ];
            actele->lm[1] = node_matrix[j+1][k  ];
            actele->lm[2] = node_matrix[j+1][k+1];
            actele->lm[3] = node_matrix[j  ][k+1];
          }



          switch (oldele->eltyp)
          {

#ifdef D_ALE
            case el_ale2:

              actele->eltyp   = el_ale2;

              /* allocate the element */
              actele->e.ale2 = (ALE2*)CCACALLOC(1,sizeof(ALE2));
              if (actele->e.ale2==NULL) dserror("Allocation of element ALE2 failed");

              /* copy element properties */
              actele->mat                = oldele->mat;
              actele->e.ale2->nGP[0]     = oldele->e.ale2->nGP[0];
              actele->e.ale2->nGP[1]     = oldele->e.ale2->nGP[1];
              actele->e.ale2->jacobi     = oldele->e.ale2->jacobi;
              break;
#endif


#ifdef D_FLUID2
            case el_fluid2:

              actele->eltyp   = el_fluid2;

              /* allocate the element */
              actele->e.f2 = (FLUID2*)CCACALLOC(1,sizeof(FLUID2));
              if (actele->e.f2==NULL) dserror("Allocation of element FLUID2 failed");

              /* copy element properties */
              actele->mat              = oldele->mat;
              actele->e.f2->nGP[0]     = oldele->e.f2->nGP[0];
              actele->e.f2->nGP[1]     = oldele->e.f2->nGP[1];
              actele->e.f2->is_ale     = oldele->e.f2->is_ale;
              actele->e.f2->create_ale = oldele->e.f2->create_ale;
              actele->e.f2->fs_on      = 0;
              break;
#endif

#ifdef D_SHELL8
            case el_shell8:

              actele->eltyp   = el_shell8;

              actele->e.s8 = (SHELL8*)CCACALLOC(1,sizeof(SHELL8));
              if (actele->e.s8==NULL) dserror("Allocation of element SHELL8 failed");

              /* copy element properties */
              amdef("intforce",&(actele->e.s8->intforce), NUMDOF_SHELL8*actele->numnp,1,"DV");
              actele->mat               = oldele->mat;
              actele->e.s8->thick       = oldele->e.s8->thick;
              actele->e.s8->nGP[0]      = oldele->e.s8->nGP[0];
              actele->e.s8->nGP[1]      = oldele->e.s8->nGP[1];
              actele->e.s8->nGP[2]      = oldele->e.s8->nGP[2];
              actele->e.s8->nGP_tri     = oldele->e.s8->nGP_tri;
              actele->e.s8->forcetyp    = oldele->e.s8->forcetyp;
              actele->e.s8->eas[0]      = oldele->e.s8->eas[0];
              actele->e.s8->eas[1]      = oldele->e.s8->eas[1];
              actele->e.s8->eas[2]      = oldele->e.s8->eas[2];
              actele->e.s8->eas[3]      = oldele->e.s8->eas[3];
              actele->e.s8->eas[4]      = oldele->e.s8->eas[4];
              actele->e.s8->ans         = oldele->e.s8->ans;
              actele->e.s8->sdc         = oldele->e.s8->sdc;


              /*--------------------- count nhyb and allocate storage for eas strains */
              for (i=0; i<5; i++) nhyb+=actele->e.s8->eas[i];
              actele->e.s8->nhyb=nhyb;
              if (nhyb>0)
              {
                amdef("alfa",&(actele->e.s8->alfa),1,nhyb,"DA");
                amzero(&(actele->e.s8->alfa));
                amdef("oalfa",&(actele->e.s8->oldalfa),1,nhyb,"DA");
                amzero(&(actele->e.s8->oldalfa));

                amdef("Dtildinv",&(actele->e.s8->Dtildinv),nhyb,nhyb,"DA");
                amzero(&(actele->e.s8->Dtildinv));

                amdef("Lt",&(actele->e.s8->Lt),nhyb,actele->numnp*NUMDOF_SHELL8,"DA");
                amzero(&(actele->e.s8->Lt));

                amdef("Rtilde",&(actele->e.s8->Rtilde),nhyb,1,"DV");
                amzero(&(actele->e.s8->Rtilde));
              }
              break;
#endif



#ifdef D_WALL1
            case el_wall1:

              actele->eltyp   = el_wall1;

              actele->e.w1 = (WALL1*)CCACALLOC(1,sizeof(WALL1));
              if (actele->e.w1==NULL) dserror("Allocation of element WALL1 failed");

              /* copy element properties */
              actele->mat               = oldele->mat;
              actele->e.w1->thick       = oldele->e.w1->thick;
              actele->e.w1->nGP[0]      = oldele->e.w1->nGP[0];
              actele->e.w1->nGP[1]      = oldele->e.w1->nGP[1];
              actele->e.w1->wtype       = oldele->e.w1->wtype;
              actele->e.w1->kintype     = oldele->e.w1->kintype;
              actele->e.w1->modeltype   = oldele->e.w1->modeltype;
              actele->e.w1->stresstyp   = oldele->e.w1->stresstyp;
              break;
#endif

            case el_none:
            default:
              dserror("Unknown element type!!");
              break;
          }



#if 0
          /* print the location vector */
          /* only for debugging */
          for (m=0; m<actele->numnp;m++)
          {
            printf(" %5i ",actele->lm[m]);
          }
          printf("\n");
#endif


          ele_counter++;
          eleid++;
        }
      }


      if (par.myrank==0)
      {
        p_counter++;
        if (p_counter == pro_fac)
        {
          printf("*");
          p_counter = 0;
        }
      }
  }  /* for (g=0; g<io_dis->ngsurf; g++) */

  if (par.myrank==0)
    printf("|\n");



  /* check the number of created nodes */
  if (node_counter != numnd )
    dserror_args(__FILE__, __LINE__,
        "Wrong number of nodes created!! numnd = %d;  node_counter = %d",
        numnd,node_counter);


  genprob.maxnode = nodeid + 1;
  genprob.nnode += numnd;




  /* check the number of created elelments */
  if (ele_counter != numele )
    dserror_args(__FILE__, __LINE__,
        "Wrong number of elements created!! numnd = %d;  ele_counter = %d",
        numnd,ele_counter);


  genprob.nele += numele;



  /* build pointers between nodes and elements */
  inp_topology(cal_dis);



  /* create all gnodes... and build the pointers between them all */
  inp_detailed_topology(cal_dis);



#if 0
  design->ndnode_fenode = (INT*)CCACALLOC(design->ndnode,sizeof(INT));
  design->dnode_fenode2 = (INT**)CCACALLOC(design->ndnode,sizeof(INT*));

  design->ndline_fenode = (INT*)CCACALLOC(design->ndline,sizeof(INT));;
  design->dline_fenode2 = (INT**)CCACALLOC(design->ndline,sizeof(INT*));;

  design->ndsurf_fenode = (INT*)CCACALLOC(design->ndsurf,sizeof(INT));;
  design->dsurf_fenode2 = (INT**)CCACALLOC(design->ndsurf,sizeof(INT*));;

  design->ndvol_fenode  = (INT*)CCACALLOC(design->ndvol,sizeof(INT));;
  design->dvol_fenode2  = (INT**)CCACALLOC(design->ndvol,sizeof(INT*));;
#endif



  for (i=0; i<design->ndnode; i++)
    design->ndnode_fenode[i] = 0;

  for (i=0; i<design->ndline; i++)
    design->ndline_fenode[i] = 0;

  for (i=0; i<design->ndsurf; i++)
    design->ndsurf_fenode[i] = 0;

  for (i=0; i<design->ndvol; i++)
    design->ndvol_fenode[i] = 0;


  /* dnodes */
  /*========*/


  /* count number of gnodes on each design node */
  for (i=0; i<io_dis->ngnode; i++)
  {
    actgnode = &(io_dis->gnode[i]);
    if (actgnode->ondesigntyp == ondnode)
    {
      id = actgnode->d.dnode->Id;
      design->ndnode_fenode[id]+=1;
    }
  }


  /* allocate dnode_fenode[i] */
  for (i=0; i<design->ndnode; i++)
  {
    if (design->ndnode_fenode[i] > 0)
      design->dnode_fenode2[i] = (INT*)CCAREALLOC(design->dnode_fenode2[i],
          design->ndnode_fenode[i]*sizeof(INT));

    if (design->ndnode_fenode[i] > 1 )
      dserror("Only one fe_node possible on each dnode !!");

  }


  /* fill dnode_fenode */
  for (i=0; i<io_dis->ngnode; i++)
  {
    actgnode = &(io_dis->gnode[i]);
    if (actgnode->ondesigntyp == ondnode)
    {
      id = actgnode->d.dnode->Id;
      design->dnode_fenode2[id][0] = actgnode->slave_node;
    }
  }



  /* dlines */
  /*========*/


  /* count number of gnodes on each design line */
  for (i=0; i<io_dis->ngline; i++)
  {
    actgline = &(io_dis->gline[i]);
    if (actgline->dline != NULL)
    {
      id = actgline->dline->Id;
      design->ndline_fenode[id] += (subdivide + 1);
    }
  }



  /* allocate dline_fenode[i] */
  for (i=0; i<design->ndline; i++)
  {
    if (design->ndline_fenode[i] > 0)
      design->dline_fenode2[i] = (INT*)CCAREALLOC(design->dline_fenode2[i],
          design->ndline_fenode[i]*sizeof(INT));
    counter = 0;


    /* fill dline_fenode */
    for (j=0; j<io_dis->ngline; j++)
    {
      actgline = &(io_dis->gline[j]);
      if (actgline->dline != NULL)
      {
        if (actgline->dline->Id != i)
          continue;
        design->dline_fenode2[i][counter++] = actgline->gnode[0]->slave_node;
        design->dline_fenode2[i][counter++] = actgline->gnode[1]->slave_node;
        /* Mittelknoten */
        for (k=0; k<subdivide-1; k++)
          design->dline_fenode2[i][counter++] = actgline->slave_node[k];
      }
    }

    if (counter != design->ndline_fenode[i])
      dserror_args(__FILE__, __LINE__,
          "Irgendetwas stimmt hier nicht !! counter = %d;  ndline_fenode[i] =%d",
          counter,design->ndline_fenode[i] );

  }




  /* dsurf */
  /*========*/


  /* count number of gnodes on each design surf */
  for (i=0; i<io_dis->ngsurf; i++)
  {
    actgsurf = &(io_dis->gsurf[i]);
    if (actgsurf->dsurf != NULL)
    {
      id = actgsurf->dsurf->Id;
      design->ndsurf_fenode[id] += (subdivide + 1) * (subdivide + 1);
    }
  }



  /* allocate dsurf_fenode[i] */
  for (i=0; i<design->ndsurf; i++)
  {
    if (design->ndsurf_fenode[i] > 0)
      design->dsurf_fenode2[i] = (INT*)CCAREALLOC(design->dsurf_fenode2[i],
          design->ndsurf_fenode[i]*sizeof(INT));
    counter = 0;


    /* fill dsurf_fenode */
    for (j=0; j<io_dis->ngsurf; j++)
    {
      actgsurf = &(io_dis->gsurf[j]);
      if (actgsurf->dsurf != NULL)
      {
        if (actgsurf->dsurf->Id != i)
          continue;
        design->dsurf_fenode2[i][counter++] = actgsurf->gnode[0]->slave_node;
        design->dsurf_fenode2[i][counter++] = actgsurf->gnode[1]->slave_node;
        design->dsurf_fenode2[i][counter++] = actgsurf->gnode[2]->slave_node;
        design->dsurf_fenode2[i][counter++] = actgsurf->gnode[3]->slave_node;
        /* Mittelknoten der Flaeche */
        for (k=0; k<subdivide-1; k++)
          for (l=0; l<subdivide-1; l++)
            design->dsurf_fenode2[i][counter++] = actgsurf->slave_node[k][l];
        /* Mittelknoten der Linien */
        for (k=0; k<actgsurf->ngline; k++)
          for (l=0; l<subdivide-1; l++)
            design->dsurf_fenode2[i][counter++] = actgsurf->gline[k]->slave_node[l];

      }
    }

    if (counter != design->ndsurf_fenode[i])
      dserror_args(__FILE__, __LINE__,
          "Irgendetwas stimmt hier nicht !! counter = %d;  ndsurf_fenode[i] =%d",
          counter,design->ndsurf_fenode[i] );

  }






  /* create the pointers between gnodes... and the design */
  inpdesign_topology_fe(cal_dis);


#ifdef DEBUG
  dstrc_exit();
#endif


  return;

}  /* END of global_subdivide_quad */







/*-----------------------------------------------------------------------*/
/*!
  \brief Brief description

  Very detailed description what the function is doing why.

  \param b      DOUBLE  (o) explanation

  \return void

  \author mn
  \date   09/05

 */
/*-----------------------------------------------------------------------*/
void global_subdivide_tri(
    FIELD        *actfield
    )

{

  DISCRET       *io_dis;
  DISCRET       *cal_dis;

  INT            subdivide;
  INT            nodeid;
  INT            eleid;

  INT            numnd;
  INT            numele;
  INT            i,j;
  INT            node_counter;
  INT            ele_counter;

  GNODE         *actgnode;
  GLINE         *actgline;
  GSURF         *actgsurf;

  ELEMENT       *actele, *oldele;

  INT            id;

  INT            counter;
  INT            k,l,g;

  INT            node_matrix[MAX_DIVIDE+1][MAX_DIVIDE+1];

  INT            pro_total=0, pro_fac=0, p=0, p_counter=0;

#ifdef DEBUG
  dstrc_enter("global_subdivide_tri");
#endif


  io_dis    = &(actfield->dis[0]);
  cal_dis   = &(actfield->dis[1]);

  subdivide = actfield->subdivide;
  nodeid    = genprob.maxnode;
  eleid     = genprob.nele;


  if (subdivide != 2)
  {
    dserror("Subdivision of tris only possible for subdivide = 2!!!");
  }


  /* calculate number of nodes on second dis */
  numnd = io_dis->ngnode + io_dis->ngline + io_dis->ngsurf ;



  /* allocate the nodes to the dis */
  cal_dis->node  = (NODE*)CCACALLOC(numnd,sizeof(NODE));
  cal_dis->numnp = numnd;
  node_counter = 0;



  /* calculate number of elements in second dis */
  numele = io_dis->ngsurf * 3;


  /* allocate elements in second dis */
  cal_dis->element = (ELEMENT*)CCACALLOC(numele,sizeof(ELEMENT));
  cal_dis->numele  = numele;
  ele_counter = 0;






  /* loop gnodes and create nodes on second dis */
  for (i=0; i<io_dis->ngnode; i++)
  {
    actgnode = &(io_dis->gnode[i]);

    /* set Id of the new node */
    cal_dis->node[node_counter].Id = nodeid;

    /* set coords of the new node */
    for (j=0; j<3; j++)
      cal_dis->node[node_counter].x[j] = actgnode->node->x[j];

    /* store slave_node Id in master gnode */
    actgnode->slave_node = nodeid;

    /* create pointer from slave to master node */
    cal_dis->node[node_counter].slave_node  = NULL;
    cal_dis->node[node_counter].master_node = actgnode->node;

    /* create pointer from master to slave node */
    actgnode->node->slave_node  = &(cal_dis->node[node_counter]);
    actgnode->node->master_node = NULL;


    node_counter++;
    nodeid++;
  }



  /* loop glines and create nodes */
  for (i=0; i<io_dis->ngline; i++)
  {

    DOUBLE x1[3],dx[3];

    actgline = &(io_dis->gline[i]);

    /* coordinates of the beginning of this gline */
    x1[0] = actgline->gnode[0]->node->x[0];
    x1[1] = actgline->gnode[0]->node->x[1];
    x1[2] = actgline->gnode[0]->node->x[2];

    /* lenght of this gline in the three directions divided by subdivide */
    dx[0] = (actgline->gnode[1]->node->x[0] - actgline->gnode[0]->node->x[0])/2;
    dx[1] = (actgline->gnode[1]->node->x[1] - actgline->gnode[0]->node->x[1])/2;
    dx[2] = (actgline->gnode[1]->node->x[2] - actgline->gnode[0]->node->x[2])/2;

    cal_dis->node[node_counter].Id = nodeid;

    /* store slave_node Id in gline */
    actgline->slave_node[0]   = nodeid;

    /* set coords of the new node */
    for (k=0; k<3; k++)
      cal_dis->node[node_counter].x[k] =  x1[k] + dx[k];

    node_counter++;
    nodeid++;

  }  /* for (i=0; i<io_dis->ngline; i++) */












  if (par.myrank==0)
  {
    printf("  creating subelements on %3i GSURFS\n",io_dis->ngsurf);

    pro_fac   = io_dis->ngsurf/50 + 1;

    pro_total = io_dis->ngsurf/pro_fac;

    printf("  0");
    for (p=0; p<pro_total; p++)
      printf(" ");
    printf("%3i\n",io_dis->ngsurf);

    printf("  |");
    for (p=0; p<pro_total; p++)
      printf("-");
    printf("|\n");


    printf("  |");
    p_counter = 0;
  }

  /* loop gsurfs and create nodes */
  for (g=0; g<io_dis->ngsurf; g++)
  {
    GNODE  *gnode0, *gnode1, *gnode2;

    actgsurf = &(io_dis->gsurf[g]);

    /* specify the nodes */
    gnode0 = actgsurf->gnode[0];
    gnode1 = actgsurf->gnode[1];
    gnode2 = actgsurf->gnode[2];


    /* generate the node */
    cal_dis->node[node_counter].Id  = nodeid;

    /* store slave_node Id in gsurf */
    actgsurf->slave_node[0][0] = nodeid;

    /* set coords of the new node */
    for (l=0; l<3; l++)
      cal_dis->node[node_counter].x[l] =
        (gnode0->node->x[l] + gnode1->node->x[l] + gnode2->node->x[l]) / 3;


    node_counter++;
    nodeid++;



    /* identify the lines */


    /* LINE 0: 0-1 */
    /*==============*/
    for (j=0; j<actgsurf->ngline; j++)
    {
      /* this line is in same direction */
      if (gnode0 == actgsurf->gline[j]->gnode[0] &&
          gnode1 == actgsurf->gline[j]->gnode[1])
      {
        actgline = actgsurf->gline[j];
        actgsurf->line_ind[0][0] = j;
        actgsurf->line_ind[1][0] = 0;
        break;
      }

      /* this line is in opposite direction */
      if (gnode1 == actgsurf->gline[j]->gnode[0] &&
          gnode0 == actgsurf->gline[j]->gnode[1])
      {
        actgline = actgsurf->gline[j];
        actgsurf->line_ind[0][0] = j;
        actgsurf->line_ind[1][0] = 1;
        break;
      }
    }  /* for (j=0; j<gsurf->ngline; j++) */




    /* LINE 1: 1-2 */
    /*==============*/
    for (j=0; j<actgsurf->ngline; j++)
    {
      /* this line is in same direction */
      if (gnode1 == actgsurf->gline[j]->gnode[0] &&
          gnode2 == actgsurf->gline[j]->gnode[1])
      {
        actgline = actgsurf->gline[j];
        actgsurf->line_ind[0][1] = j;
        actgsurf->line_ind[1][1] = 0;
        break;
      }

      /* this line is in opposite direction */
      if (gnode2 == actgsurf->gline[j]->gnode[0] &&
          gnode1 == actgsurf->gline[j]->gnode[1])
      {
        actgline = actgsurf->gline[j];
        actgsurf->line_ind[0][1] = j;
        actgsurf->line_ind[1][1] = 1;
        break;
      }
    }  /* for (j=0; j<gsurf->ngline; j++) */



    /* LINE 2: 2-0 */
    /*==============*/
    for (j=0; j<actgsurf->ngline; j++)
    {
      /* this line is in same direction */
      if (gnode2 == actgsurf->gline[j]->gnode[0] &&
          gnode0 == actgsurf->gline[j]->gnode[1])
      {
        actgline = actgsurf->gline[j];
        actgsurf->line_ind[0][2] = j;
        actgsurf->line_ind[1][2] = 0;
        break;
      }

      /* this line is in opposite direction */
      if (gnode0 == actgsurf->gline[j]->gnode[0] &&
          gnode2 == actgsurf->gline[j]->gnode[1])
      {
        actgline = actgsurf->gline[j];
        actgsurf->line_ind[0][2] = j;
        actgsurf->line_ind[1][2] = 1;
        break;
      }
    }  /* for (j=0; j<gsurf->ngline; j++) */




#if 1
    /* print the node_matrix */
    /* only for debugging */
      for (l=2; l>=0;l--)
      {
        for (j=0; j<=2;j++)
        {
          node_matrix[j][l] = -17;
        }
      }
#endif



    /* write corner nodes into node_matrix */
    node_matrix[0][0] = gnode0->slave_node;
    node_matrix[2][0] = gnode1->slave_node;
    node_matrix[0][2] = gnode2->slave_node;




    /* write line nodeid into node_matrix */
    /* line 0 */
    node_matrix[1][0] =
      actgsurf->gline[actgsurf->line_ind[0][0]]->slave_node[0];

    /* line 1 */
    node_matrix[2][1] =
      actgsurf->gline[actgsurf->line_ind[0][1]]->slave_node[0];

    /* line 2 */
    node_matrix[0][1] =
      actgsurf->gline[actgsurf->line_ind[0][2]]->slave_node[0];





    /* write the surface node in node_matrix */
    node_matrix[1][1] = actgsurf->slave_node[0][0];








#if 1
    /* print the node_matrix */
    /* only for debugging */
      for (l=2; l>=0;l--)
      {
        for (j=0; j<=2;j++)
        {
          printf(" %5i ",node_matrix[j][l]);
        }
        printf("\n");
      }
#endif




    /* ==========================================*
     *                                           *
     * generate the new elements                 *
     *                                           *
     * ==========================================*/


      oldele = actgsurf->element;
      oldele->master_ele = NULL;
      oldele->slave_ele  = NULL;


      if (oldele->distyp != tri3)
        dserror("Subdivision only possible for uniform discretizations!!");


      for (l=0; l<3; l++)
      {
        INT nhyb=0;

        actele = &(cal_dis->element[ele_counter]);


        /* set general properties */
        actele->Id      = eleid;
        actele->numnp   = 4;
        actele->distyp  = quad4;

        actele->master_ele = oldele;
        actele->slave_ele  = NULL;


        /* allocate lm */
        actele->lm = (INT*)CCACALLOC(actele->numnp,sizeof(INT));
        if (actele->lm==NULL) dserror("Allocation of lm in ELEMENT failed\n");


        /* fill lm */
        switch (l)
        {
          case 0:
            actele->lm[0] = node_matrix[0][0];
            actele->lm[1] = node_matrix[1][0];
            actele->lm[2] = node_matrix[1][1];
            actele->lm[3] = node_matrix[0][1];
            break;

          case 1:
            actele->lm[0] = node_matrix[2][0];
            actele->lm[1] = node_matrix[2][1];
            actele->lm[2] = node_matrix[1][1];
            actele->lm[3] = node_matrix[1][0];
            break;

          case 2:
            actele->lm[0] = node_matrix[0][2];
            actele->lm[1] = node_matrix[0][1];
            actele->lm[2] = node_matrix[1][1];
            actele->lm[3] = node_matrix[2][1];
            break;


          default:
            dserror("This is NOT possible!!!");
            break;
        }


        switch (oldele->eltyp)
        {

#ifdef D_ALE
          case el_ale2:

            actele->eltyp   = el_ale2;

            /* allocate the element */
            actele->e.ale2 = (ALE2*)CCACALLOC(1,sizeof(ALE2));
            if (actele->e.ale2==NULL) dserror("Allocation of element ALE2 failed");

            /* copy element properties */
            actele->mat                = oldele->mat;
            actele->e.ale2->jacobi     = oldele->e.ale2->jacobi;

            switch (oldele->e.ale2->nGP[0])
            {
              case 1:
                actele->e.ale2->nGP[0]     = 1;
                actele->e.ale2->nGP[1]     = 1;
                actele->e.ale2->nGP[2]     = 1;
                break;

              case 3:
                actele->e.ale2->nGP[0]     = 2;
                actele->e.ale2->nGP[1]     = 2;
                actele->e.ale2->nGP[2]     = 2;
                break;

              default:
                dserror("Unknown number of Gauss Points for tri element!!");
                break;
            }
            break;
#endif


#ifdef D_FLUID2
          case el_fluid2:

            actele->eltyp   = el_fluid2;

            /* allocate the element */
            actele->e.f2 = (FLUID2*)CCACALLOC(1,sizeof(FLUID2));
            if (actele->e.f2==NULL) dserror("Allocation of element FLUID2 failed");

            /* copy element properties */
            actele->mat              = oldele->mat;
            actele->e.f2->is_ale     = oldele->e.f2->is_ale;
            actele->e.f2->create_ale = oldele->e.f2->create_ale;
            actele->e.f2->fs_on      = 0;

            switch (oldele->e.f2->nGP[0])
            {
              case 1:
                actele->e.f2->nGP[0]     = 1;
                actele->e.f2->nGP[1]     = 1;
                actele->e.f2->nGP[2]     = 1;
                break;

              case 3:
                actele->e.f2->nGP[0]     = 2;
                actele->e.f2->nGP[1]     = 2;
                actele->e.f2->nGP[2]     = 2;
                break;

              default:
                dserror("Unknown number of Gauss Points for tri element!!");
                break;
            }
            break;
#endif

#ifdef D_SHELL8
            case el_shell8:

              actele->eltyp   = el_shell8;

              actele->e.s8 = (SHELL8*)CCACALLOC(1,sizeof(SHELL8));
              if (actele->e.s8==NULL) dserror("Allocation of element SHELL8 failed");

              /* copy element properties */
              amdef("intforce",&(actele->e.s8->intforce), NUMDOF_SHELL8*actele->numnp,1,"DV");
              actele->mat               = oldele->mat;
              actele->e.s8->thick       = oldele->e.s8->thick;
              actele->e.s8->nGP[2]      = oldele->e.s8->nGP[2];
              actele->e.s8->nGP_tri     = oldele->e.s8->nGP_tri;
              actele->e.s8->forcetyp    = oldele->e.s8->forcetyp;
              actele->e.s8->eas[0]      = oldele->e.s8->eas[0];
              actele->e.s8->eas[1]      = oldele->e.s8->eas[1];
              actele->e.s8->eas[2]      = oldele->e.s8->eas[2];
              actele->e.s8->eas[3]      = oldele->e.s8->eas[3];
              actele->e.s8->eas[4]      = oldele->e.s8->eas[4];
              actele->e.s8->ans         = oldele->e.s8->ans;
              actele->e.s8->sdc         = oldele->e.s8->sdc;

              switch (oldele->e.s8->nGP[0])
              {
                case 1:
                  actele->e.s8->nGP[0]     = 1;
                  actele->e.s8->nGP[1]     = 1;
                  break;

                case 3:
                  actele->e.s8->nGP[0]     = 2;
                  actele->e.s8->nGP[1]     = 2;
                  break;

                default:
                  dserror("Unknown number of Gauss Points for tri element!!");
                  break;
              }

              /*--------------------- count nhyb and allocate storage for eas strains */
              for (i=0; i<5; i++) nhyb+=actele->e.s8->eas[i];
              actele->e.s8->nhyb=nhyb;
              if (nhyb>0)
              {
                amdef("alfa",&(actele->e.s8->alfa),1,nhyb,"DA");
                amzero(&(actele->e.s8->alfa));
                amdef("oalfa",&(actele->e.s8->oldalfa),1,nhyb,"DA");
                amzero(&(actele->e.s8->oldalfa));

                amdef("Dtildinv",&(actele->e.s8->Dtildinv),nhyb,nhyb,"DA");
                amzero(&(actele->e.s8->Dtildinv));

                amdef("Lt",&(actele->e.s8->Lt),nhyb,actele->numnp*NUMDOF_SHELL8,"DA");
                amzero(&(actele->e.s8->Lt));

                amdef("Rtilde",&(actele->e.s8->Rtilde),nhyb,1,"DV");
                amzero(&(actele->e.s8->Rtilde));
              }
              break;
#endif

          case el_none:
          default:
            dserror("Unknown element type!!");
            break;
        }



#if 0
        /* print the location vector */
        /* only for debugging */
        for (m=0; m<actele->numnp;m++)
        {
          printf(" %5i ",actele->lm[m]);
        }
        printf("\n");
#endif


        ele_counter++;
        eleid++;
      }

      if (par.myrank==0)
      {
        p_counter++;
        if (p_counter == pro_fac)
        {
          printf("*");
          p_counter = 0;
        }
      }

  }  /* for (g=0; g<io_dis->ngsurf; g++) */

  if (par.myrank==0)
    printf("|\n");



  /* check the number of created nodes */
  if (node_counter != numnd )
    dserror_args(__FILE__, __LINE__,
        "Wrong number of nodes created!! numnd = %d;  node_counter = %d",
        numnd,node_counter);


  genprob.maxnode = nodeid + 1;
  genprob.nnode += numnd;




  /* check the number of created elelments */
  if (ele_counter != numele )
    dserror_args(__FILE__, __LINE__,
        "Wrong number of elements created!! numnd = %d;  ele_counter = %d",
        numnd,ele_counter);


  genprob.nele += numele;



  /* build pointers between nodes and elements */
  inp_topology(cal_dis);



  /* create all gnodes... and build the pointers between them all */
  inp_detailed_topology(cal_dis);



  design->ndnode_fenode = (INT*)CCACALLOC(design->ndnode,sizeof(INT));
  design->dnode_fenode2 = (INT**)CCACALLOC(design->ndnode,sizeof(INT*));

  design->ndline_fenode = (INT*)CCACALLOC(design->ndline,sizeof(INT));;
  design->dline_fenode2 = (INT**)CCACALLOC(design->ndline,sizeof(INT*));;

  design->ndsurf_fenode = (INT*)CCACALLOC(design->ndsurf,sizeof(INT));;
  design->dsurf_fenode2 = (INT**)CCACALLOC(design->ndsurf,sizeof(INT*));;

  design->ndvol_fenode  = (INT*)CCACALLOC(design->ndvol,sizeof(INT));;
  design->dvol_fenode2  = (INT**)CCACALLOC(design->ndvol,sizeof(INT*));;



  for (i=0; i<design->ndnode; i++)
    design->ndnode_fenode[i] = 0;

  for (i=0; i<design->ndline; i++)
    design->ndline_fenode[i] = 0;

  for (i=0; i<design->ndsurf; i++)
    design->ndsurf_fenode[i] = 0;

  for (i=0; i<design->ndvol; i++)
    design->ndvol_fenode[i] = 0;


  /* dnodes */
  /*========*/


  /* count number of gnodes on each design node */
  for (i=0; i<io_dis->ngnode; i++)
  {
    actgnode = &(io_dis->gnode[i]);
    if (actgnode->ondesigntyp == ondnode)
    {
      id = actgnode->d.dnode->Id;
      design->ndnode_fenode[id]+=1;
    }
  }


  /* allocate dnode_fenode[i] */
  for (i=0; i<design->ndnode; i++)
  {
    if (design->ndnode_fenode[i] > 0)
      design->dnode_fenode2[i] = (INT*)CCAREALLOC(design->dnode_fenode2[i],
          design->ndnode_fenode[i]*sizeof(INT));

    if (design->ndnode_fenode[i] > 1 )
      dserror("Only one fe_node possible on each dnode !!");

  }


  /* fill dnode_fenode */
  for (i=0; i<io_dis->ngnode; i++)
  {
    actgnode = &(io_dis->gnode[i]);
    if (actgnode->ondesigntyp == ondnode)
    {
      id = actgnode->d.dnode->Id;
      design->dnode_fenode2[id][0] = actgnode->slave_node;
    }
  }



  /* dlines */
  /*========*/


  /* count number of gnodes on each design line */
  for (i=0; i<io_dis->ngline; i++)
  {
    actgline = &(io_dis->gline[i]);
    if (actgline->dline != NULL)
    {
      id = actgline->dline->Id;
      design->ndline_fenode[id] += 3;
    }
  }



  /* allocate dline_fenode[i] */
  for (i=0; i<design->ndline; i++)
  {
    if (design->ndline_fenode[i] > 0)
      design->dline_fenode2[i] = (INT*)CCAREALLOC(design->dline_fenode2[i],
          design->ndline_fenode[i]*sizeof(INT));
    counter = 0;


    /* fill dline_fenode */
    for (j=0; j<io_dis->ngline; j++)
    {
      actgline = &(io_dis->gline[j]);
      if (actgline->dline != NULL)
      {
        if (actgline->dline->Id != i)
          continue;
        design->dline_fenode2[i][counter++] = actgline->gnode[0]->slave_node;
        design->dline_fenode2[i][counter++] = actgline->gnode[1]->slave_node;
        /* Mittelknoten */
        design->dline_fenode2[i][counter++] = actgline->slave_node[0];
      }
    }

    if (counter != design->ndline_fenode[i])
      dserror_args(__FILE__, __LINE__,
          "Irgendetwas stimmt hier nicht !! counter = %d;  ndline_fenode[i] =%d",
          counter,design->ndline_fenode[i] );

  }




  /* dsurf */
  /*========*/


  /* count number of gnodes on each design surf */
  for (i=0; i<io_dis->ngsurf; i++)
  {
    actgsurf = &(io_dis->gsurf[i]);
    if (actgsurf->dsurf != NULL)
    {
      id = actgsurf->dsurf->Id;
      design->ndsurf_fenode[id] += 7;
    }
  }



  /* allocate dsurf_fenode[i] */
  for (i=0; i<design->ndsurf; i++)
  {
    if (design->ndsurf_fenode[i] > 0)
      design->dsurf_fenode2[i] = (INT*)CCAREALLOC(design->dsurf_fenode2[i],
          design->ndsurf_fenode[i]*sizeof(INT));
    counter = 0;


    /* fill dsurf_fenode */
    for (j=0; j<io_dis->ngsurf; j++)
    {
      actgsurf = &(io_dis->gsurf[j]);
      if (actgsurf->dsurf != NULL)
      {
        if (actgsurf->dsurf->Id != i)
          continue;
        design->dsurf_fenode2[i][counter++] = actgsurf->gnode[0]->slave_node;
        design->dsurf_fenode2[i][counter++] = actgsurf->gnode[1]->slave_node;
        design->dsurf_fenode2[i][counter++] = actgsurf->gnode[2]->slave_node;
        /* Mittelknoten der Flaeche */
        design->dsurf_fenode2[i][counter++] = actgsurf->slave_node[0][0];
        /* Mittelknoten der Linien */
        for (k=0; k<actgsurf->ngline; k++)
          design->dsurf_fenode2[i][counter++] = actgsurf->gline[k]->slave_node[0];

      }
    }

    if (counter != design->ndsurf_fenode[i])
      dserror_args(__FILE__, __LINE__,
          "Irgendetwas stimmt hier nicht !! counter = %d;  ndsurf_fenode[i] =%d",
          counter,design->ndsurf_fenode[i] );

  }






  /* create the pointers between gnodes... and the design */
  inpdesign_topology_fe(cal_dis);


#ifdef DEBUG
  dstrc_exit();
#endif


  return;

}  /* END of global_subdivide_tri */


#endif /* ifdef SUBDIV */

/*! @} (documentation module close)*/



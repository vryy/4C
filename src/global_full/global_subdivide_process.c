/*-----------------------------------------------------------------------*/
/*!
\file
\brief Service functions for element subdivision


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


#ifdef NURBS
#include "../nurbs/nurbs_wrappers.h"
#endif



/*-----------------------------------------------------------------------*/
/*!
  \brief process all gnodes on the given dis

  Generate gnodes on the second dis (calc_dis) for all gnodes on
  the first dis (io_dis)

  \param io_dis        DISCRET  (i/o) original dis
  \param calc_dis      DISCRET  (i/o) new dis
  \param node_counter  INT      (i/o) total number of nodes
  \param nodeid        INT      (i/o) current number for the new node

  \return void

  \author mn
  \date   09/05

*/
/*-----------------------------------------------------------------------*/
void process_gnodes(
    DISCRET            *io_dis,
    DISCRET            *cal_dis,
    INT                *node_counter,
    INT                *nodeid
    )

{

  INT                   i,j;
  GNODE                *actgnode;
  INT                   node_counter_old = *node_counter;


#ifdef DEBUG
  dstrc_enter("process_gnodes");
#endif


  printf("  Generating gnodes on the second dis on %7i gnodes...",io_dis->ngnode);
  fflush(stdout);

  for (i=0; i<io_dis->ngnode; i++)
  {
    actgnode = &(io_dis->gnode[i]);

    /* set Id of the new node */
    cal_dis->node[*node_counter].Id = *nodeid;

    /* set coords of the new node */
    for (j=0; j<3; j++)
      cal_dis->node[*node_counter].x[j] = actgnode->node->x[j];

    /* store slave_node Id in master gnode */
    actgnode->slave_node = *nodeid;

    /* create pointer from slave to master node */
    cal_dis->node[*node_counter].slave_node  = NULL;
    cal_dis->node[*node_counter].master_node = actgnode->node;

    /* create pointer from master to slave node */
    actgnode->node->slave_node  = &(cal_dis->node[*node_counter]);
    actgnode->node->master_node = NULL;


    (*node_counter)++;
    (*nodeid)++;
  }

  printf(" %7i gnodes created.\n",*node_counter-node_counter_old);

#ifdef DEBUG
  dstrc_exit();
#endif

  return;

}  /* END of process_gnodes */






/*-----------------------------------------------------------------------*/
/*!
  \brief process all glines on the given dis

  Generate gnodes on the second dis (calc_dis) for all glines on
  the first dis (io_dis)

  \param io_dis        DISCRET  (i/o) original dis
  \param calc_dis      DISCRET  (i/o) new dis
  \param subdivide     INT      (i)   number of subdivisions
  \param node_counter  INT      (i/o) total number of nodes
  \param nodeid        INT      (i/o) current number for the new node

  \return void

  \author mn
  \date   09/05

*/
/*-----------------------------------------------------------------------*/
void process_glines(
    DISCRET           *io_dis,
    DISCRET           *cal_dis,
    INT                subdivide,
    INT               *node_counter,
    INT               *nodeid
    )

{

  INT                   node_counter_old = *node_counter;
  INT                   i,j,k;
  GLINE                 *actgline;
  DSURF                 *actdsurf=NULL;

  DOUBLE                x1[3],dx[3];

  INT                   gtyp = 0;           /* 1- gline on nurbline
                                               2- gline on stline or no dline
                                               3- gline in dsurface
                                             */

#ifdef DEBUG
  dstrc_enter("process_glines");
#endif


  printf("  Generating gnodes on the second dis on %7i glines...",io_dis->ngline);
  fflush(stdout);

  /* loop glines and create nodes */
  for (i=0; i<io_dis->ngline; i++)
  {

    actgline = &(io_dis->gline[i]);


    /* identify typ of gline */
    if (actgline->dline != NULL)
    {
      /* this element edge is on a design line */
      switch (actgline->dline->typ)
      {
        case stline:
          gtyp = 2;
          break;

        case nurbline:
#ifdef NURBS
          gtyp = 1;
#else
          gtyp = 2;
          dswarning(1,15);
#endif
          break;

        case arcline:
          gtyp = 2;
          dswarning(1,14);
          /*dserror("ARCLINES not supported for SUBDIVION!!");*/
          break;

        default:
          dserror("Unknown linetype not supported for SUBDIVION!!");
          break;
      }

    }


    /* this element edge is NOT on a design line */
    else
    {

      INT   same_dsurf = 0;
      DSURF *list0[25];
      INT   num0 = 0;
      DSURF *list1[25];
      INT   num1 = 0;
      INT   d,e;

      /* find all dsurfs the two gnodes are on */


      /* the gnode is on a design surface */
      if (actgline->gnode[0]->ondesigntyp == ondsurf)
        list0[num0++] = actgline->gnode[0]->d.dsurf;

      if (actgline->gnode[1]->ondesigntyp == ondsurf)
        list1[num1++] = actgline->gnode[1]->d.dsurf;


      /* the gnode is on a design line */
      if (actgline->gnode[0]->ondesigntyp == ondline)
      {
        DLINE *actdline = actgline->gnode[0]->d.dline;
        for (d=0;d<actdline->ndsurf; d++)
          list0[num0++] = actdline->dsurf[d];
      }

      if (actgline->gnode[1]->ondesigntyp == ondline)
      {
        DLINE *actdline = actgline->gnode[1]->d.dline;
        for (d=0;d<actdline->ndsurf; d++)
          list1[num1++] = actdline->dsurf[d];
      }


      if(num0 > 25 || num1 > 25)
        dserror("Too many dsurf for this gnode!!");


      for (d=0;d<num0;d++)
        for (e=0;e<num1;e++)
          if (list0[d] == list1[e])
          {
            same_dsurf = 1;
            actdsurf   = list0[d];
          }



      if (same_dsurf==1)
      {
        /* gline is on ONE dsurf */
#ifdef NURBS
        gtyp = 3;
#else
        gtyp = 2;
        dswarning(1,16);
#endif
      }

      else
      {
        /* gline is NOT on a dsurf */
        gtyp = 2;
      }

    }






    switch (gtyp)
    {

#ifdef NURBS
      case 3:
        {
          /* this gline is on a dsurface */


          INT n,m;
          INT ierr;
          DOUBLE x2[3];
          NURBSURF *actnurbsurf;


          if (actdsurf == NULL)
            dserror("No dsurf for this gline!!");

          actnurbsurf = actdsurf->props.nurbsurf;
          actdsurf == NULL;



          /* create NURBS object */
          NURBS_OBJECT_PTR cp;
          NURBS_OBJECT_PTR weights;
          NURBS_OBJECT_PTR knots_u;
          NURBS_OBJECT_PTR knots_v;
          NURBS_OBJECT_PTR nurbs_surf;

          NURBS_OBJECT_PTR p0,p1,pt;
          DOUBLE u0,u1;
          DOUBLE v0,v1;

          NURBS_OBJECT_PTR p;
          DOUBLE u,v;


          cp = nurbs_matrix_p_create(actnurbsurf->num_cp[0],actnurbsurf->num_cp[1]);

          weights = nurbs_matrix_double_create(actnurbsurf->num_cp[0],actnurbsurf->num_cp[1]);

          knots_u = nurbs_vector_double_create(actnurbsurf->num_knots[0]);
          knots_v = nurbs_vector_double_create(actnurbsurf->num_knots[1]);


          for (n=0; n<actnurbsurf->num_cp[0]; n++)
            for (m=0; m<actnurbsurf->num_cp[1]; m++)
              nurbs_matrix_p_setvalue(cp,
                  n,
                  m,
                  actnurbsurf->cp[n][m][0],
                  actnurbsurf->cp[n][m][1],
                  actnurbsurf->cp[n][m][2]
                  );


          for (n=0; n<actnurbsurf->num_cp[0]; n++)
            for (m=0; m<actnurbsurf->num_cp[1]; m++)
              nurbs_matrix_double_setvalue(weights,
                  n,
                  m,
                  actnurbsurf->weights[n][m]);


          for (n=0; n<actnurbsurf->num_knots[0]; n++)
            nurbs_vector_double_setvalue(knots_u,
                n,
                actnurbsurf->knots[0][n]);


          for (n=0; n<actnurbsurf->num_knots[1]; n++)
            nurbs_vector_double_setvalue(knots_v,
                n,
                actnurbsurf->knots[1][n]);



          nurbs_surf = nurbs_surf_create(cp, weights, knots_u, knots_v,
              actnurbsurf->degree[0], actnurbsurf->degree[1]);


          /* set the nodal coordinates */
          x1[0] = (actgline->gnode[0]->node->x[0]);
          x1[1] = (actgline->gnode[0]->node->x[1]);
          x1[2] = (actgline->gnode[0]->node->x[2]);
          p0 = nurbs_p_create(x1[0], x1[1], x1[2]);

          ierr = nurbs_surf_project(nurbs_surf,p0,&u0,&v0);

          if (ierr == 0)
          {
            /* singularitaet in der projektion */
            pt = nurbs_surf_p(nurbs_surf,u0,v0);
            if(nurbs_p_dist(p0,pt) > 1e-4)
              dserror("singularity in projektion of point onto nurbssurface and resulting point not exact!!");
          }


          x2[0] = (actgline->gnode[1]->node->x[0]);
          x2[1] = (actgline->gnode[1]->node->x[1]);
          x2[2] = (actgline->gnode[1]->node->x[2]);
          p1 = nurbs_p_create(x2[0], x2[1], x2[2]);

          ierr = nurbs_surf_project(nurbs_surf,p1,&u1,&v1);

          if (ierr == 0)
          {
            /* singularitaet in der projektion */
            pt = nurbs_surf_p(nurbs_surf,u1,v1);
            if(nurbs_p_dist(p1,pt) > 1e-4)
              dserror("singularity in projektion of point onto nurbssurface and resulting point not exact!!");
          }




          for (j=0; j<subdivide-1; j++)
          {
            cal_dis->node[*node_counter].Id = *nodeid;

            /* store slave_node Id in gline */
            actgline->slave_node[j]   = *nodeid;

            /* set local coord of the new node */
            u =  u0 + (j+1)*(u1-u0)/subdivide;
            v =  v0 + (j+1)*(v1-v0)/subdivide;

            p = nurbs_surf_p(nurbs_surf,u,v);

            cal_dis->node[*node_counter].x[0] = nurbs_p_getx(p);
            cal_dis->node[*node_counter].x[1] = nurbs_p_gety(p);
            cal_dis->node[*node_counter].x[2] = nurbs_p_getz(p);


            (*node_counter)++;
            (*nodeid)++;

          }  /* for (j=0; j<subdivide-1; j++) */

        }
        break;
#endif


      case 2:
        /* stline */


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
          cal_dis->node[*node_counter].Id = *nodeid;

          /* store slave_node Id in gline */
          actgline->slave_node[j]   = *nodeid;

          /* set coords of the new node */
          for (k=0; k<3; k++)
            cal_dis->node[*node_counter].x[k] =  x1[k] + (j+1)*dx[k];


          (*node_counter)++;
          (*nodeid)++;
        }

        break;



#ifdef NURBS
      case 1:
        /* gline on a nurbline */

        {
        INT n;
        DOUBLE x2[3];
        NURBLINE *actnurbline;


        actnurbline = actgline->dline->props.nurbline;

        /* create NURBS object */
        NURBS_OBJECT_PTR cp;
        NURBS_OBJECT_PTR weights;
        NURBS_OBJECT_PTR knots;
        NURBS_OBJECT_PTR nurbs_curve;

        NURBS_OBJECT_PTR p1,p2;
        DOUBLE u1,u2;

        NURBS_OBJECT_PTR p;
        DOUBLE u;


        cp = nurbs_vector_p_create(actnurbline->num_cp);

        weights = nurbs_vector_double_create(actnurbline->num_cp);

        knots = nurbs_vector_double_create(actnurbline->num_knots);


        for (n=0; n<actnurbline->num_cp; n++)
          nurbs_vector_p_setvalue(cp,
              n,
              actnurbline->cp[n][0],
              actnurbline->cp[n][1],
              actnurbline->cp[n][2]
              );


        for (n=0; n<actnurbline->num_knots; n++)
          nurbs_vector_double_setvalue(knots,
              n,
              actnurbline->knots[n]);


        for (n=0; n<actnurbline->num_cp; n++)
          nurbs_vector_double_setvalue(weights,
              n,
              actnurbline->weights[n]);


        nurbs_curve = nurbs_curve_create(cp, weights, knots, actnurbline->degree);

        nurbs_curve_print_ps(nurbs_curve);


        /* coordinates of the beginning of this gline */
        x1[0] = actgline->gnode[0]->node->x[0];
        x1[1] = actgline->gnode[0]->node->x[1];
        x1[2] = actgline->gnode[0]->node->x[2];

        /* local coord of the beginning */
        p1 = nurbs_p_create(x1[0], x1[1], x1[2]);
        u1 = nurbs_find_u(nurbs_curve,p1);


        /* coordinates of the end of this gline */
        x2[0] = actgline->gnode[1]->node->x[0];
        x2[1] = actgline->gnode[1]->node->x[1];
        x2[2] = actgline->gnode[1]->node->x[2];

        /* local ccord of the beginning */
        p2 = nurbs_p_create(x2[0], x2[1], x2[2]);
        u2 = nurbs_find_u(nurbs_curve,p2);



        for (j=0; j<subdivide-1; j++)
        {
          cal_dis->node[*node_counter].Id = *nodeid;

          /* store slave_node Id in gline */
          actgline->slave_node[j]   = *nodeid;

          /* set local coord of the new node */
          u =  u1 + (j+1)*(u2-u1)/subdivide;

          p = nurbs_curve_p(nurbs_curve,u);

          cal_dis->node[*node_counter].x[0] = nurbs_p_getx(p);
          cal_dis->node[*node_counter].x[1] = nurbs_p_gety(p);
          cal_dis->node[*node_counter].x[2] = nurbs_p_getz(p);


          (*node_counter)++;
          (*nodeid)++;

        }  /* for (j=0; j<subdivide-1; j++) */


        }
        break;
#endif


      default:
        dserror("Unknown gtyp for SUBDIVION!!");
        break;

    }  /* switch (gtyp) */


  }  /* for (i=0; i<io_dis->ngline; i++) */

  printf(" %7i gnodes created.\n",*node_counter-node_counter_old);

#ifdef DEBUG
  dstrc_exit();
#endif

  return;

}  /* END of process_glines */





/*-----------------------------------------------------------------------*/
/*!
  \brief process the given gsurf

  Generate gnodes on the second dis (calc_dis) for the given gsurf


  \param actgsurf      GSURF    (i)   the current gsurf
  \param cal_dis       DISCRET  (i/o) new dis
  \param subdivide     INT      (i)   number of subdivisions
  \param gnode[4]      GNODE    (o)   the four gnodes of the surface
  \param node_counter  INT      (i/o) total number of nodes
  \param nodeid        INT      (i/o) current number for the new node

  \return void

  \author mn
  \date   09/05

*/
/*-----------------------------------------------------------------------*/
void process_gsurf(
    GSURF              *actgsurf,
    DISCRET            *cal_dis,
    INT                 subdivide,
    GNODE              *gnode[4],
    INT                *node_counter,
    INT                *nodeid
    )

{

  INT                   j,k,l;

  DOUBLE                x0[3],x1[3],x2[3],x3[3];

  GLINE                *actgline = NULL;
  GNODE                *actgnode = NULL;

  DSURF_TYP             typ;


#ifdef DEBUG
  dstrc_enter("process_gsurfs");
#endif


    /* specify node 0 */
    gnode[0] = actgsurf->gnode[0];
    actgsurf->node_ind[0] = 0;



    /* find the line 0 */
    for (j=0; j<actgsurf->ngline; j++)
    {
      /* this line is in same direction */
      if (gnode[0] == actgsurf->gline[j]->gnode[0])
      {
        actgline = actgsurf->gline[j];
        actgsurf->line_ind[0][0] = j;
        actgsurf->line_ind[1][0] = 0;

        gnode[1] = actgline->gnode[1];

        /* find the node 1 */
        for (k=0; k<actgsurf->ngnode; k++)
        {
          if (gnode[1] == actgsurf->gnode[k])
          {
            actgsurf->node_ind[1] = k;
            break;
          }
        }

        break;
      }

      /* this line is in opposite direction */
      if (gnode[0] == actgsurf->gline[j]->gnode[1])
      {
        actgline = actgsurf->gline[j];
        actgsurf->line_ind[0][0] = j;
        actgsurf->line_ind[1][0] = 1;

        gnode[1] = actgline->gnode[0];

        /* find the node 1 */
        for (k=0; k<actgsurf->ngnode; k++)
        {
          if (gnode[1] == actgsurf->gnode[k])
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
      if (gnode[1] == actgsurf->gline[j]->gnode[0])
      {
        actgline = actgsurf->gline[j];
        actgsurf->line_ind[0][1] = j;
        actgsurf->line_ind[1][1] = 0;

        gnode[2] = actgline->gnode[1];

        /* find the node 2 */
        for (k=0; k<actgsurf->ngnode; k++)
        {
          if (gnode[2] == actgsurf->gnode[k])
          {
            actgsurf->node_ind[2] = k;
            break;
          }
        }

        break;
      }

      /* this line is in opposite direction */
      if (gnode[1] == actgsurf->gline[j]->gnode[1])
      {
        actgline = actgsurf->gline[j];
        actgsurf->line_ind[0][1] = j;
        actgsurf->line_ind[1][1] = 1;

        gnode[2] = actgline->gnode[0];

        /* find the node 2 */
        for (k=0; k<actgsurf->ngnode; k++)
        {
          if (gnode[2] == actgsurf->gnode[k])
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
      if (gnode[2] == actgsurf->gline[j]->gnode[1])
      {
        actgline = actgsurf->gline[j];
        actgsurf->line_ind[0][2] = j;
        actgsurf->line_ind[1][2] = 0;

        gnode[3] = actgline->gnode[0];

        /* find the node 3 */
        for (k=0; k<actgsurf->ngnode; k++)
        {
          if (gnode[3] == actgsurf->gnode[k])
          {
            actgsurf->node_ind[3] = k;
            break;
          }
        }

        break;
      }

      /* this line is in opposite direction */
      if (gnode[2] == actgsurf->gline[j]->gnode[0])
      {
        actgline = actgsurf->gline[j];
        actgsurf->line_ind[0][2] = j;
        actgsurf->line_ind[1][2] = 1;

        gnode[3] = actgline->gnode[1];

        /* find the node 3 */
        for (k=0; k<actgsurf->ngnode; k++)
        {
          if (gnode[3] == actgsurf->gnode[k])
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
      if (gnode[3] == actgsurf->gline[j]->gnode[1])
      {
        actgline = actgsurf->gline[j];
        actgsurf->line_ind[0][3] = j;
        actgsurf->line_ind[1][3] = 0;

        actgnode = actgline->gnode[0];

        if (actgnode != gnode[0])
          dserror("NODE MIXUP!!\n");

        break;
      }

      /* this line is in opposite direction */
      if (gnode[3] == actgsurf->gline[j]->gnode[0])
      {
        actgline = actgsurf->gline[j];
        actgsurf->line_ind[0][3] = j;
        actgsurf->line_ind[1][3] = 1;

        actgnode = actgline->gnode[1];

        if (actgnode != gnode[0])
          dserror("NODE MIXUP!!\n");

        break;
      }
    }  /* for (j=1; j<gsurf->ngline; j++) */




    if (actgsurf->dsurf == NULL)
      /* this element face is NOT on a design surf */
      typ = flatsurf;
    else
      /* this element face is on a design surf */
#ifdef NURBS
      typ = actgsurf->dsurf->typ;
#else
      typ = flatsurf;
      dswarning(1,16);
#endif


    switch (typ)
    {
      case flatsurf:

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
            cal_dis->node[*node_counter].Id  = *nodeid;

            /* store slave_node Id in gsurf */
            actgsurf->slave_node[j][k] = *nodeid;

            /* set coords of the new node */
            for (l=0; l<3; l++)
              cal_dis->node[*node_counter].x[l] =
                (x1[l]*(j+1) + x0[l]*(subdivide-j-1)) * (subdivide-k-1) +
                (x2[l]*(j+1) + x3[l]*(subdivide-j-1)) * (k+1);

            (*node_counter)++;
            (*nodeid)++;

          }  /* for (k=0; k<subdivide-1; k++) */

        }  /* for (j=0; j<subdivide-1; j++) */

        break;



#ifdef NURBS
      case nurbsurf:
        {
        INT n,m;
        INT ierr;
        DOUBLE x2[3];
        NURBSURF *actnurbsurf;


        actnurbsurf = actgsurf->dsurf->props.nurbsurf;

        /* create NURBS object */
        NURBS_OBJECT_PTR cp;
        NURBS_OBJECT_PTR weights;
        NURBS_OBJECT_PTR knots_u;
        NURBS_OBJECT_PTR knots_v;
        NURBS_OBJECT_PTR nurbs_surf;

        NURBS_OBJECT_PTR p0,p1,p2,p3,pt;
        DOUBLE u0,u1,u2,u3;
        DOUBLE v0,v1,v2,v3;

        NURBS_OBJECT_PTR p;
        DOUBLE u,v;


        cp = nurbs_matrix_p_create(actnurbsurf->num_cp[0],actnurbsurf->num_cp[1]);

        weights = nurbs_matrix_double_create(actnurbsurf->num_cp[0],actnurbsurf->num_cp[1]);

        knots_u = nurbs_vector_double_create(actnurbsurf->num_knots[0]);
        knots_v = nurbs_vector_double_create(actnurbsurf->num_knots[1]);


        for (n=0; n<actnurbsurf->num_cp[0]; n++)
          for (m=0; m<actnurbsurf->num_cp[1]; m++)
            nurbs_matrix_p_setvalue(cp,
                n,
                m,
                actnurbsurf->cp[n][m][0],
                actnurbsurf->cp[n][m][1],
                actnurbsurf->cp[n][m][2]
                );


        for (n=0; n<actnurbsurf->num_cp[0]; n++)
          for (m=0; m<actnurbsurf->num_cp[1]; m++)
            nurbs_matrix_double_setvalue(weights,
                n,
                m,
                actnurbsurf->weights[n][m]);


        for (n=0; n<actnurbsurf->num_knots[0]; n++)
          nurbs_vector_double_setvalue(knots_u,
              n,
              actnurbsurf->knots[0][n]);


        for (n=0; n<actnurbsurf->num_knots[1]; n++)
          nurbs_vector_double_setvalue(knots_v,
              n,
              actnurbsurf->knots[1][n]);



        nurbs_surf = nurbs_surf_create(cp, weights, knots_u, knots_v,
            actnurbsurf->degree[0], actnurbsurf->degree[1]);


        /* set the nodal coordinates */
        x0[0] = (actgsurf->gnode[actgsurf->node_ind[0]]->node->x[0]);
        x0[1] = (actgsurf->gnode[actgsurf->node_ind[0]]->node->x[1]);
        x0[2] = (actgsurf->gnode[actgsurf->node_ind[0]]->node->x[2]);
        p0 = nurbs_p_create(x0[0], x0[1], x0[2]);

        ierr = nurbs_surf_project(nurbs_surf,p0,&u0,&v0);

        if (ierr == 0)
        {
          /* singularitaet in der projektion */
          pt = nurbs_surf_p(nurbs_surf,u0,v0);
          if(nurbs_p_dist(p0,pt) > 1e-4)
            dserror("singularity in projektion of point onto nurbssurface and resulting point not exact!!");
        }



        x1[0] = (actgsurf->gnode[actgsurf->node_ind[1]]->node->x[0]);
        x1[1] = (actgsurf->gnode[actgsurf->node_ind[1]]->node->x[1]);
        x1[2] = (actgsurf->gnode[actgsurf->node_ind[1]]->node->x[2]);
        p1 = nurbs_p_create(x1[0], x1[1], x1[2]);

        ierr = nurbs_surf_project(nurbs_surf,p1,&u1,&v1);

        if (ierr == 0)
        {
          /* singularitaet in der projektion */
          pt = nurbs_surf_p(nurbs_surf,u1,v1);
          if(nurbs_p_dist(p1,pt) > 1e-4)
            dserror("singularity in projektion of point onto nurbssurface and resulting point not exact!!");
        }


        x2[0] = (actgsurf->gnode[actgsurf->node_ind[2]]->node->x[0]);
        x2[1] = (actgsurf->gnode[actgsurf->node_ind[2]]->node->x[1]);
        x2[2] = (actgsurf->gnode[actgsurf->node_ind[2]]->node->x[2]);
        p2 = nurbs_p_create(x2[0], x2[1], x2[2]);

        ierr = nurbs_surf_project(nurbs_surf,p2,&u2,&v2);

        if (ierr == 0)
        {
          /* singularitaet in der projektion */
          pt = nurbs_surf_p(nurbs_surf,u2,v2);
          if(nurbs_p_dist(p2,pt) > 1e-4)
            dserror("singularity in projektion of point onto nurbssurface and resulting point not exact!!");
        }


        x3[0] = (actgsurf->gnode[actgsurf->node_ind[3]]->node->x[0]);
        x3[1] = (actgsurf->gnode[actgsurf->node_ind[3]]->node->x[1]);
        x3[2] = (actgsurf->gnode[actgsurf->node_ind[3]]->node->x[2]);
        p3 = nurbs_p_create(x3[0], x3[1], x3[2]);

        ierr = nurbs_surf_project(nurbs_surf,p3,&u3,&v3);

        if (ierr == 0)
        {
          /* singularitaet in der projektion */
          pt = nurbs_surf_p(nurbs_surf,u3,v3);
          if(nurbs_p_dist(p3,pt) > 1e-4)
            dserror("singularity in projektion of point onto nurbssurface and resulting point not exact!!");
        }


        u0 = u0/(subdivide*subdivide);
        u1 = u1/(subdivide*subdivide);
        u2 = u2/(subdivide*subdivide);
        u3 = u3/(subdivide*subdivide);

        v0 = v0/(subdivide*subdivide);
        v1 = v1/(subdivide*subdivide);
        v2 = v2/(subdivide*subdivide);
        v3 = v3/(subdivide*subdivide);

        /* generate the nodes */
        for (j=0; j<subdivide-1; j++)
        {
          for (k=0; k<subdivide-1; k++)
          {
            cal_dis->node[*node_counter].Id  = *nodeid;

            /* store slave_node Id in gsurf */
            actgsurf->slave_node[j][k] = *nodeid;


            /* set local coord of the new node */

            u = (u1*(j+1) + u0*(subdivide-j-1)) * (subdivide-k-1) +
              (u2*(j+1) + u3*(subdivide-j-1)) * (k+1);

            v = (v1*(j+1) + v0*(subdivide-j-1)) * (subdivide-k-1) +
              (v2*(j+1) + v3*(subdivide-j-1)) * (k+1);

            p = nurbs_surf_p(nurbs_surf,u,v);

            /* set coords of the new node */
            cal_dis->node[*node_counter].x[0] = nurbs_p_getx(p);
            cal_dis->node[*node_counter].x[1] = nurbs_p_gety(p);
            cal_dis->node[*node_counter].x[2] = nurbs_p_getz(p);

            (*node_counter)++;
            (*nodeid)++;

          }  /* for (k=0; k<subdivide-1; k++) */

        }  /* for (j=0; j<subdivide-1; j++) */


        }
        break;
#endif


      default:
        dserror("Unknown surftype not supported for SUBDIVION!!");
        break;

    }


#ifdef DEBUG
  dstrc_exit();
#endif

  return;

}  /* END of process_gsurfs */



/*-----------------------------------------------------------------------*/
/*!
  \brief process the given gsurf

  Generate gnodes on the second dis (calc_dis) for the given triangular gsurf


  \param actgsurf      GSURF    (i)   the current gsurf
  \param cal_dis       DISCRET  (i/o) new dis
  \param gnode[3]      GNODE    (o)   the three gnodes of the surface
  \param node_counter  INT      (i/o) total number of nodes
  \param nodeid        INT      (i/o) current number for the new node

  \return void

  \author mn
  \date   09/05

*/
/*-----------------------------------------------------------------------*/
void process_gsurf_tri(
    GSURF              *actgsurf,
    DISCRET            *cal_dis,
    GNODE              *gnode[3],
    INT                *node_counter,
    INT                *nodeid
    )

{

  INT                   j,k,l;

  DOUBLE                x0[3],x1[3],x2[3],x3[3];

  GLINE                *actgline = NULL;
  GNODE                *actgnode = NULL;

  DSURF_TYP             typ;


#ifdef DEBUG
  dstrc_enter("process_gsurfs_tri");
#endif


    /* specify the nodes */
    gnode[0] = actgsurf->gnode[0];
    gnode[1] = actgsurf->gnode[1];
    gnode[2] = actgsurf->gnode[2];


    /* generate the node */
    cal_dis->node[*node_counter].Id  = *nodeid;

    /* store slave_node Id in gsurf */
    actgsurf->slave_node[0][0] = *nodeid;


    if (actgsurf->dsurf == NULL)
      /* this element face is NOT on a design surf */
      typ = flatsurf;
    else
      /* this element face is on a design surf */
#ifdef NURBS
      typ = actgsurf->dsurf->typ;
#else
      typ = flatsurf;
#endif

    switch (typ)
    {
      case flatsurf:

        /* set coords of the new node */
        for (l=0; l<3; l++)
          cal_dis->node[*node_counter].x[l] =
            (gnode[0]->node->x[l] + gnode[1]->node->x[l] + gnode[2]->node->x[l]) / 3;

        (*node_counter)++;
        (*nodeid)++;
        break;


#ifdef NURBS
      case nurbsurf:
        dserror("NURBSURF not yet supported for SUBDIVISION");
        break;
#endif


      default:
        dserror("Unknown surftype not supported for SUBDIVION!!");
        break;

    }

#ifdef DEBUG
  dstrc_exit();
#endif

  return;

}  /* END of process_gsurfs_tri */




#endif /* ifdef SUBDIV */

/*! @} (documentation module close)*/



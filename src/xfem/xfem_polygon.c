/*!----------------------------------------------------------------------
\file
\brief subroutines for polygonization

<pre>
Maintainer: Baris Irhan
            irhan@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/irhan/
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef D_XFEM
#include "../headers/standardtypes.h"
#include "../ls/ls_prototypes.h"
#include "xfem_prototypes.h"
/*! 
\addtogroup XFEM 
*//*! @{ (documentation module open)*/




static  LS_INT_DATA  *intdata;
static  LS_POLY_DATA *polydata;
static  DOUBLE       *pstart;
static  DOUBLE       *pend;
static  DOUBLE       *pdiag;
static  INT          *polycon[2];
static  INT           sedge;
static  INT           eedge;

static  INT          *ind;
static  DOUBLE        R[2];
static  DOUBLE        A[2][2];
static  DOUBLE        GP[2];
static  DOUBLE        xtarget[2];
static  DOUBLE        xyze[2][4];
static  DOUBLE        xyze_int[2][3];
static  DOUBLE        xyze_subp[2][3];
static  DOUBLE        polygonxy[2][2][7];
static  INT          *polygonmat[2];
static  DOUBLE       *polygonwgt[2];
static  DOUBLE       *polygonGP[2][2];
static  DOUBLE        N[4];
static  DOUBLE        gradN[2][4];

static  INT           edge_to_triangle[4] = {1,2,2,1};
static  INT           node_to_edge_tri[2][3] = {
                                                {3,1,2},
                                                {1,2,3}
                                               };

static  INT           node_to_node_tri[2][3] = {
                                                {2,3,1},
                                                {3,1,2}
                                               };

static  DOUBLE        rsI[2][4] =
                                  {
                                   {-1.0, 1.0, 1.0,-1.0},
                                   {-1.0,-1.0, 1.0, 1.0}
                                  };

static  INT           conr[2][3] =
                                   {
                                    {2,4,1},
                                    {2,3,4}
                                   };

static  INT           polygoncon[3][7] =
                                         {
                                          {1,4,5,4,2,3,5},
                                          {4,5,1,2,3,5,4},
                                          {6,6,6,7,7,7,7}
                                         };

static FILE          *f01;
static FILE          *f02;
static FILE          *f03;
static FILE          *f04;
static FILE          *f05; 


static DOUBLE         atri; /* area of triangle */
static DOUBLE         arect; /* area of rectangle */
static DOUBLE         asubp; /* area of subpolygon */
static const DOUBLE   HALF = ONE/TWO;



/*!----------------------------------------------------------------------
\brief control subroutine for polygonization

<pre>                                                            irhan 05/04
control subroutine for polygonization.
</pre>

*----------------------------------------------------------------------*/
void xfem_polygon(
  XFEMPOLYFLAG xfempolyflag,
  ELEMENT* myls2   
  )      
{
#ifdef DEBUG 
  dstrc_enter("xfem_polygon");
#endif
/*----------------------------------------------------------------------*/
  
/*------------------------------------------------- switch to do option */
  switch (xfempolyflag)
  {
/*---------------------------------------------------------- initialize */
      case xfem_poly_initialize:
        xfem_polygon_init(myls2);
        break;
/*-------------------------------------------- construct material index */
      case xfem_poly_material:
        xfem_polygon_mat(myls2);
        break;
/*----------------------------------------- construct polygon structure */
      case xfem_poly_construct:
        xfem_polygon_cons(myls2);
        break;
/*------------------------------------------------ compute Gauss points */
/*
 *  NOTE =>
 *  at the moment one Gauss point is assigned to
 *  mid-point of each subpolygon
 */
      case xfem_poly_computeGP:
        xfem_polygon_GP(myls2);
        break;
/*---------------------------------- write polygon information to files */
      case xfem_poly_write:
        xfem_polygon_write(myls2);
        break;
/*---------------------------------------------------------- open files */
      case xfem_poly_open:
        xfem_polygon_open();
        break;
/*--------------------------------------------------------- close files */
      case xfem_poly_close:
        xfem_polygon_close();
        break;
      default:
        dserror("action unknown\n");
  } /* end swtich (*action) */
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of xfem_polygon */



/*!----------------------------------------------------------------------
\brief initialize

<pre>                                                            irhan 05/04
perform initialization.
</pre>

*----------------------------------------------------------------------*/
void xfem_polygon_init(
  ELEMENT *myls2
  )
{
  INT     i;
  
#ifdef DEBUG 
  dstrc_enter("xfem_poly_init");
#endif
/*----------------------------------------------------------------------*/

  /* set intdata */
  intdata = &(myls2->e.ls2->intdata[0]);
  /* set polydata */
  polydata = &(myls2->e.ls2->polydata[0]);
  /* set start point */
  pstart = intdata->p[0];
  /* set end point */
  pend = intdata->p[1];
  /* set start edge */
  sedge = intdata->edge[0];
  /* set end edge */  
  eedge = intdata->edge[1];
  /* set polycon */
  polycon[0] = intdata->polycon[0];
  polycon[1] = intdata->polycon[1];
  /* set ind */
  ind = polydata->ind;
  /* set polygonmat */
  polygonmat[0] = polydata->polygonmat[0];
  polygonmat[1] = polydata->polygonmat[1];
  /* set polygonwgt */
  polygonwgt[0] = polydata->polygonwgt[0];
  polygonwgt[1] = polydata->polygonwgt[1];
  /* set polygonGP */
  polygonGP[0][0] = polydata->polygonGP[0][0];
  polygonGP[0][1] = polydata->polygonGP[0][1];
  polygonGP[1][0] = polydata->polygonGP[1][0];
  polygonGP[1][1] = polydata->polygonGP[1][1];
  /* set diagonal point */
  pdiag = intdata->pd;
  if (intdata->is_diagcut==0)
  {
    pdiag[0] = (pstart[0]+pend[0])/2.0;
    pdiag[1] = (pstart[1]+pend[1])/2.0;
  }
  /* set coordinate data for element */
  for(i=0; i<myls2->numnp; i++)
  {
    xyze[0][i] = myls2->node[i]->x[0];
    xyze[1][i] = myls2->node[i]->x[1];
  }
  /* set coordinate data for interface */
  xyze_int[0][0] = pstart[0];
  xyze_int[1][0] = pstart[1];
  xyze_int[0][1] = pdiag[0];
  xyze_int[1][1] = pdiag[1];
  xyze_int[0][2] = pend[0];
  xyze_int[1][2] = pend[1];
 
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of xfem_poly_init */



/*!----------------------------------------------------------------------
\brief construct the polygon coordinate array.

<pre>                                                            irhan 05/04
construct the polygon coordinate array. note that each triangular domain cut
by the interface is divided into seven subtriangles.
</pre>

*----------------------------------------------------------------------*/
void xfem_polygon_cons(
  ELEMENT *myls2
  )      
{
  INT     i;
  INT     nd1,nd2,nd3,nd4,nd5;
  INT     numtri,numtristart,numtriend;
  INT     ind1[2][2];
  
#ifdef DEBUG 
  dstrc_enter("xfem_poly_cons");
#endif
/*----------------------------------------------------------------------*/

  switch (myls2->distyp)
  {
      case tri3:   /* => 3 noded triangular element */
        ind[0] = 1;
        ind[1] = -1;

        nd1 = polycon[0][0];
        nd2 = polycon[0][1];
        nd3 = polycon[0][2];

        if (node_to_edge_tri[0][nd1-1]==eedge)
        {
          nd4 = 1;
          nd5 = 3;
        }
        else
        {
          nd4 = 3;
          nd5 = 1;          
        }

        polygonxy[0][0][0] = xyze[0][nd1-1];
        polygonxy[0][1][0] = xyze[1][nd1-1];
        polygonxy[0][0][1] = xyze[0][nd2-1];
        polygonxy[0][1][1] = xyze[1][nd2-1];
        polygonxy[0][0][2] = xyze[0][nd3-1];
        polygonxy[0][1][2] = xyze[1][nd3-1];
        
        polygonxy[0][0][3] = xyze_int[0][nd4-1];
        polygonxy[0][1][3] = xyze_int[1][nd4-1];
        polygonxy[0][0][4] = xyze_int[0][nd5-1];
        polygonxy[0][1][4] = xyze_int[1][nd5-1];
        
        polygonxy[0][0][5] = (xyze[0][nd1-1]+xyze_int[0][nd4-1]+
                              xyze_int[0][nd5-1])/3.0;
        polygonxy[0][1][5] = (xyze[1][nd1-1]+xyze_int[1][nd4-1]+
                              xyze_int[1][nd5-1])/3.0;
        polygonxy[0][0][6] = (xyze[0][nd2-1]+xyze[0][nd3-1]+
                              xyze_int[0][nd4-1]+xyze_int[0][nd5-1])/4.0;
        polygonxy[0][1][6] = (xyze[1][nd2-1]+xyze[1][nd3-1]+
                              xyze_int[1][nd4-1]+xyze_int[1][nd5-1])/4.0;
        break;
      case quad4:  /* => 4 noded quadrilateral element */
        /* check to which triangle start and end edges belong to  */
        numtristart = edge_to_triangle[sedge-1];
        numtriend = edge_to_triangle[eedge-1];
        if (numtristart==numtriend)
        {
          numtri = numtristart;
          if (numtri==1)
          {
            ind[0] = 1;
            ind[1] = 0;
          }
          else
          {
            ind[0] = 0;
            ind[1] = 1;
          }
        }
        else
        {
          ind[0] = 1;
          ind[1] = 1;
          if (numtristart==1)
          {
            ind1[0][0] = 1;
            ind1[1][0] = 2;
            ind1[0][1] = 3;
            ind1[1][1] = 2;
          }
          else
          {
            ind1[0][1] = 1;
            ind1[1][1] = 2;
            ind1[0][0] = 3;
            ind1[1][0] = 2;
          }
        }
        
        /* construct polygon coordinate data */
        for (i=0; i<2; i++)
        {
          if (ind[i]==0) continue;
          
          nd1 = polycon[i][0];
          nd2 = polycon[i][1];
          nd3 = polycon[i][2];
          
          if ((nd1==2 && i==0) || (nd1==4 && i==1))
          {
            nd4 =  ind1[1][i];
            nd5 =  ind1[0][i];
          }
          
          if ((nd1==4 && i==0) || (nd1==2 && i==1))
          {
            nd4 =  ind1[0][i];
            nd5 =  ind1[1][i];        
          }
          
          if (nd1==1 || nd1==3)
          {
            if (sedge==1 || sedge==3)
            {
              nd4 =  1;
              nd5 =  3;                  
            }
            else
            {
              nd4 =  3;
              nd5 =  1;                            
            }
          }
          
          polygonxy[i][0][0] = xyze[0][nd1-1];
          polygonxy[i][1][0] = xyze[1][nd1-1];
          polygonxy[i][0][1] = xyze[0][nd2-1];
          polygonxy[i][1][1] = xyze[1][nd2-1];
          polygonxy[i][0][2] = xyze[0][nd3-1];
          polygonxy[i][1][2] = xyze[1][nd3-1];
          
          polygonxy[i][0][3] = xyze_int[0][nd4-1];
          polygonxy[i][1][3] = xyze_int[1][nd4-1];
          polygonxy[i][0][4] = xyze_int[0][nd5-1];
          polygonxy[i][1][4] = xyze_int[1][nd5-1];
          
          polygonxy[i][0][5] = (xyze[0][nd1-1]+xyze_int[0][nd4-1]+
                                xyze_int[0][nd5-1])/3.0;
          polygonxy[i][1][5] = (xyze[1][nd1-1]+xyze_int[1][nd4-1]+
                                xyze_int[1][nd5-1])/3.0;
          polygonxy[i][0][6] = (xyze[0][nd2-1]+xyze[0][nd3-1]+
                                xyze_int[0][nd4-1]+xyze_int[0][nd5-1])/4.0;
          polygonxy[i][1][6] = (xyze[1][nd2-1]+xyze[1][nd3-1]+
                                xyze_int[1][nd4-1]+xyze_int[1][nd5-1])/4.0;
        }
        break;
      default:
        dserror("typ unknown!");
  } /* end switch(typ) */
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of xfem_poly_cons */



/*!----------------------------------------------------------------------
\brief compute Gauss points and weights for refined integration

<pre>                                                            irhan 05/04
in this subroutine Gauss points and weights are computed for refined integration.
each triangular domain cut by the interface is divided into seven subtriangles.
one Gauss point is assigned to center of each subtriangle.
</pre>

*----------------------------------------------------------------------*/
void xfem_polygon_GP(
  ELEMENT *myls2)
{
  INT              i,j;
  DOUBLE           sum;
  const DOUBLE     tol=1.0E-05;
  
#ifdef DEBUG 
  dstrc_enter("xfem_poly_GP");
#endif
/*----------------------------------------------------------------------*/


  switch (myls2->distyp)
  {
      case tri3:   /* => 3 noded triangular element */
        /* compute the area enclosed by subpolygon */
        atri = xfem_polygon_area_subtri(-1);
        /* loop over the subpolygons */
        for (j=0; j<7; j++)
        {
          /* access to coordinate information of the subpolygon */
          xfem_polygon_getsubp(0,j);
          /* compute the area of subpolygon */
          asubp = xfem_polygon_area_tri();
          /* compute the weight */
          polygonwgt[0][j] = HALF*asubp/atri;
          /* compute global coordinates of xtarget */
          xfem_polygon_target_tri();
          /* compute area coordinates of xtarget */
          xfem_polygon_compGP(myls2->distyp);
          /* set Gauss point */
          polygonGP[0][0][j] = GP[0]; 
          polygonGP[0][1][j] = GP[1];
        }
        /* simple check for weights */
        sum = 0.0;
        for (j=0; j<7; j++)
          sum += polygonwgt[0][j];
        if (abs(sum-HALF)>tol) dserror("weights are wrong!\n");
        break;
      case quad4:  /* => 4 noded quadrilateral element */
        /* compute area of the rectangle */
        xfem_polygon_area_rect();
        /* loop over the triangles */
        for (i=0; i<2; i++)
        {
          if (ind[i]==0)
          {
            /* compute the area enclosed by subpolygon */
            asubp = xfem_polygon_area_subtri(i);
            /* compute weight */
            polygonwgt[i][0] = FOUR*asubp/arect;
            /* compute global coordinates of xtarget */
            xfem_polygon_target_subtri(i);
            /* compute local coordinates of xtarget */
            xfem_polygon_compGP(myls2->distyp);
            /* set Gauss point */
            polygonGP[i][0][0] = GP[0];
            polygonGP[i][1][0] = GP[1];
            continue;      
          }
          else
          {
            /* loop over the subpolygons */
            for (j=0; j<7; j++)
            {
              /* access to coordinate information of the subpolygon */
              xfem_polygon_getsubp(i,j);
              /* compute the area of subpolygon */
              asubp = xfem_polygon_area_tri();
              /* compute the weight */
              polygonwgt[i][j] = FOUR*asubp/arect;
              /* compute global coordinates of xtarget */
              xfem_polygon_target_tri();
              /* compute local coordinates of xtarget */
              xfem_polygon_compGP(myls2->distyp);
              /* set Gauss point */
              polygonGP[i][0][j] = GP[0];
              polygonGP[i][1][j] = GP[1];
            }
          }
        }
        /* simple check for weights */
        sum = 0.0;
        for (i=0; i<2; i++)
        {
          if (ind[i]==0) sum += polygonwgt[i][0];
          else
          {
            for (j=0; j<7; j++)
              sum += polygonwgt[i][j];
          }
        }
        if (abs(sum-FOUR)>tol) dserror("weights are wrong!\n");
        break;
      default:
        dserror("typ unknown!");
  } /* end switch(typ) */
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of xfem_poly_GP */



/*!----------------------------------------------------------------------
\brief construct the coordinate array for the corresponding subtriangle

<pre>                                                            irhan 05/04
construct the coordinate array for the corresponding subtriangle.
</pre>

*----------------------------------------------------------------------*/
void xfem_polygon_getsubp(
  INT ntri,
  INT nsubp
  )
{
  INT     i,j;
  INT     nd;
  
#ifdef DEBUG 
  dstrc_enter("xfem_poly_getsubp");
#endif
/*----------------------------------------------------------------------*/

  /* initialize */  
  for (i=0; i<2; i++)
    for (j=0; j<3; j++)
      xyze_subp[i][j] = 0.0;

  /* loop */
  for (i=0; i<3; i++)
  {
    nd = polygoncon[i][nsubp]-1;
    xyze_subp[0][i] = polygonxy[ntri][0][nd];
    xyze_subp[1][i] = polygonxy[ntri][1][nd];
  }

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of xfem_poly_getsubp */



/*!----------------------------------------------------------------------
\brief compute the coordinates of the Gauss point corresponding to the
subtriangle

<pre>                                                            irhan 05/04
compute the coordinates of the Gauss point corresponding to the subtriangle.
</pre>

*----------------------------------------------------------------------*/
void xfem_polygon_compGP(
  DIS_TYP typ
  )
{
  INT              ii,j;
  INT              nd1,nd2;
  DOUBLE           det;
  DOUBLE           area;
  const DOUBLE     tol = 1.0E-06;
  const DOUBLE     nitnmax = 10;
  
#ifdef DEBUG 
  dstrc_enter("xfem_polygon_compGP");
#endif
/*----------------------------------------------------------------------*/

  switch (typ)
  {
      case tri3:
        /* initialize */
        GP[0] = 0.0;
        GP[1] = 0.0;
        /* compute */
        xyze_subp[0][0] = xtarget[0];
        xyze_subp[1][0] = xtarget[1];
        for (j=0; j<2; j++)
        {
          nd1 = node_to_node_tri[0][j];
          nd2 = node_to_node_tri[1][j]; 
          xyze_subp[0][1] = xyze[0][nd1-1];
          xyze_subp[1][1] = xyze[1][nd1-1];
          xyze_subp[0][2] = xyze[0][nd2-1];
          xyze_subp[1][2] = xyze[1][nd2-1];
          area = xfem_polygon_area_tri();
          GP[j] = area/atri;
        }
        break;
      case quad4:
        /* initial guess */
        GP[0] = 0.0;
        GP[1] = 0.0;
        /* START Newton iteration */
        for (ii=0; ii<nitnmax; ii++)
        {
          /* evaluate shape functions at GP */
          xfem_polygon_funct();
          /* evaluate derivative of shape functions at GP */
          xfem_polygon_deriv();
          /* compute residuum */
          xfem_polygon_resNewton();
          /* check norm of the residuum */
          if (sqrt(R[0]*R[0]+R[1]*R[1])<tol) return;
          /* compute tangent */
          xfem_polygon_tanNewton();
          /* compute determinant of A */
          det = A[0][0]*A[1][1] - A[1][0]*A[0][1];
          /* update GP */
          GP[0] = GP[0] - ( A[1][1]*R[0] - A[0][1]*R[1])/det;
          GP[1] = GP[1] - (-A[1][0]*R[0] + A[0][0]*R[1])/det;     
        }
        /* END Newton iteration */
        break;        
      default:
        dserror("typ unknown!");
  }

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of xfem_polygon_compGP */



/*!----------------------------------------------------------------------
\brief evaluate shape functions at the trial Gauss point

<pre>                                                            irhan 05/04
evaluate shape functions at the trial Gauss point.
</pre>

*----------------------------------------------------------------------*/
void xfem_polygon_funct()
{
  INT     i;
  
#ifdef DEBUG 
  dstrc_enter("xfem_polygon_funct");
#endif
/*----------------------------------------------------------------------*/

  /* initialize */
  for (i=0; i<4; i++)
    N[i] = 0.0;

  /* loop */
  for (i=0; i<4; i++)
    N[i] = (1 + GP[0]*rsI[0][i])*(1 + GP[1]*rsI[1][i])/4.0;

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of xfem_polygon_funct */



/*!----------------------------------------------------------------------
\brief evaluate derivative of the shape functions at the trial Gauss point

<pre>                                                            irhan 05/04
evaluate derivative of the shape functions at the trial Gauss point.
</pre>

*----------------------------------------------------------------------*/
void xfem_polygon_deriv()
{
  INT     i,j;
  INT     ind1,ind2;
  INT     ctrl[2][2] = {
                        {0,1},
                        {1,0}
                       };
#ifdef DEBUG 
  dstrc_enter("xfem_polygon_deriv");
#endif
/*----------------------------------------------------------------------*/

  /* initialize */
  for (i=0; i<2; i++)
    for (j=0; j<4; j++)
      gradN[i][j] = 0.0;

  /* loop */
  for (i=0; i<4; i++)
  {
    for (j=0; j<2; j++)
    {
      ind1 = ctrl[0][j];
      ind2 = ctrl[1][j];
      gradN[j][i] = rsI[ind1][i]*(1 + GP[ind2]*rsI[ind2][i])/4.0;
    }
  }
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of xfem_polygon_deriv */



/*!----------------------------------------------------------------------
\brief evaluate the residuum

<pre>                                                            irhan 05/04
evaluate the residuum.
</pre>

*----------------------------------------------------------------------*/
void xfem_polygon_resNewton()
{
  INT     i,j;

#ifdef DEBUG 
  dstrc_enter("xfem_polygon_resNewton");
#endif
/*----------------------------------------------------------------------*/

  /* initialize */
  for (i=0; i<2; i++)
    R[i] = 0.0;

  /* compute residuum */
  for (i=0; i<2; i++)
    {
      R[i] = xtarget[i];
      for (j=0; j<4; j++)
      {
        R[i] -= xyze[i][j]*N[j];
      }
    }
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of xfem_polygon_resNewton */



/*!----------------------------------------------------------------------
\brief evaluate the tangent

<pre>                                                            irhan 05/04
evaluate the tangent.
</pre>

*----------------------------------------------------------------------*/
void xfem_polygon_tanNewton()
{
  INT     i,j,k;

#ifdef DEBUG 
  dstrc_enter("xfem_polygon_tanNewton");
#endif
/*----------------------------------------------------------------------*/

  /* initialize */
  for (i=0; i<2; i++)
    for (j=0; j<2; j++)
      A[i][j] = 0.0;

  /* compute tangent */
  for (i=0; i<2; i++)
    for (j=0; j<2; j++)
      for (k=0; k<4; k++)
        A[i][j] -= xyze[i][k]*gradN[j][k];
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of xfem_polygon_tanNewton */



/*!----------------------------------------------------------------------
\brief write polygon data into files

<pre>                                                            irhan 05/04
write polygon data into files.
</pre>

*----------------------------------------------------------------------*/
void xfem_polygon_write(
  ELEMENT *myls2
  )
{
  INT    i,j,k;
  
#ifdef DEBUG 
  dstrc_enter("xfem_polygon_write");
#endif
/*----------------------------------------------------------------------*/

  /* write element nodal coordinate data */
  for (i=0; i<2; i++)
  {
    for(j=0; j<myls2->numnp; j++)
    {
      fprintf(f01,"%10.5E  ",xyze[i][j]);
    }
    fprintf(f01,"\n");    
  }
  fprintf(f01,"\n");

  /* write interface coordinate data */  
  for (i=0; i<2; i++)
  {
    for(j=0; j<3; j++)
    {
      fprintf(f02,"%10.5E  ",xyze_int[i][j]);
    }
    fprintf(f02,"\n");    
  }
  fprintf(f02,"\n");

  /* write polygon coordinate data */
  for (i=0; i<2; i++)
  {
    for (j=0; j<2; j++)
    {
      for(k=0; k<7; k++)
      {
        fprintf(f03,"%10.5E  ",polygonxy[i][j][k]);
      }
      fprintf(f03,"\n");    
    }
    fprintf(f03,"\n");
  }
  fprintf(f03,"\n");

  /* write polygon GP coordinate data */
  for (i=0; i<2; i++)
  {
    for (j=0; j<2; j++)
    {
      for(k=0; k<7; k++)
      {
        fprintf(f04,"%10.5E  ",polygonGP[i][j][k]);
      }
      fprintf(f04,"\n");    
    }
    fprintf(f04,"\n");
  }
  fprintf(f04,"\n");

  /* write ind */
  for (i=0; i<2; i++)
  {
    fprintf(f05,"%2d  ",ind[i]);
  }
  fprintf(f05,"\n");    
    
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of xfem_polygon_write */



/*!----------------------------------------------------------------------
\brief open files to write polygon information

<pre>                                                            irhan 05/04
open files to write polygon information.
</pre>

*----------------------------------------------------------------------*/
void xfem_polygon_open()
{
#ifdef DEBUG 
  dstrc_enter("xfem_polygon_open");
#endif
/*----------------------------------------------------------------------*/

  /* open files to write polygon information */
  f01 = fopen("src/ls/to_matlab/xfem_to_matlab/xfem_to_matlab_xyze","w");
  f02 = fopen("src/ls/to_matlab/xfem_to_matlab/xfem_to_matlab_xyze_int","w");
  f03 = fopen("src/ls/to_matlab/xfem_to_matlab/xfem_to_matlab_polygonxy","w");
  f04 = fopen("src/ls/to_matlab/xfem_to_matlab/xfem_to_matlab_polygonGP","w");
  f05 = fopen("src/ls/to_matlab/xfem_to_matlab/xfem_to_matlab_ind","w");
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of xfem_polygon_open */



/*!----------------------------------------------------------------------
\brief close files

<pre>                                                            irhan 05/04
close files.
</pre>

*----------------------------------------------------------------------*/
void xfem_polygon_close()
{
#ifdef DEBUG 
  dstrc_enter("xfem_polygon_close");
#endif
/*----------------------------------------------------------------------*/

  /* close the files */
  fclose(f01);
  fclose(f02);
  fclose(f03);
  fclose(f04);
  fclose(f05);
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of xfem_polygon_finalize */



/*!----------------------------------------------------------------------
\brief compute area of the quadrilateral element

<pre>                                                            irhan 05/04
compute area of the quadrilateral element.
</pre>

*----------------------------------------------------------------------*/
void xfem_polygon_area_rect()
{
  INT     i;
  
#ifdef DEBUG 
  dstrc_enter("xfem_polygon_area_rect");
#endif
/*----------------------------------------------------------------------*/

  /* initialize */
  arect = 0.0;
  /* add up area of the triangles */
  for (i=0; i<2; i++)
    arect += xfem_polygon_area_subtri(i);
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of xfem_polygon_area_rect */



/*!----------------------------------------------------------------------
\brief compute area of the triangular domain

<pre>                                                            irhan 05/04
compute area of the triangular domain.
</pre>

*----------------------------------------------------------------------*/
DOUBLE xfem_polygon_area_subtri(
  INT i
  )
{
  INT     j;
  INT     nd;
  
#ifdef DEBUG 
  dstrc_enter("xfem_polygon_area_subtri");
#endif
/*----------------------------------------------------------------------*/

  if (i==-1)
  {
    for (j=0; j<3; j++)
    {
      xyze_subp[0][j] = xyze[0][j];
      xyze_subp[1][j] = xyze[1][j];
    }
  }
  else
  {
    for (j=0; j<3; j++)
    {
      nd = conr[i][j] - 1;
      xyze_subp[0][j] = xyze[0][nd];
      xyze_subp[1][j] = xyze[1][nd];
    }
  }

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return xfem_polygon_area_tri();
} /* end of xfem_polygon_area_subtri */



/*!----------------------------------------------------------------------
\brief compute area of the subtriangle

<pre>                                                            irhan 05/04
compute area of the subtriangle.
</pre>

*----------------------------------------------------------------------*/
DOUBLE xfem_polygon_area_tri()
{
  DOUBLE     x21,y21,x31,y31;
  
#ifdef DEBUG 
  dstrc_enter("xfem_polygon_area_tri");
#endif
/*----------------------------------------------------------------------*/

  /* compute area of the triangle */
  x21 = xyze_subp[0][1] - xyze_subp[0][0];
  y21 = xyze_subp[1][1] - xyze_subp[1][0];
  x31 = xyze_subp[0][2] - xyze_subp[0][0];
  y31 = xyze_subp[1][2] - xyze_subp[1][0];

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return 0.5*(x21*y31 - y21*x31);
} /* end of xfem_polygon_area_tri */



/*!----------------------------------------------------------------------
\brief compute the global coordinates of the target point (mid-point of the
triangular domain)

<pre>                                                            irhan 05/04
compute the global coordinates of the target point (mid-point of the
triangular domain).
</pre>

*----------------------------------------------------------------------*/
void xfem_polygon_target_tri()
{
  INT     i;
  
#ifdef DEBUG 
  dstrc_enter("xfem_polygon_target_tri");
#endif
/*----------------------------------------------------------------------*/

  /* initialize xtarget */
  xtarget[0] = 0.0;
  xtarget[1] = 0.0;
  for (i=0; i<3; i++)
  {
    xtarget[0] += xyze_subp[0][i]/3.0;
    xtarget[1] += xyze_subp[1][i]/3.0;        
  }

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of xfem_polygon_target_tri */



/*!----------------------------------------------------------------------
\brief compute the global coordinates of the target point (mid-point of the
subtriangle)

<pre>                                                            irhan 05/04
compute the global coordinates of the target point (mid-point of the
subtriangle).
</pre>

*----------------------------------------------------------------------*/
void xfem_polygon_target_subtri(
  INT i
  )
{
  INT     j;
  INT     nd;
  
#ifdef DEBUG 
  dstrc_enter("xfem_polygon_target_subtri");
#endif
/*----------------------------------------------------------------------*/

  /* initialize xtarget */
  xtarget[0] = 0.0;
  xtarget[1] = 0.0;
  for (j=0; j<3; j++)
  {
    nd = conr[i][j]-1;
    xtarget[0] += xyze[0][nd]/3.0;
    xtarget[1] += xyze[1][nd]/3.0;        
  }

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of xfem_polygon_target_subtri */



/*!----------------------------------------------------------------------
\brief construct the material connectivity array

<pre>                                                            irhan 05/04
construct the material connectivity array. this array stores the information
about the material type to be used at each Gauss point.
</pre>

*----------------------------------------------------------------------*/
void xfem_polygon_mat(
  ELEMENT* myls2
  )
{
  INT        i,j;
  INT        sn,cn;
  DOUBLE     lset01[4];
  INT        indmat[2][7] =
                            {
                              {2,2,2,1,1,1,1},
                              {1,1,1,2,2,2,2}
                            };
  
#ifdef DEBUG 
  dstrc_enter("xfem_polygon_mat");
#endif
/*----------------------------------------------------------------------*/

  /* access to the nodal values of level set profile */
  ls2_calset1(myls2,1,lset01);
  switch (myls2->distyp)
  {
      case tri3 :
        /* check solitary node */
        sn = polycon[0][0]-1;
        if (lset01[sn]<0.0)
        {
          for (j=0; j<7; j++)
            polygonmat[0][j] = indmat[0][j];
        }
        else
        {
          for (j=0; j<7; j++)
            polygonmat[0][j] = indmat[1][j];      
        }
        break;
      case quad4 :
        /* loop over the triangles */
        for (i=0; i<2; i++)
        {
          if (ind[i]==0)
          {
            /* set proper corner node */
            cn = 2*i;
            if (lset01[cn]>0.0) polygonmat[i][0] = 1;
            else polygonmat[i][0] = 2;
            continue;
          }
          /* check solitary node */
          sn = polycon[i][0]-1;
          if (lset01[sn]<0.0)
          {
            for (j=0; j<7; j++)
              polygonmat[i][j] = indmat[0][j];
          }
          else
          {
            for (j=0; j<7; j++)
              polygonmat[i][j] = indmat[1][j];      
          }
        }
        break;
      default:
        dserror("typ unknown!");        
  }
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of xfem_polygon_mat */
/*! @} (documentation module close)*/
#endif

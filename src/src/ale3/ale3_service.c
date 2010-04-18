/*-----------------------------------------------------------------------*/
/*!
\file
\brief a lot of functions for the 3D ale element


<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

 */
/*-----------------------------------------------------------------------*/

#ifndef CCADISCRET
#ifdef D_ALE
#include "../headers/standardtypes.h"
#include "ale3.h"


/*!
\addtogroup Ale
*//*! @{ (documentation module open)*/


/*----------------------------------------------------------------------*
 |                                                         mn 06/02     |
 | vector of material laws                                              |
 | defined in global_control.c
*----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;


static DOUBLE Q12 = ONE/TWO;
static DOUBLE Q14 = ONE/FOUR;
static DOUBLE Q18 = ONE/EIGHT;


/*-----------------------------------------------------------------------*/
/*!
  \brief calculates the shape functions and derivatives

  This routine calcuates the shape functions and their derivatives at a
  point r,s,t for 3D hex and tet ale-elements.

  \param funct   DOUBLE[]   (o)   shape functions
  \param deriv   DOUBLE[][] (o)   the derivatives of the shape functions
  \param r       DOUBLE     (i)   r coordinate
  \param s       DOUBLE     (i)   s coordinate
  \param t       DOUBLE     (i)   t coordinate
  \param typ     DIS_TYP    (i)   type of dicretization
  \param option  INT        (i)   option == 0 : only functions,
                                  option == 1 : also derivatives

  \return void
  \sa calling: ---; called by: ale3_static_ke

  \author mn
  \date   06/02

 */
/*-----------------------------------------------------------------------*/
void ale3_funct_deriv(
    DOUBLE      funct[MAXNOD_ALE3],
    DOUBLE      deriv[3][MAXNOD_ALE3],
    DOUBLE      r,
    DOUBLE      s,
    DOUBLE      t,
    DIS_TYP     typ,
    INT         option
    )
{

  const DOUBLE   q18 = 1.0/8.0;
  DOUBLE         rp,sp,tp,rm,sm,tm,rrm,ssm,ttm;
  DOUBLE         t1,t2,t3,t4;

#ifdef DEBUG
  dstrc_enter("ale3_funct_deriv");
#endif

  /* if option ==0 only funtion evaluation, if option==1 also derivatives */
  rp  = 1.0+r;
  rm  = 1.0-r;
  sp  = 1.0+s;
  sm  = 1.0-s;
  tp  = 1.0+t;
  tm  = 1.0-t;
  rrm = 1.0-r*r;
  ssm = 1.0-s*s;
  ttm = 1.0-t*t;

  switch(typ)
  {
    /*  LINEAR SHAPE FUNCTIONS  SERENDIPITY (8-NODED ELEMENT) */
    case hex8:
      funct[0]=q18*rp*sm*tm;
      funct[1]=q18*rp*sp*tm;
      funct[2]=q18*rm*sp*tm;
      funct[3]=q18*rm*sm*tm;
      funct[4]=q18*rp*sm*tp;
      funct[5]=q18*rp*sp*tp;
      funct[6]=q18*rm*sp*tp;
      funct[7]=q18*rm*sm*tp;

      if (option==1)  /* check for derivative evaluation */
      {
        deriv[0][0]= q18*sm*tm  ;
        deriv[0][1]= q18*sp*tm  ;
        deriv[0][2]=-deriv[0][1];
        deriv[0][3]=-deriv[0][0];
        deriv[0][4]= q18*sm*tp  ;
        deriv[0][5]= q18*sp*tp  ;
        deriv[0][6]=-deriv[0][5];
        deriv[0][7]=-deriv[0][4];
        deriv[1][0]=-q18*tm*rp  ;
        deriv[1][1]=-deriv[1][0];
        deriv[1][2]= q18*tm*rm  ;
        deriv[1][3]=-deriv[1][2];
        deriv[1][4]=-q18*tp*rp  ;
        deriv[1][5]=-deriv[1][4];
        deriv[1][6]= q18*tp*rm  ;
        deriv[1][7]=-deriv[1][6];
        deriv[2][0]=-q18*rp*sm  ;
        deriv[2][1]=-q18*rp*sp  ;
        deriv[2][2]=-q18*rm*sp  ;
        deriv[2][3]=-q18*rm*sm  ;
        deriv[2][4]=-deriv[2][0];
        deriv[2][5]=-deriv[2][1];
        deriv[2][6]=-deriv[2][2];
        deriv[2][7]=-deriv[2][3];
      }
      break;


    case hex20: /* QUADRATIC shape functions and their natural derivatives
                   without central nodes */


      /*   Shape functions and their derivatives for a 20 noded hexaedron
       *   ==============================================================
       *
       *   Numbering of the nodes:
       *   -----------------------
       *   - this is the numbering used in GiD!!
       *   - the numbering of the brick1 element is different!!
       *
       *
       *
       *                          ^ t          / s
       *                          |           /
       *                          |          /
       *                    8     |   19    /   7
       *                    o-----|---o---------o
       *                   /|     |       /    /|
       *                  / |     |      /    / |
       *                 /  |     |     /    /  |
       *              20o   |     |    /    o18 |
       *               /  16o     |   /    /    o15
       *              /     |     |  /    /     |
       *             /      |  17 | /  6 /      |
       *          5 o---------o---------o       |
       *            |       |     *-----|---------------->
       *            |       o---------o-|-------o         r
       *            |      / 4       11 |      /3
       *            |     /             |     /
       *          13o    /              o14  /
       *            | 12o               |   o10
       *            |  /                |  /
       *            | /                 | /
       *            |/                  |/
       *            o---------o---------o
       *           1         9         2
       */



      /* form basic values */
      rp  = ONE+r;
      rm  = ONE-r;
      sp  = ONE+s;
      sm  = ONE-s;
      tp  = ONE+t;
      tm  = ONE-t;
      rrm = ONE-r*r;
      ssm = ONE-s*s;
      ttm = ONE-t*t;

      funct[ 0] = -Q18*rm*sm*tm*(TWO+r+s+t);
      funct[ 1] = -Q18*rp*sm*tm*(TWO-r+s+t);
      funct[ 2] = -Q18*rp*sp*tm*(TWO-r-s+t);
      funct[ 3] = -Q18*rm*sp*tm*(TWO+r-s+t);
      funct[ 4] = -Q18*rm*sm*tp*(TWO+r+s-t);
      funct[ 5] = -Q18*rp*sm*tp*(TWO-r+s-t);
      funct[ 6] = -Q18*rp*sp*tp*(TWO-r-s-t);
      funct[ 7] = -Q18*rm*sp*tp*(TWO+r-s-t);

      funct[ 8] =  Q14*rrm*sm*tm;
      funct[ 9] =  Q14*rp*ssm*tm;
      funct[10] =  Q14*rrm*sp*tm;
      funct[11] =  Q14*rm*ssm*tm;

      funct[12] =  Q14*rm*sm*ttm;
      funct[13] =  Q14*rp*sm*ttm;
      funct[14] =  Q14*rp*sp*ttm;
      funct[15] =  Q14*rm*sp*ttm;

      funct[16] =  Q14*rrm*sm*tp;
      funct[17] =  Q14*rp*ssm*tp;
      funct[18] =  Q14*rrm*sp*tp;
      funct[19] =  Q14*rm*ssm*tp;


      /* first derivative evaluation */
      deriv[0][ 0] =  Q18*   sm*tm*(TWO*r+s+t+ONE);
      deriv[1][ 0] =  Q18*rm*   tm*(r+TWO*s+t+ONE);
      deriv[2][ 0] =  Q18*rm*sm*   (r+s+TWO*t+ONE);

      deriv[0][ 1] =  Q18*   sm*tm*(TWO*r-s-t-ONE);
      deriv[1][ 1] = -Q18*rp*   tm*(r-TWO*s-t-ONE);
      deriv[2][ 1] = -Q18*rp*sm*   (r-s-TWO*t-ONE);

      deriv[0][ 2] =  Q18*   sp*tm*(TWO*r+s-t-ONE);
      deriv[1][ 2] =  Q18*rp*   tm*(r+TWO*s-t-ONE);
      deriv[2][ 2] = -Q18*rp*sp*   (r+s-TWO*t-ONE);

      deriv[0][ 3] =  Q18*   sp*tm*(TWO*r-s+t+ONE);
      deriv[1][ 3] = -Q18*rm*   tm*(r-TWO*s+t+ONE);
      deriv[2][ 3] =  Q18*rm*sp*   (r-s+TWO*t+ONE);

      deriv[0][ 4] =  Q18*   sm*tp*(TWO*r+s-t+ONE);
      deriv[1][ 4] =  Q18*rm*   tp*(r+TWO*s-t+ONE);
      deriv[2][ 4] = -Q18*rm*sm*   (r+s-TWO*t+ONE);

      deriv[0][ 5] =  Q18*   sm*tp*(TWO*r-s+t-ONE);
      deriv[1][ 5] = -Q18*rp*   tp*(r-TWO*s+t-ONE);
      deriv[2][ 5] =  Q18*rp*sm*   (r-s+TWO*t-ONE);

      deriv[0][ 6] =  Q18*   sp*tp*(TWO*r+s+t-ONE);
      deriv[1][ 6] =  Q18*rp*   tp*(r+TWO*s+t-ONE);
      deriv[2][ 6] =  Q18*rp*sp*   (r+s+TWO*t-ONE);

      deriv[0][ 7] =  Q18*   sp*tp*(TWO*r-s-t+ONE);
      deriv[1][ 7] = -Q18*rm*   tp*(r-TWO*s-t+ONE);
      deriv[2][ 7] = -Q18*rm*sp*   (r-s-TWO*t+ONE);


      deriv[0][ 8] = -Q12*r*sm*tm;
      deriv[1][ 8] = -Q14*rm*rp*tm;
      deriv[2][ 8] = -Q14*rm*rp*sm;

      deriv[0][ 9] =  Q14*sm*sp*tm;
      deriv[1][ 9] = -Q12*rp*s*tm;
      deriv[2][ 9] = -Q14*rp*sm*sp;

      deriv[0][10] = -Q12*r*sp*tm;
      deriv[1][10] =  Q14*rm*rp*tm;
      deriv[2][10] = -Q14*rm*rp*sp;

      deriv[0][11] = -Q14*sm*sp*tm;
      deriv[1][11] = -Q12*s*tm*rm;
      deriv[2][11] = -Q14*sm*sp*rm;


      deriv[0][12] = -Q14*sm*tm*tp;
      deriv[1][12] = -Q14*rm*tm*tp;
      deriv[2][12] = -Q12*t*rm*sm;

      deriv[0][13] =  Q14*sm*tm*tp;
      deriv[1][13] = -Q14*rp*tm*tp;
      deriv[2][13] = -Q12*t*rp*sm;

      deriv[0][14] =  Q14*sp*tm*tp;
      deriv[1][14] =  Q14*rp*tm*tp;
      deriv[2][14] = -Q12*t*rp*sp;

      deriv[0][15] = -Q14*sp*tm*tp;
      deriv[1][15] =  Q14*rm*tm*tp;
      deriv[2][15] = -Q12*t*rm*sp;


      deriv[0][16] = -Q12*r*sm*tp;
      deriv[1][16] = -Q14*rm*rp*tp;
      deriv[2][16] =  Q14*rm*rp*sm;

      deriv[0][17] =  Q14*sm*sp*tp;
      deriv[1][17] = -Q12*s*tp*rp;
      deriv[2][17] =  Q14*sm*sp*rp;

      deriv[0][18] = -Q12*r*sp*tp;
      deriv[1][18] =  Q14*rm*rp*tp;
      deriv[2][18] =  Q14*rm*rp*sp;

      deriv[0][19] = -Q14*sm*sp*tp;
      deriv[1][19] = -Q12*s*tp*rm;
      deriv[2][19] =  Q14*sm*sp*rm;
      break;

  case hex27: /* QUADRATIC shape functions and their natural derivatives
                 with central nodes                         ----*/
  {
    /* form basic values */
    DOUBLE drm1,dr00,drp1,dsm1,ds00,dsp1,dtm1,dt00,dtp1;
    DOUBLE rm1,r00,rp1,sm1,s00,sp1,tm1,t00,tp1;

    rm1=Q12*r*(r - ONE);
    r00=(ONE - r*r);
    rp1=Q12*r*(r + ONE);
    sm1=Q12*s*(s - ONE);
    s00=(ONE - s*s);
    sp1=Q12*s*(s + ONE);
    tm1=Q12*t*(t - ONE);
    t00=(ONE - t*t);
    tp1=Q12*t*(t + ONE);

    drm1 = r - Q12;
    dr00 = -TWO * r;
    drp1 = r + Q12;
    dsm1 = s - Q12;
    ds00 = -TWO * s;
    dsp1 = s + Q12;
    dtm1 = t - Q12;
    dt00 = -TWO * t;
    dtp1 = t + Q12;

    funct[0] = rp1*sp1*tp1;
    funct[1] = sm1*rp1*tp1;
    funct[2] = rm1*sm1*tp1;
    funct[3] = rm1*sp1*tp1;
    funct[4] = tm1*rp1*sp1;
    funct[5] = sm1*tm1*rp1;
    funct[6] = rm1*sm1*tm1;
    funct[7] = rm1*tm1*sp1;
    funct[8] = s00*rp1*tp1;
    funct[9] = r00*sm1*tp1;
    funct[10] = s00*rm1*tp1;
    funct[11] = r00*sp1*tp1;
    funct[12] = t00*rp1*sp1;
    funct[13] = t00*sm1*rp1;
    funct[14] = t00*rm1*sm1;
    funct[15] = t00*rm1*sp1;
    funct[16] = s00*tm1*rp1;
    funct[17] = r00*sm1*tm1;
    funct[18] = s00*rm1*tm1;
    funct[19] = r00*tm1*sp1;
    funct[20] = r00*s00*tp1;
    funct[21] = s00*t00*rp1;
    funct[22] = r00*t00*sm1;
    funct[23] = s00*t00*rm1;
    funct[24] = r00*t00*sp1;
    funct[25] = r00*s00*tm1;
    funct[26] = r00*s00*t00;

    if (option==1) /* --> first derivative evaluation */
    {
      deriv[0][0] = sp1*tp1*drp1;
      deriv[0][1] = sm1*tp1*drp1;
      deriv[0][2] = sm1*tp1*drm1;
      deriv[0][3] = sp1*tp1*drm1;
      deriv[0][4] = tm1*sp1*drp1;
      deriv[0][5] = sm1*tm1*drp1;
      deriv[0][6] = sm1*tm1*drm1;
      deriv[0][7] = tm1*sp1*drm1;
      deriv[0][8] = s00*tp1*drp1;
      deriv[0][9] = sm1*tp1*dr00;
      deriv[0][10] = s00*tp1*drm1;
      deriv[0][11] = sp1*tp1*dr00;
      deriv[0][12] = t00*sp1*drp1;
      deriv[0][13] = t00*sm1*drp1;
      deriv[0][14] = t00*sm1*drm1;
      deriv[0][15] = t00*sp1*drm1;
      deriv[0][16] = s00*tm1*drp1;
      deriv[0][17] = sm1*tm1*dr00;
      deriv[0][18] = s00*tm1*drm1;
      deriv[0][19] = tm1*sp1*dr00;
      deriv[0][20] = s00*tp1*dr00;
      deriv[0][21] = s00*t00*drp1;
      deriv[0][22] = t00*sm1*dr00;
      deriv[0][23] = s00*t00*drm1;
      deriv[0][24] = t00*sp1*dr00;
      deriv[0][25] = s00*tm1*dr00;
      deriv[0][26] = s00*t00*dr00;

      deriv[1][0] = rp1*tp1*dsp1;
      deriv[1][1] = rp1*tp1*dsm1;
      deriv[1][2] = rm1*tp1*dsm1;
      deriv[1][3] = rm1*tp1*dsp1;
      deriv[1][4] = tm1*rp1*dsp1;
      deriv[1][5] = tm1*rp1*dsm1;
      deriv[1][6] = rm1*tm1*dsm1;
      deriv[1][7] = rm1*tm1*dsp1;
      deriv[1][8] = rp1*tp1*ds00;
      deriv[1][9] = r00*tp1*dsm1;
      deriv[1][10] = rm1*tp1*ds00;
      deriv[1][11] = r00*tp1*dsp1;
      deriv[1][12] = t00*rp1*dsp1;
      deriv[1][13] = t00*rp1*dsm1;
      deriv[1][14] = t00*rm1*dsm1;
      deriv[1][15] = t00*rm1*dsp1;
      deriv[1][16] = tm1*rp1*ds00;
      deriv[1][17] = r00*tm1*dsm1;
      deriv[1][18] = rm1*tm1*ds00;
      deriv[1][19] = r00*tm1*dsp1;
      deriv[1][20] = r00*tp1*ds00;
      deriv[1][21] = t00*rp1*ds00;
      deriv[1][22] = r00*t00*dsm1;
      deriv[1][23] = t00*rm1*ds00;
      deriv[1][24] = r00*t00*dsp1;
      deriv[1][25] = r00*tm1*ds00;
      deriv[1][26] = r00*t00*ds00;

      deriv[2][0] = rp1*sp1*dtp1;
      deriv[2][1] = sm1*rp1*dtp1;
      deriv[2][2] = rm1*sm1*dtp1;
      deriv[2][3] = rm1*sp1*dtp1;
      deriv[2][4] = rp1*sp1*dtm1;
      deriv[2][5] = sm1*rp1*dtm1;
      deriv[2][6] = rm1*sm1*dtm1;
      deriv[2][7] = rm1*sp1*dtm1;
      deriv[2][8] = s00*rp1*dtp1;
      deriv[2][9] = r00*sm1*dtp1;
      deriv[2][10] = s00*rm1*dtp1;
      deriv[2][11] = r00*sp1*dtp1;
      deriv[2][12] = rp1*sp1*dt00;
      deriv[2][13] = sm1*rp1*dt00;
      deriv[2][14] = rm1*sm1*dt00;
      deriv[2][15] = rm1*sp1*dt00;
      deriv[2][16] = s00*rp1*dtm1;
      deriv[2][17] = r00*sm1*dtm1;
      deriv[2][18] = s00*rm1*dtm1;
      deriv[2][19] = r00*sp1*dtm1;
      deriv[2][20] = r00*s00*dtp1;
      deriv[2][21] = s00*rp1*dt00;
      deriv[2][22] = r00*sm1*dt00;
      deriv[2][23] = s00*rm1*dt00;
      deriv[2][24] = r00*sp1*dt00;
      deriv[2][25] = r00*s00*dtm1;
      deriv[2][26] = r00*s00*dt00;
    }
    break;
  }

    case tet4: /* LINEAR SHAPE FUNCTIONS */
/*
      t1=r;
      t2=s;
      t3=t;
      t4=ONE-r-s-t;
*/
   t1=ONE-r-s-t;
   t2=r;
   t3=s;
   t4=t;

      funct[0]= t1;
      funct[1]= t2;
      funct[2]= t3;
      funct[3]= t4;

      if(option==1) /* first derivative evaluation */
      {
        /*
        deriv[0][0]= ONE;
        deriv[0][1]= ZERO;
        deriv[0][2]= ZERO;
        deriv[0][3]=-ONE;

        deriv[1][0]= ZERO;
        deriv[1][1]= ONE;
        deriv[1][2]= ZERO;
        deriv[1][3]=-ONE;

        deriv[2][0]= ZERO;
        deriv[2][1]= ZERO;
        deriv[2][2]= ONE;
        deriv[2][3]=-ONE;
        */
      deriv[0][0]=-ONE;
      deriv[0][1]= ONE;
      deriv[0][2]= ZERO;
      deriv[0][3]= ZERO;

      deriv[1][0]=-ONE;
      deriv[1][1]= ZERO;
      deriv[1][2]= ONE;
      deriv[1][3]= ZERO;

      deriv[2][0]=-ONE;
      deriv[2][1]= ZERO;
      deriv[2][2]= ZERO;
      deriv[2][3]= ONE;

      } /* endif (option==1) */
      break;

    case tet10: /*  QUADRATIC SHAPE FUNCTIONS */

      dserror("shape functions for tet10 not yet implemented \n");
      /* form basic values */
#if 0
      t1=r;
      t2=s;
      t3=t;
      t4=ONE-r-s-t;

      funct[0] =  ;
      funct[1] =  ;
      funct[2] = ;
      funct[3] = ;
      funct[4] = ;
      funct[5] = ;
      funct[6] = ;
      funct[7] = ;
      funct[8] = ;
      funct[9] = ;


      if(option==1) /* first derivative evaluation */
      {
        deriv[0][0] = ;
        deriv[1][0] = ;
        deriv[2][0] = ;

        deriv[0][1] = ;
        deriv[1][1] = ;
        deriv[2][1] = ;

        deriv[0][2] = ;
        deriv[1][2] = ;
        deriv[2][2] = ;

        deriv[0][3] = ;
        deriv[1][3] = ;
        deriv[2][3] = ;

        deriv[0][4] = ;
        deriv[1][4] = ;
        deriv[2][4] = ;

        deriv[0][5] = ;
        deriv[1][5] = ;
        deriv[2][5] = ;

        deriv[0][6] = ;
        deriv[1][6] = ;
        deriv[2][6] = ;

        deriv[0][7] = ;
        deriv[1][7] = ;
        deriv[2][7] = ;

        deriv[0][8] = ;
        deriv[1][8] = ;
        deriv[2][8] = ;

        deriv[0][9] = ;
        deriv[1][9] = ;
        deriv[2][9] = ;
      }
      break;
#endif

    default:
      dserror("unknown typ of interpolation");
      break;
  } /* end of switch typ */

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of ale3_funct_deriv */




/*-----------------------------------------------------------------------*/
/*!
  \brief calculate operator matrix at point r,s,t

  This routine calcuates the operator matrix b at the given point r,s,t
  for an 3D-ale-element.

  \param b         DOUBLE[][]  (o)   the calculated operator matrix
  \param deriv     DOUBLE[][]  (i)   the derivatives of the shape functions
  \param xjm       DOUBLE[][]  (i)   the Jacobian matrix
  \param det       DOUBLE      (i)   the determinant of the Jacobian matrix
  \param iel       INT         (i)   number of nodes per element

  \return void
  \sa calling: ---; called by: ale3_static_ke()

  \author mn
  \date   06/02

 */
/*-----------------------------------------------------------------------*/
void ale3_bop(
    DOUBLE      b[6][3*MAXNOD_ALE3],
    DOUBLE      deriv[3][MAXNOD_ALE3],
    DOUBLE      xjm[3][3],
    DOUBLE      det,
    INT         iel
    )
{

  INT    i,node_start;
  DOUBLE dum;
  DOUBLE x1r, x2r, x3r, x1s, x2s, x3s, x1t, x2t, x3t;
  DOUBLE xi11, xi12, xi13, xi21, xi22, xi23, xi31, xi32, xi33;
  DOUBLE hr, hs, ht;
  DOUBLE h1, h2, h3;

#ifdef DEBUG
  dstrc_enter("ale3_bop");
#endif

  /* inverse of jacobian */
  x1r = xjm[0][0];
  x2r = xjm[0][1];
  x3r = xjm[0][2];
  x1s = xjm[1][0];
  x2s = xjm[1][1];
  x3s = xjm[1][2];
  x1t = xjm[2][0];
  x2t = xjm[2][1];
  x3t = xjm[2][2];

  dum=1.0/det;

  xi11=dum*(x2s*x3t - x2t*x3s);
  xi12=dum*(x3r*x2t - x2r*x3t);
  xi13=dum*(x2r*x3s - x3r*x2s);
  xi21=dum*(x3s*x1t - x3t*x1s);
  xi22=dum*(x1r*x3t - x3r*x1t);
  xi23=dum*(x3r*x1s - x1r*x3s);
  xi31=dum*(x1s*x2t - x1t*x2s);
  xi32=dum*(x2r*x1t - x1r*x2t);
  xi33=dum*(x1r*x2s - x2r*x1s);

  /* get operator b of global derivatives */
  for (i=0; i<iel; i++)
  {
    node_start = i*3;

    hr   = deriv[0][i];
    hs   = deriv[1][i];
    ht   = deriv[2][i];

    h1 = xi11*hr + xi12*hs + xi13*ht;
    h2 = xi21*hr + xi22*hs + xi23*ht;
    h3 = xi31*hr + xi32*hs + xi33*ht;

    b[0][node_start+0] = h1 ;
    b[0][node_start+1] = 0.0;
    b[0][node_start+2] = 0.0;
    b[1][node_start+0] = 0.0;
    b[1][node_start+1] = h2 ;
    b[1][node_start+2] = 0.0;
    b[2][node_start+0] = 0.0;
    b[2][node_start+1] = 0.0;
    b[2][node_start+2] = h3 ;
    b[3][node_start+0] = h2 ;
    b[3][node_start+1] = h1 ;
    b[3][node_start+2] = 0.0;
    b[4][node_start+0] = 0.0;
    b[4][node_start+1] = h3 ;
    b[4][node_start+2] = h2 ;
    b[5][node_start+0] = h3 ;
    b[5][node_start+1] = 0.0;
    b[5][node_start+2] = h1 ;
  } /* end of loop over nodes */

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of ale3_bop */




/*-----------------------------------------------------------------------*/
/*!
  \brief calculate the Jacobian matrix

  This routine calculates the Jacobian matrix  at a point r,s,t for
  a 3D ale element.

  \param deriv    DOUBLE[][] (i)   derivatives of the shape functions
  \param xjm      DOUBLE[][] (o)   the Jacobian matrix
  \param *det     DOUBLE     (i)   determinant of the Jacobian matrix
  \param *ele     ELEMENT    (i)   the element
  \param iel      INT        (i)   number of nodes of the element

  \return void
  \sa calling: ---; called by: ale3_static_ke

  \author mn
  \date   06/02

 */
/*-----------------------------------------------------------------------*/
void ale3_jaco(
    DOUBLE      deriv[3][MAXNOD_ALE3],
    DOUBLE      xjm[3][3],
    DOUBLE     *det,
    ELEMENT    *ele,
    INT         iel
    )
{

  INT i,j,l;
  DOUBLE dum;

#ifdef DEBUG
  dstrc_enter("ale3_jaco");
#endif

  /* determine jacobian at point r,s,t */
  for (i=0; i<3; i++)
  {
    for (j=0; j<3; j++)
    {
      dum=0.0;
      for (l=0; l<iel; l++)
      {
        dum += deriv[i][l]*ele->node[l]->x[j];
      }
      xjm[i][j]=dum;
    }
  }

  /* determinant of jacobian */
  *det = xjm[0][0]*xjm[1][1]*xjm[2][2]+
    xjm[0][1]*xjm[1][2]*xjm[2][0]+
    xjm[0][2]*xjm[1][0]*xjm[2][1]-
    xjm[0][2]*xjm[1][1]*xjm[2][0]-
    xjm[0][0]*xjm[1][2]*xjm[2][1]-
    xjm[0][1]*xjm[1][0]*xjm[2][2];

  if (*det<0.0) dserror("Negative Jacobian = %12.4e in element %3i",*det,ele->Id);

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of ale3_jaco */




/*-----------------------------------------------------------------------*/
/*!
  \brief coordinates and weight factors for numerical integration

  This routine  gives the coordinates and weight factors for numerical
  integration of a 3D ale element.

  \param *ele    ELEMENT    (i)   the element
  \param *data   ALE3_DATA  (o)   structure containing the coordinates and
                                  weighting factors

  \return void
  \sa calling: ---; called by: ale3_static_ke

  \author mn
  \date   06/02

 */
/*-----------------------------------------------------------------------*/
void ale3_intg(
    const ELEMENT   *ele,
    ALE3_DATA       *data
    )
{

  DOUBLE  q14, q16, q124;
  DOUBLE  palpha,pbeta;

#ifdef DEBUG
  dstrc_enter("ale3_intg");
#endif

  switch(ele->distyp)
  {
    case hex8:
    case hex20:
    case hex27:
      /*----------------------------------------------------------------------*
        |     INTEGRATION PARAMETERS FOR    H E X A H E D R A L     ELEMENTS   |
        |     GAUSS SAMPLING POINTS  AT     R/S/T-COORDINATES   RESPECTIVELY   |
        |                            AND    CORRESPONDING WEIGHTING  FACTORS   |
       *----------------------------------------------------------------------*/
      switch(ele->e.ale3->nGP[0])/* direction r */
      {
        case 3:
          data->xgpr[0] = -0.7745966692415;
          data->xgpr[1] =  0.0;
          data->xgpr[2] =  0.7745966692415;

          data->wgtr[0] =  0.5555555555556;
          data->wgtr[1] =  0.8888888888889;
          data->wgtr[2] =  0.5555555555556;
          break;
        case 2:
          data->xgpr[0] = -0.5773502691896;
          data->xgpr[1] =  0.5773502691896;

          data->wgtr[0] = 1.0            ;
          data->wgtr[1] = 1.0            ;
          break;
        case 1:
          data->xgpr[0] = 0.0;

          data->wgtr[0] = 2.0;
          break;
        default:
          dserror("unknown number of gaussian points in ale3_intg");
          break;
      }
      switch(ele->e.ale3->nGP[1])/* direction s */
      {
        case 3:
          data->xgps[0] = -0.7745966692415;
          data->xgps[1] =  0.0;
          data->xgps[2] =  0.7745966692415;

          data->wgts[0] =  0.5555555555556;
          data->wgts[1] =  0.8888888888889;
          data->wgts[2] =  0.5555555555556;
          break;
        case 2:
          data->xgps[0] = -0.5773502691896;
          data->xgps[1] =  0.5773502691896;

          data->wgts[0] = 1.0            ;
          data->wgts[1] = 1.0            ;
          break;
        case 1:
          data->xgps[0] = 0.0;

          data->wgts[0] = 2.0;
          break;
        default:
          dserror("unknown number of gaussian points in ale3_intg");
          break;
      }
      switch(ele->e.ale3->nGP[2])/* direction t */
      {
        case 3:
          data->xgpt[0] = -0.7745966692415;
          data->xgpt[1] =  0.0;
          data->xgpt[2] =  0.7745966692415;

          data->wgtt[0] =  0.5555555555556;
          data->wgtt[1] =  0.8888888888889;
          data->wgtt[2] =  0.5555555555556;
          break;
        case 2:
          data->xgpt[0] = -0.5773502691896;
          data->xgpt[1] =  0.5773502691896;

          data->wgtt[0] = 1.0            ;
          data->wgtt[1] = 1.0            ;
          break;
        case 1:
          data->xgpt[0] = 0.0;

          data->wgtt[0] = 2.0;
          break;
        default:
          dserror("unknown number of gaussian points in ale3_intg");
          break;
      }
      break; /* end case hex820 */
    case tet4:
    case tet10:
      q14 = 1.0/4.0;
      q16 = 1.0/6.0;
      q124= 1.0/24.0;
      palpha = (5.0+3.0*sqrt(5.0))/20.0;
      pbeta  = (5.0-sqrt(5.0))/20.0;

      /*----------------------------------------------------------------------*
        |     INTEGRATION PARAMETERS FOR    T E T R A H E D R A L   ELEMENTS   |
        |     GAUSS SAMPLING POINTS  AT     R/S/T-COORDINATES   RESPECTIVELY   |
        |                            AND    CORRESPONDING WEIGHTING  FACTORS   |
       *----------------------------------------------------------------------*/

      /*----------------------------------------------------------------------*
        |    GAUSS INTEGRATION         1 SAMPLING POINT, DEG.OF PRECISION 1    |
        |                              CASE 0                                  |
       *----------------------------------------------------------------------*/
      switch(ele->e.ale3->nGP[0])/* direction r */
      {
        case 1:
          data->xgpr[0]    =  q14 ;
          data->xgps[0]    =  q14 ;
          data->xgpt[0]    =  q14 ;
          data->wgtr[0]    =  q16 ;
          break;
          /*----------------------------------------------------------------*
            | GAUSS INTEGRATION    4 SAMPLING POINTS, DEG.OF PRECISION 2    |
            |                      CASE 1                                   |
           *----------------------------------------------------------------*/
        case 4:
          data->xgpr[0]    =    pbeta ;
          data->xgpr[1]    =    palpha;
          data->xgpr[2]    =    pbeta ;
          data->xgpr[3]    =    pbeta ;
          data->xgps[0]    =    pbeta ;
          data->xgps[1]    =    pbeta ;
          data->xgps[2]    =    palpha;
          data->xgps[3]    =    pbeta ;
          data->xgpt[0]    =    pbeta ;
          data->xgpt[1]    =    pbeta ;
          data->xgpt[2]    =    pbeta ;
          data->xgpt[3]    =    palpha;
          data->wgtr[0]    =    q124  ;
          data->wgtr[1]    =    q124  ;
          data->wgtr[2]    =    q124  ;
          data->wgtr[3]    =    q124  ;
          break;
        default:
          dserror("unknown number of gausian points");
          break;
      }

      break; /* end if tet410 */
    default:
      dserror("unknown typ of discretisation");
      break;
  } /* end switch */

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of ale3_intg */




/*-----------------------------------------------------------------------*/
/*!
  \brief calculates the linear elastic constitutive matrix

  This routine calculates the linear elastic isotropic constitutive
  matrix for a 3D ale element

  \param *mat   STVENANT   (i)   my material
  \param  d     DOUBLE[][] (o)   the constitutive matrix

  \return void
  \sa calling: ---; called by: ale3_static_ke

  \author mn
  \date   06/02

 */
/*-----------------------------------------------------------------------*/
void ale3_mat_linel(
    STVENANT   *mat,
    DOUBLE      d[6][6]
    )
{

  DOUBLE d1,d2,d3;
  DOUBLE ym,pv;       /* mat constants */

#ifdef DEBUG
  dstrc_enter("ale3_mat_linel");
#endif

  ym  = mat->youngs;
  pv  = mat->possionratio;

  /* evaluate basic material values */
  d1=ym*(1.0 - pv)/((1.0 + pv)*(1.0 - 2.0*pv));
  d2=ym*pv/((1.0 + pv)*(1.0 - 2.0*pv));
  d3=ym/((1.0 + pv)*2.0);

  /* set values in material-matrix */
  d[0][0]=d1;
  d[0][1]=d2;
  d[0][2]=d2;
  d[0][3]=0.0;
  d[0][4]=0.0;
  d[0][5]=0.0;

  d[1][0]=d2;
  d[1][1]=d1;
  d[1][2]=d2;
  d[1][3]=0.0;
  d[1][4]=0.0;
  d[1][5]=0.0;

  d[2][0]=d2;
  d[2][1]=d2;
  d[2][2]=d1;
  d[2][3]=0.0;
  d[2][4]=0.0;
  d[2][5]=0.0;

  d[3][0]=0.0;
  d[3][1]=0.0;
  d[3][2]=0.0;
  d[3][3]=d3;
  d[3][4]=0.0;
  d[3][5]=0.0;

  d[4][0]=0.0;
  d[4][1]=0.0;
  d[4][2]=0.0;
  d[4][3]=0.0;
  d[4][4]=d3;
  d[4][5]=0.0;

  d[5][0]=0.0;
  d[5][1]=0.0;
  d[5][2]=0.0;
  d[5][3]=0.0;
  d[5][4]=0.0;
  d[5][5]=d3;

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of ale3_mat_linel */





/*!----------------------------------------------------------------------
\brief calculates the additional stiffness matrix for hourglass stabilization

<pre>                                                              mn 06/02
This routine calcuates the additional stiffness matrix for hourglass
stabilization for a 3D element.

see also:
   (1) T. Belytschko and L.P. Bindeman:
       Assumed strain stabilization of the 8-node hexahedral element
       Comp. Meth. Appl. Mech. Eng.: 105 (1993) p. 225-260.
   (2) D.P. Flanagan and T. Belytschko:
       A uniform strain hexahedron and quadrilateral with orthogonal
       hourglass control
       Int. J. Num. Meth. Ing.: Vol. 17 (1981) p. 679-706.

</pre>
\param *ele  ELEMENT  (i)   the element
\param **s   DOUBLE   (i/o) (i) the one point quadrature matrix
                            (o) the complete, stabilized stiffness matrix
\param vol   DOUBLE   (i)   the volume of the element

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: ale3_static_ke

*----------------------------------------------------------------------*/
void ale3_hourglass(
    ELEMENT  *ele,
    DOUBLE  **s,
    DOUBLE    vol
    )
{

  MATERIAL      *actmat;

  INT            i,j,k,l;
  INT            l0=0,l1=0,l2=0,l3=0,l4=0,l5=0,l6=0,l7=0;

  DOUBLE         x[3][8];
  DOUBLE         xc[3][8];
  DOUBLE         b[3][8];

  DOUBLE         a[3][3];
  DOUBLE         ba[3];
  DOUBLE         r[3][3];

  DOUBLE         h[4][8] = {{1,1,-1,-1,-1,-1,1,1},{1,-1,-1,1,-1,1,1,-1},
                            {1,-1,1,-1,1,-1,1,-1},{-1,1,-1,1,1,-1,1,-1}};
  DOUBLE         lam[3][8] = {{-1,1,1,-1,-1,1,1,-1},
                              {-1,-1,1,1,-1,-1,1,1},
                              {-1,-1,-1,-1,1,1,1,1}};
  DOUBLE         gam[4][8];
  DOUBLE         lx[3];
  DOUBLE         hh[3][3];

  DOUBLE         gg00[8][8], gg11[8][8], gg22[8][8], gg33[8][8];
  DOUBLE         gg01[8][8], gg10[8][8], gg02[8][8], gg20[8][8];
  DOUBLE         gg12[8][8], gg21[8][8];

  DOUBLE         kstab[24][24];

  DOUBLE         c1,c2,c3;
  DOUBLE         dum;
  DOUBLE         ee,nu,mu;


#ifdef DEBUG
  dstrc_enter("ale3_hourglass");
#endif
  /* material data */
  actmat = &(mat[ele->mat-1]);
  ee = actmat->m.stvenant->youngs;
  nu = actmat->m.stvenant->possionratio;
  mu = ee / (2*(1+nu));


  /* Constants for the stabilization matrix accord. to (1) Table*/
  /* ADS */
#if 0
  c1 = 2.0/3.0;
  c2 = 2.0/9.0;
  c3 = -1.0/3.0;
#endif

  /*ASQBI */
  c1 = 1.0/(1.0 - nu);
  c2 = (1.0 + nu)/3;
  c3 = 1.0/(1.0 - nu);


  /* nodal coordinates */
  for(i=0; i<3; i++)
  {
    for(j=0; j<8; j++)
    {
      x[i][j] = ele->node[j]->x[i];
    }
  }


  /* corotational coordinate system: rotation tensor r[3][3] */
  /* accord. to (1) Appendix A, Table 9 */
  for (i=0; i<2; i++)
  {
    for (j=0; j<3; j++)
    {
      a[i][j] = 0,0;
      for (k=0; k<8; k++)
      {
        a[i][j] += lam[i][k]*x[j][k];
      }
    }
  }

  dum =(a[0][0]*a[1][0]+a[0][1]*a[1][1]+a[0][2]*a[1][2])/
    (a[0][0]*a[0][0]+a[0][1]*a[0][1]+a[0][2]*a[0][2]);

  a[1][0] = a[1][0] - dum * a[0][0];
  a[1][1] = a[1][1] - dum * a[0][1];
  a[1][2] = a[1][2] - dum * a[0][2];

  a[2][0] = a[0][1]*a[1][2] - a[0][2]*a[1][1];
  a[2][1] = a[0][2]*a[1][0] - a[0][0]*a[1][2];
  a[2][2] = a[0][0]*a[1][1] - a[0][1]*a[1][0];

  for(i = 0; i<3; i++)
  {
    ba[i] = sqrt(a[i][0]*a[i][0]+a[i][1]*a[i][1]+a[i][2]*a[i][2]);
    r[i][0] = a[i][0] / ba[i];
    r[i][1] = a[i][1] / ba[i];
    r[i][2] = a[i][2] / ba[i];
  }


  /* transforming nodal coordinates to the corotational system */
  for(i=0; i<8; i++)
  {
    for(j=0; j<3; j++)
    {
      xc[j][i] = r[j][0]*x[0][i]+r[j][1]*x[1][i]+r[j][2]*x[2][i];
    }
  }


  /* B-matrix b[3][8] accord. to (2) Appendix I, eqn (79) */
  for (i=0; i<3; i++)
  {
    j = (i+1)%3;
    k = (j+1)%3;
    for(l=0; l<8;l++)
    {
      switch(l)
      {
        case 0:
          l0=0;l1=1;l2=2;l3=3;l4=4;l5=5;l6=6;l7=7;
          break;
        case 1:
          l0=1;l1=2;l2=3;l3=0;l4=5;l5=6;l6=7;l7=4;
          break;
        case 2:
          l0=2;l1=3;l2=0;l3=1;l4=6;l5=7;l6=4;l7=5;
          break;
        case 3:
          l0=3;l1=0;l2=1;l3=2;l4=7;l5=4;l6=5;l7=6;
          break;
        case 4:
          l0=4;l1=7;l2=6;l3=5;l4=0;l5=3;l6=2;l7=1;
          break;
        case 5:
          l0=5;l1=4;l2=7;l3=6;l4=1;l5=0;l6=3;l7=2;
          break;
        case 6:
          l0=6;l1=5;l2=4;l3=7;l4=2;l5=1;l6=0;l7=3;
          break;
        case 7:
          l0=7;l1=6;l2=5;l3=4;l4=3;l5=2;l6=1;l7=0;
          break;
      }
      b[i][l0] =1/(12 * vol) * (
          xc[j][l1]*((xc[k][l5] - xc[k][l2])-(xc[k][l3] - xc[k][l4]))
          + xc[j][l2]* (xc[k][l1] - xc[k][l3])
          + xc[j][l3]*((xc[k][l2] - xc[k][l7])-(xc[k][l4] - xc[k][l1]))
          + xc[j][l4]*((xc[k][l7] - xc[k][l5])-(xc[k][l1] - xc[k][l3]))
          + xc[j][l5]* (xc[k][l4] - xc[k][l1])
          + xc[j][l7]* (xc[k][l3] - xc[k][l4]) );
    }
  }


  /* gamma vectors, accord. to (1) eqn (2.12b) */
  for(i=0; i<4; i++)
  {
    for(j=0; j<8; j++)
    {
      gam[i][j] = 0.125 * h[i][j];
      for(k=0; k<3; k++)
      {
        dum = h[i][0]*xc[k][0]+h[i][1]*xc[k][1]+h[i][2]*xc[k][2]+
          h[i][3]*xc[k][3]+h[i][4]*xc[k][4]+h[i][5]*xc[k][5]+
          h[i][6]*xc[k][6]+h[i][7]*xc[k][7];
        gam[i][j] -= 0.125 * dum * b[k][j];
      }
    }
  }


  /* lambda * x (auxiliary vector) */
  lx[0] = lam[0][0]*xc[0][0]+lam[0][1]*xc[0][1]+lam[0][2]*xc[0][2]+
    lam[0][3]*xc[0][3]+lam[0][4]*xc[0][4]+lam[0][5]*xc[0][5]+
    lam[0][6]*xc[0][6]+lam[0][7]*xc[0][7];
  lx[1] = lam[1][0]*xc[1][0]+lam[1][1]*xc[1][1]+lam[1][2]*xc[1][2]+
    lam[1][3]*xc[1][3]+lam[1][4]*xc[1][4]+lam[1][5]*xc[1][5]+
    lam[1][6]*xc[1][6]+lam[1][7]*xc[1][7];
  lx[2] = lam[2][0]*xc[2][0]+lam[2][1]*xc[2][1]+lam[2][2]*xc[2][2]+
    lam[2][3]*xc[2][3]+lam[2][4]*xc[2][4]+lam[2][5]*xc[2][5]+
    lam[2][6]*xc[2][6]+lam[2][7]*xc[2][7];


  /* H_ij, accord. to (1) eqns. (3.15d) and (3.15e) */
  hh[0][0] = 1.0/3.0 * (lx[1]*lx[2])/lx[0];
  hh[1][1] = 1.0/3.0 * (lx[2]*lx[0])/lx[1];
  hh[2][2] = 1.0/3.0 * (lx[0]*lx[1])/lx[2];
  hh[0][1] = 1.0/3.0 * lx[2];
  hh[1][0] = 1.0/3.0 * lx[2];
  hh[0][2] = 1.0/3.0 * lx[1];
  hh[2][0] = 1.0/3.0 * lx[1];
  hh[1][2] = 1.0/3.0 * lx[0];
  hh[2][1] = 1.0/3.0 * lx[0];


  /* stabalization matrix with respect to the corotational ccord. system. */
  /* rearranging orders of dofs, accord. to (1) eqns. (3.15a) to (3.15c) */
  for(i=0; i<8; i++)
  {
    for(j=0; j<8; j++)
    {
      gg00[i][j] = gam[0][i] * gam[0][j];
      gg11[i][j] = gam[1][i] * gam[1][j];
      gg22[i][j] = gam[2][i] * gam[2][j];
      gg33[i][j] = gam[3][i] * gam[3][j];
      gg01[i][j] = gam[0][i] * gam[1][j];
      gg10[i][j] = gam[1][i] * gam[0][j];
      gg02[i][j] = gam[0][i] * gam[2][j];
      gg20[i][j] = gam[2][i] * gam[0][j];
      gg12[i][j] = gam[1][i] * gam[2][j];
      gg21[i][j] = gam[2][i] * gam[1][j];

      /* kstab 00 */
      kstab[i*3][j*3]     = 2*mu* (hh[0][0]*(c1*(gg11[i][j] + gg22[i][j])
            + c2*gg33[i][j]) + 0.5 * (hh[1][1] + hh[2][2]) * gg00[i][j]);

      /* kstab 11 */
      kstab[i*3+1][j*3+1] = 2*mu* (hh[1][1]*(c1*(gg22[i][j] + gg00[i][j])
            + c2*gg33[i][j]) + 0.5 * (hh[2][2] + hh[0][0]) * gg11[i][j]);

      /* kstab 22 */
      kstab[i*3+2][j*3+2] = 2*mu* (hh[2][2]*(c1*(gg00[i][j] + gg11[i][j])
            + c2*gg33[i][j]) + 0.5 * (hh[0][0] + hh[1][1]) * gg22[i][j]);

      /* kstab 01 */
      kstab[i*3][j*3+1]   = 2*mu* (hh[0][1]*(c3*gg10[i][j]+0.5*gg01[i][j]));

      /* kstab 10 */
      kstab[i*3+1][j*3]   = 2*mu* (hh[1][0]*(c3*gg01[i][j]+0.5*gg10[i][j]));

      /* kstab 02 */
      kstab[i*3][j*3+2]   = 2*mu* (hh[0][2]*(c3*gg20[i][j]+0.5*gg02[i][j]));

      /* kstab 20 */
      kstab[i*3+2][j*3]   = 2*mu* (hh[2][0]*(c3*gg02[i][j]+0.5*gg20[i][j]));

      /* kstab 12 */
      kstab[i*3+1][j*3+2] = 2*mu* (hh[1][2]*(c3*gg21[i][j]+0.5*gg12[i][j]));

      /* kstab 21 */
      kstab[i*3+2][j*3+1] = 2*mu* (hh[2][1]*(c3*gg12[i][j]+0.5*gg21[i][j]));
    }
  }


  /* transforming kstab to the global coordinate system and */
  /* and adding to the one point quadrature matrix */
  for(i=0; i<8; i++)
  {
    for(j=0;j<8;j++)
    {
      for(k=0;k<3;k++)
      {
        for(l=0;l<3;l++)
        {
          s[i*3+k][j*3+l] = s[i*3+k][j*3+l]
            + (r[0][k]*kstab[i*3+0][j*3+0] + r[1][k]*kstab[i*3+1][j*3+0]
                + r[2][k]*kstab[i*3+2][j*3+0]) * r[0][l]
            + (r[0][k]*kstab[i*3+0][j*3+1] + r[1][k]*kstab[i*3+1][j*3+1]
                + r[2][k]*kstab[i*3+2][j*3+1]) * r[1][l]
            + (r[0][k]*kstab[i*3+0][j*3+2] + r[1][k]*kstab[i*3+1][j*3+2]
                + r[2][k]*kstab[i*3+2][j*3+2]) * r[2][l];
        }
      }
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of ale3_hourglass */



/*! @} (documentation module close)*/


#endif /* ifdef D_ALE */
#endif


/*!----------------------------------------------------------------------
\file
\brief shape functions and their natural derivatives for fluid3 element

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

------------------------------------------------------------------------*/
#ifdef D_FLUID3
#include "../headers/standardtypes.h"
#include "fluid3_prototypes.h"
#include "../fluid_full/fluid_prototypes.h"
#include "../fsi_full/fsi_prototypes.h"
static DOUBLE Q12 = ONE/TWO;
static DOUBLE Q14 = ONE/FOUR;
static DOUBLE Q18 = ONE/EIGHT;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;
/*!---------------------------------------------------------------------
\brief shape functions and their natural derivatives for hexaeder

<pre>                                                         genk 05/02

In this routine the shape functions and their natural first and second
derivatives with respect to r/s/t are evaluated for
 H E X A H E D E R

   Numbering of the nodes:

                           ^ t
                           |
                           |
                           |
                    8      |  15        7
                    o---------o---------o
                   /|                  /|
                  / |                 / |
                 /  |                /	|
              16o   |     o       14o   |
               /    o20       o    /	o19
              /     |             /     |
             /      |  13      6 /	|
          5 o---------o---------o	|
            |   o   |     o   	|   o   |  ---------->
            |       o---------o-|-------o           s
            |      / 4       11 |      /3
            |     /             |     /
          17o    /    o         o18  /
            | 12o         o     |   o10
            |  /                |  /
            | /                 | /
            |/	                |/
            o---------o---------o
	    1	/     9         2
	       /
	      /
	     /
	    r

   GiD:

                           ^ t
                           |
                           |
                           |
                    8      |  19        7
                    o---------o---------o
                   /|                  /|
                  / |                 / |
                 /  |                /	|
              20o   |   26o       18o   |
               /    o16     24o    /	o15
              /     |             /     |
             /      |  17      6 /	|
          5 o---------o---------o   23	|
            |   o   |   27o   	|   o   |  ---------->
            |  25   o---------o-|-------o           s
            |      / 4       11 |      /3
            |     /             |     /
          13o    /  22o         o14  /
            | 12o         o     |   o10
            |  /         21     |  /
            | /                 | /
            |/	                |/
            o---------o---------o
	    1	/     9         2
	       /
	      /
	     /
	    r



   PROBLEM: GID has a different numbering of the element nodes than this one.
            So either the shape functions for hex20 and hex27 (see drawing)
	    has to be adapted or during the input phase the numbering has to
	    be adapted to the shape functions.
	    This is all in progress and should be done for fluid3 and
	    brick1 the same way!!!!

   There are no HEX27 Elements in brick1 so we just go ahead here and
   use the GiD numbering for HEX27.

</pre>
\param  *funct    DOUBLE   (o)    shape functions
\param **deriv    DOUBLE   (o)    1st natural deriv. of shape funct.
\param **deriv2   DOUBLE   (o)    2nd natural deriv. of shape funct.
\param   r        DOUBLE   (i)    coordinate
\param   s        DOUBLE   (i)    coordinate
\param   t        DOUBLE   (i)    coordinate
\param   typ      DIS_TYP  (i)    element type
\param   icode    INT      (i)    evaluation flag
\return void
\warning shape functions for hex20/hex27/tet10 not implemented yet!!!

------------------------------------------------------------------------*/
void f3_hex(
               DOUBLE     *funct,
               DOUBLE    **deriv,
               DOUBLE    **deriv2,
               DOUBLE      r,
               DOUBLE      s,
               DOUBLE      t,
               DIS_TYP     typ,
               INT         icode
            )
{
DOUBLE rp,rm,sp,sm,tp,tm;
DOUBLE rrm,ssm,ttm;

#ifdef DEBUG
dstrc_enter("f3_rec");
#endif

/*------------------------------- selection of polynomial interpolation */
switch (typ)
{
case hex8: /* LINEAR shape functions and their natural derivatives ----*/
/*--------------------------------------------------- form basic values */
   rp=ONE+r;
   rm=ONE-r;
   sp=ONE+s;
   sm=ONE-s;
   tp=ONE+t;
   tm=ONE-t;

   funct[0]=Q18*rp*sm*tm;
   funct[1]=Q18*rp*sp*tm;
   funct[2]=Q18*rm*sp*tm;
   funct[3]=Q18*rm*sm*tm;
   funct[4]=Q18*rp*sm*tp;
   funct[5]=Q18*rp*sp*tp;
   funct[6]=Q18*rm*sp*tp;
   funct[7]=Q18*rm*sm*tp;

   if(icode>1) /* --> first derivative evaluation */
   {
      deriv[0][0]= Q18*sm*tm  ;
      deriv[0][1]= Q18*sp*tm  ;
      deriv[0][2]=-deriv[0][1];
      deriv[0][3]=-deriv[0][0];
      deriv[0][4]= Q18*sm*tp  ;
      deriv[0][5]= Q18*sp*tp  ;
      deriv[0][6]=-deriv[0][5];
      deriv[0][7]=-deriv[0][4];

      deriv[1][0]=-Q18*tm*rp  ;
      deriv[1][1]=-deriv[1][0];
      deriv[1][2]= Q18*tm*rm  ;
      deriv[1][3]=-deriv[1][2];
      deriv[1][4]=-Q18*tp*rp  ;
      deriv[1][5]=-deriv[1][4];
      deriv[1][6]= Q18*tp*rm  ;
      deriv[1][7]=-deriv[1][6];

      deriv[2][0]=-Q18*rp*sm  ;
      deriv[2][1]=-Q18*rp*sp  ;
      deriv[2][2]=-Q18*rm*sp  ;
      deriv[2][3]=-Q18*rm*sm  ;
      deriv[2][4]=-deriv[2][0];
      deriv[2][5]=-deriv[2][1];
      deriv[2][6]=-deriv[2][2];
      deriv[2][7]=-deriv[2][3];
   } /* endif (icode>1) */
   if(icode==3) /* --> second derivative evaluation */
   {
      deriv2[0][0] =  ZERO;
      deriv2[1][0] =  ZERO;
      deriv2[2][0] =  ZERO;
      deriv2[3][0] = -Q18*tm;
      deriv2[4][0] = -Q18*sm;
      deriv2[5][0] =  Q18*rp;

      deriv2[0][1] =  ZERO;
      deriv2[1][1] =  ZERO;
      deriv2[2][1] =  ZERO;
      deriv2[3][1] = -deriv2[3][0];
      deriv2[4][1] = -Q18*sp;
      deriv2[5][1] = -deriv2[5][0];

      deriv2[0][2] =  ZERO;
      deriv2[1][2] =  ZERO;
      deriv2[2][2] =  ZERO;
      deriv2[3][2] =  deriv2[3][0];
      deriv2[4][2] = -deriv2[4][1];
      deriv2[5][2] = -Q18*rm;

      deriv2[0][3] =  ZERO;
      deriv2[1][3] =  ZERO;
      deriv2[2][3] =  ZERO;
      deriv2[3][3] = -deriv2[3][0];
      deriv2[4][3] = -deriv2[4][0];
      deriv2[5][3] = -deriv2[5][2];

      deriv2[0][4] =  ZERO;
      deriv2[1][4] =  ZERO;
      deriv2[2][4] =  ZERO;
      deriv2[3][4] = -Q18*tp;
      deriv2[4][4] = -deriv2[4][0];
      deriv2[5][4] = -deriv2[5][0];

      deriv2[0][5] =  ZERO;
      deriv2[1][5] =  ZERO;
      deriv2[2][5] =  ZERO;
      deriv2[3][5] = -deriv2[3][4];
      deriv2[4][5] = -deriv2[4][1];
      deriv2[5][5] =  deriv2[5][0];

      deriv2[0][6] =  ZERO;
      deriv2[1][6] =  ZERO;
      deriv2[2][6] =  ZERO;
      deriv2[3][6] =  deriv2[3][4];
      deriv2[4][6] =  deriv2[4][1];
      deriv2[5][6] = -deriv2[5][2];

      deriv2[0][7] =  ZERO;
      deriv2[1][7] =  ZERO;
      deriv2[2][7] =  ZERO;
      deriv2[3][7] = -deriv2[3][4];
      deriv2[4][7] =  deriv2[4][0];
      deriv2[5][7] =  deriv2[5][2];
   } /* endif icode==3) */
break;

case hex20: /* QUADRATIC shape functions and their natural derivatives
                         without central nodes                      ----*/

   dserror("shape functions for hex20 are not correct!!! Compare with f3f_calfunctderiv.f \n");
/*--------------------------------------------------- form basic values */
   rp=ONE+r;
   rm=ONE-r;
   sp=ONE+s;
   sm=ONE-s;
   tp=ONE+t;
   tm=ONE-t;
   rrm=ONE-r*r;
   ssm=ONE-s*s;
   ttm=ONE-t*t;

   funct[0] =Q18*rp*sm*tm*(rp+sm+tm-FIVE);
   funct[1] =Q18*rp*sp*tm*(rp+sp+tm-FIVE);
   funct[2] =Q18*rm*sp*tm*(rm+sp+tm-FIVE);
   funct[3] =Q18*rm*sm*tm*(rm+sm+tm-FIVE);
   funct[4] =Q18*rp*sm*tp*(rp+sm+tp-FIVE);
   funct[5] =Q18*rp*sp*tp*(rp+sp+tp-FIVE);
   funct[6] =Q18*rm*sp*tp*(rm+sp+tp-FIVE);
   funct[7] =Q18*rm*sm*tp*(rm+sm+tp-FIVE);
   funct[8] =Q14*rp*ssm*tm;
   funct[9] =Q14*rrm*sp*tm;
   funct[10]=Q14*rm*ssm*tm;
   funct[11]=Q14*rrm*sm*tm;
   funct[12]=Q14*rp*ssm*tp;
   funct[13]=Q14*rrm*sp*tp;
   funct[14]=Q14*rm*ssm*tp;
   funct[15]=Q14*rrm*sm*tp;
   funct[16]=Q14*rp*sm*ttm;
   funct[17]=Q14*rp*sp*ttm;
   funct[18]=Q14*rm*sp*ttm;
   funct[19]=Q14*rm*sm*ttm;  /* analytisch gecheckt und fuer OK erklaert!!! */

   if(icode>1) /* --> first derivative evaluation */
   {
      deriv[0][0] = Q18*sm*tm*(TWO*rp+sm+tm-FIVE);
      deriv[1][0] =-Q18*tm*rp*(TWO*sm+tm+rp-FIVE);
      deriv[2][0] =-Q18*rp*sm*(TWO*tm+rp+sm-FIVE);

      deriv[0][1] = Q18*sp*tm*(TWO*rp+sp+tm-FIVE);
      deriv[1][1] = Q18*tm*rp*(TWO*sp+tm+rp-FIVE);
      deriv[2][1] =-Q18*rp*sp*(TWO*tm+rp+sp-FIVE);

      deriv[0][2] =-Q18*sp*tm*(TWO*rm+sp+tm-FIVE);
      deriv[1][2] = Q18*tm*rm*(TWO*sp+tm+rm-FIVE);
      deriv[2][2] =-Q18*rm*sp*(TWO*tm+rm+sp-FIVE);

      deriv[0][3] =-Q18*sm*tm*(TWO*rm+sm+tm-FIVE);
      deriv[1][3] =-Q18*tm*rm*(TWO*sm+tm+rm-FIVE);
      deriv[2][3] =-Q18*rm*sm*(TWO*tm+rm+sm-FIVE);

      deriv[0][4] = Q18*sm*tp*(TWO*rp+sm+tp-FIVE);
      deriv[1][4] =-Q18*tp*rp*(TWO*sm+tp+rp-FIVE);
      deriv[2][4] = Q18*rp*sm*(TWO*tp+rp+sm-FIVE);

      deriv[0][5] = Q18*sp*tp*(TWO*rp+sp+tp-FIVE);
      deriv[1][5] = Q18*tp*rp*(TWO*sp+tp+rp-FIVE);
      deriv[2][5] = Q18*rp*sp*(TWO*tp+rp+sp-FIVE);

      deriv[0][6] =-Q18*sp*tp*(TWO*rm+sp+tp-FIVE);
      deriv[1][6] = Q18*tp*rm*(TWO*sp+tp+rm-FIVE);
      deriv[2][6] = Q18*rm*sp*(TWO*tp+rm+sp-FIVE);

      deriv[0][7] =-Q18*sm*tp*(TWO*rm+sm+tp-FIVE);
      deriv[1][7] =-Q18*tp*rm*(TWO*sm+tp+rm-FIVE);
      deriv[2][7] = Q18*rm*sm*(TWO*tp+rm+sm-FIVE);

      deriv[0][8] = Q14*ssm*tm;
      deriv[1][8] =-Q12*s*tm*rp;
      deriv[2][8] =-Q14*ssm*rp;

      deriv[0][9] =-Q12*r*sp*tm;
      deriv[1][9] = Q14*rrm*tm;
      deriv[2][9] =-Q14*rrm*sp;

      deriv[0][10]=-deriv[0][8];
      deriv[1][10]=-Q12*s*tm*rm;
      deriv[2][10]=-Q14*ssm*rm;

      deriv[0][11]=-Q12*r*sm*tm;
      deriv[1][11]=-deriv[1][9];
      deriv[2][11]=-Q14*rrm*sm;

      deriv[0][12]= Q14*ssm*tp;
      deriv[1][12]=-Q12*s*tp*rp;
      deriv[2][12]=-deriv[2][8];

      deriv[0][13]=-Q12*r*sp*tp;
      deriv[1][13]= Q14*rrm*tp;
      deriv[2][13]=-deriv[2][8];

      deriv[0][14]=-deriv[0][12];
      deriv[1][14]=-Q12*s*tp*rm;
      deriv[2][14]=-deriv[2][10];

      deriv[0][15]=-Q12*r*sm*tp;
      deriv[1][15]=-deriv[1][13];
      deriv[2][15]=-deriv[2][11];

      deriv[0][16]= Q14*sm*ttm;
      deriv[1][16]=-Q14*ttm*rp;
      deriv[2][16]=-Q12*t*rp*sm;

      deriv[0][17]= Q14*sp*ttm;
      deriv[1][17]=-deriv[1][16];
      deriv[2][17]=-Q12*t*rp*sp;

      deriv[0][18]=-deriv[0][17];
      deriv[1][18]= Q14*ttm*rm;
      deriv[2][18]=-Q12*t*rm*sp;

      deriv[0][19]=-deriv[0][16];
      deriv[1][19]=-deriv[1][18];
      deriv[2][19]=-Q12*t*rm*sm;
   } /* endif (icode>1) */
   if(icode==3) /* --> second derivative evaluation  */
   {
      deriv2[0][0] = Q14*sm*tm;
      deriv2[1][0] = Q14*tm*rp;
      deriv2[2][0] = Q14*rp*sm;
      deriv2[3][0] =-Q18*(tm*(2*rp+sm+tm-FIVE+sm*tm));
      deriv2[4][0] =-Q18*(sm*(2*rp+sm+tm-FIVE+sm*tm));
      deriv2[5][0] = Q18*(rp*(2*sm+tm+rp-FIVE+tm*rp));

      deriv2[0][1] = Q14*sp*tm;
      deriv2[1][1] = deriv2[2][1];
      deriv2[2][1] = Q14*rp*sp;
      deriv2[3][1] =-Q18*(tm*(2*rp+sp+tm-FIVE+sp*tm));
      deriv2[4][1] =-Q18*(sp*(2*rp+sp+tm-FIVE+sp*tm));
      deriv2[5][1] =-Q18*(rp*(2*sp+tm+rp-FIVE+tm*rp));

      deriv2[0][2] =-deriv2[1][2];
      deriv2[1][2] = Q14*tm*rm;
      deriv2[2][2] = Q14*rm*sp;
      deriv2[3][2] =-Q18*(tm*(2*rm+sp+tm-FIVE+sp*tm));
      deriv2[4][2] = Q18*(sp*(2*rm+sp+tm-FIVE+sp*tm));
      deriv2[5][2] =-Q18*(rm*(2*sp+tm+rm-FIVE+tm*rm));

      deriv2[0][3] =-deriv2[1][1];
      deriv2[1][3] = deriv2[2][3];
      deriv2[2][3] = Q14*rm*sm;
      deriv2[3][3] =-Q18*(tm*(2*rm+sm+tm-FIVE+sm*tm));
      deriv2[4][3] = Q18*(sm*(2*rm+sm+tm-FIVE+sm*tm));
      deriv2[5][3] = Q18*(rm*(2*sm+tm+rm-FIVE+tm*rm));

      deriv2[0][4] = Q14*sm*tp;
      deriv2[1][4] = Q14*tp*rp;
      deriv2[2][4] = deriv2[3][1];
      deriv2[3][4] =-Q18*(tp*(2*rp+sm+tp-FIVE+sm*tp));
      deriv2[4][4] = Q18*(sm*(2*rp+sm+tp-FIVE+sm*tp));
      deriv2[5][4] =-Q18*(rp*(2*sm+tp+rp-FIVE+tp*rp));

      deriv2[0][5] = Q14*sp*tp;
      deriv2[1][5] = deriv2[2][5];
      deriv2[2][5] = deriv2[3][2];
      deriv2[3][5] =-Q18*(tp*(2*rp+sp+tp-FIVE+sp*tp));
      deriv2[4][5] = Q18*(sp*(2*rp+sp+tp-FIVE+sp*tp));
      deriv2[5][5] = Q18*(rp*(2*sp+tp+rp-FIVE+tp*rp));

      deriv2[0][6] =-deriv2[1][6];
      deriv2[1][6] = Q14*tp*rm;
      deriv2[2][6] = deriv2[3][3];
      deriv2[3][6] =-Q18*(tp*(2*rm+sp+tp-FIVE+sp*tp));
      deriv2[4][6] =-Q18*(sp*(2*rm+sp+tp-FIVE+sp*tp));
      deriv2[5][6] = Q18*(rm*(2*sp+tp+rm-FIVE+tp*rm));

      deriv2[0][7] =-deriv2[1][5];
      deriv2[1][7] = deriv2[2][7];
      deriv2[2][7] = deriv2[3][4];
      deriv2[3][7] =-Q18*(tp*(2*rm+sm+tp-FIVE+sm*tp));
      deriv2[4][7] =-Q18*(sm*(2*rm+sm+tp-FIVE+sm*tp));
      deriv2[5][7] =-Q18*(rm*(2*sm+tp+rm-FIVE+tp*rm));

      deriv2[0][8] = ZERO;
      deriv2[1][8] = -Q12*tm*rp;
      deriv2[2][8] = ZERO;
      deriv2[3][8] =-Q12*s*tm;
      deriv2[4][8] =-Q14*ssm;
      deriv2[5][8] = Q12*s*rp;

      deriv2[0][9]=-Q12*sp*tm;
      deriv2[1][9]= ZERO;
      deriv2[2][9]= ZERO;
      deriv2[3][9]=-Q12*r*tm;
      deriv2[4][9]= Q12*r*sp;
      deriv2[5][9]=-Q14*rrm ;

      deriv2[0][10]= ZERO;
      deriv2[1][10]= -Q12*tm*rm;
      deriv2[2][10]= ZERO;
      deriv2[3][10]= Q12*s*tm;
      deriv2[4][10]=-deriv2[4][8];
      deriv2[5][10]= Q12*s*rm;

      deriv2[0][11]=-Q12*sm*tm;
      deriv2[1][11]= ZERO;
      deriv2[2][11]= ZERO;
      deriv2[3][11]= Q12*r*tm;
      deriv2[4][11]= Q12*r*sm;
      deriv2[5][11]=-deriv2[5][9];

      deriv2[0][12]= ZERO;
      deriv2[1][12]= -Q12*tp*rp;
      deriv2[2][12]= ZERO;
      deriv2[3][12]=-Q12*s*tp;
      deriv2[4][12]=-deriv2[4][8];
      deriv2[5][12]=-deriv2[5][8];

      deriv2[0][13]=-Q12*sp*tp;
      deriv2[1][13]= ZERO;
      deriv2[2][13]= ZERO;
      deriv2[3][13]=-Q12*r*tp;
      deriv2[4][13]=-deriv2[4][9];
      deriv2[5][13]=-deriv2[5][9];

      deriv2[0][14]= ZERO;
      deriv2[1][14]= -Q12*tp*rm;
      deriv2[2][14]= ZERO;
      deriv2[3][14]= Q12*s*tp;
      deriv2[4][14]= deriv2[4][8];
      deriv2[5][14]=-deriv2[5][10];

      deriv2[0][15]=-Q12*sm*tp;
      deriv2[1][15]= ZERO;
      deriv2[2][15]= ZERO;
      deriv2[3][15]= Q12*r*tp;
      deriv2[4][15]=-deriv2[4][11];
      deriv2[5][15]= deriv2[5][9];

      deriv2[0][16]= ZERO;
      deriv2[1][16]= ZERO;
      deriv2[2][16]= ZERO;
      deriv2[3][16]=-Q14*ttm;
      deriv2[4][16]=-Q12*t*sm;
      deriv2[5][16]= Q12*t*rp;

      deriv2[0][17]= ZERO;
      deriv2[1][17]= ZERO;
      deriv2[2][17]= ZERO;
      deriv2[3][17]= Q14*ttm;
      deriv2[4][17]=-Q12*t*sp;
      deriv2[5][17]=-deriv2[5][16];

      deriv2[0][18]= ZERO;
      deriv2[1][18]= ZERO;
      deriv2[2][18]= ZERO;
      deriv2[3][18]= deriv2[3][16];
      deriv2[4][18]= Q12*t*sp;
      deriv2[5][18]= Q12*t*rm;

      deriv2[0][19]= ZERO;
      deriv2[1][19]= ZERO;
      deriv2[2][19]= ZERO;
      deriv2[3][19]= deriv2[3][17];
      deriv2[4][19]= Q12*t*sm;
      deriv2[5][19]=-deriv2[5][18];
   } /* endif (icode==3) */
break;

case hex27: /* QUADRATIC shape functions and their natural derivatives
               with central nodes                         ----*/
/*--------------------------------------------------- form basic values */
{
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

  if (icode>1) /* --> first derivative evaluation */
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
  if (icode==3) /* --> second derivative evaluation */
  {
    deriv2[0][0] = sp1*tp1;
    deriv2[0][1] = sm1*tp1;
    deriv2[0][2] = sm1*tp1;
    deriv2[0][3] = sp1*tp1;
    deriv2[0][4] = tm1*sp1;
    deriv2[0][5] = sm1*tm1;
    deriv2[0][6] = sm1*tm1;
    deriv2[0][7] = tm1*sp1;
    deriv2[0][8] = s00*tp1;
    deriv2[0][9] = -2*sm1*tp1;
    deriv2[0][10] = s00*tp1;
    deriv2[0][11] = -2*sp1*tp1;
    deriv2[0][12] = t00*sp1;
    deriv2[0][13] = t00*sm1;
    deriv2[0][14] = t00*sm1;
    deriv2[0][15] = t00*sp1;
    deriv2[0][16] = s00*tm1;
    deriv2[0][17] = -2*sm1*tm1;
    deriv2[0][18] = s00*tm1;
    deriv2[0][19] = -2*tm1*sp1;
    deriv2[0][20] = -2*s00*tp1;
    deriv2[0][21] = s00*t00;
    deriv2[0][22] = -2*t00*sm1;
    deriv2[0][23] = s00*t00;
    deriv2[0][24] = -2*t00*sp1;
    deriv2[0][25] = -2*s00*tm1;
    deriv2[0][26] = -2*s00*t00;

    deriv2[1][0] = rp1*tp1;
    deriv2[1][1] = rp1*tp1;
    deriv2[1][2] = rm1*tp1;
    deriv2[1][3] = rm1*tp1;
    deriv2[1][4] = tm1*rp1;
    deriv2[1][5] = tm1*rp1;
    deriv2[1][6] = rm1*tm1;
    deriv2[1][7] = rm1*tm1;
    deriv2[1][8] = -2*rp1*tp1;
    deriv2[1][9] = r00*tp1;
    deriv2[1][10] = -2*rm1*tp1;
    deriv2[1][11] = r00*tp1;
    deriv2[1][12] = t00*rp1;
    deriv2[1][13] = t00*rp1;
    deriv2[1][14] = t00*rm1;
    deriv2[1][15] = t00*rm1;
    deriv2[1][16] = -2*tm1*rp1;
    deriv2[1][17] = r00*tm1;
    deriv2[1][18] = -2*rm1*tm1;
    deriv2[1][19] = r00*tm1;
    deriv2[1][20] = -2*r00*tp1;
    deriv2[1][21] = -2*t00*rp1;
    deriv2[1][22] = r00*t00;
    deriv2[1][23] = -2*t00*rm1;
    deriv2[1][24] = r00*t00;
    deriv2[1][25] = -2*r00*tm1;
    deriv2[1][26] = -2*r00*t00;

    deriv2[2][0] = rp1*sp1;
    deriv2[2][1] = sm1*rp1;
    deriv2[2][2] = rm1*sm1;
    deriv2[2][3] = rm1*sp1;
    deriv2[2][4] = rp1*sp1;
    deriv2[2][5] = sm1*rp1;
    deriv2[2][6] = rm1*sm1;
    deriv2[2][7] = rm1*sp1;
    deriv2[2][8] = s00*rp1;
    deriv2[2][9] = r00*sm1;
    deriv2[2][10] = s00*rm1;
    deriv2[2][11] = r00*sp1;
    deriv2[2][12] = -2*rp1*sp1;
    deriv2[2][13] = -2*sm1*rp1;
    deriv2[2][14] = -2*rm1*sm1;
    deriv2[2][15] = -2*rm1*sp1;
    deriv2[2][16] = s00*rp1;
    deriv2[2][17] = r00*sm1;
    deriv2[2][18] = s00*rm1;
    deriv2[2][19] = r00*sp1;
    deriv2[2][20] = r00*s00;
    deriv2[2][21] = -2*s00*rp1;
    deriv2[2][22] = -2*r00*sm1;
    deriv2[2][23] = -2*s00*rm1;
    deriv2[2][24] = -2*r00*sp1;
    deriv2[2][25] = r00*s00;
    deriv2[2][26] = -2*r00*s00;

    deriv2[3][0] = tp1*drp1*dsp1;
    deriv2[3][1] = tp1*dsm1*drp1;
    deriv2[3][2] = tp1*drm1*dsm1;
    deriv2[3][3] = tp1*drm1*dsp1;
    deriv2[3][4] = tm1*drp1*dsp1;
    deriv2[3][5] = tm1*dsm1*drp1;
    deriv2[3][6] = tm1*drm1*dsm1;
    deriv2[3][7] = tm1*drm1*dsp1;
    deriv2[3][8] = tp1*ds00*drp1;
    deriv2[3][9] = tp1*dr00*dsm1;
    deriv2[3][10] = tp1*ds00*drm1;
    deriv2[3][11] = tp1*dr00*dsp1;
    deriv2[3][12] = t00*drp1*dsp1;
    deriv2[3][13] = t00*dsm1*drp1;
    deriv2[3][14] = t00*drm1*dsm1;
    deriv2[3][15] = t00*drm1*dsp1;
    deriv2[3][16] = tm1*ds00*drp1;
    deriv2[3][17] = tm1*dr00*dsm1;
    deriv2[3][18] = tm1*ds00*drm1;
    deriv2[3][19] = tm1*dr00*dsp1;
    deriv2[3][20] = 4*r*s*tp1;
    deriv2[3][21] = t00*ds00*drp1;
    deriv2[3][22] = t00*dr00*dsm1;
    deriv2[3][23] = t00*ds00*drm1;
    deriv2[3][24] = t00*dr00*dsp1;
    deriv2[3][25] = 4*r*s*tm1;
    deriv2[3][26] = 4*r*s*t00;

    deriv2[4][0] = sp1*drp1*dtp1;
    deriv2[4][1] = sm1*drp1*dtp1;
    deriv2[4][2] = sm1*drm1*dtp1;
    deriv2[4][3] = sp1*drm1*dtp1;
    deriv2[4][4] = sp1*dtm1*drp1;
    deriv2[4][5] = sm1*dtm1*drp1;
    deriv2[4][6] = sm1*drm1*dtm1;
    deriv2[4][7] = sp1*drm1*dtm1;
    deriv2[4][8] = s00*drp1*dtp1;
    deriv2[4][9] = sm1*dr00*dtp1;
    deriv2[4][10] = s00*drm1*dtp1;
    deriv2[4][11] = sp1*dr00*dtp1;
    deriv2[4][12] = sp1*dt00*drp1;
    deriv2[4][13] = sm1*dt00*drp1;
    deriv2[4][14] = sm1*dt00*drm1;
    deriv2[4][15] = sp1*dt00*drm1;
    deriv2[4][16] = s00*dtm1*drp1;
    deriv2[4][17] = sm1*dr00*dtm1;
    deriv2[4][18] = s00*drm1*dtm1;
    deriv2[4][19] = sp1*dr00*dtm1;
    deriv2[4][20] = s00*dr00*dtp1;
    deriv2[4][21] = s00*dt00*drp1;
    deriv2[4][22] = 4*r*t*sm1;
    deriv2[4][23] = s00*dt00*drm1;
    deriv2[4][24] = 4*r*t*sp1;
    deriv2[4][25] = s00*dr00*dtm1;
    deriv2[4][26] = 4*r*t*s00;

    deriv2[5][0] = rp1*dsp1*dtp1;
    deriv2[5][1] = rp1*dsm1*dtp1;
    deriv2[5][2] = rm1*dsm1*dtp1;
    deriv2[5][3] = rm1*dsp1*dtp1;
    deriv2[5][4] = rp1*dtm1*dsp1;
    deriv2[5][5] = rp1*dsm1*dtm1;
    deriv2[5][6] = rm1*dsm1*dtm1;
    deriv2[5][7] = rm1*dtm1*dsp1;
    deriv2[5][8] = rp1*ds00*dtp1;
    deriv2[5][9] = r00*dsm1*dtp1;
    deriv2[5][10] = rm1*ds00*dtp1;
    deriv2[5][11] = r00*dsp1*dtp1;
    deriv2[5][12] = rp1*dt00*dsp1;
    deriv2[5][13] = rp1*dt00*dsm1;
    deriv2[5][14] = rm1*dt00*dsm1;
    deriv2[5][15] = rm1*dt00*dsp1;
    deriv2[5][16] = rp1*ds00*dtm1;
    deriv2[5][17] = r00*dsm1*dtm1;
    deriv2[5][18] = rm1*ds00*dtm1;
    deriv2[5][19] = r00*dtm1*dsp1;
    deriv2[5][20] = r00*ds00*dtp1;
    deriv2[5][21] = 4*s*t*rp1;
    deriv2[5][22] = r00*dt00*dsm1;
    deriv2[5][23] = 4*s*t*rm1;
    deriv2[5][24] = r00*dt00*dsp1;
    deriv2[5][25] = r00*ds00*dtm1;
    deriv2[5][26] = 4*s*t*r00;
  }
  break;
}

/*----------------------------------------------------------------------*/
default:
  dserror("distyp unknown\n");
} /* end switch (typ) */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f3_hex */

/*!---------------------------------------------------------------------
\brief shape functions and their natural derivatives for rectangles

<pre>                                                         genk 02/04

In this routine the shape functions and their natural first and second
derivatives with respect to r/s are evaluated for
 R E C T A N G L E S

Numbering of the nodes:


                    ^ s
                    |
              1     |4    0
              o-----o-----o
              |           |
              |           |7
            5 o     o     o -----> r
              |     8     |
              |           |
              o-----o-----o
              2     6     3

</pre>
\param  *funct    DOUBLE   (o)    shape functions
\param **deriv    DOUBLE   (o)    1st natural deriv. of shape funct.
\param **deriv2   DOUBLE   (o)    2nd natural deriv. of shape funct.
\param   r        DOUBLE   (i)    coordinate
\param   s        DOUBLE   (i)    coordinate
\param   typ      DIS_TYP  (i)    element type
\param   icode    INT      (i)    evaluation flag
\return void

------------------------------------------------------------------------*/
void f3_rec(
            DOUBLE     *funct,
            DOUBLE    **deriv,
            DOUBLE    **deriv2,
            DOUBLE      r,
            DOUBLE      s,
            DIS_TYP     typ,
            INT         icode
            )
{
INT    i,ii;
DOUBLE rp,rm,sp,sm,r2,s2;
DOUBLE rh,sh,rs;
DOUBLE rhm,shm,rhp,shp;

#ifdef DEBUG
dstrc_enter("f3_rec");
#endif

/*------------------------------- selection of polynomial interpolation */
switch (typ)
{
/*----------------------------------------------------------------------*/
case quad4: /* LINEAR shape functions and their natural derivatives ----*/
/*--------------------------------------------------- form basic values */
   rp=ONE+r;
   rm=ONE-r;
   sp=ONE+s;
   sm=ONE-s;

   dsassert(funct!=NULL,"cannot write to NULL-pointer!\n");
   funct[0]=Q14*rp*sp;
   funct[1]=Q14*rm*sp;
   funct[2]=Q14*rm*sm;
   funct[3]=Q14*rp*sm;

   if(icode>1) /* --> first derivative evaluation */
   {
      dsassert(deriv!=NULL,"cannot write to NULL-pointer!\n");
      deriv[0][0]= Q14*sp;
      deriv[1][0]= Q14*rp;

      deriv[0][1]=-Q14*sp;
      deriv[1][1]= Q14*rm;

      deriv[0][2]=-Q14*sm;
      deriv[1][2]=-Q14*rm;

      deriv[0][3]= Q14*sm;
      deriv[1][3]=-Q14*rp;
   } /* endif (icode>1) */

   if(icode==3) /* --> second derivative evaluation */
   {
      dsassert(deriv2!=NULL,"cannot write to NULL-pointer!\n");
      deriv2[0][0]= ZERO;
      deriv2[1][0]= ZERO;
      deriv2[2][0]= Q14;

      deriv2[0][1]= ZERO;
      deriv2[1][1]= ZERO;
      deriv2[2][1]=-Q14;

      deriv2[0][2]= ZERO;
      deriv2[1][2]= ZERO;
      deriv2[2][2]= Q14;

      deriv2[0][3]= ZERO;
      deriv2[1][3]=ZERO;
      deriv2[2][3]=-Q14;
   } /* endif (icode==3) */
break;

/*----------------------------------------------------------------------*/
case quad8: /* QUADRATIC shape functions and their natural derivatives
               without central node ------------------------------------*/
/*--------------------------------------------------- form basic values */
   rp=ONE+r;
   rm=ONE-r;
   sp=ONE+s;
   sm=ONE-s;
   r2=ONE-r*r;
   s2=ONE-s*s;

   funct[4]=Q12*r2*sp;
   funct[5]=Q12*rm*s2;
   funct[6]=Q12*r2*sm;
   funct[7]=Q12*rp*s2;
   funct[0]=Q14*rp*sp-Q12*(funct[4]+funct[7]);
   funct[1]=Q14*rm*sp-Q12*(funct[4]+funct[5]);
   funct[2]=Q14*rm*sm-Q12*(funct[5]+funct[6]);
   funct[3]=Q14*rp*sm-Q12*(funct[6]+funct[7]);

   if(icode>1) /* --> first derivative evaluation */
   {
      deriv[0][0]= Q14*sp;
      deriv[1][0]= Q14*rp;

      deriv[0][1]=-Q14*sp;
      deriv[1][1]= Q14*rm;

      deriv[0][2]=-Q14*sm;
      deriv[1][2]=-Q14*rm;

      deriv[0][3]= Q14*sm;
      deriv[1][3]=-Q14*rp;

      deriv[0][4]=-ONE*r*sp;
      deriv[1][4]= Q12*r2;

      deriv[0][5]=-Q12*  s2;
      deriv[1][5]=-ONE*rm*s;

      deriv[0][6]=-ONE*r*sm;
      deriv[1][6]=-Q12*r2;

      deriv[0][7]= Q12*  s2;
      deriv[1][7]=-ONE*rp*s;

      deriv[0][0]-= Q12*(deriv[0][4]+deriv[0][7]);
      deriv[1][0]-= Q12*(deriv[1][4]+deriv[1][7]);

      for(i=1;i<4;i++)
      {
         ii=i+3;
         deriv[0][i] -= Q12*(deriv[0][ii]+deriv[0][ii+1]);
	 deriv[1][i] -= Q12*(deriv[1][ii]+deriv[1][ii+1]);
      } /* end loop over i */
   } /* endif (icode>1) */

   if(icode==3) /* --> second derivative evaluation */
   {
      deriv2[0][0]= ZERO;
      deriv2[1][0]= ZERO;
      deriv2[2][0]= Q14;

      deriv2[0][1]= ZERO;
      deriv2[1][1]= ZERO;
      deriv2[2][1]=-Q14;

      deriv2[0][2]= ZERO;
      deriv2[1][2]= ZERO;
      deriv2[2][2]= Q14;

      deriv2[0][3]= ZERO;
      deriv2[1][3]= ZERO;
      deriv2[2][3]=-Q14;

      deriv2[0][4]=-(ONE+s);
      deriv2[1][4]= ZERO;
      deriv2[2][4]=-r;

      deriv2[0][5]= ZERO;
      deriv2[1][5]=-(ONE-r);
      deriv2[2][5]= s;

      deriv2[0][6]=-(ONE-s);
      deriv2[1][6]= ZERO;
      deriv2[2][6]= r;

      deriv2[0][7]= ZERO;
      deriv2[1][7]=-(ONE+r);
      deriv2[2][7]=-s;

      deriv2[0][0] -= Q12*(deriv2[0][4]+deriv2[0][7]);
      deriv2[1][0] -= Q12*(deriv2[1][4]+deriv2[1][7]);
      deriv2[2][0] -= Q12*(deriv2[2][4]+deriv2[2][7]);

      for(i=1;i<4;i++)
      {
         ii=i+3;
         deriv2[0][i] -= Q12*(deriv2[0][ii]+deriv2[0][ii+1]);
	 deriv2[1][i] -= Q12*(deriv2[1][ii]+deriv2[1][ii+1]);
	 deriv2[2][i] -= Q12*(deriv2[2][ii]+deriv2[2][ii+1]);
      } /* end loop over i */
   } /* endif (icode==3) */
break;

/*----------------------------------------------------------------------*/
case quad9: /* QUADRATIC  shape functions and their natural derivatives
               with central node ---------------------------------------*/
/*--------------------------------------------------- form basic values */
   rp=ONE+r;
   rm=ONE-r;
   sp=ONE+s;
   sm=ONE-s;
   r2=ONE-r*r;
   s2=ONE-s*s;
   rh=Q12*r;
   sh=Q12*s;
   rs=rh*sh;
   rhp=r+Q12;
   rhm=r-Q12;
   shp=s+Q12;
   shm=s-Q12;

   funct[0]= rs*rp*sp;
   funct[1]=-rs*rm*sp;
   funct[2]= rs*rm*sm;
   funct[3]=-rs*rp*sm;
   funct[4]= sh*sp*r2;
   funct[5]=-rh*rm*s2;
   funct[6]=-sh*sm*r2;
   funct[7]= rh*rp*s2;
   funct[8]= r2*s2;

   if(icode>1) /* --> first derivative evaluation */
   {
      deriv[0][0]= rhp*sh*sp;
      deriv[1][0]= shp*rh*rp;

      deriv[0][1]= rhm*sh*sp;
      deriv[1][1]=-shp*rh*rm;

      deriv[0][2]=-rhm*sh*sm;
      deriv[1][2]=-shm*rh*rm;

      deriv[0][3]=-rhp*sh*sm;
      deriv[1][3]= shm*rh*rp;

      deriv[0][4]=-TWO*r*sh*sp;
      deriv[1][4]= shp*r2;

      deriv[0][5]= rhm*s2;
      deriv[1][5]= TWO*s*rh*rm;

      deriv[0][6]= TWO*r*sh*sm;
      deriv[1][6]= shm*r2;

      deriv[0][7]= rhp*s2;
      deriv[1][7]=-TWO*s*rh*rp;

      deriv[0][8]=-TWO*r*s2;
      deriv[1][8]=-TWO*s*r2;
   } /* endif (icode>1) */

   if(icode==3) /* --> second derivative evaluation */
   {
      deriv2[0][0]= sh*sp;
      deriv2[1][0]= rh*rp;
      deriv2[2][0]= shp*rhp;

      deriv2[0][1]= sh*sp;
      deriv2[1][1]=-rh*rm;
      deriv2[2][1]= shp*rhm;

      deriv2[0][2]=-sh*sm;
      deriv2[1][2]=-rh*rm;
      deriv2[2][2]= shm*rhm;

      deriv2[0][3]=-sh*sm;
      deriv2[1][3]= rh*rp;
      deriv2[2][3]= shm*rhp;

      deriv2[0][4]=-TWO*sh*sp;
      deriv2[1][4]= r2;
      deriv2[2][4]=-TWO*r*shp;

      deriv2[0][5]= s2;
      deriv2[1][5]= TWO*rh*rm;
      deriv2[2][5]=-TWO*s*rhm;

      deriv2[0][6]= TWO*sh*sm;
      deriv2[1][6]= r2;
      deriv2[2][6]=-TWO*r*shm;

      deriv2[0][7]= s2;
      deriv2[1][7]=-TWO*rh*rp;
      deriv2[2][7]=-TWO*s*rhp;

      deriv2[0][8]=-TWO*s2;
      deriv2[1][8]=-TWO*r2;
      deriv2[2][8]= TWO*s*TWO*r;
   } /* endif (icode==3) */
break;

/*----------------------------------------------------------------------*/
default:
   dserror("distyp unknown");
} /* end switch(typ) */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f3_rec */

/*!---------------------------------------------------------------------
\brief shape functions and their natural derivatives for tetraeder

<pre>                                                         genk 08/02

In this routine the shape functions and their natural first and second
derivatives with respect to r/s/t are evaluated for
T E T R A E D E R

</pre>
\param  *funct    DOUBLE   (o)    shape functions
\param **deriv    DOUBLE   (o)    1st natural deriv. of shape funct.
\param **deriv2   DOUBLE   (o)    2nd natural deriv. of shape funct.
\param   r        DOUBLE   (i)    coordinate
\param   s        DOUBLE   (i)    coordinate
\param   t        DOUBLE   (i)    coordinate
\param   typ      DIS_TYP  (i)    element type
\param   icode    INT      (i)    evaluation flag
\return void
\warning shape functions for TET10 not implemented yet!!!

------------------------------------------------------------------------*/
void f3_tet(
            DOUBLE     *funct,
            DOUBLE    **deriv,
            DOUBLE    **deriv2,
            DOUBLE      r,
            DOUBLE      s,
            DOUBLE      t,
            DIS_TYP     typ,
            INT         icode
            )
{
DOUBLE t1,t2,t3,t4;


#ifdef DEBUG
dstrc_enter("f3_tet");
#endif
/*------------------------------- selection of polynomial interpolation */
switch (typ)
{
case tet4: /* LINEAR shape functions and their natural derivatives -----*/
/*--------------------------------------------------- form basic values */

  /*
   Numbering of the nodes:
   -----------------------
   - this is the numbering used in GiD!!


          4 o---
            |\  ---
            |  \   -o3
            |   \  / \
            |     \   \
            |    / \   \
            |   /    \  \
            |  /      \  \
            | /         \ \
            |/            \\
            o---------------o
           1                2
   */


   t1=ONE-r-s-t;
   t2=r;
   t3=s;
   t4=t;

   funct[0]= t1;
   funct[1]= t2;
   funct[2]= t3;
   funct[3]= t4;

   if(icode>1) /* --> first derivative evaluation */
   {
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
   } /* endif (icode>1) */
break;

case tet10: /*  QUADRATIC shape functions and their natural derivatives */

   dserror("shape functions for tet10 not implemented yet!\n");
/*--------------------------------------------------- form basic values */
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


   if(icode>1) /* --> first derivative evaluation */
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
   if(icode==3) /* --> second derivative evaluation  */
   {
      deriv2[0][0] =  ;
      deriv2[1][0] =  ;
      deriv2[2][0] =  ;
      deriv2[3][0] = ;
      deriv2[4][0] = ;
      deriv2[5][0] = ;

      deriv2[0][1] =  ;
      deriv2[1][1] =  ;
      deriv2[2][1] =  ;
      deriv2[3][1] = ;
      deriv2[4][1] = ;
      deriv2[5][1] = ;

      deriv2[0][2] =  ;
      deriv2[1][2] = ;
      deriv2[2][2] =  ;
      deriv2[3][2] = ;
      deriv2[4][2] = ;
      deriv2[5][2] = ;

      deriv2[0][3] = ;
      deriv2[1][3] =  ;
      deriv2[2][3] =  ;
      deriv2[3][3] = ;
      deriv2[4][3] = ;
      deriv2[5][3] = ;

      deriv2[0][4] =  ;
      deriv2[1][4] =  ;
      deriv2[2][4] =  ;
      deriv2[3][4] = ;
      deriv2[4][4] = ;
      deriv2[5][4] = ;

      deriv2[0][5] =  ;
      deriv2[1][5] =  ;
      deriv2[2][5] =  ;
      deriv2[3][5] = ;
      deriv2[4][5] = ;
      deriv2[5][5] = ;

      deriv2[0][6] = ;
      deriv2[1][6] =  ;
      deriv2[2][6] =  ;
      deriv2[3][6] = ;
      deriv2[4][6] = ;
      deriv2[5][6] = ;

      deriv2[0][7] = ;
      deriv2[1][7] = ;
      deriv2[2][7] = ;
      deriv2[3][7] = ;
      deriv2[4][7] = ;
      deriv2[5][7] = ;

      deriv2[0][8] =  ;
      deriv2[1][8] =  ;
      deriv2[2][8] =  ;
      deriv2[3][8] = ;
      deriv2[4][8] = ;
      deriv2[5][8] =  ;

      deriv2[0][9]= ;
      deriv2[1][9]=  ;
      deriv2[2][9]=  ;
      deriv2[3][9]= ;
      deriv2[4][9]= ;
      deriv2[5][9]=  ;
   }
#endif

break;

/*----------------------------------------------------------------------*/
default:
   dserror("distyp unknown\n");
} /* end switch(typ) */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f3_tet */

/*!---------------------------------------------------------------------
\brief shape functions and their natural derivatives for triangles

<pre>                                                         genk 02/04

In this routine the shape functions and their natural first and second
derivatives with respect to r/s are evaluated for
T R I A N G L E S

</pre>
\param  *funct    DOUBLE   (o)    shape functions
\param **deriv    DOUBLE   (o)    1st natural deriv. of shape funct.
\param **deriv2   DOUBLE   (o)    2nd natural deriv. of shape funct.
\param   r        DOUBLE   (i)    coordinate
\param   s        DOUBLE   (i)    coordinate
\param   typ      DIS_TYP  (i)    element type
\param   icode    INT      (i)    evaluation flag
\return void

------------------------------------------------------------------------*/
void f3_tri(
            DOUBLE     *funct,
            DOUBLE    **deriv,
            DOUBLE    **deriv2,
            DOUBLE      r,
            DOUBLE      s,
            DIS_TYP     typ,
            INT         icode
	   )
{

DOUBLE rr,rs,ss;

#ifdef DEBUG
dstrc_enter("f3_tri");
#endif

/*------------------------------- selection of polynomial interpolation */
switch (typ)
{
/*----------------------------------------------------------------------*/
case tri3: /* LINEAR shape functions and their natural derivatives -----*/
/*----------------------------------------------------------------------*/

   funct[0]=ONE-r-s;
   funct[1]=r;
   funct[2]=s;

   if(icode>1) /* --> first derivative evaluation */
   {
      deriv[0][0]=-ONE;
      deriv[1][0]=-ONE;
      deriv[0][1]= ONE;
      deriv[1][1]=ZERO;
      deriv[0][2]=ZERO;
      deriv[1][2]= ONE;
   } /* endif (icode>1) */
break;

/*----------------------------------------------------------------------*/
case tri6: /* QUADRATIC shape functions and ther natural derivatives ---*/
/*--------------------------------------------------- form basic values */
   rr=r*r;
   ss=s*s;
   rs=r*s;

   funct[0]=(ONE-TWO*r-TWO*s)*(ONE-r-s);
   funct[1]=TWO*rr-r;
   funct[2]=TWO*ss-s;
   funct[3]=FOUR*(r-rr-rs);
   funct[4]=FOUR*rs;
   funct[5]=FOUR*(s-rs-ss);

   if(icode>1) /* --> first derivative evaluation */
   {
      deriv[0][0]=-THREE+FOUR*(r+s);
      deriv[1][0]= deriv[0][0];

      deriv[0][1]= FOUR*r-ONE;
      deriv[1][1]= ZERO;

      deriv[0][2]= ZERO;
      deriv[1][2]= FOUR*s-ONE;

      deriv[0][3]= FOUR*(ONE-TWO*r-s);
      deriv[1][3]=-FOUR*r;

      deriv[0][4]= FOUR*s;
      deriv[1][4]= FOUR*r;

      deriv[0][5]=-FOUR*s;
      deriv[1][5]= FOUR*(ONE-r-TWO*s);
   } /* endif (icode>1) */

   if(icode==3) /* --> second derivative evaluation */
   {
      deriv2[0][0]= FOUR;
      deriv2[1][0]= FOUR;
      deriv2[2][0]= FOUR;

      deriv2[0][1]= FOUR;
      deriv2[1][1]= ZERO;
      deriv2[2][1]= ZERO;

      deriv2[0][2]= ZERO;
      deriv2[1][2]= FOUR;
      deriv2[2][2]= ZERO;

      deriv2[0][3]=-EIGHT;
      deriv2[1][3]= ZERO;
      deriv2[2][3]=-FOUR;

      deriv2[0][4]= ZERO;
      deriv2[1][4]= ZERO;
      deriv2[2][4]= FOUR;

      deriv2[0][5]= ZERO;
      deriv2[1][5]=-EIGHT;
      deriv2[2][5]=-FOUR;
   } /* endif (icode==3) */
break;

/*----------------------------------------------------------------------*/
default:
   dserror("distyp unknown");
} /* end swtich(typ) */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f3_tri */

/*!---------------------------------------------------------------------
\brief jacobian matrix

<pre>                                                         genk 05/02

In this routine the jacobian matrix and its determinant is calculated

</pre>
\param **xyze      DOUBLE   (i)    nodal coordinates
\param **deriv     DOUBLE   (i)    natural deriv. of shape funcs
\param **xjm       DOUBLE   (o)    jacobian matrix
\param  *det       DOUBLE   (o)    determinant of jacobian matrix
\param  *ele       ELEMENT  (i)    actual element
\param   iel       INT      (i)    num. of nodes of act. ele
\return void

------------------------------------------------------------------------*/
void f3_jaco(DOUBLE    **xyze,
             DOUBLE    **deriv,
             DOUBLE    **xjm,
             DOUBLE     *det,
             ELEMENT    *ele,
             INT         iel)
{
INT i,j,l;
DOUBLE dum;

#ifdef DEBUG
dstrc_enter("f3_jaco");
#endif

/*-------------------------------- determine jacobian at point r,s,t ---*/
for (i=0; i<3; i++)
{
   for (j=0; j<3; j++)
   {
      dum=ZERO;
      for (l=0; l<iel; l++)
      {
         dum += deriv[i][l]*xyze[j][l];
      }
      xjm[i][j]=dum;
   } /* end loop j */
} /* end loop i */
/*------------------------------------------ determinant of jacobian ---*/
*det = xjm[0][0]*xjm[1][1]*xjm[2][2]+
       xjm[0][1]*xjm[1][2]*xjm[2][0]+
       xjm[0][2]*xjm[1][0]*xjm[2][1]-
       xjm[0][2]*xjm[1][1]*xjm[2][0]-
       xjm[0][0]*xjm[1][2]*xjm[2][1]-
       xjm[0][1]*xjm[1][0]*xjm[2][2];

if(*det<ZERO)
{
   printf("\n");
   printf("GLOBAL ELEMENT %i\n",ele->Id);
   printf("NEGATIVE JACOBIAN DETERMINANT: %lf\n",*det);
#ifdef PARALLEL
   dserror("Stopped not regulary!\n");
#else
#ifdef D_FSI
   if (genprob.probtyp==prb_fluid && genprob.numfld==2)
   {
      fluid_mf(99);
      dserror("Stopped regulary!\n");
   }
   else if (genprob.probtyp==prb_fsi)
   {
      dyn_fsi(99);
      dserror("Stopped regulary!\n");
   }
   else dserror("Stopped not regulary!\n");
#else
   dserror("Stopped not regulary!\n");
#endif /* endif D_FSI */
#endif /* endif PARALLEL */
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f3_jaco */

/*!---------------------------------------------------------------------
\brief jacobian matrix

<pre>                                                         genk 02/04

In this routine the jacobian matrix and its determinant is calculated

</pre>
\param **xyze      DOUBLE   (i)    nodal coordinates
\param **deriv     DOUBLE   (i)    natural deriv. of shape funcs
\param **xjm       DOUBLE   (o)    jacobian matrix
\param  *det       DOUBLE   (o)    determinant of jacobian matrix
\param  *iedgnod   INT      (i)    edgenodes
\param   iel       INT      (i)    num. of nodes of act. ele
\param  *ele       ELEMENT  (i)    actual element
\return void

------------------------------------------------------------------------*/
void f3_edgejaco( DOUBLE    **xyze,
                  DOUBLE    **deriv,
                  DOUBLE    **xjm,
                  DOUBLE     *det,
                  INT        *iedgnod,
                  INT         iel,
                  ELEMENT    *ele
               )

{
INT k;
INT node;

#ifdef DEBUG
dstrc_enter("f3_edgejaco");
#endif

/*---------------------------------- determine jacobian at point r,s ---*/
xjm[0][0] = ZERO ;
xjm[0][1] = ZERO ;
xjm[1][0] = ZERO ;
xjm[1][1] = ZERO ;

for (k=0; k<iel; k++) /* loop all nodes of the element */
{
     node = iedgnod[k];
     xjm[0][0] += deriv[0][k] * xyze[0][node] ;
     xjm[0][1] += deriv[0][k] * xyze[1][node] ;
     xjm[1][0] += deriv[1][k] * xyze[0][node] ;
     xjm[1][1] += deriv[1][k] * xyze[1][node] ;
} /* end loop over iel */

/*------------------------------------------ determinant of jacobian ---*/
*det = xjm[0][0]* xjm[1][1] - xjm[1][0]* xjm[0][1];

if(*det<ZERO)
{
   printf("\n");
   printf("GLOBAL ELEMENT %i\n",ele->Id);
   printf("NEGATIVE JACOBIAN DETERMINANT: %lf\n",*det);
#ifdef PARALLEL
   dserror("Stopped not regulary!\n");
#else
#ifdef D_FSI
   if (genprob.probtyp==prb_fluid && genprob.numfld==2)
   {
      fluid_mf(99);
      dserror("Stopped regulary!\n");
   }
   else if (genprob.probtyp==prb_fsi)
   {
      dyn_fsi(99);
      dserror("Stopped regulary!\n");
   }
   else dserror("Stopped not regulary!\n");
#else
   dserror("Stopped not regulary!\n");
#endif /* endif D_FSI */
#endif /* endif PARALLEL */
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f3_edgejaco */

/*!---------------------------------------------------------------------
\brief calculate metric of element edge

<pre>                                                         genk 04/04

In this routine the metric of the element edge is calculated

</pre>
\return void

------------------------------------------------------------------------*/
void f3_tvmr(DOUBLE   **x,
                DOUBLE     akov[3][3],
                DOUBLE    *funct,
                DOUBLE   **deriv,
                INT       *iedgnod,
                INT        iel)
{
INT    idim,ialpha,inode,node;

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("f3_tvmr");
#endif

/*------------------------------------ interpolation of kovariant a1,a2 */
for (ialpha=0; ialpha<2; ialpha++)
{
   for (idim=0; idim<3; idim++)
   {
      akov[idim][ialpha]=ZERO;
      for (inode=0; inode<iel; inode++)
      {
         node=iedgnod[inode];
         akov[idim][ialpha] +=
            deriv[ialpha][inode] * x[idim][node];
      }
   }
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of f3_tvmr */

/*!---------------------------------------------------------------------
\brief global derivates

<pre>                                                         genk 05/02

In this routine the global derivatives w.r.t. x,y,z at point r,s,t are
calculated.

</pre>
\param **derxy     DOUBLE   (o)    global derivatives wrt. x/y/z
\param **deriv     DOUBLE   (i)    derivatives of shape functions
\param **xjm       DOUBLE   (i)    jacobian matrix
\param **xji       DOUBLE   (-)    inverse of jacobian
\param   det       DOUBLE   (i)    jacobian determinant
\param   iel       INT      (i)    number of nodes in actual element
\return void

------------------------------------------------------------------------*/
void f3_gder(
               DOUBLE   **derxy,
               DOUBLE   **deriv,
               DOUBLE   **xjm,
               DOUBLE   **xji,
               DOUBLE     det,
               INT        iel
            )
{
INT    k;

#ifdef DEBUG
dstrc_enter("f3_gder");
#endif

/*------------------------------------------------------- initialistion */
for(k=0;k<iel;k++)
{
   derxy[0][k]=ZERO;
   derxy[1][k]=ZERO;
   derxy[2][k]=ZERO;
} /* end of loop over k */

/*------------------------------------------------- inverse of jacobian */
xji[0][0] = (  xjm[1][1]*xjm[2][2] - xjm[2][1]*xjm[1][2])/det;
xji[1][0] = (- xjm[1][0]*xjm[2][2] + xjm[2][0]*xjm[1][2])/det;
xji[2][0] = (  xjm[1][0]*xjm[2][1] - xjm[2][0]*xjm[1][1])/det;
xji[0][1] = (- xjm[0][1]*xjm[2][2] + xjm[2][1]*xjm[0][2])/det;
xji[1][1] = (  xjm[0][0]*xjm[2][2] - xjm[2][0]*xjm[0][2])/det;
xji[2][1] = (- xjm[0][0]*xjm[2][1] + xjm[2][0]*xjm[0][1])/det;
xji[0][2] = (  xjm[0][1]*xjm[1][2] - xjm[1][1]*xjm[0][2])/det;
xji[1][2] = (- xjm[0][0]*xjm[1][2] + xjm[1][0]*xjm[0][2])/det;
xji[2][2] = (  xjm[0][0]*xjm[1][1] - xjm[1][0]*xjm[0][1])/det;

/*---------------------------------------- calculate global derivatives */
for (k=0;k<iel;k++)
{
   derxy[0][k] +=   xji[0][0] * deriv[0][k] \
                  + xji[0][1] * deriv[1][k] \
                  + xji[0][2] * deriv[2][k] ;
   derxy[1][k] +=   xji[1][0] * deriv[0][k] \
                  + xji[1][1] * deriv[1][k] \
                  + xji[1][2] * deriv[2][k] ;
   derxy[2][k] +=   xji[2][0] * deriv[0][k] \
                  + xji[2][1] * deriv[1][k] \
                  + xji[2][2] * deriv[2][k] ;
} /* end of loop over k */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f3_gder */

/*!---------------------------------------------------------------------
\brief global derivates

<pre>                                                         genk 04/02

In this routine the global derivatives w.r.t. x,y at point r,s are
calculated.

</pre>
\param **derxy     DOUBLE   (o)    global derivatives wrt. x/y
\param **deriv     DOUBLE   (i)    derivatives of shape functions
\param **xjm       DOUBLE   (i)    jacobian matrix
\param   det       DOUBLE   (i)    jacobian determinant
\param   iel       INT      (i)    number of nodes in actual element
\return void

------------------------------------------------------------------------*/
void f3_edgegder(
               DOUBLE   **derxy,
               DOUBLE   **deriv,
               DOUBLE   **xjm,
               DOUBLE     det,
               INT        iel
	    )
{
INT    k;

#ifdef DEBUG
dstrc_enter("f3_edgegder");
#endif

for (k=0;k<iel;k++) /* loop all nodes of the element */
{
   derxy[0][k] =(  xjm[1][1]*deriv[0][k]  + (-xjm[0][1]*deriv[1][k]))/det;
   derxy[1][k] =((-xjm[1][0]*deriv[0][k]) +   xjm[0][0]*deriv[1][k] )/det;
} /* end of loop over iel */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_gder */

/*!---------------------------------------------------------------------
\brief global coordinates

<pre>                                                         genk 05/02

In this routine the global coordinates for given shape function values
are set.

</pre>
\param *funct      DOUBLE   (i)    shape functions
\param *ele        DOUBLE   (i)    actual element
\param  iel        DOUBLE   (i)    number of nodes in act. element
\param *gcoor      DOUBLE   (o)    global coordinates
\return void

------------------------------------------------------------------------*/
void f3_gcoor(
               DOUBLE     *funct,
               ELEMENT    *ele,
               INT         iel,
               DOUBLE     *gcoor
             )
{
INT i;

#ifdef DEBUG
dstrc_enter("f3_gcoor");
#endif

gcoor[0]=ZERO;
gcoor[1]=ZERO;
gcoor[2]=ZERO;

for(i=0;i<iel;i++)
{
   gcoor[0] += funct[i] * ele->node[i]->x[0];
   gcoor[1] += funct[i] * ele->node[i]->x[1];
   gcoor[2] += funct[i] * ele->node[i]->x[2];
} /* end of loop over i */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f3_gcoor */

/*!---------------------------------------------------------------------
\brief second global derivatives

<pre>                                                         genk 05/02

In this routine the second global derivatives w.r.t x/y/z at point r,s,t
are calculated.

</pre>
\param **xyze     DOUBLE   (-)    nodal coordinates
\param **xjm      DOUBLE   (i)    jacobian matrix
\param **bm       DOUBLE   (-)    working array
\param **xder2    DOUBLE   (-)    working array
\param **derxy    DOUBLE   (i)    glob. coord.. deriv.
\param **derxy2   DOUBLE   (o)    2nd. glob. coord. deriv.
\param **deriv2   DOUBLE   (i)    2nd. nat. deriv. of shape funcs
\param   iel      INT      (i)    number of nodes of actual ele
\return void

------------------------------------------------------------------------*/
void f3_gder2(
               DOUBLE     **xyze,
               DOUBLE     **xjm,
               DOUBLE     **bm,
               DOUBLE     **xder2,
               DOUBLE     **derxy,
               DOUBLE     **derxy2,
               DOUBLE     **deriv2,
               INT          iel
	     )
{
INT i,j;
DOUBLE r0,r1,r2,r3,r4,r5;

#ifdef DEBUG
dstrc_enter("f3_gder2");
#endif

/*--------------------------- calculate elements of jacobian_bar matrix */
bm[0][0] = xjm[0][0]*xjm[0][0];
bm[1][0] = xjm[1][0]*xjm[1][0];
bm[2][0] = xjm[2][0]*xjm[2][0];
bm[3][0] = xjm[0][0]*xjm[1][0];
bm[4][0] = xjm[0][0]*xjm[2][0];
bm[5][0] = xjm[1][0]*xjm[2][0];

bm[0][1] = xjm[0][1]*xjm[0][1];
bm[1][1] = xjm[1][1]*xjm[1][1];
bm[2][1] = xjm[2][1]*xjm[2][1];
bm[3][1] = xjm[0][1]*xjm[1][1];
bm[4][1] = xjm[0][1]*xjm[2][1];
bm[5][1] = xjm[1][1]*xjm[2][1];

bm[0][2] = xjm[0][2]*xjm[0][2];
bm[1][2] = xjm[1][2]*xjm[1][2];
bm[2][2] = xjm[2][2]*xjm[2][2];
bm[3][2] = xjm[0][2]*xjm[1][2];
bm[4][2] = xjm[0][2]*xjm[2][2];
bm[5][2] = xjm[1][2]*xjm[2][2];

bm[0][3] = TWO*xjm[0][0]*xjm[0][1];
bm[1][3] = TWO*xjm[1][0]*xjm[1][1];
bm[2][3] = TWO*xjm[2][0]*xjm[2][1];
bm[3][3] = xjm[0][0]*xjm[1][1]+xjm[1][0]*xjm[0][1];
bm[4][3] = xjm[0][0]*xjm[2][1]+xjm[2][0]*xjm[0][1];
bm[5][3] = xjm[1][0]*xjm[2][1]+xjm[2][0]*xjm[1][1];

bm[0][4] = TWO*xjm[0][0]*xjm[0][2];
bm[1][4] = TWO*xjm[1][0]*xjm[1][2];
bm[2][4] = TWO*xjm[2][0]*xjm[2][2];
bm[3][4] = xjm[0][0]*xjm[1][2]+xjm[1][0]*xjm[0][2];
bm[4][4] = xjm[0][0]*xjm[2][2]+xjm[2][0]*xjm[0][2];
bm[5][4] = xjm[1][0]*xjm[2][2]+xjm[2][0]*xjm[1][2];

bm[0][5] = TWO*xjm[0][1]*xjm[0][2];
bm[1][5] = TWO*xjm[1][1]*xjm[1][2];
bm[2][5] = TWO*xjm[2][1]*xjm[2][2];
bm[3][5] = xjm[0][1]*xjm[1][2]+xjm[1][1]*xjm[0][2];
bm[4][5] = xjm[0][1]*xjm[2][2]+xjm[2][1]*xjm[0][2];
bm[5][5] = xjm[1][1]*xjm[2][2]+xjm[2][1]*xjm[1][2];

/*-------------------------------------- inverse of jacobian_bar matrix */
math_unsym_inv6x6(bm);

/*----------------------------------------------------------- initialise*/
for (i=0;i<3;i++)
{
   for (j=0;j<6;j++) xder2[j][i]=ZERO;
}
for (i=0;i<iel;i++)
{
   for (j=0;j<6;j++) derxy2[j][i]=ZERO;
}

/*----------------------- determine 2nd derivatives of coord.-functions */
for (i=0;i<iel;i++)
{
   xder2[0][0] += deriv2[0][i] * xyze[0][i];
   xder2[1][0] += deriv2[1][i] * xyze[0][i];
   xder2[2][0] += deriv2[2][i] * xyze[0][i];
   xder2[3][0] += deriv2[3][i] * xyze[0][i];
   xder2[4][0] += deriv2[4][i] * xyze[0][i];
   xder2[5][0] += deriv2[5][i] * xyze[0][i];

   xder2[0][1] += deriv2[0][i] * xyze[1][i];
   xder2[1][1] += deriv2[1][i] * xyze[1][i];
   xder2[2][1] += deriv2[2][i] * xyze[1][i];
   xder2[3][1] += deriv2[3][i] * xyze[1][i];
   xder2[4][1] += deriv2[4][i] * xyze[1][i];
   xder2[5][1] += deriv2[5][i] * xyze[1][i];

   xder2[0][2] += deriv2[0][i] * xyze[2][i];
   xder2[1][2] += deriv2[1][i] * xyze[2][i];
   xder2[2][2] += deriv2[2][i] * xyze[2][i];
   xder2[3][2] += deriv2[3][i] * xyze[2][i];
   xder2[4][2] += deriv2[4][i] * xyze[2][i];
   xder2[5][2] += deriv2[5][i] * xyze[2][i];
} /* end of loop over i */

/*--------------------------------- calculate second global derivatives */
for (i=0;i<iel;i++)
{
   r0 = deriv2[0][i] - xder2[0][0]*derxy[0][i] - xder2[0][1]*derxy[1][i] \
                     - xder2[0][2]*derxy[2][i];
   r1 = deriv2[1][i] - xder2[1][0]*derxy[0][i] - xder2[1][1]*derxy[1][i] \
                     - xder2[1][2]*derxy[2][i];
   r2 = deriv2[2][i] - xder2[2][0]*derxy[0][i] - xder2[2][1]*derxy[1][i] \
                     - xder2[2][2]*derxy[2][i];
   r3 = deriv2[3][i] - xder2[3][0]*derxy[0][i] - xder2[3][1]*derxy[1][i] \
                     - xder2[3][2]*derxy[2][i];
   r4 = deriv2[4][i] - xder2[4][0]*derxy[0][i] - xder2[4][1]*derxy[1][i] \
                     - xder2[4][2]*derxy[2][i];
   r5 = deriv2[5][i] - xder2[5][0]*derxy[0][i] - xder2[5][1]*derxy[1][i] \
                     - xder2[5][2]*derxy[2][i];

   derxy2[0][i] += bm[0][0]*r0 + bm[0][1]*r1 + bm[0][2]*r2 \
                +  bm[0][3]*r3 + bm[0][4]*r4 + bm[0][5]*r5;
   derxy2[1][i] += bm[1][0]*r0 + bm[1][1]*r1 + bm[1][2]*r2 \
                +  bm[1][3]*r3 + bm[1][4]*r4 + bm[1][5]*r5;
   derxy2[2][i] += bm[2][0]*r0 + bm[2][1]*r1 + bm[2][2]*r2 \
                +  bm[2][3]*r3 + bm[2][4]*r4 + bm[2][5]*r5;
   derxy2[3][i] += bm[3][0]*r0 + bm[3][1]*r1 + bm[3][2]*r2 \
                +  bm[3][3]*r3 + bm[3][4]*r4 + bm[3][5]*r5;
   derxy2[4][i] += bm[4][0]*r0 + bm[4][1]*r1 + bm[4][2]*r2 \
                +  bm[4][3]*r3 + bm[4][4]*r4 + bm[4][5]*r5;
   derxy2[5][i] += bm[5][0]*r0 + bm[5][1]*r1 + bm[5][2]*r2 \
                +  bm[5][3]*r3 + bm[5][4]*r4 + bm[5][5]*r5;
} /* end of loop over i */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f3_gder2 */

#endif

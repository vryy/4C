/*!----------------------------------------------------------------------
\file
\brief contains the routine
 - s9_xint:  superposition of weighted function values
 and the functions
 - DOUBLE s9con: which calculates the continuity coefficients of thickness
                 coordinate according to Dis. Braun (p.116)
 - DOUBLE s9notr:decide whether jlay lies on trajectory of reference layer
                 to ilay  -> very similar to s9con
 - DOUBLE s9ksi: independent part of ksi
                 (dependent part of ksi is integrated over the thickness
                 in the material tensor)


<pre>
Maintainer: Stefan Hartmann
            hartmann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hartmann/
            0711 - 685-6120
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SHELL9
#include "../headers/standardtypes.h"
#include "shell9.h"

/*!
\addtogroup SHELL9
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief superposition of weighted function values

<pre>                     m.gee 6/01              modified by    sh 10/02
This routine calculates the superposition of weighted function values
</pre>
\param  DOUBLE *result (o)  result of superposition
\param  DOUBLE *values (i)  values at nodal points
\param  DOUBLE *funct  (i)  shape functions at GP
\param  INT     iel    (1)  number of nodes to this element

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: s9static_keug() [s9_static_keug.c]
                             s9_stress()     [s9_stress.c]

*----------------------------------------------------------------------*/
void s9_xint(DOUBLE *result, DOUBLE *values, DOUBLE *funct, INT iel)
{
INT              i;
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("s9_xint");
#endif
/*----------------------------------------------------------------------*/
*result=0.0;
for (i=0; i<iel; i++) *result += (funct[i]*values[i]);
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of s9_xint */



/*!----------------------------------------------------------------------
\brief calculate the continuity coeficients

<pre>                                                            sh 10/02
This function calculates calculate the continuity coefficients
of thickness coordinate
.........................................................................
.      -> pic. 7.4 on p. 116 of Dis. Braun                              .
.                                                                       .
.    if ilay==jlay:                                                     .
.             zeta = (0.5/condfac)*(e3 +- (1-kronecker(jlay,ref.lay)))  .
.    if ilay!=jlay:                                                     .
.             zeta = +- (1/condfac)*(1-0.5*kronecker(jlay,ref.lay))     .
.........................................................................
</pre>
\param  DOUBLE    e3     (i)  zeta_kl
\param  INT       numlay (i)  number of kinematic layers
\param  INT       ilay   (i)  kin layer where point is located
\param  INT       jlay   (1)  kin layer over which is looped

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: s9_tmtr()       [s9_tmtr.c]
                             s9_tvhe()       [s9_tvhe.c]
                             s9_tvma()       [s9_tvma.c]
                             s9_vthv()       [s9_vthv.c]
                             s9_ans_tvhe_q() [s9_ans.c]
                             s9out_nodal_dis(); s9_out_gid_allcoords(); s9_out_gid_sol_dis() [s9_out.c]

*----------------------------------------------------------------------*/
DOUBLE s9con(DOUBLE    e3,       /* zeta_kl */
             INT       numlay,   /* number of kinematic layers */
             INT       ilay,     /* kin layer where point is located */
             INT       jlay,     /* kin layer over which is looped */
             DOUBLE    condfac)
{
DOUBLE        e33;
DOUBLE        s9con;
DOUBLE        delta;
DOUBLE        fac;
INT           mid;    /* reference layer or the one under the reference area */
INT           mod;
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("s9con");
#endif
/*----------------------------------------------------------------------*/
e33   = 1.0;
delta = 0.0;
fac   = 0.0;
/*--------set ilay, jlay +1 --------------------------------------------*/
ilay += 1;
jlay += 1;
/*---------- check if coefficient is zero ------------------------------*/
mod = numlay % 2;    /*divide modulo 2*/

if (mod == 0)
{
   mid = numlay/2;
   if (ilay <= mid)       /* actual layer is situated under the ref. area */
   {
      fac = -1.0;
      if (jlay > mid ) e33 = 0.0;
      if (jlay < ilay) e33 = 0.0;
   }
   if (ilay > mid)        /* actual layer is situated above the ref. area */
   {
      fac = 1.0;
      if (jlay <= mid) e33 = 0.0;
      if (jlay > ilay) e33 = 0.0;
   }
}
else                    /* uneven number of kinematic layers */
{
   mid = (numlay+1)/2;
   if (jlay == mid) delta = 1.0; /*jlay equals the ref. layer */

   if (ilay < mid)       /* actual layer is situated under the ref. layer */
   {
      fac = -1.0;
      if (jlay > mid ) e33 = 0.0;
      if (jlay < ilay) e33 = 0.0;
   }
   if (ilay == mid)       /* actual layer is the reference layer */
   {
      if (jlay != ilay) e33 = 0.0;
   }
   if (ilay > mid)        /* actual layer is situated above the ref. layer */
   {
      fac = 1.0;
      if (jlay < mid ) e33 = 0.0;
      if (jlay > ilay) e33 = 0.0;
   }
}

if (e33 == 0.0)
{
   s9con = e33;
   goto end;
}
/*----------------------------------------------------------------------*/
if (ilay == jlay)
{
   e3 = (0.5/condfac) * ( e3 + fac * (1.0 - delta) );
}
else
{
   e3 = (fac/condfac) * ( 1.0 - 0.5 * delta );
}
s9con = e3;
/*----------------------------------------------------------------------*/
end:
/*A3_IST_EINHALB*/ /*s9con = 2.0 * s9con;*/
s9con = (1.0/A3FAC_SHELL9) * s9con;
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return s9con;
} /* end of s9con */




/*!----------------------------------------------------------------------
\brief decide whether jlay lies on trajectory of reference layer to ilay

<pre>                                                            sh 10/02
This function decides whether jlay lies on trajectory of reference layer
to ilay  -> very similar to s9con
</pre>
\param  INT       numlay (i)  number of kinematic layers
\param  INT       ilay   (i)  kin layer where point is located
\param  INT       jlay   (1)  kin layer over which is looped

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: s9_tvbo()       [s9_tvbo.c]
                             s9_tvkg()       [s9_tvkg.c]
                             s9_ans_bbar_q() [s9_ans.c]

*----------------------------------------------------------------------*/
DOUBLE s9notr(INT       numlay,   /* number of kinematic layers */
              INT       ilay,     /* kin layer where point is located */
              INT       jlay)     /* kin layer over which is looped */
{
DOUBLE        e33;
DOUBLE        s9notr;
INT           mid;    /* reference layer or the one under the reference area */
INT           mod;
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("s9notr");
#endif
/*----------------------------------------------------------------------*/
e33   = 1.0;
/*--------set ilay, jlay +1 --------------------------------------------*/
ilay += 1;
jlay += 1;
/*---------- check if coefficient is zero ------------------------------*/
mod = numlay % 2;    /*divide modulo 2*/

if (mod == 0)
{
   mid = numlay/2;
   if (ilay <= mid)       /* actual layer is situated under the ref. area */
   {
      if (jlay > mid ) e33 = 0.0;
      if (jlay < ilay) e33 = 0.0;
   }
   if (ilay > mid)        /* actual layer is situated above the ref. area */
   {
      if (jlay <= mid) e33 = 0.0;
      if (jlay > ilay) e33 = 0.0;
   }
}
else                    /* uneven number of kinematic layers */
{
   mid = (numlay+1)/2;
   if (ilay < mid)       /* actual layer is situated under the ref. layer */
   {
      if (jlay > mid ) e33 = 0.0;
      if (jlay < ilay) e33 = 0.0;
   }
   if (ilay == mid)       /* actual layer is the reference layer */
   {
      if (jlay != ilay) e33 = 0.0;
   }
   if (ilay > mid)        /* actual layer is situated above the ref. layer */
   {
      if (jlay < mid ) e33 = 0.0;
      if (jlay > ilay) e33 = 0.0;
   }
}

s9notr = e33;

/*A3_IST_EINHALB*/ /*s9notr = 2.0 * s9notr;*/
s9notr = (1.0/A3FAC_SHELL9) * s9notr;
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return s9notr;
} /* end of s9notr */


/*!----------------------------------------------------------------------
\brief independent part of ksi

<pre>                                                            sh 10/02
This function gets the independent part of ksi (dependent part of ksi is
integrated over the thickness in the material tensor)
decides whether jlay lies on trajectory of reference layer  to ilay
-> very similar to s9con
.........................................................................
.                                                                       .
. ilay == jlay:                                                         .
.        zeta = +- (0.5/condfac)*(1-Kronecker(jlay,ref.lay) )           .
. ilay != jlay:                                                         .
.        zeta = +- (1/condfac)*(1- 0.5*Kronecker(jlay,ref.lay) )        .
.                                                                       .
.........................................................................
</pre>
\param  INT       numlay (i)  number of kinematic layers
\param  INT       ilay   (i)  kin layer where point is located
\param  INT       jlay   (1)  kin layer over which is looped

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: s9_tvbo()       [s9_tvbo.c]
                             s9_ans_bbar_q() [s9_ans.c]

*----------------------------------------------------------------------*/
DOUBLE s9ksi(INT       numlay,   /* number of kinematic layers */
             INT       ilay,     /* kin layer where point is located */
             INT       jlay,     /* kin layer over which is looped */
             DOUBLE    condfac)
{
DOUBLE        s9ksi;
DOUBLE        delta;
DOUBLE        fac;
INT           mid;    /* reference layer or the one under the reference area */
INT           mod;
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("s9ksi");
#endif
/*--------set ilay, jlay +1 --------------------------------------------*/
ilay += 1;
jlay += 1;
/*---------- check if coefficient is zero ------------------------------*/
mod = numlay % 2;    /*divide modulo 2*/

if (mod == 0)
{
   delta = 0.0;
}
else                    /* uneven number of kinematic layers */
{
   mid = (numlay+1)/2;
   if (jlay == mid)       /*  jlay is equal to the reference layer */
   {
      delta = 1.0;
   }
   else
   {
      delta = 0.0;
   }
}
/*--------------- lower or upper half ---------------------------------*/
if (ilay <= numlay/2 ) fac = -1.0;
else                   fac =  1.0;
/*----------------------------------------------------------------------*/
if (ilay == jlay)
{
   s9ksi = fac * ( 0.5/condfac) * (1.0 - delta);
}
else
{
   s9ksi = (fac/condfac) * (1.0 - 0.5*delta);
}
/*----------------------------------------------------------------------*/
/*A3_IST_EINHALB*/ /*s9ksi = 2.0 * s9ksi;*/
s9ksi = (1.0/A3FAC_SHELL9) * s9ksi;
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return s9ksi;
} /* end of s9ksi */
/*----------------------------------------------------------------------*/
#endif /*D_SHELL9*/
/*! @} (documentation module close)*/

/*!----------------------------------------------------------------------
\file
\brief contains the routine
 - s9_cdia: which calculates the equivalent element length of one element

*----------------------------------------------------------------------*/
#ifndef CCADISCRET
#ifdef D_SHELL9
#include "../headers/standardtypes.h"
#include "shell9.h"

/*!
\addtogroup SHELL9
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief calculates the equivalent internal length

<pre>                                                            sh 10/03
This routine calculates the equivalent internal length of one element
                   ___      ___
                  \        \
dia = chi* sqrt ( /___  *  /___ * det J * w(lr) * w(ls)  )  (Dis. Menrath p.68)
                   lr       ls

with chi = sqrt(2) for Quad4
           1       for Quad8/9

Note: for this calculation, the JACOBIAN is calculated in the middle surface
      -> no need to account for different layers
        -> numklay == 1, nummlay = 1;
</pre>
\param  ELEMENT   *ele     (i)  the element structure
\param  S9_DATA   *data    (i)  element integration data
\param  DOUBLE    *funct   (-)  shape functions at GP
\param  DOUBLE   **deriv   (-)  shape function derivatives at GP
\param  DOUBLE   **xjm     (-)  jacobian matrix
\param  DOUBLE   **x       (-)  coordinates of nodal points (global coordinate system)
\param  DOUBLE  ***a3r     (-)  normed (not shared) director in ref/cur config. (->s9a3ref_extern)
\
\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: s9init() [s9_init.c]

*----------------------------------------------------------------------*/
void s9_cdia(ELEMENT   *ele,
             S9_DATA   *data,
             DOUBLE    *funct,
             DOUBLE   **deriv,
             DOUBLE   **xjm,
             DOUBLE   **x,
             DOUBLE  ***a3r)
{
INT      k;
INT      iel;                 /* numnp to this element */
INT      nir,nis;             /* num GP in r/s direction */
INT      lr,ls;               /* loopers over GP */
DOUBLE   facr,facs;           /* weights at GP */
DOUBLE   e1,e2;               /*GP-coords*/

DOUBLE **a3ref;               /* elements directors (lenght 1) */
DOUBLE   hte[MAXNOD_SHELL9];  /* element thickness at nodal points */
DOUBLE   h2;
DOUBLE   condfac;             /* sdc conditioning factor */
DOUBLE   klayhgt_dummy[1];

DOUBLE   det;                 /*jacobi-determinant*/
DOUBLE   dia;                 /*equivalent element length*/

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("s9_cdia");
#endif
/*----------------------------------------------------------------------*/
iel     = ele->numnp;
nir     = ele->e.s9->nGP[0];
nis     = ele->e.s9->nGP[1];
condfac = ele->e.s9->sdc;

a3ref   = ele->e.s9->a3ref.a.da;
klayhgt_dummy[0] = 100.0;

/*-- initialize ---*/
dia = 0.0;

for (k=0; k<iel; k++)           /*loop over all nodes per layer*/
{
   hte[k] = ele->e.s9->thick_node.a.dv[k];
   /*if (ele->e.s9->dfield == 0)*/      /*Layerthicknes, norm(a3L) = HL */
   h2 = ele->e.s9->thick_node.a.dv[k] * condfac;
   /*h2 = 0.5*h2;*/ /*A3_IST_EINHALB halber Direktor*/
   h2 = A3FAC_SHELL9 * h2;
   /*else if (ele->e.s9->dfield == 1)*/ /*half of shell thickness, norm(a3) = H/2*/
   /*  h2 = ele->e.s9->thick_node.a.dv[k]/2. * condfac;*/

   a3r[0][k][0] = a3ref[0][k] * h2;
   a3r[1][k][0] = a3ref[1][k] * h2;
   a3r[2][k][0] = a3ref[2][k] * h2;

   x[0][k] = ele->node[k]->x[0];
   x[1][k] = ele->node[k]->x[1];
   x[2][k] = ele->node[k]->x[2];

} /*end loop over nodes*/

for (lr=0; lr<nir; lr++)   /* loop in r-direction */
{
   /*================================== gaussian point and weight at it */
   e1   = data->xgpr[lr];
   facr = data->wgtr[lr];
   for (ls=0; ls<nis; ls++)  /* loop in s-direction*/
   {
      /*=============================== gaussian point and weight at it */
      e2   = data->xgps[ls];
      facs = data->wgts[ls];
      /*-------------------- shape functions at gp e1,e2 on mid surface */
      s9_funct_deriv(funct,deriv,e1,e2,ele->distyp,1);
      /*------------------------------------ compute jacobian matrix ---*/
      s9jaco(funct,deriv,x,xjm,hte,a3r,0.0,iel,&det,0,1,klayhgt_dummy);

      dia += det * facr*facs;
    }/*============================================= end of loop over ls */
}/*================================================ end of loop over lr */

if      (iel == 4)             dia = sqrt(2.) * sqrt(dia);
else if (iel == 8 || iel == 9) dia =            sqrt(dia);
else dserror("only Quad4/8/9 implemented -> s9_cdia.c");

/*-- write dia to element ----------------------------------------------*/
ele->e.s9->dia = dia;

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of s9_cdia */
/*----------------------------------------------------------------------*/
#endif /*D_SHELL9*/
/*! @} (documentation module close)*/

#endif

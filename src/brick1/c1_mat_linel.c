/*!----------------------------------------------------------------------
\file
\brief contains the routine 'c1_mat_linel' to establish local material law
       stress-strain law for isotropic material for a 3D hex element

<pre>
Maintainer: Andreas Lipka
            lipka@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/lipka/
            0711 - 685-6575
</pre>

*----------------------------------------------------------------------*/
#ifndef CCADISCRET
#ifdef D_BRICK1

#include "../headers/standardtypes.h"
#include "brick1.h"
#include "brick1_prototypes.h"

/*!
\addtogroup BRICK1
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief establish local material

<pre>                                                              al 06/02
This routine to establish local material law
       stress-strain law for isotropic material for for a 3D-hex-element.

</pre>
\param       youngs   DOUBLE  (i)   young's modulus
\param possionratio   DOUBLE  (i)   poisson's ratio
\param            d   DOUBLE**(o)   constitutive matrix

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: c1_cint()

*----------------------------------------------------------------------*/
void c1_mat_linel(DOUBLE youngs,
                  DOUBLE possionratio,
                  DOUBLE **d)
{
DOUBLE d1,d2,d3;
DOUBLE ym,pv;/*------------------------------------------ mat constants */
#ifdef DEBUG
dstrc_enter("c1_mat_linel");
#endif
/*----------------------------------------------------------------------*/
ym  = youngs;
pv  = possionratio;
/*----------------------------------- evaluate basic material values ---*/
d1=ym*(1.0 - pv)/((1.0 + pv)*(1.0 - 2.0*pv));
d2=ym*pv/((1.0 + pv)*(1.0 - 2.0*pv));
d3=ym/((1.0 + pv)*2.0);
/*------------------------------------ set values in material-matrix ---*/
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
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of c1_mat_linel */

/*!----------------------------------------------------------------------
\brief establish local porous material

<pre>                                                              al 06/02
This routine to establish local porous material law
       stress-strain law for isotropic material for for a 3D-hex-element.

</pre>
\param       youngs   DOUBLE    (i)   young's modulus
\param      matdata MATERIAL*   (i)   material data
\param            d   DOUBLE**  (i)   constitutive matrix

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: c1_call_mat()

*----------------------------------------------------------------------*/
void c1_mat_stvpor(MATERIAL  *mat,
                   DOUBLE *matdata,
                   DOUBLE **d)
{
DOUBLE d1,d2,d3;
/*------------------------------------------ mat constants */
DOUBLE ym, pv, dn, rd, ex;
#ifdef DEBUG
dstrc_enter("c1_mat_stvpor");
#endif
/*----------------------------------------------------------------------*/
  /* input value: dn = mat->m.stvenpor->density     ;*/
  /* current values */
  dn = matdata[0];


  ym = mat->m.stvenpor->youngs      ;
  pv = mat->m.stvenpor->possionratio;
  rd = mat->m.stvenpor->refdens     ;
  ex = mat->m.stvenpor->exponent    ;

  ym = pow((dn/rd),ex)*ym;
/*----------------------------------- evaluate basic material values ---*/
d1=ym*(1.0 - pv)/((1.0 + pv)*(1.0 - 2.0*pv));
d2=ym*pv/((1.0 + pv)*(1.0 - 2.0*pv));
d3=ym/((1.0 + pv)*2.0);
/*------------------------------------ set values in material-matrix ---*/
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
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of c1_mat_stvpor */

/*!----------------------------------------------------------------------
\brief establish derivatives of local porous material

<pre>                                                              al 06/02
This routine to establish derivatives of llocal porous material law
       stress-strain law for isotropic material for for a 3D-hex-element.

</pre>
\param       youngs   DOUBLE    (i)   young's modulus
\param      matdata MATERIAL*   (i)   material data
\param            d   DOUBLE**  (i)   constitutive matrix

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: c1_call_matd()

*----------------------------------------------------------------------*/
void c1_matd_stvpor(MATERIAL  *mat,
                   DOUBLE *matdata,
                   DOUBLE **d)
{
DOUBLE d1,d2,d3;
/*------------------------------------------ mat constants */
DOUBLE dym, ym, pv, dn, rd, ex;
#ifdef DEBUG
dstrc_enter("c1_matd_stvpor");
#endif
/*----------------------------------------------------------------------*/
  /* input value: dn = mat->m.stvenpor->density     ;*/
  /* current values */
  dn = matdata[0];


  ym = mat->m.stvenpor->youngs      ;
  pv = mat->m.stvenpor->possionratio;
  rd = mat->m.stvenpor->refdens     ;
  ex = mat->m.stvenpor->exponent    ;

  dym=pow((dn/rd),(ex-1.));
  dym=(ym*ex*dym)/rd;

/* derivative for material law: */
/*  ym = pow((dn/rd),ex)*ym; */
/*----------------------------------- evaluate basic material values ---*/
d1=dym*(1.0 - pv)/((1.0 + pv)*(1.0 - 2.0*pv));
d2=dym*pv/((1.0 + pv)*(1.0 - 2.0*pv));
d3=dym/((1.0 + pv)*2.0);
/*------------------------------------ set values in material-matrix ---*/
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
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of c1_matd_stvpor */

/*!----------------------------------------------------------------------
\brief establish derivatives of local porous material

<pre>                                                              al 06/02
This routine to establish derivatives of llocal porous material law
       stress-strain law for isotropic material for for a 3D-hex-element.

</pre>
\param       youngs   DOUBLE    (i)   young's modulus
\param      matdata MATERIAL*   (i)   material data
\param            d   DOUBLE**  (i)   constitutive matrix

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: c1_call_matd()

*----------------------------------------------------------------------*/
void c1_mat_mfoc(  MATERIAL  *mat,
                   DOUBLE *matdata,
                   DOUBLE **d)
{
DOUBLE d1,d2,d3;
/*------------------------------------------ mat constants */
DOUBLE es, ym, pr, dens, denss, oce, ocf, denmin, denmax;
#ifdef DEBUG
dstrc_enter("c1_mat_mfoc");
#endif
/*----------------------------------------------------------------------*/
  /* input value: dn = mat->m.mfoc->dens     ;*/
  /* current values */
  dens = matdata[0];


  es = mat->m.mfoc->es      ;
  pr = mat->m.mfoc->pr;
/*dn = mat->m.mfoc->dens     ;*/
  denss = mat->m.mfoc->denss     ;
  oce     =  mat->m.mfoc->oce;
  ocf     =  mat->m.mfoc->ocf;
  denmin  =  mat->m.mfoc->denmin;
  denmax  =  mat->m.mfoc->denmax;

  ym = ocf * pow((dens/denss),oce)*es;
/*----------------------------------- evaluate basic material values ---*/
d1=ym*(1.0 - pr)/((1.0 + pr)*(1.0 - 2.0*pr));
d2=ym*pr/((1.0 + pr)*(1.0 - 2.0*pr));
d3=ym/((1.0 + pr)*2.0);
/*------------------------------------ set values in material-matrix ---*/
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
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of c1_mat_mfoc */

/*!----------------------------------------------------------------------
\brief establish derivatives of local porous material

<pre>                                                              al 06/02
This routine to establish derivatives of llocal porous material law
       stress-strain law for isotropic material for for a 3D-hex-element.

</pre>
\param       youngs   DOUBLE    (i)   young's modulus
\param      matdata MATERIAL*   (i)   material data
\param            d   DOUBLE**  (i)   constitutive matrix

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: c1_call_matd()

*----------------------------------------------------------------------*/
void c1_mat_mfcc(  MATERIAL  *mat,
                   DOUBLE *matdata,
                   DOUBLE **d)
{
DOUBLE d1,d2,d3;
/*------------------------------------------ mat constants */
DOUBLE es, ym, pr, dens, denss, cce, ccf, denmin, denmax;
#ifdef DEBUG
dstrc_enter("c1_mat_mfcc");
#endif
/*----------------------------------------------------------------------*/
  /* input value: dn = mat->m.mfcc->dens     ;*/
  /* current values */
  dens = matdata[0];


  es = mat->m.mfcc->es      ;
  pr = mat->m.mfcc->pr;
/*dn = mat->m.mfcc->dens     ;*/
  denss = mat->m.mfcc->denss     ;
  cce     =  mat->m.mfcc->cce;
  ccf     =  mat->m.mfcc->ccf;
  denmin  =  mat->m.mfcc->denmin;
  denmax  =  mat->m.mfcc->denmax;

  ym = ccf *( 0.5 * pow((dens/denss),cce) + 0.3 * dens/denss)*es;
/*----------------------------------- evaluate basic material values ---*/
d1=ym*(1.0 - pr)/((1.0 + pr)*(1.0 - 2.0*pr));
d2=ym*pr/((1.0 + pr)*(1.0 - 2.0*pr));
d3=ym/((1.0 + pr)*2.0);
/*------------------------------------ set values in material-matrix ---*/
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
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of c1_mat_mfcc */

/*!----------------------------------------------------------------------
\brief establish derivatives of local porous material

<pre>                                                              al 06/02
This routine to establish derivatives of llocal porous material law
       stress-strain law for isotropic material for for a 3D-hex-element.

</pre>
\param       youngs   DOUBLE    (i)   young's modulus
\param      matdata MATERIAL*   (i)   material data
\param            d   DOUBLE**  (i)   constitutive matrix

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: c1_call_matd()

*----------------------------------------------------------------------*/
void c1_matd_mfoc( MATERIAL  *mat,
                   DOUBLE *matdata,
                   DOUBLE **d)
{
/*----------------------------------------------------------------------*/
DOUBLE d1,d2,d3;
DOUBLE es, dym, pr, dens, denss, oce, ocf, denmin, denmax;
#ifdef DEBUG
dstrc_enter("c1_matd_mfoc");
#endif
/*----------------------------------------------------------------------*/
  /* input value: dn = mat->m.mfoc->dens     ;*/
  /* current values */
  dens = matdata[0];


  es = mat->m.mfoc->es      ;
  pr = mat->m.mfoc->pr;
/*dn = mat->m.mfoc->dens     ;*/
  denss = mat->m.mfoc->denss     ;
  oce     =  mat->m.mfoc->oce;
  ocf     =  mat->m.mfoc->ocf;
  denmin  =  mat->m.mfoc->denmin;
  denmax  =  mat->m.mfoc->denmax;

/*----------------------------------- evaluate basic material values ---*/
/* derivative for material law: */
dym= ocf * oce/dens * pow((dens/denss),(oce)) * es;

d1=dym*(1.0 - pr)/((1.0 + pr)*(1.0 - 2.0*pr));
d2=dym*pr/((1.0 + pr)*(1.0 - 2.0*pr));
d3=dym/((1.0 + pr)*2.0);
/*------------------------------------ set values in material-matrix ---*/
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
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of c1_matd_mfoc */

/*!----------------------------------------------------------------------
\brief establish derivatives of local porous material

<pre>                                                              al 06/02
This routine to establish derivatives of llocal porous material law
       stress-strain law for isotropic material for for a 3D-hex-element.

</pre>
\param       youngs   DOUBLE    (i)   young's modulus
\param      matdata MATERIAL*   (i)   material data
\param            d   DOUBLE**  (i)   constitutive matrix

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: c1_call_matd()

*----------------------------------------------------------------------*/
void c1_matd_mfcc( MATERIAL  *mat,
                   DOUBLE *matdata,
                   DOUBLE **d)
{
/*----------------------------------------------------------------------*/
DOUBLE d1,d2,d3;
DOUBLE es, dym, pr, dens, denss, cce, ccf, denmin, denmax;
#ifdef DEBUG
dstrc_enter("c1_matd_mfcc");
#endif
/*----------------------------------------------------------------------*/
  /* input value: dn = mat->m.mfcc->dens     ;*/
  /* current values */
  dens = matdata[0];


  es = mat->m.mfcc->es      ;
  pr = mat->m.mfcc->pr;
/*dn = mat->m.mfcc->dens     ;*/
  denss = mat->m.mfcc->denss     ;
  cce     =  mat->m.mfcc->cce;
  ccf     =  mat->m.mfcc->ccf;
  denmin  =  mat->m.mfcc->denmin;
  denmax  =  mat->m.mfcc->denmax;

/*----------------------------------- evaluate basic material values ---*/
/* derivative for material law: */
dym= ccf * (dens/pow(denss,cce)+ 0.3/denss) * es;

d1=dym*(1.0 - pr)/((1.0 + pr)*(1.0 - 2.0*pr));
d2=dym*pr/((1.0 + pr)*(1.0 - 2.0*pr));
d3=dym/((1.0 + pr)*2.0);
/*------------------------------------ set values in material-matrix ---*/
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
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of c1_matd_mfcc */

/*!----------------------------------------------------------------------
\brief evaluate stresses for elastic material

<pre>                                                              al 06/02
This routine to evaluate stresses for elastic material for a 3D-hex-element.

</pre>
\param      strain   DOUBLE*    (i)   global strains
\param           d   MATERIAL** (i)   material matrices
\param      stress   DOUBLE*    (o)   forces/moments/additional terms

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: c1_call_matd()

*----------------------------------------------------------------------*/
void c1mefm(DOUBLE *strain, /* global strains                           */
            DOUBLE **d,     /* material matrices                        */
            DOUBLE *stress) /* forces/moments/additional terms          */
{
/*----------------------------------------------------------------------*/
INT    i, j;
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("c1mefm");
#endif
/*-------------------------------------------------- global stresses ---*/
  for (i=0; i<6; i++) stress[i] = 0.0;
  for (i=0; i<6; i++) for (j=0; j<6; j++) stress[i] += d[i][j]*strain[j];
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of c1mefm */
/*----------------------------------------------------------------------*/

/*!----------------------------------------------------------------------
\brief establish of local porous material

<pre>                                                              al 06/02
This routine to establish of llocal porous material law
       stress-strain law for orthotropic material for for a 3D-hex-element.

</pre>
\param     emod1   DOUBLE  (i)   young's modulus
\param     emod2   DOUBLE  (i)   young's modulus
\param     emod3   DOUBLE  (i)   young's modulus
\param    xnue23   DOUBLE  (i)   poisson's ratio
\param    xnue13   DOUBLE  (i)   poisson's ratio
\param    xnue12   DOUBLE  (i)   poisson's ratio
\param    gmod12   DOUBLE  (i)   shear modulus
\param    gmod23   DOUBLE  (i)   shear modulus
\param    gmod13   DOUBLE  (i)   shear modulus
\param         c   DOUBLE**(o)   constitutive matrix

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: c1_call_mat()

*----------------------------------------------------------------------*/
void c1_mat_elorth(DOUBLE   emod1 ,
                   DOUBLE   emod2 ,
                   DOUBLE   emod3 ,
                   DOUBLE   xnue23,
                   DOUBLE   xnue13,
                   DOUBLE   xnue12,
                   DOUBLE   gmod12,
                   DOUBLE   gmod23,
                   DOUBLE   gmod13,
                   DOUBLE      **c)
{
INT i,j;
DOUBLE xnue31        = 0.;
DOUBLE xnue32        = 0.;
DOUBLE xnue21        = 0.;
DOUBLE emod          = 0.;
DOUBLE delta         = 0.;
#ifdef DEBUG
dstrc_enter("c1_mat_elorth");
#endif
/*----------------------------------------------------------------------*/
  for (i=0; i<6; i++) for (j=0; j<6; j++) c[i][j] = 0.;
/*----------------------------------------------------------------------*/
  xnue31=xnue13*emod3/emod1;
  xnue32=xnue23*emod3/emod2;
  xnue21=xnue12*emod2/emod1;
  emod=emod1*emod2*emod3;
  delta=1.-xnue13*xnue31-xnue23*xnue32-xnue12*xnue21;
  delta=(delta-2.*xnue31*xnue23*xnue12)/emod;
  c[0][0]=(1.-xnue23*xnue32)/(emod2*emod3*delta);
  c[1][1]=(1.-xnue13*xnue31)/(emod1*emod3*delta);
  c[2][2]=(1.-xnue12*xnue21)/(emod1*emod2*delta);
  c[1][0]=(xnue12+xnue13*xnue32)/(emod1*emod3*delta);
  c[2][0]=(xnue13+xnue12*xnue23)/(emod1*emod2*delta);
  c[2][1]=(xnue23+xnue21*xnue13)/(emod1*emod2*delta);
  c[3][3]=gmod12;
  c[5][5]=gmod13;
  c[4][4]=gmod23;
  c[0][1]=c[1][0];
  c[0][2]=c[2][0];
  c[1][2]=c[2][1];
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of c1_mat_elorth */
#endif
/*! @} (documentation module close)*/
#endif

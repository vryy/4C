/*!----------------------------------------------------------------------
\file
\brief contains the routine 'c1_mat_neohook' to establish local material law
       stress-strain law for isotropic material for a 3D hex element

<pre>
Maintainer: Andreas Lipka
            lipka@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/lipka/
            0711 - 685-6575
</pre>

*----------------------------------------------------------------------*/
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
\param       youngs   DOUBLE    (i)   young's modulus
\param possionratio   DOUBLE    (i)   poisson's ratio
\param         disd   DOUBLE*   (i)   displacement derivatives
\param       stress   DOUBLE*   (o)   ele stress (-resultant) vector
\param            d   DOUBLE**  (o)   constitutive matrix

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: c1_cint()

*----------------------------------------------------------------------*/
void c1_mat_neohook(DOUBLE youngs,     /* young's modulus              */
                    DOUBLE possionratio, /* poisson's ratio            */
                    DOUBLE *disd,      /* displacement derivatives     */
                    DOUBLE *stress,  /* ele stress (-resultant) vector */
                    DOUBLE **d)      /* material matrix                */
{
/*----------------------------------------------------------------------*/
INT i,j,k,l;
/*----------------------------------------------------------------------*/
DOUBLE ym,pv;/*------------------------------------------ mat constants */
/*----------------------------------------------------------------------*/
DOUBLE xl, g, xj, J, f1;
DOUBLE sp[3][3];
DOUBLE c[3][3][3][3];
DOUBLE     gk[3][3] = {1.,0.,0.,0.,1.,0.,0.,0.,1.};
/*----------------------------------------------------------------------*/
INT    cci;
DOUBLE invFN[9]; /* inverse deformation gradient            */
DOUBLE F[3][3];  /* deformation gradient                    */
DOUBLE FTF[9];
DOUBLE FFT[9];
DOUBLE C[9];
DOUBLE Ct[3][3];
DOUBLE b[9];
DOUBLE bt[3][3];
DOUBLE Wene,Ib;
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("c1_mat_neohook");
#endif
/*----------------------------------------------------------------------*/
ym  = youngs;
pv  = possionratio;
/*----------------------------------------------------------------------*/
/* determine def. gradient and inverse */
    disd[0] += 1.;
    disd[1] += 1.;
    disd[2] += 1.;
    c1invf(disd,invFN,&J);
    J=1./J;

  F[0][0] = disd[0]; /* [1,1] */
  F[1][1] = disd[1]; /* [2,2] */
  F[2][2] = disd[2]; /* [3,3] */
  F[0][1] = disd[3]; /* [2,1] */
  F[1][0] = disd[4]; /* [1,2] */
  F[1][2] = disd[5]; /* [3,2] */
  F[2][1] = disd[6]; /* [2,3] */
  F[0][2] = disd[7]; /* [3,1] */
  F[2][0] = disd[8]; /* [1,3] */

  /* right cauchy green C: C = FT F = U^2 */
  for (i=0; i<9; i++)   FTF[i] = 0.0;
  cci=0;
  for (i=0; i<3; i++)
  {
    for (j=0; j<3; j++)
    {
      for (k=0; k<3; k++)
      {
        FTF[cci] += F[k][j]*F[k][i];
      }
      cci++;
    }
  }
  for (i=0; i<9; i++)   C[i] = FTF[i];

  c1inv3(C);
  cci = 0; /* copy back to ansiC style */
  for (i=0; i<3; i++) for (j=0; j<3; j++)  Ct[j][i] = C[cci++];


  /* left cauchy green b: b = F FT = V^2 */
  for (i=0; i<9; i++)   FFT[i] = 0.0;
  cci=0;
  for (i=0; i<3; i++)
  {
    for (j=0; j<3; j++)
    {
      for (k=0; k<3; k++)
      {
        FFT[cci] += F[i][k]*F[j][k];
      }
      cci++;
    }
  }
  for (i=0; i<9; i++)   b[i] = FFT[i];
  bt[0][0] = b[0];
  bt[1][0] = b[1];
  bt[2][0] = b[2];
  bt[0][1] = b[3];
  bt[1][1] = b[4];
  bt[2][1] = b[5];
  bt[0][2] = b[6];
  bt[1][0] = b[7];
  bt[2][2] = b[8];
/*----------------------------------------------------------------------*/
  /* first invariant */
  Ib = bt[0][0]+bt[1][1]+bt[2][2];
/*----------------------------------------------------------------------*/
  xl=ym*pv/(1.0+pv)/(1.0-2.0*pv); /* lambda */
  g=ym/2.0/(1.0+pv);              /* mue =G */
  xj=J;
  if(xj<1.0E-6)xj=1.0E-6;
  f1=xl*log(xj)-g; /* ansiC(log)==ln */
/*----------------------------------------------------------------------*/
  Wene = 0.5*xl*(log(xj))*(log(xj))-g*log(xj)+0.5*g*(Ib-3.0);
/*--------------------------------------------------- pk2-spannungen ---*/
  for (i=0; i<3; i++)
  {
    for (j=0; j<3; j++)
    {
       sp[i][j]=f1*Ct[i][j]+g*gk[i][j];
  }}

  stress[0]=sp[0][0];
  stress[1]=sp[1][1];
  stress[2]=sp[2][2];
  stress[3]=sp[0][1];
  stress[4]=sp[1][2];
  stress[5]=sp[0][2];
/*-------------------------------------------------- werkstofftensor ---*/
  for (k=0; k<3; k++)
  {
    for (l=0; l<3; l++)
    {
      c[0][0][k][l]= xl*Ct[0][0]*Ct[k][l] -
                     f1*(Ct[0][k]*Ct[0][l]+Ct[0][l]*Ct[0][k]) ;

      c[1][0][k][l]= xl*Ct[1][0]*Ct[k][l] -
                     f1*(Ct[1][k]*Ct[0][l]+Ct[1][l]*Ct[0][k]) ;

      c[1][1][k][l]= xl*Ct[1][1]*Ct[k][l] -
                     f1*(Ct[1][k]*Ct[1][l]+Ct[1][l]*Ct[1][k]) ;

      c[2][0][k][l]= xl*Ct[2][0]*Ct[k][l] -
                     f1*(Ct[2][k]*Ct[0][l]+Ct[2][l]*Ct[0][k]) ;

      c[2][1][k][l]= xl*Ct[2][1]*Ct[k][l] -
                     f1*(Ct[2][k]*Ct[1][l]+Ct[2][l]*Ct[1][k]) ;

      c[2][2][k][l]= xl*Ct[2][2]*Ct[k][l] -
                     f1*(Ct[2][k]*Ct[2][l]+Ct[2][l]*Ct[2][k]) ;
  }}

 d[0][0]=c[0][0][0][0];
 d[0][1]=c[0][0][1][1];
 d[0][2]=c[0][0][2][2];
 d[0][3]=c[0][0][1][0];
 d[0][4]=c[0][0][2][1];
 d[0][5]=c[0][0][2][0];

 d[1][0]=c[1][1][0][0];
 d[1][1]=c[1][1][1][1];
 d[1][2]=c[1][1][2][2];
 d[1][3]=c[1][1][1][0];
 d[1][4]=c[1][1][2][1];
 d[1][5]=c[1][1][2][0];

 d[2][0]=c[2][2][0][0];
 d[2][1]=c[2][2][1][1];
 d[2][2]=c[2][2][2][2];
 d[2][3]=c[2][2][1][0];
 d[2][4]=c[2][2][2][1];
 d[2][5]=c[2][2][2][0];

 d[3][0]=c[1][0][0][0];
 d[3][1]=c[1][0][1][1];
 d[3][2]=c[1][0][2][2];
 d[3][3]=c[1][0][1][0];
 d[3][4]=c[1][0][2][1];
 d[3][5]=c[1][0][2][0];

 d[4][0]=c[2][1][0][0];
 d[4][1]=c[2][1][1][1];
 d[4][2]=c[2][1][2][2];
 d[4][3]=c[2][1][1][0];
 d[4][4]=c[2][1][2][1];
 d[4][5]=c[2][1][2][0];

 d[5][0]=c[2][0][0][0];
 d[5][1]=c[2][0][1][1];
 d[5][2]=c[2][0][2][2];
 d[5][3]=c[2][0][1][0];
 d[5][4]=c[2][0][2][1];
 d[5][5]=c[2][0][2][0];
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of c1_mat_neohook */

/*!----------------------------------------------------------------------
\brief establish local porous NeoHook-material

<pre>                                                              al 06/02
This routine to establish local material law
       stress-strain law for isotropic material for for a 3D-hex-element.

</pre>
\param          mat   MATERIAL* (i)   material data ...
\param      matdata   DOUBLE*   (i)   variable mat. data, opti.
\param         disd   DOUBLE*   (i)   displacement derivatives
\param       stress   DOUBLE*   (o)   ele stress (-resultant) vector
\param            d   DOUBLE**  (o)   constitutive matrix
\param          ste   DOUBLE*   (o)   energy

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: c1_cint()

*----------------------------------------------------------------------*/
void c1_mat_nhmfcc( MATERIAL  *mat,     /* material data ...           */
                    DOUBLE *matdata,    /* variable mat. data, opti.   */
                    DOUBLE *disd,      /* displacement derivatives     */
                    DOUBLE *stress,  /* ele stress (-resultant) vector */
                    DOUBLE **d,         /* material matrix             */
                    DOUBLE *ste)
{
/*----------------------------------------------------------------------*/
INT i,j,k,l;
/*----------------------------------------------------------------------*/
DOUBLE es, ym, pv, dens, denss, cce, ccf, denmin, denmax;
/*----------------------------------------------------------------------*/
DOUBLE xl, g, xj, J, f1;
DOUBLE sp[3][3];
DOUBLE c[3][3][3][3];
DOUBLE gk[3][3] = {1.,0.,0.,0.,1.,0.,0.,0.,1.};
/*----------------------------------------------------------------------*/
INT    cci;
DOUBLE invFN[9]; /* inverse deformation gradient            */
DOUBLE F[3][3];  /* deformation gradient                    */
DOUBLE FTF[9];
DOUBLE FFT[9];
DOUBLE C[9];
DOUBLE Ct[3][3];
DOUBLE b[9];
DOUBLE bt[3][3];
DOUBLE Wene,Ib,Ic;
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("c1_mat_nhmfcc");
#endif
/*----------------------------------------------------------------------*/
  /* input value: dn = mat->m.mfcc->dens     ;*/
  /* current values */
  dens    = matdata[0];
  es      = mat->m.nhmfcc->es      ;
  pv      = mat->m.nhmfcc->pr;
  denss   = mat->m.nhmfcc->denss     ;
  cce     =  mat->m.nhmfcc->cce;
  ccf     =  mat->m.nhmfcc->ccf;
  denmin  =  mat->m.nhmfcc->denmin;
  denmax  =  mat->m.nhmfcc->denmax;

  ym = ccf *( 0.5 * pow((dens/denss),cce) + 0.3 * dens/denss)*es;
/*----------------------------------------------------------------------*/
/* determine def. gradient and inverse */
    disd[0] += 1.;
    disd[1] += 1.;
    disd[2] += 1.;
    c1invf(disd,invFN,&J);
    J=1./J;

  for (i=0; i<3; i++) for (j=0; j<3; j++)   F[i][j] = 0.0;
  F[0][0] = disd[0]; /* [1,1] */
  F[1][1] = disd[1]; /* [2,2] */
  F[2][2] = disd[2]; /* [3,3] */
  F[0][1] = disd[3]; /* [2,1] */
  F[1][0] = disd[4]; /* [1,2] */
  F[1][2] = disd[5]; /* [3,2] */
  F[2][1] = disd[6]; /* [2,3] */
  F[0][2] = disd[7]; /* [3,1] */
  F[2][0] = disd[8]; /* [1,3] */

  /* right cauchy green C: C = FT F = U^2 */
  for (i=0; i<9; i++)   FTF[i] = 0.0;
  cci=0;
  for (i=0; i<3; i++)
  {
    for (j=0; j<3; j++)
    {
      for (k=0; k<3; k++)
      {
        FTF[cci] += F[k][j]*F[k][i];
      }
      cci++;
    }
  }
  for (i=0; i<9; i++)   C[i] = FTF[i];

  /* first invariant */
  Ic = C[0]+C[4]+C[8];

  c1inv3(C);
  cci = 0; /* copy to ansiC style */
  for (i=0; i<3; i++) for (j=0; j<3; j++)  Ct[j][i] = C[cci++];


  /* left cauchy green b: b = F FT = V^2 */
  for (i=0; i<9; i++)   FFT[i] = 0.0;
  cci=0;
  for (i=0; i<3; i++)
  {
    for (j=0; j<3; j++)
    {
      for (k=0; k<3; k++)
      {
        FFT[cci] += F[i][k]*F[j][k];
      }
      cci++;
    }
  }
  for (i=0; i<9; i++)   b[i] = FFT[i];
  bt[0][0] = b[0];
  bt[1][0] = b[1];
  bt[2][0] = b[2];
  bt[0][1] = b[3];
  bt[1][1] = b[4];
  bt[2][1] = b[5];
  bt[0][2] = b[6];
  bt[1][0] = b[7];
  bt[2][2] = b[8];
/*----------------------------------------------------------------------*/
  /* first invariant */
  Ib = bt[0][0]+bt[1][1]+bt[2][2];
/*----------------------------------------------------------------------*/
  xl=ym*pv/(1.0+pv)/(1.0-2.0*pv); /* lambda */
  g=ym/2.0/(1.0+pv);              /* mue =G */
  xj=J;
  if(xj<1.0E-6)xj=1.0E-6;
  f1=xl*log(xj)-g; /* ansiC(log)==ln */
/*----------------------------------------------------------------------*/
  Wene = 0.5*xl*(log(xj))*(log(xj))-g*log(xj)+0.5*g*(Ib-3.0);
  *ste = Wene;
/*--------------------------------------------------- pk2-spannungen ---*/
  for (i=0; i<3; i++)
  {
    for (j=0; j<3; j++)
    {
       sp[i][j]=f1*Ct[i][j]+g*gk[i][j];
  }}
  /* */
  stress[0]=sp[0][0];
  stress[1]=sp[1][1];
  stress[2]=sp[2][2];
  stress[3]=sp[0][1];
  stress[4]=sp[1][2];
  stress[5]=sp[0][2];
/*-------------------------------------------------- werkstofftensor ---*/
  for (k=0; k<3; k++)
  {
    for (l=0; l<3; l++)
    {
      c[0][0][k][l]= xl*Ct[0][0]*Ct[k][l] -
                     f1*(Ct[0][k]*Ct[0][l]+Ct[0][l]*Ct[0][k]) ;

      c[1][0][k][l]= xl*Ct[1][0]*Ct[k][l] -
                     f1*(Ct[1][k]*Ct[0][l]+Ct[1][l]*Ct[0][k]) ;

      c[1][1][k][l]= xl*Ct[1][1]*Ct[k][l] -
                     f1*(Ct[1][k]*Ct[1][l]+Ct[1][l]*Ct[1][k]) ;

      c[2][0][k][l]= xl*Ct[2][0]*Ct[k][l] -
                     f1*(Ct[2][k]*Ct[0][l]+Ct[2][l]*Ct[0][k]) ;

      c[2][1][k][l]= xl*Ct[2][1]*Ct[k][l] -
                     f1*(Ct[2][k]*Ct[1][l]+Ct[2][l]*Ct[1][k]) ;

      c[2][2][k][l]= xl*Ct[2][2]*Ct[k][l] -
                     f1*(Ct[2][k]*Ct[2][l]+Ct[2][l]*Ct[2][k]) ;
  }}
/**/
 d[0][0]=c[0][0][0][0];
 d[0][1]=c[0][0][1][1];
 d[0][2]=c[0][0][2][2];
 d[0][3]=c[0][0][1][0];
 d[0][4]=c[0][0][2][1];
 d[0][5]=c[0][0][2][0];

 d[1][0]=c[1][1][0][0];
 d[1][1]=c[1][1][1][1];
 d[1][2]=c[1][1][2][2];
 d[1][3]=c[1][1][1][0];
 d[1][4]=c[1][1][2][1];
 d[1][5]=c[1][1][2][0];

 d[2][0]=c[2][2][0][0];
 d[2][1]=c[2][2][1][1];
 d[2][2]=c[2][2][2][2];
 d[2][3]=c[2][2][1][0];
 d[2][4]=c[2][2][2][1];
 d[2][5]=c[2][2][2][0];

 d[3][0]=c[1][0][0][0];
 d[3][1]=c[1][0][1][1];
 d[3][2]=c[1][0][2][2];
 d[3][3]=c[1][0][1][0];
 d[3][4]=c[1][0][2][1];
 d[3][5]=c[1][0][2][0];

 d[4][0]=c[2][1][0][0];
 d[4][1]=c[2][1][1][1];
 d[4][2]=c[2][1][2][2];
 d[4][3]=c[2][1][1][0];
 d[4][4]=c[2][1][2][1];
 d[4][5]=c[2][1][2][0];

 d[5][0]=c[2][0][0][0];
 d[5][1]=c[2][0][1][1];
 d[5][2]=c[2][0][2][2];
 d[5][3]=c[2][0][1][0];
 d[5][4]=c[2][0][2][1];
 d[5][5]=c[2][0][2][0];
/*----------------------------------------------------------------------*/
 disd[0] -= 1.;
 disd[1] -= 1.;
 disd[2] -= 1.;
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of c1_mat_nhmfcc */

/*!----------------------------------------------------------------------
\brief establish derivatives of  local foam NeoHook - material

<pre>                                                              al 06/02
This routine to establish derivatives of  local material law
       stress-strain law for isotropic material for for a 3D-hex-element.

</pre>
\param          mat   MATERIAL* (i)   material data ...
\param      matdata   DOUBLE*   (i)   variable mat. data, opti.
\param         disd   DOUBLE*   (i)   displacement derivatives
\param       stress   DOUBLE*   (o)   ele stress (-resultant) vector
\param            d   DOUBLE**  (o)   derivative of constitutive matrix
\param          ste   DOUBLE*   (o)   energy

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: c1_cint()

*----------------------------------------------------------------------*/
void c1_dmat_nhmfcc( MATERIAL  *mat,
                    DOUBLE *matdata,
                    DOUBLE *disd,
                    DOUBLE *stress,
                    DOUBLE **d,
                    DOUBLE *ste)
{
/*----------------------------------------------------------------------*/
INT i,j,k,l;
/*----------------------------------------------------------------------*/
DOUBLE es, ym, dym, pv, dens, denss, cce, ccf, denmin, denmax;
/*----------------------------------------------------------------------*/
DOUBLE xl, g, xj, J, f1;
DOUBLE sp[3][3];
DOUBLE c[3][3][3][3];
DOUBLE gk[3][3] = {1.,0.,0.,0.,1.,0.,0.,0.,1.};
/*----------------------------------------------------------------------*/
INT    cci;
DOUBLE invFN[9]; /* inverse deformation gradient            */
DOUBLE F[3][3];  /* deformation gradient                    */
DOUBLE FTF[9];
DOUBLE FFT[9];
DOUBLE C[9];
DOUBLE Ct[3][3];
DOUBLE b[9];
DOUBLE bt[3][3];
DOUBLE Wene,Ib,Ic;
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("c1_dmat_nhmfcc");
#endif
/*----------------------------------------------------------------------*/
  /* input value: dn = mat->m.mfcc->dens     ;*/
  /* current values */
  dens = matdata[0];
  es = mat->m.nhmfcc->es      ;
  pv = mat->m.nhmfcc->pr;
  denss = mat->m.nhmfcc->denss     ;
  cce     =  mat->m.nhmfcc->cce;
  ccf     =  mat->m.nhmfcc->ccf;
  denmin  =  mat->m.nhmfcc->denmin;
  denmax  =  mat->m.nhmfcc->denmax;
/*----------------------------------------------------------------------*/
/*  ym = ccf *( 0.5 * pow((dens/denss),cce) + 0.3 * dens/denss)*es; */
/*----------------------------------------------------------------------*/
  ym = ccf *( 0.5 * pow((dens/denss),cce) + 0.3 * dens/denss)*es;
/*----------------------------------- evaluate basic material values ---*/
/* derivative for material law: */
  dym= ccf * (0.5*cce*pow((dens/denss),(cce-1.0))+ 0.3/denss) * es;
  /*dym= ccf * (dens/pow(denss,cce)+ 0.3/denss) * es;*/
  ym = dym;
/*----------------------------------------------------------------------*/
/* determine def. gradient and inverse */
    disd[0] += 1.;
    disd[1] += 1.;
    disd[2] += 1.;
    c1invf(disd,invFN,&J);
    J=1./J;

  for (i=0; i<3; i++) for (j=0; j<3; j++)   F[i][j] = 0.0;
  F[0][0] = disd[0]; /* [1,1] */
  F[1][1] = disd[1]; /* [2,2] */
  F[2][2] = disd[2]; /* [3,3] */
  F[0][1] = disd[3]; /* [2,1] */
  F[1][0] = disd[4]; /* [1,2] */
  F[1][2] = disd[5]; /* [3,2] */
  F[2][1] = disd[6]; /* [2,3] */
  F[0][2] = disd[7]; /* [3,1] */
  F[2][0] = disd[8]; /* [1,3] */

  /* right cauchy green C: C = FT F = U^2 */
  for (i=0; i<9; i++)   FTF[i] = 0.0;
  cci=0;
  for (i=0; i<3; i++)
  {
    for (j=0; j<3; j++)
    {
      for (k=0; k<3; k++)
      {
        FTF[cci] += F[k][j]*F[k][i];
      }
      cci++;
    }
  }
  for (i=0; i<9; i++)   C[i] = FTF[i];


  /* first invariant */
  Ic = C[0]+C[4]+C[8];


  c1inv3(C);
  cci = 0; /* copy to ansiC style */
  for (i=0; i<3; i++) for (j=0; j<3; j++)  Ct[j][i] = C[cci++];


  /* left cauchy green b: b = F FT = V^2 */
  for (i=0; i<9; i++)   FFT[i] = 0.0;
  cci=0;
  for (i=0; i<3; i++)
  {
    for (j=0; j<3; j++)
    {
      for (k=0; k<3; k++)
      {
        FFT[cci] += F[i][k]*F[j][k];
      }
      cci++;
    }
  }
  for (i=0; i<9; i++)   b[i] = FFT[i];
  bt[0][0] = b[0];
  bt[1][0] = b[1];
  bt[2][0] = b[2];
  bt[0][1] = b[3];
  bt[1][1] = b[4];
  bt[2][1] = b[5];
  bt[0][2] = b[6];
  bt[1][0] = b[7];
  bt[2][2] = b[8];
/*----------------------------------------------------------------------*/
  /* first invariant */
  Ib = bt[0][0]+bt[1][1]+bt[2][2];
/*----------------------------------------------------------------------*/
  xl=ym*pv/(1.0+pv)/(1.0-2.0*pv); /* lambda */
  g=ym/2.0/(1.0+pv);              /* mue =G */
  xj=J;
  if(xj<1.0E-6)xj=1.0E-6;
  f1=xl*log(xj)-g; /* ansiC(log)==ln */
/*----------------------------------------------------------------------*/
  Wene = 0.5*xl*(log(xj))*(log(xj))-g*log(xj)+0.5*g*(Ib-3.0);
  *ste = Wene;
/*--------------------------------------------------- pk2-spannungen ---*/
  /* lala*/
  for (i=0; i<3; i++)
  {
    for (j=0; j<3; j++)
    {
       sp[i][j]=f1*Ct[i][j]+g*gk[i][j];
  }}
  /* */
  stress[0]=sp[0][0]; /* ? */
  stress[1]=sp[1][1];
  stress[2]=sp[2][2];
  stress[3]=sp[0][1];
  stress[4]=sp[1][2];
  stress[5]=sp[0][2];
/*-------------------------------------------------- werkstofftensor ---*/
  for (k=0; k<3; k++)
  {
    for (l=0; l<3; l++)
    {
      c[0][0][k][l]= xl*Ct[0][0]*Ct[k][l] -
                     f1*(Ct[0][k]*Ct[0][l]+Ct[0][l]*Ct[0][k]) ;

      c[1][0][k][l]= xl*Ct[1][0]*Ct[k][l] -
                     f1*(Ct[1][k]*Ct[0][l]+Ct[1][l]*Ct[0][k]) ;

      c[1][1][k][l]= xl*Ct[1][1]*Ct[k][l] -
                     f1*(Ct[1][k]*Ct[1][l]+Ct[1][l]*Ct[1][k]) ;

      c[2][0][k][l]= xl*Ct[2][0]*Ct[k][l] -
                     f1*(Ct[2][k]*Ct[0][l]+Ct[2][l]*Ct[0][k]) ;

      c[2][1][k][l]= xl*Ct[2][1]*Ct[k][l] -
                     f1*(Ct[2][k]*Ct[1][l]+Ct[2][l]*Ct[1][k]) ;

      c[2][2][k][l]= xl*Ct[2][2]*Ct[k][l] -
                     f1*(Ct[2][k]*Ct[2][l]+Ct[2][l]*Ct[2][k]) ;
  }}
/**/
 d[0][0]=c[0][0][0][0];
 d[0][1]=c[0][0][1][1];
 d[0][2]=c[0][0][2][2];
 d[0][3]=c[0][0][1][0];
 d[0][4]=c[0][0][2][1];
 d[0][5]=c[0][0][2][0];

 d[1][0]=c[1][1][0][0];
 d[1][1]=c[1][1][1][1];
 d[1][2]=c[1][1][2][2];
 d[1][3]=c[1][1][1][0];
 d[1][4]=c[1][1][2][1];
 d[1][5]=c[1][1][2][0];

 d[2][0]=c[2][2][0][0];
 d[2][1]=c[2][2][1][1];
 d[2][2]=c[2][2][2][2];
 d[2][3]=c[2][2][1][0];
 d[2][4]=c[2][2][2][1];
 d[2][5]=c[2][2][2][0];

 d[3][0]=c[1][0][0][0];
 d[3][1]=c[1][0][1][1];
 d[3][2]=c[1][0][2][2];
 d[3][3]=c[1][0][1][0];
 d[3][4]=c[1][0][2][1];
 d[3][5]=c[1][0][2][0];

 d[4][0]=c[2][1][0][0];
 d[4][1]=c[2][1][1][1];
 d[4][2]=c[2][1][2][2];
 d[4][3]=c[2][1][1][0];
 d[4][4]=c[2][1][2][1];
 d[4][5]=c[2][1][2][0];

 d[5][0]=c[2][0][0][0];
 d[5][1]=c[2][0][1][1];
 d[5][2]=c[2][0][2][2];
 d[5][3]=c[2][0][1][0];
 d[5][4]=c[2][0][2][1];
 d[5][5]=c[2][0][2][0];
/*----------------------------------------------------------------------*/
 disd[0] -= 1.;
 disd[1] -= 1.;
 disd[2] -= 1.;
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of c1_dmat_nhmfcc */
#endif
/*! @} (documentation module close)*/

/*!----------------------------------------------------------------------
\file
\brief contains routine for calculating contributions of surface
       phenomena to element stiffness matrix of brick1 (hex8 only)

<pre>
Maintainer: Lena Wiechert
            wiechert@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/wiechert
            089/28915303
</pre>

*----------------------------------------------------------------------*/
#ifndef CCADISCRET

#include "../headers/standardtypes.h"
#include "brick1.h"
#include "brick1_prototypes.h"

static void getnodesgp(
                       INT     gsurfno,
                       INT    *nodevec,
                       INT    *nxgp,
                       INT    *nygp,
                       INT    *nzgp,
                       DOUBLE *xgp,
                       DOUBLE *ygp,
                       DOUBLE *zgp,
                       DOUBLE *wgtx,
                       DOUBLE *wgty,
                       DOUBLE *wgtz,
                       INT    *rstgeo
                      );

static void gamma_calc(
                       DOUBLE *gamma,
                       DOUBLE *dgamma,
                       DOUBLE *con_quot,
                       DOUBLE *A_old,
                       DOUBLE A_new,
                       DOUBLE dt,
                       INT    step,
                       DOUBLE k1,
                       DOUBLE k2,
                       DOUBLE C,
                       DOUBLE m1,
                       DOUBLE m2,
                       DOUBLE gamma_0,
                       DOUBLE gamma_min,
                       DOUBLE gamma_min_eq,
                       DOUBLE con_quot_max
                      );
/*!
\addtogroup BRICK1
*//*! @{ (documentation module open)*/


/*!----------------------------------------------------------------------
\brief routine for calculating stiffness due to constant surface tension
or surfactant

<pre>                                                            lw 04/06
This routine calculates the contribution of interfacial energy to the
element stiffness matrix for the surface nodes of a BRICK HEX8
element.

K_ij=(dgamma/dA)*(dA/du_i)*(dA/du_j)+gamma*[d^2A/(du_i*du_j)]
(cf. Kowe et al. 1986)

Since the applied surfactant model only covers the steady state,
arbitrary starting values are assumed.

</pre>
\param surface_flag  INT      (i)   =0 for surfactant,
                                    =1 for const. surface tension
\param           k1  DOUBLE   (i)   adsorption coefficient (surfactant)
\param           k2  DOUBLE   (i)   desorption coefficient (surfactant)
\param            C  DOUBLE   (i)   bulk surfactant concentration
\param           m1  DOUBLE   (i)   first isotherm slope (surfactant)
\param           m2  DOUBLE   (i)   second isotherm slope (surfactant)
\param      gamma_0  DOUBLE   (i)   gamma_0 (surfactant)
\param    gamma_min  DOUBLE   (i)   minimal surface stress (surfactant)
\param    gamma_min  DOUBLE   (i)   minimal equilibrium surface stress
                                    (surfactant)
\param     xyz_curr  DOUBLE   (i)   current coordinates of element nodes
\param       K_surf  DOUBLE** (o)   stiffness matrix due to surface energy
\param     fie_surf  DOUBLE*  (o)   internal force due to surface energy
\param     node_vec  INT*     (o)   node numbers of gsurf
\param           dt  DOUBLE   (i)   time step size
\param            r  INT      (i)   current gsurf number
\param         step  INT      (i)   current step number
\param        A_old  DOUBLE* (i/o)  interfacial area
\param     con_quot  DOUBLE* (i/o)  surfactant concentration

\warning There is nothing special to this routine
\return void
\sa calling: area_calc(), gamma_calc(); called by: c1_cint.c
*----------------------------------------------------------------------*/

void c1_surf(
             INT     surface_flag,
             DOUBLE  k1,
             DOUBLE  k2,
             DOUBLE  C,
             DOUBLE  m1,
             DOUBLE  m2,
             DOUBLE  gamma_0,
             DOUBLE  gamma_min,
             DOUBLE  gamma_min_eq,
             DOUBLE  const_gamma,
             DOUBLE *xyz_curr,
             DOUBLE **K_surf,
             DOUBLE *fie_surf,
             INT    *nodevec,
             DOUBLE  dt,
             INT     r,
             INT     step,
             DOUBLE *A_old,
             DOUBLE *con_quot
            )
{
  DOUBLE gamma, dgamma;
  DOUBLE A, Adiff[12];
  DOUBLE **Adiff2;
  DOUBLE con_quot_eq;
  DOUBLE con_quot_max;
  INT i, j;
  DOUBLE con, a_old;

  /*-------------------------------------------------------------------*/
  #ifdef DEBUG
  dstrc_enter("c1_surf");
  #endif
  /*-------------------------------------------------------------------*/

  Adiff2=(double**)CCACALLOC(12, sizeof(double*));
  for (i=0;i<12;i++)
  {
    Adiff2[i]=(double*)CCACALLOC(12, sizeof(double));
  }

  /*------------ calculation of derivatives of current interfacial area
   * -------------------------------- with respect to the displacements*/
  area_calc(&A, Adiff, Adiff2, nodevec, xyz_curr, r);

  if (surface_flag==0)                                   /* SURFACTANT */
  {
    con_quot_max=(gamma_min_eq-gamma_min)/m2+1.;
    con_quot_eq=(k1*C)/(k1*C+k2);

    /*------------calculation of current surface stress and its partial
     *-----------------derivative with respect to the interfacial area */
    a_old=*A_old;

    if (step<100)             /* gradual application of surface stress */
    {
      gamma=gamma_0-m1*con_quot_eq;
      con=con_quot_eq;
      *con_quot=con;
      dgamma=0.;
    }
    else
    {
      con=*con_quot;
      gamma_calc(&gamma, &dgamma, &con, &a_old, A, dt, step, k1, k2, C, m1,
               m2, gamma_0, gamma_min, gamma_min_eq, con_quot_max);
      *con_quot=con;
    }
  }

  if (surface_flag==1)                             /* SURFACE TENSION */
  {
    gamma=const_gamma;
    dgamma=0.;
  }

  for (i=0;i<12;i++)
  {
    for (j=0;j<12;j++)
    {
      K_surf[i][j]=dgamma*Adiff[i]*Adiff[j]+gamma*Adiff2[i][j];

      if (step<100)         /* gradual application of surface stress */
      {
        if (step<99)
          K_surf[i][j]*=-2*pow(step/98., 3)+3*pow(step/98., 2);
        else
          K_surf[i][j]*=1.0;
      }
    }
  }

  /*------calculation of current internal force due to surface energy*/
  for (i=0;i<12;i++)
  {
    fie_surf[i]=gamma*Adiff[i];

    if (step<100)           /* gradual application of surface stress */
    {
      if (step<99)
        fie_surf[i]*=-2*pow(step/98., 3)+3*pow(step/98., 2);
      else
        fie_surf[i]*=1.0;
    }
  }

  /*-------------------------------------------------------------update*/
  if (surface_flag==0)
  {
    *A_old=A;
  }

  for (i=0;i<12;i++)
  {
    CCAFREE(Adiff2[i]);
  }
  CCAFREE(Adiff2);

  /*-------------------------------------------------------------------*/
  #ifdef DEBUG
  dstrc_exit();
  #endif
  return;
} /* end of c1_surf */


/*!---------------------------------------------------------------------
\brief routine for calculating surface area and its derivatives

<pre>                                                           lw 04/06
This routine calculates the current interfacial area and its first and
second derivatives with respect to the displacements.

</pre>
\param           A DOUBLE*  (o)   current surface area
\param       Adiff DOUBLE*  (o)   first partial derivatives of area
\param      Adiff2 DOUBLE** (o)   second partial derivative of area
\param     nodevec DOUBLE*  (o)   node numbers of gsurf
\param    xyz_curr DOUBLE*  (i)   current coordinates of nodes
\param    gsurf_no INT      (i)   current gsurf number

\warning There is nothing special to this routine
\return void
\sa calling: getnodesgp(); called by: c1_surf()

*----------------------------------------------------------------------*/

void area_calc(
               DOUBLE *A,
               DOUBLE *Adiff,
               DOUBLE **Adiff2,
               INT    *nodevec,
               DOUBLE *xyz_curr,
               INT     gsurfno
              )
{
  /*-------------------------------------------------------------------*/
  #ifdef DEBUG
  dstrc_enter("area_calc");
  #endif
  /*-------------------------------------------------------------------*/

  /* calculation of area currently only for 4 Gauss points implemented */

  DOUBLE xgp[2], ygp[2], zgp[2];        /* coordinates of Gauss points */
  DOUBLE wgtx[2], wgty[2], wgtz[2];     /* Gaussian weights */
  DOUBLE wgtr, wgts, wgtt;
  INT i, j, l, m, n, o, p, gpr, gps, gpt;
  INT nxgp, nygp, nzgp;                 /* number of Gauss points */
  DOUBLE e1, e2, e3;
  INT rstgeo;                           /* 0=r,1=s,2=t as normal */
  DOUBLE det[3], ddet[3][12], ddet2[3][12][12], jacobi_deriv[12];
  INT var1, var2;
  DOUBLE Jac, temp;
  static INT var=0;
  INT                 iel=8;            /* numnp to this element */
  static ARRAY    funct_a;              /* shape functions */
  static DOUBLE  *funct;
  static ARRAY    deriv_a;              /* derivatives of shape functions */
  static DOUBLE **deriv;
  static ARRAY    xjm_a;                /* jacobian matrix */
  static DOUBLE **xjm;
  DOUBLE det0;

  if (var==0)
  {
    funct     = amdef("funct"  ,&funct_a,8,1 ,"DV");
    deriv     = amdef("deriv"  ,&deriv_a,3,8,"DA");
    xjm       = amdef("xjm"    ,&xjm_a ,3,3,"DA");
  }

  /*------------determination of related surface nodes and Gauss points*/
  getnodesgp(gsurfno, nodevec, &nxgp, &nygp, &nzgp, xgp, ygp, zgp, wgtx,
             wgty, wgtz, &rstgeo);

  /*-----------------------------------------------------initialization*/
  temp=0.;

  for (i=0;i<12;i++)
  {
    Adiff[i]=0.;
    for (j=0;j<12;j++)
    {
      Adiff2[i][j]=0.;
    }
  }

  /*------------------------------- start loop over integration points */
  for (gpr=0; gpr<nxgp; gpr++) {
     e1   = xgp[gpr];                  /* gp-coordinate in r direction */
     wgtr=wgtx[gpr];
  for (gps=0; gps<nygp; gps++) {
     e2   = ygp[gps];                  /* gp-coordinate in s direction */
     wgts=wgty[gps];
  for (gpt=0; gpt<nzgp; gpt++) {
     e3   = zgp[gpt];                  /* gp-coordinate in t direction */
     wgtt=wgtz[gpt];

     /*------------------ derivatives of shape functions at this point */
     c1_funct_deriv(funct,deriv,e1,e2,e3,iel,1);

     /*------------------------------ determination of jacobian matrix */
     c1_jaco (deriv,xjm,&det0,xyz_curr,iel);

     /*-- calculation of appropriate minor determiants of the Jacobian */
     switch (rstgeo)
     {
       case 0:
         det[0] = xjm[1][1]*xjm[2][2] - xjm[2][1]*xjm[1][2];
         det[1] = xjm[1][2]*xjm[2][0] - xjm[2][2]*xjm[1][0];
         det[2] = xjm[1][0]*xjm[2][1] - xjm[2][0]*xjm[1][1];
         Jac = sqrt( det[0]*det[0] + det[1]*det[1] + det[2]*det[2] );
         l=1;
         m=2;
       break;
       case 1:
         det[0] = xjm[0][1]*xjm[2][2] - xjm[2][1]*xjm[0][2];
         det[1] = xjm[0][2]*xjm[2][0] - xjm[2][2]*xjm[0][0];
         det[2] = xjm[0][0]*xjm[2][1] - xjm[2][0]*xjm[0][1];
         Jac = sqrt( det[0]*det[0] + det[1]*det[1] + det[2]*det[2] );
         l=0;
         m=2;
       break;
       case 2:
         det[0] = xjm[0][1]*xjm[1][2] - xjm[1][1]*xjm[0][2];
         det[1] = xjm[0][2]*xjm[1][0] - xjm[1][2]*xjm[0][0];
         det[2] = xjm[0][0]*xjm[1][1] - xjm[1][0]*xjm[0][1];
         Jac = sqrt( det[0]*det[0] + det[1]*det[1] + det[2]*det[2] );
         l=0;
         m=1;
       break;
     }

     temp+=Jac*wgtr*wgts*wgtt;

     /*--------------- derivation of minor determiants of the Jacobian
      *----------------------------- with respect to the displacements */
     for (n=0;n<4;n++)
     {
       i=nodevec[n];

       ddet[0][3*n]=0.;
       ddet[0][3*n+1]=deriv[l][i]*xjm[m][2]-deriv[m][i]*xjm[l][2];
       ddet[0][3*n+2]=deriv[m][i]*xjm[l][1]-deriv[l][i]*xjm[m][1];

       ddet[1][3*n]=deriv[m][i]*xjm[l][2]-deriv[l][i]*xjm[m][2];
       ddet[1][3*n+1]=0.;
       ddet[1][3*n+2]=deriv[l][i]*xjm[m][0]-deriv[m][i]*xjm[l][0];

       ddet[2][3*n]=deriv[l][i]*xjm[m][1]-deriv[m][i]*xjm[l][1];
       ddet[2][3*n+1]=deriv[m][i]*xjm[l][0]-deriv[l][i]*xjm[m][0];
       ddet[2][3*n+2]=0.;

       jacobi_deriv[n*3]=1/Jac*(det[2]*ddet[2][3*n]+det[1]*ddet[1][3*n]);
       jacobi_deriv[n*3+1]=1/Jac*(det[2]*ddet[2][3*n+1]+det[0]*ddet[0][3*n+1]);
       jacobi_deriv[n*3+2]=1/Jac*(det[0]*ddet[0][3*n+2]+det[1]*ddet[1][3*n+2]);

     }

     /*--- calculation of first derivative of current interfacial area
      *----------------------------- with respect to the displacements */
     for (i=0;i<12;i++)
     {
       Adiff[i]+=jacobi_deriv[i]*wgtr*wgts*wgtt;
     }

     /*-------- second derivation of minor determiants of the Jacobian
      *----------------------------- with respect to the displacements */
     for (n=0;n<3;n++)
     {
       for (o=0;o<12;o++)
       {
         for (p=0;p<12;p++)
           ddet2[n][o][p]=0.;
       }
     }

     for (n=0;n<4;n++)
     {
       for (o=0;o<4;o++)
       {
         i=nodevec[n];
         j=nodevec[o];

         ddet2[0][n*3+1][o*3+2]=deriv[l][i]*deriv[m][j]-deriv[m][i]*deriv[l][j];
         ddet2[0][n*3+2][o*3+1]=-ddet2[0][n*3+1][o*3+2];

         ddet2[1][n*3][o*3+2]=deriv[m][i]*deriv[l][j]-deriv[l][i]*deriv[m][j];
         ddet2[1][n*3+2][o*3]=-ddet2[1][n*3][o*3+2];

         ddet2[2][n*3][o*3+1]=deriv[l][i]*deriv[m][j]-deriv[m][i]*deriv[l][j];
         ddet2[2][n*3+1][o*3]=-ddet2[2][n*3][o*3+1];
       }
     }

     /*- calculation of second derivatives of current interfacial area
      *----------------------------- with respect to the displacements */
     for (i=0;i<12;i++)
     {
       if (i==0 || i==3 || i==6 || i==9)
       {
         var1=1;
         var2=2;
       }
       if (i==1 || i==4 || i==7 || i==10)
       {
         var1=0;
         var2=2;
       }
       if (i==2 || i==5 || i==8 || i==11)
       {
         var1=0;
         var2=1;
       }

       for (j=0;j<12;j++)
       {
         Adiff2[i][j]+=(-1/Jac*jacobi_deriv[j]*jacobi_deriv[i]+1/Jac*(ddet[var1][i]*ddet[var1][j]
          +det[var1]*ddet2[var1][i][j]+ddet[var2][i]*ddet[var2][j]+det[var2]*ddet2[var2][i][j]))*wgtr*wgts*wgtt;
       }
     }
  }}}

  *A=temp;
  var++;

  /*-------------------------------------------------------------------*/
  #ifdef DEBUG
  dstrc_exit();
  #endif
  return;
} /* end of area_calc */

/*!---------------------------------------------------------------------
\brief routine for calculating current surface stress

<pre>                                                           lw 04/06
This routine calculates the current surface stress for a finite
element according to Otis'adsorption-limited SURFACTANT model
(Otis et al. 1994). This model only covers the steady state.

</pre>
\param          gamma DOUBLE* (o)  generalized surface energy [N/mm]
\param         dgamma DOUBLE* (o)  derivative of gamma with resp. to A
\param       con_quot DOUBLE*(i/o) non-dimensionalized concentration
\param     d_con_quot DOUBLE*(i/o) derivative of non-dimensionalized
                                   surfactant mass
\param          A_old DOUBLE  (i)  old surface area
\param          A_new DOUBLE  (i)  surface area
\param             dt DOUBLE  (i)  time step size
\param             k1 DOUBLE  (i)  adsorption coefficient [mm^3/(g*ms)]
\param             k2 DOUBLE  (i)  desorption coefficient [1/ms]
\param              C DOUBLE  (i)  bulk concentration [g/mm^3]
\param             m1 DOUBLE  (i)  isotherm slope first regime [N/mm]
\param             m2 DOUBLE  (i)  isotherm slope second regime [N/mm]
\param        gamma_0 DOUBLE  (i)  surface tension of water [N/mm]
\param      gamma_min DOUBLE  (i)  minimum surface energy [N/mm]
\param   gamma_min_eq DOUBLE  (i)  minimum equilibrium surface energy [N/mm]
\param   con_quot_max DOUBLE  (i)  max. non-dimensionalized concentration

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: c1_surf()

*----------------------------------------------------------------------*/

static void gamma_calc(
                       DOUBLE *gamma,
                       DOUBLE *dgamma,
                       DOUBLE *con_quot,
                       DOUBLE *A_old,
                       DOUBLE A_new,
                       DOUBLE dt,
                       INT    step,
                       DOUBLE k1,
                       DOUBLE k2,
                       DOUBLE C,
                       DOUBLE m1,
                       DOUBLE m2,
                       DOUBLE gamma_0,
                       DOUBLE gamma_min,
                       DOUBLE gamma_min_eq,
                       DOUBLE con_quot_max
                      )
{
  /* use of non-dimensionalized interfacial surfactant concentration:
   * con_quot=Gamma/Gamma_min_eq
   * mass of surfactant molecules M_old is also divided
   * by minimum equilibrium interfacial surfactant concentration!*/

  DOUBLE con_quot_new;          /*new non-dimensionalized concentration*/
  DOUBLE M_old;                 /*old mass of surfactant*/
  DOUBLE con;
  DOUBLE a_old;

  /*-------------------------------------------------------------------*/
  #ifdef DEBUG
  dstrc_enter("gamma_calc");
  #endif
  /*-------------------------------------------------------------------*/

  con=*con_quot;
  a_old=*A_old;
  M_old=con*a_old;

  /*-------------------------------------------------------------------*
   |     DETERMINATION OF CURRENT NON-DIMENSIONALIZED INTERFACIAL      |
   |             SURFACTANT CONCENTRATION CON_QUOT_NEW                 |
   *-------------------------------------------------------------------*/

  /*----------------Regime 1: Langmuir kinetics (adsorption/desorption)*/

  if (con<1.)
  {
    /* assumed continuous drop of adsorption/desorption coefficients
     * when approaching insoluble regime (cf. Morris et al. 2001,
     * p108) is only needed for high bulk concentrations C! To
     * simplify matters, drop is left out */

    con_quot_new=(1.0/dt*M_old+k1*C*A_new)/(A_new*(1.0/dt+k1*C+k2));
  }
  else
  {
    /*------------------------------------Regime 2: Insoluble monolayer*/
    if (con<con_quot_max)
    {
      con_quot_new=M_old/A_new;
    }
    /*----------------------------Regime 3: "Squeeze out"/Film collapse*/
    else
    {
      if (A_new<a_old)
      {
        con_quot_new=con_quot_max;
      }
      else
      {
        con_quot_new=M_old/A_new;
      }
    }
  }

  /* interfacial surfactant concentration must not be larger than the
   * maximum interfacial concentration */

  if (con_quot_new>con_quot_max)
  {
    con_quot_new=con_quot_max;
  }

  *con_quot=con_quot_new;

  /*-------------------------------------------------------------------*
   |        DETERMINATION OF CURRENT GENERALIZED SURFACE ENERGY        |
   *-------------------------------------------------------------------*/

  /*----------------Regime 1: Langmuir kinetics (adsorption/desorption)*/

  if (con_quot_new<1.)
  {
    *gamma=gamma_0-m1*con_quot_new;
    *dgamma=m1/dt*M_old/(A_new*A_new*(1/dt+k1*C+k2));
  }
  else
  {
    /*------------------------------------Regime 2: Insoluble monolayer*/
    if (con_quot_new<con_quot_max)
    {
      *gamma=gamma_min_eq-m2*(con_quot_new-1.);
      *dgamma=m2*con_quot_new/A_new;
    }
    /*----------------------------Regime 3: "Squeeze out"/Film collapse*/
    else
    {
      *gamma=gamma_min;
      *dgamma=0.;
    }
  }
 /*-------------------------------------------------------------------*/
  #ifdef DEBUG
  dstrc_exit();
  #endif
  return;
} /* end of gamma_calc */


/*!---------------------------------------------------------------------
\brief routine for determining nodes and gp of current gsurf

\param        gsurfno INT     (i)  current gsurf number
\param        nodevec INT*    (o)  nodes of gsurf
\param           nxgp INT*    (i)  numbers of gp in r-direction
\param           nygp INT*    (i)  numbers of gp in s-direction
\param           nzgp INT*    (i)  numbers of gp in t-direction
\param            xgp DOUBLE* (o)  gp r-direction
\param            ygp DOUBLE* (o)  gp s-direction
\param            zgp DOUBLE* (o)  gp t-direction
\param           wgtx DOUBLE* (o)  gp weights r-direction
\param           wgty DOUBLE* (o)  gp weights s-direction
\param           wgtz DOUBLE* (o)  gp weights t-direction
\param         rstgeo INT*    (o)  direction of surface normal
                                   0=r,1=s,2=t

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: area_calc()

*----------------------------------------------------------------------*/
static void getnodesgp(
                       INT     gsurfno,
                       INT    *nodevec,
                       INT    *nxgp,
                       INT    *nygp,
                       INT    *nzgp,
                       DOUBLE *xgp,
                       DOUBLE *ygp,
                       DOUBLE *zgp,
                       DOUBLE *wgtx,
                       DOUBLE *wgty,
                       DOUBLE *wgtz,
                       INT    *rstgeo
                      )
{
  /*-------------------------------------------------------------------*/
  #ifdef DEBUG
  dstrc_enter("getnodesgp");
  #endif
  /*-------------------------------------------------------------------*/

   switch (gsurfno)
   {
   case 0: case 5:
     *nxgp=2;
     *nygp=2;
     *nzgp=1;

     if (gsurfno==5)
     {
       zgp[0]=1.;
       nodevec[0]=7;
       nodevec[1]=6;
       nodevec[2]=5;
       nodevec[3]=4;
     }
     else
     {
       zgp[0]=-1.;
       nodevec[0]=0;
       nodevec[1]=1;
       nodevec[2]=2;
       nodevec[3]=3;
     }

     xgp[0]=1/sqrt(3);
     ygp[0]=1/sqrt(3);

     xgp[1]=-1/sqrt(3);
     ygp[1]=-1/sqrt(3);

     wgtx[0]=1.0;
     wgtx[1]=1.0;
     wgty[0]=1.0;
     wgty[1]=1.0;
     wgtz[0]=1.0;
     wgtz[1]=1.0;

     *rstgeo=2;
   break;

   case 1: case 3:
     *nxgp=1;
     *nygp=2;
     *nzgp=2;

     if (gsurfno==1)
     {
       xgp[0]=1.;
       nodevec[0]=4;
       nodevec[1]=5;
       nodevec[2]=1;
       nodevec[3]=0;
     }
     else
     {
       xgp[0]=-1.;
       nodevec[0]=3;
       nodevec[1]=2;
       nodevec[2]=6;
       nodevec[3]=7;
     }

     ygp[0]=1/sqrt(3);
     zgp[0]=1/sqrt(3);

     ygp[1]=-1/sqrt(3);
     zgp[1]=-1/sqrt(3);

     wgtx[0]=1.0;
     wgtx[1]=1.0;
     wgty[0]=1.0;
     wgty[1]=1.0;
     wgtz[0]=1.0;
     wgtz[1]=1.0;

     *rstgeo=0;
   break;

   case 2: case 4:
     *nxgp=2;
     *nygp=1;
     *nzgp=2;

     if (gsurfno==2)
     {
       ygp[0]=1.;
       nodevec[0]=5;
       nodevec[1]=6;
       nodevec[2]=2;
       nodevec[3]=1;
     }
     else
     {
       ygp[0]=-1.;
       nodevec[0]=7;
       nodevec[1]=4;
       nodevec[2]=0;
       nodevec[3]=3;
     }

     xgp[0]=1/sqrt(3);
     zgp[0]=1/sqrt(3);

     xgp[1]=-1/sqrt(3);
     zgp[1]=-1/sqrt(3);

     wgtx[0]=1.0;
     wgtx[1]=1.0;
     wgty[0]=1.0;
     wgty[1]=1.0;
     wgtz[0]=1.0;
     wgtz[1]=1.0;

     *rstgeo=1;
   break;
   }
   /*-------------------------------------------------------------------*/
   #ifdef DEBUG
   dstrc_exit();
   #endif
   return;
}/* end of getnodesgp */

/*---------------------------------------------------------------------*/
/*! @} (documentation module close)*/
#endif

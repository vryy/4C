/*!----------------------------------------------------------------------
\file
\brief contains the routines 'wge_damvar' -> calculates actual damage 
       variable and derivative
       routine 'wge_epsequiv' -> calculates local equivalent strains and
       Tangente F = (partial eps_eq)/(partial epsilon)
       routine 'wge_stress' -> calculates elastic and damaged stresses
       routine 'wge_tangent' -> calculates elasto-damage-tangent
       routine 'wge_condense' condenses the 3-3 component of the tangents

<pre>
Maintainer: Andrea Hund
            hund@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hund/
            0771 - 685-6122
</pre> 
*----------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "wallge.h"
#include "wallge_prototypes.h"

/*! 
\addtogroup WALLGE 
*//*! @{ (documentation module open)*/


/*!----------------------------------------------------------------------
\brief routine 'wge_damvar' -> calculates actual damage 
       variable and derivative (partial D)/(partial eps_equiv)
       damtyp = 1 :         linear softening damage
       damtyp = 2 :         exponential damage law


\param   damtyp        INT      I: Damage law
\param   kappa         DOUBLE   I: actual nonl. eqvival. strain 
\param   kappa_0       DOUBLE   I: parameter of exp. damage law
\param  *eps_vnl       DOUBLE   I: nonlocal equiv. strain
\param   alpha         DOUBLE   I: parameter of exp. damage law
\param   beta          DOUBLE   I: parameter of exp. damage law
\param  *damage        DOUBLE   O: actual damage variable(scalar)
\param  *dam_deriv     DOUBLE   O: (partial D)/(partial eps_equi)

\return void                                               
\sa calling:   nothing; 
    called by: wge_mat_damage();
*----------------------------------------------------------------------*/
void wge_damvar(INT     damtyp,        /*  Damage law                   */
                DOUBLE  kappa,         /*  actual nonl. eqvival. strain */
                DOUBLE  kappa_0,       /* parameter of exp. damage law  */
                DOUBLE *eps_vnl,       /* nonlocal equiv. strain        */
                DOUBLE  alpha,         /* parameter of exp. damage law  */
                DOUBLE  beta,          /* parameter of exp. damage law  */
                DOUBLE *damage,        /* actual damage variable        */
                DOUBLE *dam_deriv)     /* (partial D)/(partial eps_equi)*/       
{
#ifdef D_WALLGE
/*----------------------------------------------------------------------*/

#ifdef DEBUG 
dstrc_enter("wge_damvar");
#endif
/*----------------------------------------------------------------------*/
switch(damtyp)
{
case 1:
  dserror("linear softening damage law not implemented");
break;
/*----------------------------------------------------------------------*/
case 2:
  /*-------------------------------------------------------- elastic ---*/
  if(kappa <= kappa_0 && *eps_vnl <= kappa_0)
  {
    *damage    = 0.0;
    *dam_deriv = 0.0;
     # if 0
     printf(" Elastisch: D=%le\n",*damage);
     # endif
  }
  /*-------------------------------------------------------- damaged ---*/
  else
  {
    *damage    = 1 - (kappa_0*(1-alpha+alpha*exp(beta*kappa_0-beta*kappa)))/kappa;
    *dam_deriv = (kappa_0*(1-alpha))/(kappa*kappa)
               + (alpha*kappa_0*exp(beta*kappa_0-beta*kappa))*(beta/kappa+1/(kappa*kappa));
     # if 0
     printf(" Schaedigung: D=%le\n",*damage);
     # endif
  }
break;
/*----------------------------------------------------------------------*/
default:
  dserror("damage law not specified in input file (must be 2)");
break;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

#endif /*D_WALLGE*/
return; 
} /* end of wge_damvar */

/*!----------------------------------------------------------------------
\brief routine 'wge_epsequiv' -> calculates local equivalent strains and
       Tangente F = (partial eps_eq)/(partial epsilon)
       equival = 1:      eps_v= energy release rate
       equival = 2:
       equival = 3:
       equival = 4:     de Vree

\param   equival      INT      I: definition of equiv. strains
\param   epsilon      DOUBLE   I: strain (tensor) 
\param   delta        DOUBLE   I: kroneker delta (tensor)
\param   k_f          DOUBLE   I: nonlocal equiv. strain
\param   nue          DOUBLE   I: parameter for de Vree
\param  *eps_vl       DOUBLE   O: poisson rate
\param   F_ed         DOUBLE   O: equ. strain derivative (tensor)

\return void                                               
\sa calling:   nothing; 
    called by: wge_mat_damage();
*----------------------------------------------------------------------*/
void wge_epsequiv(INT     equival,     /* definition of equiv. strains */
                  DOUBLE  epsilon[3][3], /* strain tensor              */
                  DOUBLE  delta[3][3],   /* kroneker delta             */
                  DOUBLE  k_f,           /* parameter for de Vree      */
                  DOUBLE  nue,           /* poisson rate               */
                  DOUBLE *eps_vl,        /* local equiv. strains       */
                  DOUBLE  F_ed[3][3])    /* equ. strain derivative     */         
{
#ifdef D_WALLGE
/*----------------------------------------------------------------------*/
INT    i,j;
DOUBLE I1,J2;
DOUBLE scalar,fac1,fac2,wurzel;
DOUBLE ev_I1, ev_J2;
#ifdef DEBUG 
dstrc_enter("wge_epsequiv");
#endif
/*--------------------------------------------- elasto-damge tangent ---*/
switch(equival)
{
case 4:
  I1 = epsilon[0][0]+epsilon[1][1]+epsilon[2][2];
  scalar= 0.0;
  for (i=0; i<3; i++)
  {
   for (j=0; j<3; j++)
   {
    scalar += epsilon[i][j]*epsilon[i][j];
   }
  } 
  J2     = (I1 * I1)/6.0 - scalar/2.0;
  fac1   = (k_f - 1.0)/(1.0 - 2.0* nue);
  fac2   = (12.0*k_f)/((1.0+nue)*(1.0+nue));
  wurzel = sqrt(fac1 * fac1 * I1 * I1 - fac2*J2 );
  if(wurzel == 0.0)  
  {
    wurzel = 1.0E-50;
    /*
    printf("Achtung Wurzel = 0 \n");
    */
  }
  *eps_vl = (fac1*I1)/(2.0*k_f) + wurzel/(2*k_f);
  
  ev_I1  = (fac1/(2*k_f))*(1 + (fac1*I1)/wurzel);
  ev_J2  = (-3.0)/((1.0+nue)*(1.0+nue)*wurzel);
  for (i=0; i<3; i++)
  {
   for (j=0; j<3; j++)
   {
    F_ed[i][j] = ev_I1*delta[i][j] + ev_J2*((I1*delta[i][j])/3.0-epsilon[i][j]);
   }
  } 
break;
default:
  dserror(" these euivalent strains are not implemented");
break;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

#endif /*D_WALLGE*/
return; 
} /* end of wge_epsequiv */

/*!----------------------------------------------------------------------
\brief routine 'wge_stress' -> calculates elastic and damaged stresses 

\param   youngs     DOUBLE   I: youngs modulus
\param   nue        DOUBLE   I: poisson ratio
\param   delta      DOUBLE   I: kroneker delta (tensor)
\param   epsilon    DOUBLE   I: strain (tensor)
\param   damage     DOUBLE   I: actual damage variable
\param   sigma_el   DOUBLE   O: elastic stress (tensor)
\param   sigma      DOUBLE   O: stress (tensor)

\return void                                               
\sa calling:   nothing; 
    called by: wge_mat_damage();
*----------------------------------------------------------------------*/
void wge_stress(DOUBLE   youngs,           /*  youngs modulus          */
                DOUBLE   nue,              /*  poisson ratio           */
                DOUBLE   delta[3][3],      /*  kroneker delta          */
                DOUBLE   epsilon[3][3],    /*  strain tensor           */
                DOUBLE   damage,           /*  actual damage variable  */
                DOUBLE   sigma_el[3][3],   /*  elastic stress tensor   */ 
                DOUBLE   sigma[3][3])      /*  stress tensor           */      
{
#ifdef D_WALLGE
/*----------------------------------------------------------------------*/
INT    i,j;
DOUBLE mu,lambda,trace;

#ifdef DEBUG 
dstrc_enter("wge_elastic");
#endif
/*--------------------------------------------- elasto-damge tangent ---*/

lambda = (youngs*nue)/((1.0+nue)*(1.0-2.0*nue));
mu     =  youngs/(2.0*(1.0+nue));

/*---------------------------------------- elastic and real stresses ---*/
trace = epsilon[0][0]+epsilon[1][1]+epsilon[2][2];
for (i=0; i<3; i++)
{
 for (j=0; j<3; j++)
 {
  sigma_el[i][j] = lambda*trace*delta[i][j]+2.0*mu*epsilon[i][j];
  sigma[i][j]    = (1 - damage)*sigma_el[i][j];
 }
} 
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

#endif /*D_WALLGE*/
return; 
} /* end of wge_stress */


/*!----------------------------------------------------------------------
\brief routine 'wge_tangent' -> calculates elasto-damage-tangent

\param   youngs     DOUBLE   I: youngs modulus
\param   nue        DOUBLE   I: poisson ratio
\param   delta      DOUBLE   I: kroneker delta (tensor)
\param   damage     DOUBLE   I: actual damage variable
\param   C_ed       DOUBLE   O: elasto-damage-tangent (tensor)

\return void                                               
\sa calling:   nothing; 
    called by: wge_mat_damage();
*----------------------------------------------------------------------*/
void wge_tangent(DOUBLE   youngs,           /*  youngs modulus         */
                 DOUBLE   nue,              /*  poisson ratio          */
                 DOUBLE   delta[3][3],      /*  kroneker delta         */
                 DOUBLE   damage,           /*  damaage variable       */
                 DOUBLE   C_ed[3][3][3][3]) /*  elasto-damage-tangent  */         
{
#ifdef D_WALLGE
/*----------------------------------------------------------------------*/
INT    i,j,k,l;
DOUBLE mu,lambda;

#ifdef DEBUG 
dstrc_enter("wge_tangent");
#endif
/*--------------------------------------------- elasto-damge tangent ---*/

lambda = (youngs*nue)/((1.0+nue)*(1.0-2.0*nue));
mu     =  youngs/(2.0*(1.0+nue));

for (i=0; i<3; i++)
{
 for (j=0; j<3; j++)
 {
  for (k=0; k<3; k++)
  {
   for (l=0; l<3; l++)
   {
    C_ed[i][j][k][l] = (1 - damage)*(lambda*delta[i][j]*delta[k][l]
                     + mu*(delta[i][k]*delta[j][l]+delta[i][l]*delta[j][k])); 
   }
  } 
 }
} 
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

#endif /*D_WALLGE*/
return; 
} /* end of wge_tangent */

/*!----------------------------------------------------------------------
\brief routine 'wge_condense' condenses the 3-3 component of the tangents

\param  **D      DOUBLE     I/O: 1.elasto-damage-tangent-matrix
\param   *E      DOUBLE     I/O: 2.elasto-damage-tangent-matrix
\param   *F      DOUBLE     I/O: 3.elasto-damage-tangent-matrix
\param    wtype  WALLGE_TYPE  I: plane-stress or plane strain
\param    nue    DOUBLE       I: poisson ratio
 
\return void                                               
\sa calling:   nothing; 
    called by: wge_mat_damage();

*----------------------------------------------------------------------*/
void wge_condense(DOUBLE    **D,     /* 1.elasto-damage-tangent-matrix */
                  DOUBLE     *E,     /* 2.elasto-damage-tangent-matrix */
                  DOUBLE     *F,     /* 3.elasto-damage-tangent-matrix */
                  WALLGE_TYPE wtype, /* plane-stress or plane strain   */
                  DOUBLE      nue)   /* poisson ratio                  */         
{
#ifdef D_WALLGE
/*----------------------------------------------------------------------*/
INT  i,j;
DOUBLE D_con[4][4];
DOUBLE F_con[3];
DOUBLE factor;
#ifdef DEBUG 
dstrc_enter("wge_condense");
#endif
/*----------------------------------------------------------------------*/
switch(wtype)
{
/*----------------------------------------------------------------------*/
case pl_strain:

  for (i=0; i<4; i++)
  {
   D[i][3] = 0.0;
   D[3][i] = 0.0;
  }
  E[3] = 0.0;
  F[3] = 0.0;
  
break;
/*----------------------------------------------------------------------*/
case pl_stress:

  if(D[3][3] != 0.0)
  {
    for(i=0; i<3; i++)
    {
     for (j=0; j<3; j++)
     {
       D_con[i][j] = D[i][j] - (D[i][3]*D[3][j])/D[3][3];
     }  
    }  
    for(i=0; i<3; i++)
    {
     for (j=0; j<3; j++)
     {
       D[i][j] = D_con[i][j];
     }  
    }
    for(i=0; i<4; i++)
    {
      D[i][3] = 0.0;
      D[3][i] = 0.0;
    }
  }  
  E[3] = 0.0;
  factor   = nue/(1.0 - nue); 
  F_con[0] = F[0] - factor * F[3];
  F_con[1] = F[1] - factor * F[3];
  F_con[2] = F[2];
  for(i=0; i<3; i++)
  {
    F[i] = F_con[i];
  }
  F[3] = 0.0;
  
break;
default:
  dserror("wall-type(planestrain/planestress) not specified");
break;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

#endif /*D_WALLGE*/
return; 
} /* end of wge_condense */




/*----------------------------------------------------------------------*/
/*! @} (documentation module close)*/

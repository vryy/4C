/*!-----------------------------------------------------------------------
\file
\brief contains the routine 'w1_mat_dam_mp' which is a isotropic damage
 model (mazars/pijadier-cabot -> PDH-Peerlings) in vector-matrix-format
*-----------------------------------------------------------------------*/
#ifdef D_WALL1
#include "../headers/standardtypes.h"
#include "wall1.h"
#include "wall1_prototypes.h"

/*!
\addtogroup WALL1
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief  contains a isotropic damage model

<pre>                                                              mn 05/03
This routine computes the constitutive  matrix of the interface element

</pre>
\param *mat          STVENANT   (i)   blabal

\warning There is nothing special to this routine
\return void
\sa calling:   ---;
    called by: w1_call_mat();

*----------------------------------------------------------------------*/
void w1_mat_dam_mp(DOUBLE     youngs,
                   DOUBLE     nue,
                   DOUBLE     kappa_0,
                   DOUBLE     alph,
                   DOUBLE     beta,
                   ELEMENT   *ele,
                   WALL_TYPE  wtype,
                   DOUBLE   **bop,
                   DOUBLE    *gop,
                   DOUBLE    *alpha,
                   INT        ip,
                   DOUBLE    *stress,
                   DOUBLE   **D,
                   INT        istore,
                   INT        newval)
{
INT i,j,k;
INT yip;
/*------------------------------------------------------------------------*/
INT flag;
/*------------------------------------------------------------------------*/
DOUBLE disd[5];    /* displacement derivatives*/
DOUBLE strain[4];  /* strains*/
DOUBLE meps[3];    /* main strains*/
DOUBLE eps_eq;     /* equivalent strain*/
DOUBLE kappa;      /* maximum value of equivalent strain in the past*/
DOUBLE kappa_n;    /* maximum value of equivalent strain in the past (load step n)*/
DOUBLE d;          /* damage variable*/
DOUBLE e1,e2,e3,wurzel,bruch,help1,help2,help3;   /* working values*/
DOUBLE ev_e1,ev_e2,ev_e3;/* (partial equiv. strain)/(partial main strains)*/
DOUBLE D_ev;         /* (partial D)/(partial equiv. strain)*/
DOUBLE e_eps[3][3];/* (partial main strains)/(partial strain components)*/
/*------------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("w1_mat_dam_mp");
#endif
/*------------------------------------------------------------------------*/

/*------------------------------------------------------------------------*/
/*        Routine is called by if_cal_stress -> only get stresses         */
/*------------------------------------------------------------------------*/
if(newval==1)
{
  stress[0] = ele->e.w1->elewa[0].ipwa[ip].sig[0];
  stress[1] = ele->e.w1->elewa[0].ipwa[ip].sig[1];
  stress[2] = ele->e.w1->elewa[0].ipwa[ip].sig[2];
goto end;
}
/*------------------------------------------------------------------------*/
# if 0


yip    = ele->e.w1->elewa[0].ipwa[ip].yip;
/*------------------------------------------------------------------------*/
/*    Global praedictor: yip > 0 (excaption: first load step:yip=-1)      */
/*------------------------------------------------------------------------*/


if (yip>0)
{
  kappa = ele->e.w1->elewa[0].ipwa[ip].kappa;
  if(kappa>kappa_0)
  {
   d = 1 -(kappa_0*(1 - alph + alph * exp(beta*kappa_0 - beta*kappa)))/kappa;
  }
  else
  {
   d = 0.0;
  }
  switch(wtype)
  {
  case plane_stress:
/*-------------------------------  schabernack fuer lokale Schwaechung ---*/
flag=0;
for (k=0; k<4; k++)
{
if(ele->node[k]->x[0]==3.0 )
{
flag=1;
}
}
for (i=0; i<4; i++)
{
if(ele->node[i]->x[0]==2.0  && flag==1)
{
youngs = youngs*0.95;
}
}
/*------------------------------------------------ Ende des Bloedsinns ---*/

    e1 = youngs/(1. - nue*nue);
    e2 = nue * e1;
    e3 = e1 * (1. - nue)/2.;
    D[0][0] = (1.0 - d) *e1;
    D[0][1] = (1.0 - d) *e2;
    D[0][2] = 0.;
    D[1][0] = (1.0 - d) *e2;
    D[1][1] = (1.0 - d) *e1;
    D[1][2] = 0.;
    D[2][0] = 0.;
    D[2][1] = 0.;
    D[2][2] = (1.0 - d) *e3;
  break;
  default:
    dserror("only plane-stress damage implemented");
  break;
  }
  stress[0] = ele->e.w1->elewa[0].ipwa[ip].sig[0];
  stress[1] = ele->e.w1->elewa[0].ipwa[ip].sig[1];
  stress[2] = ele->e.w1->elewa[0].ipwa[ip].sig[2];
  ele->e.w1->elewa[0].ipwa[ip].yip = -yip;
}/*end of: if (yip>0) */


# endif

/*------------------------------------------------------------------------*/
/*    Global correktor or updateing or praedictor of first loadstep       */
/*------------------------------------------------------------------------*/

# if 0

else if (yip < 0)
{
# endif
/*-------------------------------  schabernack fuer lokale Schwaechung ---*/
# if 0
flag=0;
for (k=0; k<4; k++)
{
if(ele->node[k]->x[0]==3.0 )
{
flag=1;
}
}
for (i=0; i<4; i++)
{
if(ele->node[i]->x[0]==2.0  && flag==1)
{
youngs = youngs*0.95;
}
}
# endif
/*------------------------------------------------ Ende des Bloedsinns ---*/
  /*----------------------------------------------------------------------*/
  /*-------------- maximale bisher aufgetretene equivalente Verzerrung ---*/
  kappa_n = ele->e.w1->elewa[0].ipwa[ip].kappa;
  /*--------------------------------- compute displacement derivatives ---*/
  w1_disd (ele,bop,gop,alpha,wtype,disd);
  /*------------------------------------- get actual strains -> strain ---*/
  /*------------------------------------------- strain[2]=eps12+eps21! ---*/
  w1_eps (disd,wtype,strain);
  switch(wtype)
  {
  case plane_stress:
  strain[3] = (- nue *(strain[0]+strain[1]))/(1 - nue);
  /*------------------------------------------- calculate main strains ---*/
  wurzel  = sqrt(((strain[0]-strain[1])*(strain[0]-strain[1]))/4.0 + (strain[2]*strain[2])/4.0);
  meps[0] = (strain[0]+strain[1])/2.0 + wurzel;
  meps[1] = (strain[0]+strain[1])/2.0 - wurzel;
  meps[2] = strain[3];
  /*------------------------------- calculate actual equivalent strain ---*/
  eps_eq  = (sqrt((meps[0]+FABS(meps[0]))*(meps[0]+FABS(meps[0]))+
                  (meps[1]+FABS(meps[1]))*(meps[1]+FABS(meps[1]))+
                  (meps[2]+FABS(meps[2]))*(meps[2]+FABS(meps[2]))  ))/2.0;
  /*----------------------------------------------------------------------*/
  /*--------------------------------------------- damage is increasing ---*/
  if(eps_eq>kappa_n)
  {
   kappa = eps_eq;
  }
  /*---------------------------------------- else: kappa bleibt gleich ---*/
  else
  {
   kappa = kappa_n;
  }
  /*---------------  exponential damage law (elastic for kappa<kappa_0)---*/
  if(kappa>kappa_0)
  {
   d = 1 -(kappa_0*(1 - alph + alph * exp(beta*kappa_0 - beta*kappa)))/kappa;
  }
  else
  {
   d = 0.0;
  }
  /*---------------------------------------------- constitutive matrix ---*/
    e1 = youngs/(1. - nue*nue);
    e2 = nue * e1;
    e3 = e1 * (1. - nue)/2.;
    D[0][0] = (1.0 - d) *e1;
    D[0][1] = (1.0 - d) *e2;
    D[0][2] = 0.;
    D[1][0] = (1.0 - d) *e2;
    D[1][1] = (1.0 - d) *e1;
    D[1][2] = 0.;
    D[2][0] = 0.;
    D[2][1] = 0.;
    D[2][2] = (1.0 - d) *e3;
  /*------------------------------------------------------- stresses ---*/
    stress[0]=D[0][0]*strain[0] +D[0][1]*strain[1] +D[0][2]*strain[2];
    stress[1]=D[1][0]*strain[0] +D[1][1]*strain[1] +D[1][2]*strain[2];
    stress[2]=D[2][0]*strain[0] +D[2][1]*strain[1] +D[2][2]*strain[2];
    stress[3]=0.0;
  /*----- if damage is increasing: constitutive matrix to be adjusted ---*/
  if(eps_eq>kappa_n  && eps_eq>kappa_0)
  {
    D_ev  = (kappa_0 * (1 - alph))/(kappa * kappa) +
            (alph * kappa_0 * exp(beta*kappa_0 - beta*kappa))*(beta/kappa + 1.0/(kappa*kappa));
    ev_e1 = ((meps[0]+FABS(meps[0]))*(1. + FSIGN(meps[0])))/(4.* eps_eq);
    ev_e2 = ((meps[1]+FABS(meps[1]))*(1. + FSIGN(meps[1])))/(4.* eps_eq);
    ev_e3 = ((meps[2]+FABS(meps[2]))*(1. + FSIGN(meps[2])))/(4.* eps_eq);
    bruch = (strain[0] - strain[1])/(2.0*wurzel);
    e_eps[0][0] = (1.0 + bruch)/2.0;
    e_eps[0][1] = (1.0 - bruch)/2.0;
    e_eps[0][2] = strain[2]/(2.0 * wurzel);
    e_eps[1][0] = (1.0 - bruch)/2.0;
    e_eps[1][1] = (1.0 + bruch)/2.0;
    e_eps[1][2] = ( -strain[2])/(2.0 * wurzel);
    e_eps[2][0] = - (nue)/(1-nue);
    e_eps[2][1] = - (nue)/(1-nue);
    e_eps[2][2] = 0.0;
    help1 = D_ev*(ev_e1*e_eps[0][0]+ev_e2*e_eps[1][0]+ev_e3*e_eps[2][0]);
    help2 = D_ev*(ev_e1*e_eps[0][1]+ev_e2*e_eps[1][1]+ev_e3*e_eps[2][1]);
    help3 = D_ev*(ev_e1*e_eps[0][2]+ev_e2*e_eps[1][2]);
  /*----------------------------------------------------------------------*/
    D[0][0] -= help1 * (e1*strain[0] + e2*strain[1]);
    D[0][1] -= help2 * (e1*strain[0] + e2*strain[1]);
    D[0][2] -= help3 * (e1*strain[0] + e2*strain[1]);
    D[1][0] -= help1 * (e2*strain[0] + e1*strain[1]);
    D[1][1] -= help2 * (e2*strain[0] + e1*strain[1]);
    D[1][2] -= help3 * (e2*strain[0] + e1*strain[1]);
    D[2][0] -= help1 * (e3*strain[2]);
    D[2][1] -= help2 * (e3*strain[2]);
    D[2][2] -= help3 * (e3*strain[2]);
  }
  break;
  default:
    dserror("only plane-stress damage implemented");
  break;
  }
  yip=1;
  /*----------------------------------------------------------------------*/
  /*              update the converged results of a loadstep              */
  /*----------------------------------------------------------------------*/
  if(istore==1)
  {
    ele->e.w1->elewa[0].ipwa[ip].yip    = yip;
    ele->e.w1->elewa[0].ipwa[ip].kappa  = kappa;
    ele->e.w1->elewa[0].ipwa[ip].damage = d;
    ele->e.w1->elewa[0].ipwa[ip].aequistrain = eps_eq;
    ele->e.w1->elewa[0].ipwa[ip].sig[0] = stress[0];
    ele->e.w1->elewa[0].ipwa[ip].sig[1] = stress[1];
    ele->e.w1->elewa[0].ipwa[ip].sig[2] = stress[2];
  }/*end of: if (istore==1) */

# if 0
}/*end of: if (yip<0) */
# endif

/*----------------------------------------------------------------------*/
end:
/*----------------------------------------------------------------------*/



#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of w1_mat_dam_mp */

/*----------------------------------------------------------------------*/
#endif /*D_WALL1*/
/*! @} (documentation module close)*/

/*!-----------------------------------------------------------------------
\file
\brief contains the routine 'if_mat' which is the material law
 for the interface element

<pre>
Maintainer: Andrea Hund
            hund@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hund/
            0711 - 685-6122
</pre>
*-----------------------------------------------------------------------*/
#ifndef CCADISCRET
#ifdef D_INTERF
#include "../headers/standardtypes.h"
#include "interf.h"
#include "interf_prototypes.h"
#include "../struct2_ml/s2ml.h"
#include "../wall1/wall1.h"

/*!
\addtogroup INTERF
*/
/*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief  contains a material law for the interface element

<pre>                                                              ah 05/03
This routine computes the constitutive  matrix of the interface element

</pre>
\param   *ele          ELEMENT   (I)   actual element (macro)
\param   *mat          MATERIAL  (I)   actual material
\param  **bop          DOUBLE    (I)   B-operator
\param  **D            DOUBLE    (O)   material tangent
\param   *T            DOUBLE    (O)   stresses
\param    ip           INT       (I)   ID of actual GP
\param    istore       INT       (I)   is it update?
\param    newval       INT       (I)   is it stress calculation
\param    smallscale   INT       (I)   is it call from sm-ele in multiscale
\param   *actsmele     ELEMENT   (I)   if multiscale: actual sm-element
\param   *jumpu_tot    DOUBLE    (I)   if multiscale: total displ. jump
\param   *DELTAjumpu_tot DOUBLE  (I)   if multi: tot increm. displ.jump

\warning There is nothing special to this routine
\return void
\sa calling:   ---;
    called by: if_static_ke();

*----------------------------------------------------------------------*/
void if_mat(ELEMENT   *ele,        /* actual element (macro)               */
            MATERIAL  *mat,        /* actual material                      */
            DOUBLE   **bop,        /* B-operator                           */
            DOUBLE   **D,          /* material tangent                     */
            DOUBLE    *T,          /* stresses                             */
            INT        ip,         /* ID of actual GP                      */
            INT        istore,     /* is it update?                        */
            INT        newval,     /* is it stress calculation             */
            INT        smallscale, /* is it call from sm-ele in multiscale */
            ELEMENT   *actsmele,   /* if multiscale: actual sm-element     */
            DOUBLE    *jumpu_tot,  /* if multiscale: total displ. jump     */
            DOUBLE    *DELTAjumpu_tot) /* if multi: tot increm. displ.jump */
{
INT i,j;
INT ID=0;        /* ID of actual submesh-element if it's multiscale */
INT yip;
DOUBLE E,K,G,thick,delta_n,delta_t,mu,Ynmax,Ytmax;
DOUBLE dn,Yn;
DOUBLE dn_neu=0.0;
DOUBLE dt_neu=0.0;
DOUBLE ut_pl_neu=0.0;
DOUBLE dt,dt_tr,Yt,Yt_tr,Tt_tr,ut_pl,f_tr;
DOUBLE deltalambda;
DOUBLE disjump[2];
DOUBLE Deltadisjump[2];

/*------------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("if_mat");
#endif
/*---------------------------------------------- material parameters -----*/

E       = mat->m.ifmat->emod;   /*- Normalzugsteifigkeit-*/
K       = mat->m.ifmat->kmod;   /*- Normaldrucksteifigkeit-*/
G       = mat->m.ifmat->gmod;   /*- Schubsteifigkeit,die bei Dekohaesion abnimmt-*/
thick   = mat->m.ifmat->dick;   /*- Pseudodicke des Interfaces -*/
/*Q       = mat->m.ifmat->qmod;   - konstanter Anteil der Schubsteifigkeit-*/
delta_n = mat->m.ifmat->deltan; /*- max Normalzugverschiebung -> vollst Dekohaesion erreicht-*/
delta_t = mat->m.ifmat->deltat; /*- max Schubverschiebung -> vollst Dekohaesion erreicht-*/
mu      = mat->m.ifmat->mu;     /*- Reibungskoeffizient-*/
/*------------------------------------------------------------------------*/

if(smallscale == 0)
{
  yip    = ele->e.interf->elewa[0].ipwa[ip].yip;
}
#ifdef D_MLSTRUCT
else
{
  ID  = actsmele->Id;
  yip = ele->e.w1->sm_eledata[ID].sm_GPdata[ip].yip;
}
#endif /* D_MLSTRUCT */
/*------------------------------------------------------------------------*/
/*    Global corrector or updateing or praedictor of first loadstep       */
/*------------------------------------------------------------------------*/
Ynmax  = (E * delta_n * delta_n) / (TWO * thick);
Ytmax  = (G * delta_t * delta_t) / (TWO * thick);

if(smallscale == 0)
{
  /* calc.tang. and normal displ jumps [un],[ut],[DELTAun],[DELTAut] ---*/
  if_jumpu(ele,bop,disjump,Deltadisjump);
  /*-------------------------------- Werte des letzten Lastschrittes ---*/
  dn = ele->e.interf->elewa[0].ipwa[ip].dn;
  dt = ele->e.interf->elewa[0].ipwa[ip].dt;
  ut_pl  = ele->e.interf->elewa[0].ipwa[ip].jump_ut_pl;
}
#ifdef D_MLSTRUCT
else  /*----------------- smallscale ==1, disjump already calculated ---*/
{
  disjump[0] = jumpu_tot[0];
  disjump[1] = jumpu_tot[1];
  Deltadisjump[0] = DELTAjumpu_tot[0];
  Deltadisjump[1] = DELTAjumpu_tot[1];
  /*-------------------------------- Werte des letzten Lastschrittes ---*/
  dn =  ele->e.w1->sm_eledata[ID].sm_GPdata[ip].dn;
  dt =  ele->e.w1->sm_eledata[ID].sm_GPdata[ip].dt;
  ut_pl = ele->e.w1->sm_eledata[ID].sm_GPdata[ip].utpl;
}
#endif /* D_MLSTRUCT */
/*---------------------- is the interface in tension or compression? ---*/
/*----------------------------------------------------------------------*/
/*                            tension case                              */
/*----------------------------------------------------------------------*/
/*  printf("Elem-ID: %d , GP-ID: %d \n",ele->Id,ip); */
if (disjump[1]>0.0)
{

  /*-------------------------------------------- 1. Normal-direction ---*/
  /*--------------------- Normaldekohaesionsfortschritt (Schaedigung)---*/
  if (disjump[1]*Deltadisjump[1]>=0.0)
  {
      /* printf("Normalzug - Normaldekohaesionsfortschritt\n"); */
     Yn     = (E * disjump[1] * disjump[1])/ (TWO * thick);
     dn_neu = 1.0 - (1.0 - sqrt(Yn/Ynmax)) * (1.0 - sqrt(Yn/Ynmax));
     if (dn_neu < dn)
     {
       dn_neu = dn;
     }
     if (Yn > Ynmax)
     {
       dn_neu = 1.0;
     }
     T[1]      = (E * disjump[1] * (1.0 - dn_neu)) / thick;
     D[1][0] = 0.0;
     if (dn_neu <1.0)
     {
       D[1][1] = (E * (1.0 - dn_neu)) / thick
               - (E * E * disjump[1] *disjump[1]*(1.0/sqrt(Yn/Ynmax) -1.0))
                  / (thick * thick * Ynmax);
     }
     else
     {
       D[1][1] = (E * 1.0E-8)/ thick;
# if 0
       D[1][1] = 0.0;
# endif
       }
  }
  /*----------------------------- Normaldekohaesionsstop (elastisch) ---*/
  else
  {
    /* printf("Normalzug - Normaldekohaesionsstop un: %f DELTAun: %f\n",disjump[1],Deltadisjump[1]); */
    dn_neu  = dn;
    T[1]      = (E * disjump[1] * (1.0 - dn_neu)) / thick;
    D[1][0] = 0.0;
    D[1][1] = ( E * (1.0 - dn_neu)) / thick;
  }
  /*---------------------------------------- 2. Tangential-direction ---*/
  ut_pl_neu = ut_pl;
  /*---------------- Tangentialdekohaesionsfortschritt (Schaedigung) ---*/
  if (disjump[0]*Deltadisjump[0]>=0.0)
  {
     /* printf("Normalzug - Schubdekohaesionsfortschritt\n"); */
     Yt     = (G*(disjump[0]- ut_pl)*(disjump[0]- ut_pl))/(TWO * thick);
     dt_neu = 1.0 - (1.0 - sqrt(Yt/Ytmax)) * (1.0 - sqrt(Yt/Ytmax));
     if (dt_neu < dt)
     {
       dt_neu = dt;
     }
     if (Yt > Ytmax)
     {
       dt_neu = 1.0;
     }
     T[0] = ( G *(1.0 - dt_neu) * (disjump[0]- ut_pl))/thick;
     D[0][1] = 0.0;
     if (dt_neu <1.0)
     {

       D[0][0] = (G * (1.0 - dt_neu))/thick
               - (G*G*(disjump[0]-ut_pl)*(disjump[0]-ut_pl)*(1.0/sqrt(Yt/Ytmax) -1.0))
                 /(thick * thick * Ytmax);
     }
     else
     {
       D[0][0] = (G * 1.0E-8)/ thick;
# if 0
       D[0][0] = 0;
# endif
     }
  }
  /*------------------------- Tangentialdekohaesionsstop (elastisch) ---*/
  else
  {
     /*----------- Achtung: allererster praediktor hat wohl dt=0  ---*/
     /* printf("Normalzug - Schubdekohaesionsstop ut: %f DELTAut: %f\n",disjump[0],Deltadisjump[0]); */
     dt_neu  = dt;
     T[0]    = ( G *(1.0 - dt_neu) * (disjump[0]- ut_pl))/thick;
     D[0][1] = 0.0;
     D[0][0] = (G * (1.0 - dt_neu))/thick;
  }
  yip=1;
}/* end of: if (disjump[1]>1) */
/*----------------------------------------------------------------------*/
/*                           compression case                           */
/*----------------------------------------------------------------------*/
else if (disjump[1]<0.0)
{
   /*---------------------------------- 1. Normal-direction (elastic)---*/
   T[1]      = (K * disjump[1]) / thick;
   D[1][1]   = K / thick;
   dn_neu    = dn;
   /*--------------------------------------- 2. Tangential-direction ---*/
   /*--------------------------------- is it decohesion or friction? ---*/
   Yt_tr   = (G * disjump[0] * disjump[0])/ (TWO * thick);
   dt_tr   = 1.0 - (1.0 - sqrt(Yt_tr/Ytmax)) * (1.0 - sqrt(Yt_tr/Ytmax));
   if (dt_tr < dt)
   {
      dt_tr = dt;
   }
   if (Yt_tr > Ytmax)
   {
     dt_tr = 1.0;
   }
   /*----------------------------------------------- decohesion phase---*/
   if (dt_tr <1.0)
   {
     /*------------------------------------- Dekohaesionsfortschritt ---*/
     if(disjump[0]*Deltadisjump[0]>=0.0)
     {
       /* printf("Normaldruck - Schubdekohaesionsfortschritt\n"); */
       dt_neu  = dt_tr;
       T[0]    = ((G * (1.0 - dt_neu) - (mu * K * disjump[1])/(delta_t))* disjump[0])/ thick;
       D[1][0] = 0.0;
       D[0][1] = -(mu * K * disjump[0])/(delta_t * thick);
       D[0][0] = (G * (1.0 - dt_neu) - (mu * K * disjump[1])/(delta_t)) / thick
               - (G * G * disjump[0] * disjump[0]*(1.0/sqrt(Yt_tr/Ytmax) -1))
                 /(thick * thick * Ytmax);
     }
     /*------------------------------------------- Dekohaesionsstop ---*/
     else
     {
       /* printf("Normaldruck - Schubdekohaesionsstop ut: %f DELTAut: %f\n",disjump[0],Deltadisjump[0]); */
       dt_neu  = dt;
       T[0]    = ((G * (1.0 - dt_neu) - (mu * K * disjump[1])/(delta_t))* disjump[0]) / thick;
       D[1][0] = 0.0;
       D[0][1] = -(mu * K * disjump[0])/(delta_t * thick);
       D[0][0] = (G * (1.0 - dt_neu) - (mu * K * disjump[1])/(delta_t)) / thick;
     }
     ut_pl_neu = ut_pl;
     yip = 1;
   }
   /*------------------------------------------------- friction phase---*/
   else
   {
     dt_neu = dt_tr;
     Tt_tr  = - ( mu * K * disjump[1] * (disjump[0] - ut_pl)) / (thick * delta_t);
     f_tr   =  FABS (Tt_tr) + mu * T[1];
   /*------------------------------------------------------- plastic ---*/
     if (f_tr >=0.0)
     {
        /* printf("Normaldruck - Schubreibung\n"); */
       deltalambda = FABS(disjump[0] - ut_pl) - delta_t;
       ut_pl_neu   = ut_pl + deltalambda * FSIGN(Tt_tr);
       T[0]        = Tt_tr +
                    (mu* K* disjump[1]* deltalambda* FSIGN(Tt_tr))/(thick*delta_t);
       D[1][0] = 0.0;
# if 0
       D[0][0] = 0.0;
# endif
   /*----------------- Achtung diese Steigung ist eigentlich Null!!! ---*/
       D[0][0] = (mu * K * 1.0E-4)/ thick;
   /*-------------------------------------------------------------------*/
       D[0][1] = - (mu * K)/(thick * FSIGN(T[0]));
       yip = 2;
     }
   /*------------------------------------------------------- elastic ---*/
     else
     {
       /* printf("Normaldruck - Schubelastisch\n"); */
       ut_pl_neu = ut_pl;
       T[0]        = Tt_tr;
       D[1][0]   = 0.0;
       D[0][1]   = - (mu * K * (disjump[0] - ut_pl))/(thick * delta_t);
       D[0][0]   = - (mu * K * disjump[1]) / (thick * delta_t);
       yip = 1;
     }
   }
}/* end of: if (disjump[1]<1) */
/*----------------------------------------------------------------------*/
/*         no normal displacement or first loadstep's predictor         */
/*----------------------------------------------------------------------*/
else if (disjump[1]==0.0)
{
  /*-------------------------------------------- 1. Normal-direction ---*/
  T[1]    = 0.0;
  D[1][1] = ( K + E )/(TWO*thick); /*-- dunno if tension or compression ---*/
  D[0][1] = 0.0;
  D[1][0] = 0.0;
  dn_neu  = dn;
  /*---------------------------------------- 2. Tangential-direction ---*/
  ut_pl_neu = ut_pl;
  /*---------------- Tangentialdekohaesionsfortschritt (Schaedigung) ---*/
  if (disjump[0]*Deltadisjump[0]>0.0)
  {
     Yt     = (G*(disjump[0]- ut_pl)*(disjump[0]- ut_pl))/(TWO * thick);
     dt_neu = 1.0 - (1.0 - sqrt(Yt/Ytmax)) * (1.0 - sqrt(Yt/Ytmax));
     if (dt_neu < dt)
     {
       dt_neu = dt;
     }
     if (Yt > Ytmax)
     {
       dt_neu = 1.0;
     }
     T[0] = ( G *(1.0 - dt_neu) * (disjump[0]- ut_pl))/thick;
     D[0][1] = 0.0;
     if (dt_neu <1.0)
     {

       D[0][0] = (G * (1.0 - dt_neu))/thick
               - (G*G*(disjump[0]-ut_pl)*(disjump[0]-ut_pl)*(1.0/sqrt(Yt/Ytmax) -1.0))
                 /(thick * thick * Ytmax);
     }
     else
     {
       D[0][0] = (G * 1.0E-8)/ thick;
# if 0
       D[0][0] = 0.0;
# endif
     }
  }
  /*------------------------- Tangentialdekohaesionsstop (elastisch) ---*/
  else
  {
     /*----------- Achtung: allererster praediktor hat wohl dt=0  ---*/
     dt_neu  = dt;
     T[0]      = ( G *(1.0 - dt_neu) * (disjump[0]- ut_pl))/thick;
     D[0][1] = 0.0;
     D[0][0] = (G * (1.0 - dt_neu))/thick;
  }
  yip=1;
}/* end of: if (disjump[1]==0) */
/*----------------------------------------------------------------------*/
/*              update the converged results of a loadstep              */
/*----------------------------------------------------------------------*/
if(istore==1)
{
  if(smallscale ==0)
  {
    ele->e.interf->elewa[0].ipwa[ip].Tt = T[0];
    ele->e.interf->elewa[0].ipwa[ip].Tn = T[1];
    ele->e.interf->elewa[0].ipwa[ip].dt = dt_neu;
    ele->e.interf->elewa[0].ipwa[ip].dn = dn_neu;
    ele->e.interf->elewa[0].ipwa[ip].jump_ut_pl = ut_pl_neu;
    for (i=0; i<2; i++)
    {
      for (j=0; j<2; j++)
      {
        ele->e.interf->elewa[0].ipwa[ip].Q[i][j] = D[i][j];
      }
    }
    ele->e.interf->elewa[0].ipwa[ip].yip = yip;
  }
#ifdef D_MLSTRUCT
  else  /*smallscale ==1*/
  {
    ele->e.w1->sm_eledata[ID].sm_GPdata[ip].dn = dn_neu;
    ele->e.w1->sm_eledata[ID].sm_GPdata[ip].dt = dt_neu;
    ele->e.w1->sm_eledata[ID].sm_GPdata[ip].utpl = ut_pl_neu;
    ele->e.w1->sm_eledata[ID].sm_GPdata[ip].T[0] = T[0];
    ele->e.w1->sm_eledata[ID].sm_GPdata[ip].T[1] = T[1];
    for (i=0; i<2; i++)
    {
      for (j=0; j<2; j++)
      {
        ele->e.w1->sm_eledata[ID].sm_GPdata[ip].D[i][j] = D[i][j];
      }
    }
    ele->e.w1->sm_eledata[ID].sm_GPdata[ip].yip = yip;
  }
#endif /* D_MLSTRUCT */
}/*end of: if (istore==1) */


/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of if_mat */




/*!----------------------------------------------------------------------
\brief  contains a material law for the interface element

<pre>                                                              ah 05/03
This routine computes the constitutive  matrix of the interface element

</pre>
\param   *ele          ELEMENT   (I)   actual element (macro)
\param   *mat          MATERIAL  (I)   actual material
\param  **bop          DOUBLE    (I)   B-operator
\param  **D            DOUBLE    (O)   material tangent
\param   *T            DOUBLE    (O)   stresses
\param    ip           INT       (I)   ID of actual GP
\param    istore       INT       (I)   is it update?
\param    newval       INT       (I)   is it stress calculation
\param    smallscale   INT       (I)   is it call from sm-ele in multiscale
\param   *actsmele     ELEMENT   (I)   if multiscale: actual sm-element
\param   *jumpu_tot    DOUBLE    (I)   if multiscale: total displ. jump
\param   *DELTAjumpu_tot DOUBLE  (I)   if multi: tot increm. displ.jump

\warning There is nothing special to this routine
\return void
\sa calling:   ---;
    called by: if_static_ke();

*----------------------------------------------------------------------*/
void if_mat_thermodyn(ELEMENT   *ele,        /* actual element (macro)               */
            MATERIAL  *mat,        /* actual material                      */
            DOUBLE   **bop,        /* B-operator                           */
            DOUBLE   **D,          /* material tangent                     */
            DOUBLE    *T,          /* stresses [0]: tang. [1]:normal       */
            INT        ip,         /* ID of actual GP                      */
            INT        istore,     /* is it update?                        */
            INT        newval,     /* is it stress calculation             */
            INT        smallscale, /* is it call from sm-ele in multiscale */
            ELEMENT   *actsmele,   /* if multiscale: actual sm-element     */
            DOUBLE    *jumpu_tot,  /* if multiscale: total displ. jump     */
            DOUBLE    *DELTAjumpu_tot) /* if multi: tot increm. displ.jump */
{
INT equival,damtyp; /* kind of equivalent strain measure, kind of damage evol.law */
DOUBLE fac_t,fac_n; /* fac=0 -> damaged but unloading, fac = 1 damaged and loading*/
DOUBLE E,nu,dick,kappa0_n,alpha_n,beta_n,kappa0_t,alpha_t,beta_t;/* mat.param */
DOUBLE kappa_n_last,kappa_t_last;  /* history variables */
DOUBLE epsv_n=0.0;  /* actual equivalent strain, normal and tang. */
DOUBLE epsv_t=0.0;  /* actual equivalent strain, normal and tang. */
DOUBLE epsv_t_deriv=0.0;  /* (partial eps_v)/(partial [u]) */
DOUBLE epsv_n_deriv=0.0;  /* (partial eps_v)/(partial [u]) */
DOUBLE kappa_n,kappa_t;  /* actual value of kappa_n and kappa_t */
DOUBLE damage_n,damage_t;  /* damage variables */
DOUBLE damage_n_deriv,damage_t_deriv;  /* (partial d)/(partial kappa) */
DOUBLE Q_n, Q_t;       /* Normal and tang. Stiffness*/
DOUBLE T_elast[2];       /*elastic stress vector [0]: tang. [1]:normal*/
DOUBLE disjump[2];       /*displacement jump [0]: tang. [1]:normal*/
DOUBLE Deltadisjump[2];  /*incremental displ. jump[0]: tang. [1]:normal*/
DOUBLE DeltaD[2][2];  /* change of mat.tangent w.r.t.secant due to loading*/
DOUBLE tol = 1.0E-10;  /* Tolerance for Fallunterschiedung loading-unloading*/

INT ID=0;        /* ID of actual submesh-element if it's multiscale */

/*------------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("if_mat_thermodyn");
#endif
/*------------------------------------------------ material parameters ---*/
E        = mat->m.interf_therm->emod;    /*- E-Modul-*/
nu       = mat->m.interf_therm->nu;      /*- Querdehnzahl-*/
dick     = mat->m.interf_therm->dick;    /*- Pseudodicke des Interfaces-*/
kappa0_n = mat->m.interf_therm->kappa0_n;/*- Schaedigungsschwellwert normal -*/
alpha_n  = mat->m.interf_therm->alpha_n; /*- Schaedigungsgrenzwert -*/
beta_n   = mat->m.interf_therm->beta_n;  /*- Schaedigungsevolutionsgeschwindigkeit-*/
kappa0_t = mat->m.interf_therm->kappa0_t; /*- dasselbe tangential-*/
alpha_t  = mat->m.interf_therm->alpha_t;  /*- dasselbe tangential-*/
beta_t   = mat->m.interf_therm->beta_t;   /*- dasselbe tangential-*/

equival  = mat->m.interf_therm->equival;    /*-definition of equivalent strains -*/
damtyp   = mat->m.interf_therm->damtyp;     /*-definition of damage evolution law -*/

Q_n = E*(1.0-nu)/(dick*(1.0+nu)*(1.0-2*nu));
Q_t = E*(1.0-2*nu)/(2*dick*(1.0+nu)*(1.0-2*nu));

if(smallscale == 0)
{
/*---------------------------------------------- get history variables ---*/
  kappa_n_last = ele->e.interf->elewa[0].ipwa[ip].kappa_n;
  kappa_t_last = ele->e.interf->elewa[0].ipwa[ip].kappa_t;

/*---------- calculate displacement jumps[un],[ut],[DELTAun],[DELTAut] ---*/
  if_jumpu(ele,bop,disjump,Deltadisjump);/*disjump[0] tang. disjump[1] normal*/
}
#ifdef D_MLSTRUCT
else
{
  ID  = actsmele->Id;
/*---------------------------------------------- get history variables ---*/
  kappa_n_last =  ele->e.w1->sm_eledata[ID].sm_GPdata[ip].kappa_n;
  kappa_t_last =  ele->e.w1->sm_eledata[ID].sm_GPdata[ip].kappa_t;

/*---------------- get displacement jumps[un],[ut],[DELTAun],[DELTAut] ---*/
  disjump[0] = jumpu_tot[0];
  disjump[1] = jumpu_tot[1];
}
#endif /* D_MLSTRUCT */

/*----------------------------------------- calculate elastic stresses ---*/
T_elast[0] = Q_t * disjump[0];
T_elast[1] = Q_n * disjump[1];

/*---------------- calculate actual equvalent strain epsv_n and epsv_t ---*/
/*  Equival = 1 -  ENERGY RELEASE RATE CONCEPT - THERMODYNAMICS           */
/*  Equival = 4 -  eps_v_n/t = abs([u_n/t])                               */
/*------------------------------------------------------------------------*/
switch(equival)
{
  case 1:
    epsv_t = Q_t * disjump[0] * disjump[0]/2.0;
    epsv_n = Q_n * disjump[1] * disjump[1]/2.0;
    epsv_t_deriv = T_elast[0];
    epsv_n_deriv = T_elast[1];
  break;
  case 4:
    epsv_t = sqrt(disjump[0] * disjump[0]);
    epsv_n = sqrt(disjump[1] * disjump[1]);
    epsv_t_deriv = FSIGN(disjump[0]);
    epsv_n_deriv = FSIGN(disjump[1]);
  break;
  default:
    dserror(" unknown typ of interface equivalent strains ");
  break;
}

/*------------------------------------------------------------------------*/
/*                               NORMAL DIRECTION                         */
/*------------------------------------------------------------------------*/

/*------------------------------------------------------------ ELASTIC ---*/
if (epsv_n <= kappa0_n && kappa_n_last <= kappa0_n)
{
  damage_n = 0.0;
  DeltaD[1][1] = 0.0;
  kappa_n = epsv_n;
}/*end: if (epsv_n <= kappa0_n && kappa_n_last <= kappa0_n)*/
/*------------------------------------------------------------ DAMAGED ---*/
else
{
  /*---------------------------------------------------------- loading ---*/
  if (epsv_n - kappa_n_last >= -tol)
  {
    fac_n = 1.0;
    kappa_n = epsv_n;
  }
  /*-------------------------------------------------------- unloading ---*/
  else
  {
    fac_n = 0.0;
    kappa_n = kappa_n_last;
  }
  /*-------------------------------------------------- damage variable ---*/
  damage_n = 1.0 - (kappa0_n/kappa_n)*(1.0-alpha_n+alpha_n*exp(beta_n*(kappa0_n-kappa_n)));

  /*---------------- change of mat.tangent w.r.t.secant due to loading ---*/
  damage_n_deriv = kappa0_n*(1.0-alpha_n)/(kappa_n*kappa_n)+
                   kappa0_n*alpha_n*exp(beta_n*(kappa0_n-kappa_n))*
                   (1.0/kappa_n+beta_n)/kappa_n;

  DeltaD[1][1] = fac_n * damage_n_deriv * epsv_n_deriv * T_elast[1];

}/*end: else to if (epsv_n <= kappa0_n && kappa_n_last <= kappa0_n)*/



/*------------------------------------------------------------------------*/
/*                             TANGENTIAL DIRECTION                       */
/*------------------------------------------------------------------------*/

/*------------------------------------------------------------ ELASTIC ---*/
if (epsv_t <= kappa0_t && kappa_t_last <= kappa0_t)
{
  damage_t = 0.0;
  DeltaD[0][0] = 0.0;
  kappa_t = epsv_t;
}/*end: if (epsv_t <= kappa0_t && kappa_t_last <= kappa0_t)*/
/*------------------------------------------------------------ DAMAGED ---*/
else
{
  /*---------------------------------------------------------- loading ---*/
  if (epsv_t - kappa_t_last >= -tol)
  {
    fac_t = 1.0;
    kappa_t = epsv_t;
  }
  /*-------------------------------------------------------- unloading ---*/
  else
  {
    fac_t = 0.0;
    kappa_t = kappa_t_last;
  }
  /*-------------------------------------------------- damage variable ---*/
  damage_t = 1.0 - (kappa0_t/kappa_t)*(1.0-alpha_t+alpha_t*exp(beta_t*(kappa0_t-kappa_t)));

  /*---------------- change of mat.tangent w.r.t.secant due to loading ---*/
  damage_t_deriv = kappa0_t*(1.0-alpha_t)/(kappa_t*kappa_t)+
                   kappa0_t*alpha_t*exp(beta_t*(kappa0_t-kappa_t))*
                   (1.0/kappa_t+beta_t)/kappa_t;

  DeltaD[0][0] = fac_t * damage_t_deriv * epsv_t_deriv * T_elast[0];

}/*end: else to if (epsv_t <= kappa0_t && kappa_t_last <= kappa0_t)*/

/*----------------------------------------------------------- STRESSES ---*/
T[1] = (1.0 -damage_n) * T_elast[1];
T[0] = (1.0 -damage_t) * T_elast[0];
/*--------------------------------------------------- MATERIAL TANGENT ---*/
D[1][0] = 0.0;
D[0][1] = 0.0;
D[1][1] = (1.0 -damage_n) * Q_n - DeltaD[1][1];
D[0][0] = (1.0 -damage_t) * Q_t - DeltaD[0][0];
/*-------------------------------------------- STORE HISTORY IF UPDATE ---*/
if(istore==1 && smallscale ==0)
{
   ele->e.interf->elewa[0].ipwa[ip].dt = damage_t;
   ele->e.interf->elewa[0].ipwa[ip].dn = damage_n;
   ele->e.interf->elewa[0].ipwa[ip].kappa_t = kappa_t;
   ele->e.interf->elewa[0].ipwa[ip].kappa_n = kappa_n;
   ele->e.interf->elewa[0].ipwa[ip].Tt = T[0];
   ele->e.interf->elewa[0].ipwa[ip].Tn = T[1];
}
#ifdef D_MLSTRUCT
else if(istore==1 && smallscale ==1)
{
   ele->e.w1->sm_eledata[ID].sm_GPdata[ip].T[0] = T[0];
   ele->e.w1->sm_eledata[ID].sm_GPdata[ip].T[1] = T[1];
   ele->e.w1->sm_eledata[ID].sm_GPdata[ip].dt = damage_t;
   ele->e.w1->sm_eledata[ID].sm_GPdata[ip].dn = damage_n;
   ele->e.w1->sm_eledata[ID].sm_GPdata[ip].kappa_t = kappa_t;
   ele->e.w1->sm_eledata[ID].sm_GPdata[ip].kappa_n = kappa_n;

}
#endif /* D_MLSTRUCT */
/*------------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of if_mat_thermodyn */
















/*------------------------------------------------------------------------*/
#endif /*D_INTERF*/
/*! @} (documentation module close)*/
#endif

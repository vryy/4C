/*!----------------------------------------------------------------------
\file
\brief contains the routine 'wge_mat_damage' which calculates stresses and
       material tangentes for isotropic gradient enhanced damage

*----------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "wallge.h"
#include "wallge_prototypes.h"

/*!
\addtogroup WALLGE
*//*! @{ (documentation module open)*/


/*!----------------------------------------------------------------------
\brief contains the routine 'wge_mat_damage' which calculates stresses and
       material tangentes for isotropic gradient enhanced damage

*----------------------------------------------------------------------*/
void wge_mat_damage(INT      equival,  /* flag for equivalent strains     */
                    INT      damtyp,   /* flag for Damage-Typ             */
                    DOUBLE   youngs,   /* young's modulus                 */
                    DOUBLE   nue,      /* poisson's ratio                 */
                    DOUBLE   kappa_0,  /* initial damage equivalent strain*/
                    DOUBLE   kappa_m,  /* factor for damage-law           */
                    DOUBLE   alpha,    /* factor for exp. damage-law      */
                    DOUBLE   beta,     /* factor for exp. damage-law      */
                    DOUBLE   k_fac,    /* de Vree                         */
                    ELEMENT *ele,      /* actual element                  */
                    DOUBLE **bopd,     /* B-operator for displacements    */
                    DOUBLE  *functe,   /* Ansatz-funct. for equiv. strain */
                    DOUBLE **bope,     /* B-operator for equiv. strain    */
                    INT      ip,       /* integration point Id            */
                    DOUBLE  *stress,   /* stress vector                   */
                    DOUBLE  *eps_vl,   /* local equivalent strain         */
                    DOUBLE  *eps_vnl,  /* nonlocal equivalent strain      */
                    DOUBLE  *grad_eps_vnl,/* grad of nonlocal equi.strain */
                    DOUBLE **D,        /* 1. Material tangent             */
                    DOUBLE  *E,        /* 2. Material tangent             */
                    DOUBLE  *F,        /* 3. Material tangent             */
                    INT      istore,
                    INT      newval)

{
#ifdef D_WALLGE
/*----------------------------------------------------------------------*/
INT i,j,k;
INT yip;                     /* loading yip=1, unloading yip=0          */
/*-------------------------------------------- for local degredation ---*/
INT flag;
/*------------------------------------------------------------------------*/
DOUBLE kappa_n, kappa;       /* max. of ever reached nonl equiv.strain  */
DOUBLE damage,dam_derv;      /* damage, (partial D)/(partial eps_equiv) */
DOUBLE strain[4];            /* strain vector                           */
DOUBLE epsilon[3][3];        /* strain tensor                           */
DOUBLE sigma[3][3];          /* stress tensor                           */
DOUBLE sigma_el[3][3];       /* elastic stress tensor                   */
DOUBLE delta[3][3];          /* Kroneker delta (unit tensor)            */
DOUBLE C_ed[3][3][3][3];     /* 1. tangent (d sig_ed)/(d eps)           */
DOUBLE E_ed[3][3];           /* 2. tangent (d sig_ed)/(d eps_v_nl)      */
DOUBLE F_ed[3][3];           /* 3. tangent (d eps_v)/(d eps)            */
DOUBLE tol = -1.0E-10;



#ifdef DEBUG
dstrc_enter("wge_mat_damage");
#endif
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/*        Routine is called by if_cal_stress -> only get stresses       */
/*----------------------------------------------------------------------*/
if(newval==1)
{
 for (i=0; i<4; i++)
 {
  stress[i] = ele->e.wallge->elwa[0].iptwa[ip].sig[i];
 }
goto end;
}
/*----------------------------------------------------------------------*/
/*        Routine is called by wge_static_ke -> stresses+tangents       */
/*----------------------------------------------------------------------*/

/*-------------------- local degredation to initiate localization ---*/
# if 0
flag=0;
for (k=0; k<4; k++)
{
if(ele->node[k]->x[0]<=0.05000000001 )
{
flag=1;
}
}
for (i=0; i<4; i++)
{
if(ele->node[i]->x[0]>=0.0400001  && flag==1)
{
kappa_0 = kappa_0*0.9;
}
}
# endif
/*--------- compute strain nonl.equiv.strain, grad nonl.equiv.strain ---*/
wge_strain(ele,bopd,functe,bope,strain,eps_vnl,grad_eps_vnl);
/*------------------------------------------------ compute strain_33 ---*/
if(ele->e.wallge->wgetype==pl_strain)  strain[3]=0.0;
if(ele->e.wallge->wgetype==pl_stress)  strain[3]=-(nue*(strain[0]+strain[1]))/(1.0-nue);
/*------------------------------ from strain-vector to strain-tensor ---*/
w1_4to9(strain,epsilon);
/*------------------------------ max. ever reached equivalent strain ---*/
kappa_n = ele->e.wallge->elwa[0].iptwa[ip].kappa;

/*----------------------------------------------------------------------*/
/*        Is it still ELASTIC, DAMAGE INCREASE or UN-/RELOADING         */
/*----------------------------------------------------------------------*/

/*----------------------------- loading (elastic or already damaged) ---*/
if(*eps_vnl - kappa_n >= tol)
{
 yip = 1;
 kappa = *eps_vnl;
 # if 0
 printf("Element=%d",ele->Id);
 printf(" GP=%d",ip);
 printf(" Belastung: kappa=%le",kappa);
 # endif
}
/*-------------------- un- or reloading (elastic or already damaged) ---*/
else
{
 yip = 0;
 kappa = kappa_n;
 # if 0
 printf("Element=%d",ele->Id);
 printf(" GP=%d",ip);
 printf(" Entlastung: kappa=%le",kappa);
 # endif
}
/*----------------------------------------------------------------------*/
/*             compute stresses, local equiv. strain and tangentes      */
/*----------------------------------------------------------------------*/
/*----------------------------------------- compute Damamge variable ---*/
wge_damvar(damtyp,kappa,kappa_0,eps_vnl,alpha,beta,&damage,&dam_derv);
/*------------------------------------------- produce Kroneker Delta ---*/
w1_kroneker(delta);
/*------------------- compute local equivalent strain and derivative ---*/
wge_epsequiv(equival,epsilon,delta,k_fac,nue,eps_vl,F_ed);
/*--------------------------- calculate elastic and damaged stresses ---*/
wge_stress(youngs,nue,delta,epsilon,damage,sigma_el,sigma);
/*-------------------------------- calculate elasto-damage tangente1 ---*/
wge_tangent(youngs,nue,delta,damage,C_ed);
/*-------------------------------- calculate elasto-damage tangente2 ---*/
/*- if elastic:dam_deriv=0->E=0, if damaged un-or reloading:yip=0->E=0 -*/
for (i=0; i<3; i++)
{
 for (j=0; j<3; j++)
 {
  E_ed[i][j] = - yip * dam_derv * sigma_el[i][j];
 }
}
/*----------------------------------------------------------------------*/
/*               from tensor- back to vetor/matrix-notation             */
/*----------------------------------------------------------------------*/

/*------ from stresstensor to stress vector (just another arangement)---*/
w1_9to4(sigma,stress);
/*--- from tangenttensor to tangentmatrix D (just another arangement)---*/
w1_81to16(C_ed,D);
/*--- from tangenttensor to tangentvector E (just another arangement)---*/
w1_9to4(E_ed,E);
/*--- from tangenttensor to tangentvector F (just another arangement)---*/
w1_9to4(F_ed,F);
/*--------------------------------- condensation of tangent matrixes ---*/
wge_condense(D,E,F,ele->e.wallge->wgetype,nue);

/*----------------------------------------------------------------------*/
/*           store statevariables of converged state (update)           */
/*----------------------------------------------------------------------*/
if(istore==1)
{
  ele->e.wallge->elwa[0].iptwa[ip].kappa  = kappa;
  ele->e.wallge->elwa[0].iptwa[ip].damage = damage;
  ele->e.wallge->elwa[0].iptwa[ip].aequistrain = *eps_vl;
  ele->e.wallge->elwa[0].iptwa[ip].aequistrain_nl = *eps_vnl;
  ele->e.wallge->elwa[0].iptwa[ip].sig[0] = stress[0];
  ele->e.wallge->elwa[0].iptwa[ip].sig[1] = stress[1];
  ele->e.wallge->elwa[0].iptwa[ip].sig[2] = stress[2];
  ele->e.wallge->elwa[0].iptwa[ip].sig[3] = stress[3];
}/*end of: if (istore==1) */

/*----------------------------------------------------------------------*/
end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

#endif /*D_WALLGE*/
return;
} /* end of wge_mat_damage */
/*----------------------------------------------------------------------*/
/*! @} (documentation module close)*/

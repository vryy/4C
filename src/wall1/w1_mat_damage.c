#ifdef D_WALL1
#include "../headers/standardtypes.h"
#include "wall1.h"
#include "wall1_prototypes.h"
/*----------------------------------------------------------------------*
 | constitutive matrix - forces - linear elastic- Damage - 2D   he 04/03|
 | plane stress, plane strain, rotational symmetry                      |
 *----------------------------------------------------------------------*/
void w1_mat_damage(DOUBLE ym,      /* young's modulus                   */
                   DOUBLE pv,      /* poisson's ratio                   */
                   INT    Equival, /* flag for equivalent strains       */
                   INT    Damtyp,  /* flag for Damage-Typ               */
                   DOUBLE Kappa_0, /* initial damage equivalent strain  */
                   DOUBLE Kappa_m, /* factor for damage-law             */
                   DOUBLE Alpha,   /* factor for expon. damage-function */       
                   DOUBLE Beta,    /* factor for expon. damage-function */
                   DOUBLE k_fac,   /* factor for de Vree                */       
                   ELEMENT   *ele, /* actual element                    */
                   WALL_TYPE wtype,/* plane stress/strain...            */         
                   DOUBLE **bop,   /* derivative operator               */
                   DOUBLE  *gop,
                   DOUBLE  *alpha,
                   INT ip,         /* integration point Id              */
                   DOUBLE *stress, /* vector of stresses                */
                   DOUBLE **d,     /* constitutive matrix               */
                   INT istore,     /* controls storing of stresses      */
                   INT newval)     /* controls eval. of stresses        */ 
{
/*----------------------------------------------------------------------*/
INT i,j,k,l;
/*-------------------------------------------- for local degredation ---*/
INT flag;
/*------------------------------------------------------------------------*/
DOUBLE yip;
DOUBLE disd[5];
DOUBLE sig[4];
DOUBLE sig_esz[4];
DOUBLE eps[4];
DOUBLE eps_esz[4];
DOUBLE strain[4];
DOUBLE d4_esz[4];
DOUBLE eta,kappa;
DOUBLE damage,d_mal_eps;
DOUBLE dam_deriv;
DOUBLE sigma[3][3];
DOUBLE epsilon[3][3];
DOUBLE sigma_el[3][3];
DOUBLE eta_der[3][3];
DOUBLE c_tan[3][3][3][3];
DOUBLE c_sec[3][3][3][3];
DOUBLE c_el[3][3][3][3];
DOUBLE delta_c[3][3][3][3];
DOUBLE delta_d[4][4];
DOUBLE delta[3][3];
DOUBLE tol = -1.0E-10;
DOUBLE d_con[4][4];
#ifdef DEBUG 
dstrc_enter("w1_mat_damage");
#endif

/*---------------------------------- for stresses ----------------------*/
  if(newval==1)
  {
    for (i=0; i<4; i++)  
    {
     sig[i]    = ele->e.w1->elewa[0].ipwa[ip].sig[i];
     stress[i] = sig[i];
    }
    goto end;
  }

/*-------------------- local degredation to initiate localization ---*/
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
if(ele->node[i]->x[0]==2.5  && flag==1)
{
ym = ym*0.95;
}
}
# endif
/*------------------------------- compute displacement derivatives ---*/        
  w1_disd (ele,bop,gop,alpha,wtype,disd) ;                  

/*--------------- get actual strains -> strain -----------------------*/
  w1_eps (disd,wtype,strain);

/*----------------------------- plain stress-ellen -------------------*/
# if 0
 if (wtype == plane_stress) 
 {
  for (i=0; i<4; i++)  
  {
   eps_esz[i] = ele->e.w1->elewa[0].ipwa[ip].eps_esz[i];
   sig_esz[i] = ele->e.w1->elewa[0].ipwa[ip].sig_esz[i];
   d4_esz[i]  = ele->e.w1->elewa[0].ipwa[ip].d4_esz[i];
  }
  
  d_mal_eps = 0.0;
  for (i=0; i<3; i++)  
  {
   d_mal_eps = d_mal_eps + d4_esz[i]*(strain[i]-eps_esz[i]);
  }
    
  if(d4_esz[3] != 0.0)
  {
   strain[3] = eps_esz[3] - (d_mal_eps+sig_esz[3])/d4_esz[3]; 
  }
  
  if(d4_esz[3] == 0.0)
  {
   strain[3] = eps_esz[3]; 
  }
 }
# endif 
/*----------------------------- plain stress-andrea -------------------*/
 if (wtype == plane_stress) 
 {
  strain[3] = - (pv*(strain[0]+strain[1]))/(1.0 - pv);
 }
/*-- calculate strainvector to straintensor only for strains -----------*/
  w1_4to9(strain,epsilon);
/*--------------- calculate strainvector to straintensor ---------------*/
  w1_kroneker(delta);
/*------------------------- get elast. stresses ------------------------*/
  w1_stress_ela(ym,pv,epsilon,sigma_el,delta);
/*----------------- get actual equival. strains an derivatives-> eta ---*/
  w1_equi_eps(epsilon,sigma_el,delta,Equival,&eta,eta_der,ym,pv,k_fac);
/*----------------- get history-parameter kappa ------------------------*/
  kappa  = ele->e.w1->elewa[0].ipwa[ip].kap;
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*--------Switch between ELASTIC, DAMAGED: loading, unloading ----------*/
/*----------------------------------------------------------------------*/
/*-------- ELASTIC -----------------------------------------------------*/
if (eta <= Kappa_0 && kappa <= Kappa_0)
{
 damage = 0.0;
 for(i=0; i<3; i++)
 for(j=0; j<3; j++)
 for(k=0; k<3; k++)
 for(l=0; l<3; l++)
 delta_c[i][j][k][l] = 0.0; 
 # if 0
 printf("Element=%d",ele->Id);
 printf(" GP=%d",ip);
 printf(" Elastisch, D=%le",damage);
 printf(" kappa=%le\n",kappa);
 # endif 
}

/*-------- DAMAGED ----------------------------------------------------*/
else
{
 # if 0
 printf("Element=%d",ele->Id);
 printf(" GP=%d",ip);
 # endif 
/*-------- loading ----------------------------------------------------*/
if (eta-kappa >= tol)
{
 kappa = eta;
 yip   = 1.0;
 # if 0
 printf(" schaedigung: kappa=%le",kappa);
 # endif 
}
/*-------- unloading --------------------------------------------------*/
else
{
 kappa = kappa;
 yip   = 0.0;
 # if 0
  printf(" Entlastung: kappa=%le",kappa);
 # endif 
}

/*----- get actual damage-variable and derivative -> damage,dam_deriv -*/
  w1_dam_typ(&damage,&dam_deriv,kappa,Kappa_0,Kappa_m,Alpha,Beta,Damtyp);
 # if 0
 printf(" D=%le\n",damage);
 # endif 

/*----------------------------- get delta_c for c_tan -----------------*/
 for(i=0; i<3; i++)
 for(j=0; j<3; j++)
 for(k=0; k<3; k++)
 for(l=0; l<3; l++)
 delta_c[i][j][k][l] = -dam_deriv*yip*eta_der[k][l]*sigma_el[i][j]; 
}

/*----------------- Calculate Elastizitaetstensor ----------------------*/
  w1_mat_ela(ym,pv,delta,c_el); 

/*---------------------------- get Secantentensor ---------------------*/
  w1_sec(damage,c_el,c_sec);

/*---------------------------- get stresstensor -----------------------*/
  w1_stress(c_sec,epsilon,sigma);

/*----------------- Calculate Stresstensor to Stressvector ------------*/
  w1_9to4(sigma,sig); 
  
/*----- Calculate 4-stufiger Sekantentensor to 2-stufigem Tensor ------*/
  w1_81to16(c_sec,d); 

/*----- Calculate 4-stufiger Tensor delta_c to 2-stufigem Tensor ------*/
  w1_81to16_1(delta_c,delta_d); 

/*-----------------------------calculate c_tan ------------------------*/
 for(i=0; i<4; i++)
 {
  for(j=0; j<4; j++)
  {
   d[i][j] = d[i][j] + delta_d[i][j];
  }
 }

/*----------------------------- plain strain -------------------------*/
 if (wtype == plane_strain) 
 {
  for (i=0; i<4; i++)
  {
   d[i][3] = 0.0;
   d[3][i] = 0.0;
  }
 }

/*----------------------------- plain stress-ellen  ------------------*/
# if 0
 if (wtype == plane_stress) 
 {
   for (i=0; i<4; i++)  
   {
    ele->e.w1->elewa[0].ipwa[ip].eps_esz[i] = strain[i];
    ele->e.w1->elewa[0].ipwa[ip].sig_esz[i] = sig[i];
    ele->e.w1->elewa[0].ipwa[ip].d4_esz[i]  = d[3][i];
   }
  w1_cond(sig,d);
 }
# endif 
/*----------------------------- plain stress-andrea  ------------------*/
 if (wtype == plane_stress)
 { 
  if(d[3][3] != 0.0)
  {
    for(i=0; i<3; i++)
    {
     for (j=0; j<3; j++)
     {
      d_con[i][j]=d[i][j]-d[i][3]*d[3][j]/d[3][3];
     }  
    }  
    for(i=0; i<3; i++)
    {
     for (j=0; j<3; j++)
     {
      d[i][j]=d_con[i][j];
     }  
    }
    for(i=0; i<4; i++)
    {
      d[i][3] = 0.0;
      d[3][i] = 0.0;
    }
  }
 }
/*----------------------------- put sig to stress ---------------------*/  
  for (i=0; i<4; i++) stress[i] = sig[i];
  
/*---------------------------------------------------------------------*/  
end:
/*---------------------------------------------------------------------*/
  if(istore==1)
  {
   for (i=0; i<4; i++)  ele->e.w1->elewa[0].ipwa[ip].sig[i] = stress[i];

/*----------------- put history-parameter kappa -----------------------*/
   ele->e.w1->elewa[0].ipwa[ip].kap = kappa;
/*----------------- put damage on GP for output -----------------------*/
   ele->e.w1->elewa[0].ipwa[ip].kappa  = kappa;
   ele->e.w1->elewa[0].ipwa[ip].dam = damage;
   ele->e.w1->elewa[0].ipwa[ip].damage = damage;
   ele->e.w1->elewa[0].ipwa[ip].aequistrain = eta;
  }

/*---------------------------------------------------------------------*/

#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of w1_mat_damage */
#endif

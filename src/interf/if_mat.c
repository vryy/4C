/*!-----------------------------------------------------------------------
\file
\brief contains the routine 'if_mat' which is the material law 
 for the interface element
*-----------------------------------------------------------------------*/
#ifdef D_INTERF
#include "../headers/standardtypes.h"
#include "interf.h"
#include "interf_prototypes.h"

/*! 
\addtogroup INTERF
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief  contains the material law for the interface element

<pre>                                                              mn 05/03 
This routine computes the constitutive  matrix of the interface element

</pre>
\param *mat          STVENANT   (i)   blabal

\warning There is nothing special to this routine
\return void                                               
\sa calling:   ---; 
    called by: if_static_ke();

*----------------------------------------------------------------------*/
void if_mat(ELEMENT   *ele,
            MATERIAL  *mat, 
            DOUBLE   **bop,
            DOUBLE   **D,
            DOUBLE    *T,
            INT        ip,
            DOUBLE     istore,
            DOUBLE     newval) 
{
INT i,j;
INT yip;
DOUBLE E,K,G,thick,Q,delta_n,delta_t,mu,Ynmax,Ytmax;
DOUBLE dn,dn_neu,Yn;
DOUBLE dt,dt_neu,dt_tr,Yt,Yt_tr,Tt_tr,ut_pl,ut_pl_neu,f_tr;
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

/*------------------------------------------------------------------------*/
/*        Routine is called by if_cal_stress -> only get stresses         */
/*------------------------------------------------------------------------*/
/*if(newval==1)*/
/*{*/
/*  T[0]   = ele->e.interf->elewa[0].ipwa[ip].Tt;*/
/*  T[1]   = ele->e.interf->elewa[0].ipwa[ip].Tn;*/
/*goto end;*/
/*}*/
/*------------------------------------------------------------------------*/
yip    = ele->e.interf->elewa[0].ipwa[ip].yip;
/*------------------------------------------------------------------------*/
/*    Global praedictor: yip > 0 (excaption: first load step:yip=-1)      */
/*------------------------------------------------------------------------*/

if (yip>0)
{
  T[0] = ele->e.interf->elewa[0].ipwa[ip].Tt;
  T[1] = ele->e.interf->elewa[0].ipwa[ip].Tn;
  for (i=0; i<2; i++)
  {
    for (j=0; j<2; j++)
    {
      D[i][j] = ele->e.interf->elewa[0].ipwa[ip].Q[i][j];
    }
  }
  ele->e.interf->elewa[0].ipwa[ip].yip = -yip;
}/*end of: if (yip>0) */

/*------------------------------------------------------------------------*/
/*    Global correktor or updateing or praedictor of first loadstep       */
/*------------------------------------------------------------------------*/

else if (yip < 0)
{
  /*----------------------------------------------------------------------*/
  Ynmax  = (E * delta_n * delta_n) / (TWO * thick);
  Ytmax  = (G * delta_t * delta_t) / (TWO * thick);
  /*-------------------------------------------- reinitalization to zero---*/
  /*-- calc. tang. and normal displ jumps [un],[ut],[DELTAun],[DELTAut] ---*/        
  if_jumpu(ele,bop,disjump,Deltadisjump);
  /*----------------------------------- Werte des letzten Lastschrittes ---*/
  dn = ele->e.interf->elewa[0].ipwa[ip].dn;
  dt = ele->e.interf->elewa[0].ipwa[ip].dt;
  ut_pl  = ele->e.interf->elewa[0].ipwa[ip].jump_ut_pl;
  /*---------------------- is the interface in tension or compression? ---*/
  /*----------------------------------------------------------------------*/
  /*                            tension case                              */
  /*----------------------------------------------------------------------*/
  if (disjump[1]>0.0)
  {
    /*-------------------------------------------- 1. Normal-direction ---*/
    /*--------------------- Normaldekohaesionsfortschritt (Schaedigung)---*/
    if (disjump[1]*Deltadisjump[1]>0.0)
    {
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
      dn_neu  = dn;
      T[1]      = (E * disjump[1] * (1.0 - dn_neu)) / thick;
      D[1][0] = 0.0;
      D[1][1] = ( E * (1.0 - dn_neu)) / thick;
    }
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
         D[0][0] = 0;
# endif
       } 
    }
    /*------------------------- Tangentialdekohaesionsstop (elastisch) ---*/
    else
    {
       /*----------- Achtung: allererster praediktor hat wohl dt=0  ---*/
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
       if(disjump[0]*Deltadisjump[0]>0.0)
       {
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
  }/*end of: if (istore==1) */  
}/*end of: if (yip<0) */

/*----------------------------------------------------------------------*/
end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of if_mat */

/*----------------------------------------------------------------------*/
#endif /*D_INTERF*/
/*! @} (documentation module close)*/

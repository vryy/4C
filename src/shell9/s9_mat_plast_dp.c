/*!----------------------------------------------------------------------
\file
\brief contains the routine 
 - s9_mat_plast_dp: which calculates the constitutive matrix for a
                    'Drucker Prager'-Plasticity with linear hardening
                    for a shell9 element


<pre>
Maintainer: Stefan Hartmann
            hartmann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hartmann/
            0771 - 685-6120
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SHELL9
#include "../headers/standardtypes.h"
#include "../materials/mat_prototypes.h"
#include "shell9.h"

/*! 
\addtogroup SHELL9 
*//*! @{ (documentation module open)*/ 

/*!----------------------------------------------------------------------
\brief  shell9 element: consitutive matrix for 'Drucker Prager'-Plasticity                                    

<pre>                                                            sh 09/03
This routine calculates the constitutive matrix and forces for a 'Drucker Prager'
Plasticity with a combined linear isotropic/kinematic hardening law.
</pre>
\param  DOUBLE     ym        (i)  young's modulus              
\param  DOUBLE     pv        (i)  poisson's ratio              
\param  DOUBLE     sigy      (i)  yield stress                 
\param  DOUBLE     eh        (i)  hardening modulus            
\param  DOUBLE     gf        (i)  fracture energy              
\param  DOUBLE     betah     (i)  controls the iso/kin hard.   
\param  DOUBLE     phi       (i)  angle of friction  
\param  ELEMENT   *ele       (i)  actual element               
\param  INT        ip        (i)  integration point Id         
\param  INT        actlay    (i)  actual layer                 
\param  DOUBLE     stress[6] (o)  vector of stresses [11,22,33,12,23,13]
\param  DOUBLE     strain[6] (i)  vector of strains  [11,22,33,12,23,13]
\param  DOUBLE   **d         (o)  constitutive matrix          
\param  INT        istore    (i)  controls storing of stresses 
\param  INT        newval    (i)  controls eval. of stresses   

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: s9_call_mat()     [s9_call_mat.c]

*----------------------------------------------------------------------*/
void s9_mat_plast_dp(DOUBLE     ym,        
                     DOUBLE     pv,        
                     DOUBLE     sigy,      
                     DOUBLE     eh,        
                     DOUBLE     gf,        
                     DOUBLE     betah,     
                     DOUBLE     phi,       
                     ELEMENT   *ele,       
                     INT        ip,        
                     INT        actlay,    
                     DOUBLE     stress[6], 
                     DOUBLE     strain[6], 
                     DOUBLE   **d,         
                     INT        istore,    
                     INT        newval)    
{
INT i;
DOUBLE sig[6];     /*stresses from last update -> WA [11,22,33,12,23,13]*/
DOUBLE eps[6];     /*strains from last updat -> WA [11,22,33,12,23,13]*/
DOUBLE qn[6];      /*backstress vector from last update -> WA [11,22,33,12,23,13]*/
INT iupd=0;
INT yip;
DOUBLE epstn;
DOUBLE dia;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("s9_mat_plast_dp");
#endif
/*----------------------------------------------------------------------*/
  dia     = ele->e.s9->dia;
/*----------------------------- get old values -> sig, eps,epstn,yip ---*/
  for (i=0; i<6; i++)
  {
    sig[i] = ele->e.s9->elewa[actlay].ipwa[ip].sig[i];
    eps[i] = ele->e.s9->elewa[actlay].ipwa[ip].eps[i];
    qn[i]  = ele->e.s9->elewa[actlay].ipwa[ip].qn[i];
  }
  yip   = ele->e.s9->elewa[actlay].ipwa[ip].yip;
  epstn = ele->e.s9->elewa[actlay].ipwa[ip].epstn;

  if(newval==1)
  {
    for (i=0; i<6; i++)  stress[i] = sig[i];
    goto end;
  }
  /*-------------------------------------------------------------------*/

  if(newval!=1)  /*Check of yield criteria with possible return*/
  {
   /*Aufruf der Materialroutine, allg. 3D*/
   mat_pl_dp_lin_main(ym,
                      pv,
                      sigy,
                      eh,
                      gf,
                      betah,
                      phi,
                      stress,   /*stress*/
                      strain,
                      d,        /*Material-Matrix to be calculated 3D*/
                      &iupd,    /*to be modified*/
                      &yip,     /*to be modified*/ 
                      &epstn,   /*to be modified*/
                      sig,
                      eps,
                      qn,
                      dia);
  }

/*----------------------------- put new values -> sig, eps,epstn,yip ---*/
end:
/*----------------------------------------------------------------------*/
   if(istore==1 || iupd==1)
   {
     for (i=0; i<6; i++)
     {
       ele->e.s9->elewa[actlay].ipwa[ip].sig[i] = stress[i];  
       ele->e.s9->elewa[actlay].ipwa[ip].eps[i] = strain[i];
       ele->e.s9->elewa[actlay].ipwa[ip].qn[ i] = qn[i] ;
     }
     ele->e.s9->elewa[actlay].ipwa[ip].epstn = epstn;
     ele->e.s9->elewa[actlay].ipwa[ip].yip   = yip  ;
   }
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of s9_mat_plast_dp */
/*----------------------------------------------------------------------*/
#endif /*D_SHELL9*/
/*! @} (documentation module close)*/

/*!----------------------------------------------------------------------
\file
\brief contains the routine 
 - s9_mat_plast_epc: which calculates the constitutive matrix for the
                     elastoplastic concrete material model described
                     in Dis. Menrath and Haufe
*----------------------------------------------------------------------*/
#ifdef D_SHELL9
#include "../headers/standardtypes.h"
#include "../materials/mat_prototypes.h"
#include "shell9.h"

/*! 
\addtogroup SHELL9 
*//*! @{ (documentation module open)*/ 

/*!----------------------------------------------------------------------
\brief  shell9 element: consitutive matrix for 'Elastoplastic Concrete'
        -> Menrath                                   

<pre>                                                            sh 10/03
This routine calculates the constitutive matrix and forces for for the
elastoplastic concrete material model described in Dis. Menrath and Haufe
</pre>
\param  DOUBLE     Ec        (i)  young's modulus of concrete             
\param  DOUBLE     vc        (i)  poisson's ratio of concrete             
\param  DOUBLE     ftm       (i)  tensile strength of concrete              
\param  DOUBLE     fcm       (i)  compressive strength of concrete              
\param  DOUBLE     gt        (i)  tensile fracture energie of concrete              
\param  DOUBLE     gc        (i)  compressive fracture energie of concrete              
\param  DOUBLE     dam       (i)  d_user for damage factor              
\param  DOUBLE     gamma1    (i)  fitting parameter              
\param  DOUBLE     gamma2    (i)  fitting parameter              
\param  DOUBLE     gamma3    (i)  fitting parameter              
\param  DOUBLE     gamma4    (i)  fitting parameter              
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
void s9_mat_plast_epc(DOUBLE     Ec,        
                      DOUBLE     vc,        
                      DOUBLE     ftm,       
                      DOUBLE     fcm,       
                      DOUBLE     gt,        
                      DOUBLE     gc, 
                      DOUBLE     gamma1,    
                      DOUBLE     gamma2,    
                      DOUBLE     gamma3,    
                      DOUBLE     gamma4,    
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
DOUBLE eps[6];     /*strains from last update -> WA [11,22,33,12,23,13]*/
DOUBLE kappa_t;    /*uniaxial plastic strain in tension*/
DOUBLE kappa_c;    /*uniaxial plastic strain in compression*/
DOUBLE dia;        /*internal length*/
INT iupd=0;
INT yip;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("s9_mat_plast_epc");
#endif
/*----------------------------- get internal length of this element- ---*/
  dia     = ele->e.s9->dia;
/*--------------- get old values -> sig,eps,dfds,yip,kappa_c,kappa_t ---*/
  for (i=0; i<6; i++)
  {
    sig[i] = ele->e.s9->elewa[actlay].ipwa[ip].sig[i];
    eps[i] = ele->e.s9->elewa[actlay].ipwa[ip].eps[i];
  }
  yip     = ele->e.s9->elewa[actlay].ipwa[ip].yip;
  kappa_t = ele->e.s9->elewa[actlay].ipwa[ip].kappa_t;
  kappa_c = ele->e.s9->elewa[actlay].ipwa[ip].kappa_c;

  if(newval==1)
  {
    for (i=0; i<6; i++)  stress[i] = sig[i];
    goto end;
  }
  /*-------------------------------------------------------------------*/

  if(newval!=1)  /*Check of yield criteria with possible return*/
  {
   /*Aufruf der Materialroutine, allg. 3D -> SH-Version*/
   mat_pl_epc_main(Ec,
                   vc,
                   ftm,
                   fcm,
                   gt,
                   gc,
                   gamma1,
                   gamma2,
                   gamma3,
                   gamma4,
                   dia,
                   stress,    /*stress*/
                   strain,
                   d,          /*Material-Matrix to be calculated 3D*/
                   &iupd,      /*to be modified*/
                   &yip,       /*to be modified*/ 
                   &kappa_t,   /*to be modified*/
                   &kappa_c,   /*to be modified*/
                   sig,
                   eps);
}

/*----------------------------- put new values -> sig, eps,epstn,yip ---*/
end:
/*----------------------------------------------------------------------*/
   if(istore==1 || iupd==1)
   {
     for (i=0; i<6; i++)
     {
       ele->e.s9->elewa[actlay].ipwa[ip].sig[i]  = stress[i];  
       ele->e.s9->elewa[actlay].ipwa[ip].eps[i]  = strain[i];
     }
     ele->e.s9->elewa[actlay].ipwa[ip].yip     =  yip; 
     ele->e.s9->elewa[actlay].ipwa[ip].kappa_t =  kappa_t;
     ele->e.s9->elewa[actlay].ipwa[ip].kappa_c =  kappa_c;
   }
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of s9_mat_plast_epc */
/*----------------------------------------------------------------------*/
#endif /*D_SHELL9*/
/*! @} (documentation module close)*/

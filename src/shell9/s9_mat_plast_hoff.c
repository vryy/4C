/*!----------------------------------------------------------------------
\file
\brief contains the routine 
 - s9_mat_plast_hoff: which calculates the constitutive matrix for an
                      anisotropic plastic material model, based on the
                      hoffman yield criterion 
                      -> Dis. Hoermann 

*----------------------------------------------------------------------*/
#ifdef D_SHELL9
#include "../headers/standardtypes.h"
#include "../materials/mat_prototypes.h"
#include "shell9.h"

/*! 
\addtogroup SHELL9 
*//*! @{ (documentation module open)*/ 

/*!----------------------------------------------------------------------
\brief  shell9 element: consitutive matrix for von Anisotropic-Plasticity,
based on the Hoffman-criterion (see Dis. Hoermann p.57f)                                    

<pre>                                                            sh 03/03
This routine calculates the constitutive matrix and forces for an 
Anisotropic-Plasticity model with a linear hardening law, based on the
Hoffman-criterion. 
Within this routine, everything is done in a cartesian coordinate system
with the following sorting of stresses and strains:
"brick" [11,22,33,12,23,13]
</pre>
\param  PL_HOFF   *mat       (i)  material properties for hoffman material 
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
void s9_mat_plast_hoff(
                 PL_HOFF   *mat,       /*!< material properties          */
                 ELEMENT   *ele,       /*!< actual element               */
                 INT        ip,        /*!< integration point Id         */
                 INT        actlay,    /*!< actual layer                 */
                 DOUBLE     stress[6], /*!< vector of stresses [11,22,33,12,23,13]  */
                 DOUBLE     strain[6], /*!< vector of strains  [11,22,33,12,23,13]  */
                 DOUBLE   **d,         /*!< constitutive matrix          */
                 INT        istore,    /*!< controls storing of stresses */
                 INT        newval)    /*!< controls eval. of stresses   */
{
INT i;
DOUBLE sig[6];     /*stresses from last update -> WA [11,22,33,12,23,13]*/
DOUBLE eps[6];     /*strains from last updat -> WA [11,22,33,12,23,13]*/
DOUBLE dkappa[6];  /**/
DOUBLE gamma[6];   /**/

DOUBLE rkappa[9];  /**/
DOUBLE dhard;

INT iupd=0;
INT yip;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("s9_mat_plast_hoff");
#endif
/*----------------------------------------------------------------------*/

/*----------------------------- get old values -> sig, eps, ..., yip ---*/
  for (i=0; i<6; i++)
  {
    sig[i]    = ele->e.s9->elewa[actlay].ipwa[ip].sig[i];
    eps[i]    = ele->e.s9->elewa[actlay].ipwa[ip].eps[i];
    dkappa[i] = ele->e.s9->elewa[actlay].ipwa[ip].dkappa[i];
    gamma[i]  = ele->e.s9->elewa[actlay].ipwa[ip].gamma[i];
  }

  for (i=0; i<9; i++)
  {
    rkappa[i] = ele->e.s9->elewa[actlay].ipwa[ip].rkappa[i];
  }

  yip   = ele->e.s9->elewa[actlay].ipwa[ip].yip;
  dhard = ele->e.s9->elewa[actlay].ipwa[ip].dhard;

  if(newval==1)
  {
    for (i=0; i<6; i++)  stress[i] = sig[i];
    goto end;
  }
  /*-------------------------------------------------------------------*/

  if(newval!=1)  /*Check of yield criteria with possible return*/
  {
   /*Aufruf der Materialroutine, allg. 3D -> src/materials/mat_pl_hoff.c*/
   mat_pl_hoff_main( mat,
                     stress,   /*stress*/
                     d,        /*Material-Matrix to be calculated 3D*/
                    &iupd,    /*to be modified*/
             /*zusaetzliche Parameter*/                       
                    &yip,     /*to be modified*/ 
                    &dhard,   /*to be modified*/
                     strain,
                     sig,
                     eps,
                     dkappa,
                     gamma,
                     rkappa);
  }

/*------------------------ put new values -> sig,eps,dhard,yip, ... ----*/
end:
/*----------------------------------------------------------------------*/
  if(istore==1 || iupd==1)
  {
    for (i=0; i<6; i++)
    {
      ele->e.s9->elewa[actlay].ipwa[ip].sig[i]    = stress[i]; 
      ele->e.s9->elewa[actlay].ipwa[ip].eps[i]    = strain[i];
      ele->e.s9->elewa[actlay].ipwa[ip].dkappa[i] = dkappa[i];
      ele->e.s9->elewa[actlay].ipwa[ip].gamma[i]  = gamma[i] ;
    }
    for (i=0; i<9; i++)
    {
      ele->e.s9->elewa[actlay].ipwa[ip].rkappa[i] = rkappa[i];
    }
    ele->e.s9->elewa[actlay].ipwa[ip].yip   = yip   ;
    ele->e.s9->elewa[actlay].ipwa[ip].dhard = dhard ;
  }
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of s9_mat_plast_hoff */
/*----------------------------------------------------------------------*/
#endif /*D_SHELL9*/
/*! @} (documentation module close)*/

/*!----------------------------------------------------------------------
\file
\brief contains the routine 'w1_mat_plast_mises_3D' which calculates 
 the constitutive matrix and forces using a von Mises-Plasticity 
 with a combined linear isotropic/kinematic hardening law
 the material formulation is 3D, so that the 2D-conditions from a plane
 calculation (Wall1-Element) have to be blown up ('w1mat_trans_up before) 
 calling 'mat_plast_mises_3D', which is a Element independent 3D-Routine 
 for calculating the 3D constitutive matrix. The calculated 3D-based 
 values have to be condesed back to either plane_stress or plane_strain
 conditions ('w1mat_trans_down').              
 (rotational symmetry is not implemented)   
 contains the routine 'mat_plast_mises_3D' which calcualtes the
 constitutive matrix - forces - von Mises - 3D which is element type 
 independent  
 
<pre>
Maintainer: Stefan Hartmann
            hartmann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hartmann/
            0771 - 685-6120
</pre>

*----------------------------------------------------------------------*/
#ifdef D_WALL1
#ifdef D_MAT

#include "../headers/standardtypes.h"
#include "../materials/mat_prototypes.h"
#include "wall1.h"
#include "wall1_prototypes.h"


/*! 
\addtogroup WALL1
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief  calculates the constitutive matrix for a von Mises-Plasticity

<pre>                                                             sh 08/02 
This routine calculates the constitutive matrix and forces using a 3D
von Mises-Plasticity with a combined linear isotropic/kinematic hardening 
law. Needed routines: w1mat_trans_up (2D->3D), w1mat_trans_down (3D->2D),
mat_plast_mises_3D (constitutive matrix, general 3D, element independent).
Works for plane_strain & plane_stress. Rotational symmetrie is not 
implemented yet.

</pre>
\param  DOUBLE      ym        (i)  young's modulus of concrete            
\param  DOUBLE      pv        (i)  poisson's ratio of concrete            
\param  DOUBLE      sigy      (i)  yield stress              
\param  DOUBLE      hard      (i)  hardening modulus             
\param  DOUBLE      gf        (i)  fracture energy              
\param  DOUBLE      betah     (i)  controls the iso/kin hard.             
\param  ELEMENT    *ele       (i)  actual element                  
\param  WALL_TYPE   wtype     (i)  plane stress/strain...          
\param  DOUBLE    **bop       (i)  B-Operator             
\param  DOUBLE     *gop       (i)  for incompatible modes                                  
\param  DOUBLE     *alpha     (i)  for incompatible modes   
\param  INT         ip        (i)  integration point Id        
\param  DOUBLE     *stressc   (o)  vector of stresses condensed
\param  DOUBLE    **d         (o)  constitutive matrix         
\param  INT         istore    (i)  controls storing of stresses     
\param  INT         newval    (i)  controls eval. of stresses     

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: w1_call_mat()

*----------------------------------------------------------------------*/
void w1_mat_plast_mises_3D(DOUBLE     ym,
                           DOUBLE     pv,
                           DOUBLE     sigy,
                           DOUBLE     hard,
                           DOUBLE     gf,
                           DOUBLE     betah,
                           ELEMENT   *ele,
                           WALL_TYPE  wtype,
                           DOUBLE   **bop,
                           DOUBLE    *gop,
                           DOUBLE    *alpha,
                           INT        ip,
                           DOUBLE    *stressc,
                           DOUBLE   **d,
                           INT        istore,
                           INT        newval)
{
INT i;
DOUBLE stress3D[6];  /*actual stresses*/
DOUBLE strain3D[6];  /*actual strains from displacements*/
DOUBLE strain[4];    /*actual strains from displacements [11,22,12,33]*/
DOUBLE sig3D[6];     /*stresses from last update -> WA [4]->[6]*/
DOUBLE eps3D[6];     /*strains from last updat -> WA [4]->[6]*/
DOUBLE qn3D[6];      /*backstress vector from last update -> WA [4]->[6]*/
DOUBLE eps[4];       /*strains from last update -> WA*/
INT iupd=0;
INT yip,yipc;
DOUBLE epstn;
DOUBLE dia;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("w1_mat_plast_mises_3D");
#endif
/*----------------------------------------------------------------------*/
  for (i=0; i<4; i++)
  {
    sig3D[i] = ele->e.w1->elewa[0].ipwa[ip].sig[i];  /*[11,22,12,33]  -> wall      */
    eps[i]   = ele->e.w1->elewa[0].ipwa[ip].eps[i];  /*[11,22,12,33]  -> wall      */
    qn3D[i]  = ele->e.w1->elewa[0].ipwa[ip].qn[i];   /*[11,22,12,33]  -> wall      */
  }
  for (i=4; i<6; i++)
  {
    sig3D[i]    = 0.;
    eps3D[i]    = 0.;
    strain3D[i] = 0.;
    qn3D[i]     = 0.;
  }

/*-- get additional strain e_zz --------------------------------------------------*/
    w1mat_trans_up(ym,pv,ele,wtype,bop,gop,alpha,ip,strain);
    
    /* do sorting for 3D-Material law*/
    w1_vec_switch(sig3D,2,3);      
    w1_vec_switch(qn3D,2,3);      

    /* 3D-Routine needs origianal shear strains (not doubled as for vector-matrix notation)*/
    strain3D[0] =       strain[0];
    strain3D[1] =       strain[1];
    strain3D[2] =       strain[3];
    strain3D[3] = 0.5 * strain[2];
    eps3D[0]    =       eps[0];
    eps3D[1]    =       eps[1];
    eps3D[2]    =       eps[3];
    eps3D[3]    = 0.5 * eps[2];
    
    /*Werte aus Elementinformationen -> Material soll unabhaengig von Ele sein*/
    dia   = ele->e.w1->elewa[0].dia;
    yip   = ele->e.w1->elewa[0].ipwa[ip].yip;
    yipc  = yip;       /*make a copy of yip for correct update */
    epstn = ele->e.w1->elewa[0].ipwa[ip].epstn;
    /**/
/*----------------------------------------------------------------------*/
    if(newval==1) /*sorting [11,12,22,33]*/
    {
      for (i=0; i<4; i++)  stressc[i] = ele->e.w1->elewa[0].ipwa[ip].sigc[i];
      goto end;
    }
/*----------------------------------------------------------------------*/
   if(newval!=1)  /*Check of yield criteria with possible return*/
   {
#ifdef D_MAT
    /*Aufruf der Materialroutine, allg. 3D*/
    mat_pl_mises_lin_main(
                       ym,
                       pv,
                       sigy,
                       hard,
                       gf,
                       betah,
                       stress3D, /*stress3d to be calculated (output)*/
                       strain3D, /*strain3d  (input)*/
                       d,        /*Material-Matrix to be calculated 3D*/
                       &iupd,    /*to be modified*/
                       &yip,     /*to be modified*/ 
                       &epstn,   /*to be modified*/
                       sig3D,    /*(input)*/
                       eps3D,    /*(input)*/
                       qn3D,     /*qn3d to be calculated (output)*/
                       dia);
#else
   dserror("mat-funcs needed but not defined!\n");
#endif
   }
    /* do sorting back to wall [11,22,12,33] */
    w1_vec_switch(stress3D,2,3);      
    w1_vec_switch(sig3D,2,3);      
    w1_vec_switch(qn3D,2,3);      
    w1_matrix_switch(d,2,3,6);      

    w1mat_trans_down(d,         /*Material-Matrix to be condensed 3D->2D*/
                     ele,
                     wtype,
                     ip,
                     yipc,
                     stressc,    /*condensed stresses [11,22,12,33] for calculation*/
                     sig3D,      /*[11,22,12,33]*/
                     eps,        /*[11,22,12,33]*/
                     stress3D,   /*to be condensed [11,22,12,33]*/
                     strain,     /*[11,22,12,33]*/
                     qn3D);      /*to be condensed [11,22,12,33]*/
/*----------------------------- put new values -> sig, eps,epstn,yip ---*/
end:
/*----------------------------------------------------------------------*/
   if(istore==1 || iupd==1)
   {
     for (i=0; i<4; i++)
     {
       ele->e.w1->elewa[0].ipwa[ip].sig[i]  = stress3D[i];  
       ele->e.w1->elewa[0].ipwa[ip].sigc[i] = stressc[i];   /*condensed stress*/ 
       ele->e.w1->elewa[0].ipwa[ip].eps[i]  = strain[i];
       ele->e.w1->elewa[0].ipwa[ip].qn[ i]  = qn3D[i] ;
     }
     ele->e.w1->elewa[0].ipwa[ip].epstn = epstn;
     ele->e.w1->elewa[0].ipwa[ip].yip   = yip  ;
   }
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of w1_mat_plast_mises_3D */
/*----------------------------------------------------------------------*/
#endif /* D_MAT */
#endif /*D_WALL1*/
/*! @} (documentation module close)*/


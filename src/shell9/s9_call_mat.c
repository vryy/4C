/*!----------------------------------------------------------------------
\file
\brief contains the routine 
 - s9_call_mat: which calls the different material laws implemented for this 
                element
 - s9_getdensity: get density out of material law


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
\brief call material laws                                      

<pre>                     m.gee 12/01             modified by    sh 01/03
This routine calls the different material laws implemented for this 
element.
</pre>
\param  ELEMENT   *ele       (i)  the element structure -> 'ele->e.s9->stresses.a.d3'
\param  MULTIMAT  *multimat  (i)  the material structure (shell9 -> multimat)
\param  DOUBLE     stress[6] (o)  PK_II stresses from constitutive law
\param  DOUBLE     strain[6] (o)  green-lagrange strains from metrics
\param  DOUBLE   **C         (o)  constitutive matrix
\param  DOUBLE   **gmkovc,.. (i)  all the metrics in ref/cur configuration
\param  INT        rot_axis  (i)  rotation axis of laminat (1=x-; 2=y-; 3=z-axis)
\param  DOUBLE     phi       (i)  angle of rotation about rot_axis
\param  INT        ip        (i)  the actual integration point
\param  INT        actlay    (i)  the actual layer
\param  INT        istore    (i)  controls storing of new stresses to wa
\param  INT        newval    (i)  controls evaluation of new stresses

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: s9_stress()     [s9_stress.c]
                             s9static_keug() [s9_static_keug.c]

*----------------------------------------------------------------------*/
void s9_call_mat(ELEMENT    *ele,
                 MULTIMAT   *multimat,  /*!< material of actual material layer */
                 DOUBLE      stress[6],
                 DOUBLE      strain[6],
                 DOUBLE    **C,         /*!< constitutive matrix               */
                 DOUBLE    **gmkovc,
                 DOUBLE    **gmkovr,
                 DOUBLE    **gmkonr,
                 DOUBLE    **gkovr,
                 DOUBLE    **gkonr,
                 INT         rot_axis,  /*!< rotation axis of laminat          */
                 DOUBLE      phi,       /*!< angle of rotation about rot_axis  */
                 INT         ip,        /*!< the actual integration point      */
                 INT         actlay,    /*!< the actual layer                  */
                 INT         istore,    /*!< controls storing of new stresses to wa */
                 INT         newval)    /*!< controls evaluation of new stresses    */
{
DOUBLE g1[3],g2[3],g3[3];
DOUBLE T[6][6];
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("s9_call_mat");
#endif
/*------------------------------------------------ switch material type */
switch(multimat->mattyp)
{
case m_stvenant:/*-------------------------- ST.VENANT-KIRCHHOFF-MATERIAL */
   s9_eps(strain,gmkovc,gmkovr); 
   
   /*get material matrix directly in curvilinear coordsys*/
   s9_mat_linel(multimat->m.stvenant,gmkonr,C);

/********** for TESTING issues ************************/
   /*calculate C-Matrix in orthonormal basis*/
/*   mat_el_iso(multimat->m.stvenant->youngs,
              multimat->m.stvenant->possionratio,
              C);*/
   /*sort C-Matrix from "brick" to "shell9" */ 
/*   s9_Msort_bs9(C);*/ 

   /*transform the C-Matrix from cartesian to curvilinear*/
/*   s9_Ccacu_sym(C,gkonr);*/
   
   /*do modification of matrix due to shear correction*/
/*   s9shear_cor(C,1.2); */
/********** for TESTING issues ************************/
   
   s9_mat_stress1(stress,strain,C);
break;
case m_el_orth:/*-----------------linear elastic, orthotropic material */
   s9_eps(strain,gmkovc,gmkovr); 

   /*linear elastic, orthotropic material matrix, formulated in cartesian coordinates -> [11,22,33,12,23,13]*/
   mat_el_orth(multimat->m.el_orth->emod1,
               multimat->m.el_orth->emod2,
               multimat->m.el_orth->emod3,
               multimat->m.el_orth->xnue12,
               multimat->m.el_orth->xnue23,
               multimat->m.el_orth->xnue13,
               multimat->m.el_orth->gmod12,
               multimat->m.el_orth->gmod23,
               multimat->m.el_orth->gmod13,
               C);

   /* get basis vectors of material coord. sys according to cartesian coord. sys*/
   s9_rot(phi,rot_axis,gkovr,g1,g2,g3);
   /* get transformation matrix T_sup_epsilon*/
   s9_Teps(g1,g2,g3,T);
   /*transform C from orthonormal material coord. sys. to global cartesian*/
   s9_Cmaca(C, T);

   /*change sorting from [11,22,33,12,23,13] -> [11,12,22,13,23,33]*/
   s9_Msort_bs9(C); 

   /*transform the C-Matrix from cartesian to curvilinear*/
   s9_Ccacu_sym(C,gkonr);   
   
   s9_mat_stress1(stress,strain,C);
break;
case m_neohooke:/*------------------------------ kompressible neo-hooke */
   dserror("neohooke not yet implemented");
break;
case m_pl_hoff:/*--- anisotropic plasticity model based on hoffman-criterion */
   s9_eps(strain,gmkovc,gmkovr); 

   /* transform strains from curvilinear to cartesian */
   s9_Ecuca(strain,gkonr);   
   /* do sorting from "shell9" = [11,12,22,13,23,33] -> "brick" = [11,22,33,12,23,13] */
   s9_Vsort_s9b(strain);  

   strain[3] = 2. * strain[3];     /*write strains as vector for transformation*/
   strain[4] = 2. * strain[4];     /*from cartesian coord. sys to material coord. sys*/
   strain[5] = 2. * strain[5]; 

   /* get basis vectors of material coord. sys according to cartesian coord. sys*/
   s9_rot(phi,rot_axis,gkovr,g1,g2,g3);
   /* get transformation matrix T_sup_epsilon*/
   s9_Teps(g1,g2,g3,T);
   /* transform strains from cartesian to material/laminat coord. sys.*/
   s9_Ecama(strain, T);

   /*von "Hoffman" Plasticity, formulated in cartesian coordinate system, sorting: [11,22,33,12,23,13] */
   s9_mat_plast_hoff(multimat->m.pl_hoff,  /* material properties          */
                     ele,                  /* actual element               */
                     ip,                   /* integration point Id         */
                     actlay,               /* actual layer                 */
                     stress,               /* vector of stresses [11,22,33,12,23,13]  */
                     strain,               /* vector of strains  [11,22,33,12,23,13]  */
                     C,                    /* constitutive matrix          */
                     istore,               /* controls storing of stresses */
                     newval);              /* controls eval. of stresses   */
                

   /*transform C from orthonormal material coord. sys. to global cartesian*/
   s9_Cmaca(C, T);
   /*change sorting from [11,22,33,12,23,13] -> [11,12,22,13,23,33]*/
   s9_Msort_bs9(C); 
   /*transform the C-Matrix from cartesian to curvilinear*/
   s9_Ccacu_unsym(C,gkonr);  

   /*transform stresses from orthonormal material coord. sys. to global cartesian*/
   s9_Smaca(stress, T);
   /* do sorting from "brick" = [11,22,33,12,23,13] -> "shell9" = [11,12,22,13,23,33] */
   s9_Vsort_bs9(stress);  
   /*transform  stresses from cartesian to curvilinear */
   s9_Scacu(stress,gkonr);
break;
case m_pl_mises:/*--------------------- von mises material law ---*/
   s9_eps(strain,gmkovc,gmkovr);        /* get strains: curvilinear*/

   /* transform strains from curvilinear to cartesian */
   s9_Ecuca(strain,gkonr);  
   /* do sorting from "shell9" = [11,12,22,13,23,33] -> "brick" = [11,22,33,12,23,13] */
   s9_Vsort_s9b(strain); 

   /*von Mises Plasticity, formulated in cartesian coordinate system, sorting: [11,22,33,12,23,13] */
   s9_mat_plast_mises(multimat->m.pl_mises->youngs,
                      multimat->m.pl_mises->possionratio,
                      multimat->m.pl_mises->Sigy,
                      multimat->m.pl_mises->Hard,
                      multimat->m.pl_mises->GF,
                      multimat->m.pl_mises->betah,
                      ele,
                      ip,
                      actlay,
                      stress,
                      strain,
                      C,
                      istore,
                      newval);

   /* do sorting from "brick" = [11,22,33,12,23,13] -> "shell9" = [11,12,22,13,23,33] */
   s9_Msort_bs9(C);
   /*transform the C-Matrix from cartesian to curvilinear*/
   s9_Ccacu_sym(C,gkonr);   
   /* do sorting from "brick" = [11,22,33,12,23,13] -> "shell9" = [11,12,22,13,23,33] */
   s9_Vsort_bs9(stress);  
   /*transform  stresses from cartesian to curvilinear */
   s9_Scacu(stress,gkonr);
break;
case m_pl_dp:/*------------------- drucker prager material law ---*/
   s9_eps(strain,gmkovc,gmkovr);        /* get strains: curvilinear*/

   /* transform strains from curvilinear to cartesian */
   s9_Ecuca(strain,gkonr);  
   /* do sorting from "shell9" = [11,12,22,13,23,33] -> "brick" = [11,22,33,12,23,13] */
   s9_Vsort_s9b(strain); 
   
   /*Drucker Prager Plasticity, formulated in cartesian coordinate system, sorting: [11,22,33,12,23,13] */
   s9_mat_plast_dp(multimat->m.pl_dp->youngs,       
                   multimat->m.pl_dp->possionratio,      
                   multimat->m.pl_dp->Sigy,     
                   multimat->m.pl_dp->Hard,       
                   multimat->m.pl_dp->GF,
                   multimat->m.pl_dp->betah,    
                   multimat->m.pl_dp->PHI,      
                   ele,      
                   ip,       
                   actlay,   
                   stress,
                   strain,
                   C,        
                   istore,   
                   newval);   
                 
   /* do sorting from "brick" = [11,22,33,12,23,13] -> "shell9" = [11,12,22,13,23,33] */
   s9_Msort_bs9(C);
   /*transform the C-Matrix from cartesian to curvilinear*/
   if (FABS(multimat->m.pl_dp->betah)<EPS8) s9_Ccacu_sym(C,gkonr);   /*only kinematic hardening -> ass. in apex*/
   else                                     s9_Ccacu_unsym(C,gkonr); /*isotropic hardening -> nonass. in apex*/ 
   /* do sorting from "brick" = [11,22,33,12,23,13] -> "shell9" = [11,12,22,13,23,33] */
   s9_Vsort_bs9(stress);  
   /*transform  stresses from cartesian to curvilinear */
   s9_Scacu(stress,gkonr);
break;
case m_pl_epc:/*------ Elasto Plastic Concrete (Menrath) material law ---*/
   s9_eps(strain,gmkovc,gmkovr);        /* get strains: curvilinear*/

   /* transform strains from curvilinear to cartesian */
   s9_Ecuca(strain,gkonr);  
   /* do sorting from "shell9" = [11,12,22,13,23,33] -> "brick" = [11,22,33,12,23,13] */
   s9_Vsort_s9b(strain); 
   
   /*Elasto Plastic Concrete (Menrath), formulated in cartesian coordinate system, sorting: [11,22,33,12,23,13] */
   s9_mat_plast_epc(multimat->m.pl_epc->youngs,       /* young's modulus              */
                    multimat->m.pl_epc->possionratio, /* poisson's ratio              */
                    multimat->m.pl_epc->ftm,          /* tensile strength             */
                    multimat->m.pl_epc->fcm,          /* compressive strength         */
                    multimat->m.pl_epc->gt,           /* tensile fracture energie     */
                    multimat->m.pl_epc->gc,           /* compressive fracture energie */
                    multimat->m.pl_epc->gamma1,       /* fitting parameter            */
                    multimat->m.pl_epc->gamma2,       /* fitting parameter            */
                    multimat->m.pl_epc->gamma3,       /* fitting parameter            */
                    multimat->m.pl_epc->gamma4,       /* fitting parameter            */
                    ele,                              /* actual element               */
                    ip,                               /* integration point Id         */
                    actlay,                           /* actual layer                 */
                    stress,                           /* vector of stresses [11,22,33,12,23,13]  */
                    strain,                           /* vector of strains  [11,22,33,12,23,13]  */
                    C,                                /* constitutive matrix          */
                    istore,                           /* controls storing of stresses */
                    newval);                          /* controls eval. of stresses   */
                 
   /* do sorting from "brick" = [11,22,33,12,23,13] -> "shell9" = [11,12,22,13,23,33] */
   s9_Msort_bs9(C);


   /*transform the C-Matrix from cartesian to curvilinear*/
/*TEST ===> check when to use sym/unsym!!!!!!!!!!1*/
   s9_Ccacu_unsym(C,gkonr);
/*TEST*/


   /* do sorting from "brick" = [11,22,33,12,23,13] -> "shell9" = [11,12,22,13,23,33] */
   s9_Vsort_bs9(stress);  
   /*transform  stresses from cartesian to curvilinear */
   s9_Scacu(stress,gkonr);
break;


default:
   dserror("Ilegal typ of material for this element");
break;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of s9_call_mat */



/*!----------------------------------------------------------------------
\brief get density out of material law                                      

<pre>                                                         m.gee 12/01             
This routine gets density out of material law
</pre>
\param  MATERIAL  *mat      (i)  the material structure
\param  DOUBLE    *density  (o)  density

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: -----

*----------------------------------------------------------------------*/
void s9_getdensity(MATERIAL   *mat, DOUBLE *density)
{
#ifdef DEBUG 
dstrc_enter("s9_getdensity");
#endif
/*----------------------------------------------------------------------*/
/*------------------------------------------------ switch material type */
switch(mat->mattyp)
{
case m_stvenant:/*-------------------------- ST.VENANT-KIRCHHOFF-MATERIAL */
   *density = mat->m.stvenant->density;
break;
case m_neohooke:/*------------------------------ kompressible neo-hooke */
   *density = mat->m.neohooke->density;
break;
default:
   dserror("Ilegal typ of material for this element");
break;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of s9_getdensity */


/*----------------------------------------------------------------------*/
#endif /*D_SHELL9*/
/*! @} (documentation module close)*/

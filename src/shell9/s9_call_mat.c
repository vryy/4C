/*!----------------------------------------------------------------------
\file
\brief contains the routine 
 - s9_call_mat: which calls the different material laws implemented for this 
                element
 - s9_getdensity: get density out of material law

*----------------------------------------------------------------------*/
#ifdef D_SHELL9
#include "../headers/standardtypes.h"
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
\param  ELEMENT   *ele      (i)  the element structure -> 'ele->e.s9->stresses.a.d3'
\param  MULTIMAT  *multimat (i)  the material structure (shell9 -> multimat)
\param  double    *stress   (o)  PK_II stresses from constitutive law
\param  double    *strain   (o)  green-lagrange strains from metrics
\param  double   **C        (o)  constitutive matrix
\param  double   **gmkovc,..(i)  all the metrics in ref/cur configuration
\param  int        rot_axis (i)  rotation axis of laminat (1=x-; 2=y-; 3=z-axis)
\param  double     phi      (i)  angle of rotation about rot_axis

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: s9_stress()     [s9_stress.c]
                             s9static_keug() [s9_static_keug.c]

*----------------------------------------------------------------------*/
void s9_call_mat(ELEMENT    *ele,
                 MULTIMAT   *multimat,    /* material of actual material layer */
                 double     *stress,
                 double     *strain,
                 double    **C,
                 double    **gmkovc,
                 double    **gmkonc,
                 double    **gmkovr,
                 double    **gmkonr,
                 double    **gkovc,
                 double    **gkonc,
                 double    **gkovr,
                 double    **gkonr,
                 int         rot_axis,
                 double      phi)
{
#ifdef DEBUG 
dstrc_enter("s9_call_mat");
#endif
/*------------------------------------------------ switch material type */
switch(multimat->mattyp)
{
case m_stvenant:/*-------------------------- ST.VENANT-KIRCHHOFF-MATERIAL */
   s9_eps(strain,gmkovc,gmkovr); 
   s9_mat_linel(multimat->m.stvenant,gmkonr,C);
   s9_mat_stress1(stress,strain,C);

   /********** for TESTING issues ************************/
   /*calculate C-Matrix and stresses in orthonormal basis*/
/*   s9_mat_linel3D(multimat->m.stvenant->youngs,
                  multimat->m.stvenant->possionratio,
                  C);
/*   s9_matrix_sort(C,"brick","s9");   

   /*transform strains from curvilinear to cartesian*/ 
/*   s9tcuca (strain,gkovr,"CUCA","E");

   s9_mat_stress1(stress,strain,C);

   /*transform  strains and stresses from cartesian to curvilinear */
/*   s9tcuca(strain,gkovr,"CACU","E");
   s9tcuca(stress,gkovr,"CACU","S");

   /*transform the C-Matrix from cartesian to curvilinear*/
/*   s9T4sym(C,gkonr);   
   s9shear_cor(C,1.2); /*do modification of matrix due to shear correction*/

   /*transform the C-Matrix from cartesian to curvilinear*/
/*   s9_mat_linel3D(multimat->m.stvenant->youngs,
                  multimat->m.stvenant->possionratio,
                  C);
   s9_Tcacu(gkonr,C);   
   s9_matrix_sort(C,"brick","s9");   
/*   s9shear_cor(C,1.2); /*do modification of matrix due to shear correction*/

/*   s9_mat_stress1(stress,strain,C);
 */
break;
case m_orthotropic:/*-----------------linear elastic, orthotropic material */
   s9_eps(strain,gmkovc,gmkovr); 
   s9_mat_orth3D(multimat->m.orthotropic,C);

   /*transform C from orthonormal material coord. sys. to global cartesian*/
   s9_Tmaca (gkovr,phi,rot_axis,C);
   /*change sorting from [11,22,33,12,23,13] -> [11,12,22,13,23,33]*/
   s9_matrix_sort(C,"brick","s9");  
   /*transform the C-Matrix from cartesian to curvilinear*/
   s9T4sym(C,gkonr);   
   
   s9_mat_stress1(stress,strain,C);
break;
case m_neohooke:/*------------------------------ kompressible neo-hooke */
   dserror("neohooke not yet implemented");
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
\param  double    *density  (o)  density

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: -----

*----------------------------------------------------------------------*/
void s9_getdensity(MATERIAL   *mat, double *density)
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

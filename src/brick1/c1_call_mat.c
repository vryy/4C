/*!----------------------------------------------------------------------
\file
\brief contains the routines 'c1_call_mat' and  'c1_call_matd'
       control programs for formulation of material law and its derivatives
       select proper material law 
       and evaluation of element stresses 

*----------------------------------------------------------------------*/
#ifdef D_BRICK1

#include "../headers/standardtypes.h"
#include "brick1.h"
#include "brick1_prototypes.h"

/*! 
\addtogroup BRICK1 
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief control program for formulation of material law

<pre>                                                              al 06/02
This routine selects proper material law evaluates element stresses.

</pre>
\param      ele ELEMENT * (i)   element data
\param      mat MATERIAL* (i)   material data
\param       ip int       (i)   integration point Id      
\param   stress double*   (o)   stress vector
\param   strain double*   (o)   strain vector
\param        d double**  (o)   constitutive matrix
\param     disd double*   (i)   displacement derivatives 
\param  g[6][6] double    (i)   transformation matrix s(glob)=g*s(loc)
\param gi[6][6] double    (i)   inverse of g          s(loc) =gi*s(glob)
\param   istore int       (i)   controls storing of new stresses to wa
\param   newval int       (i)   controls evaluation of new stresses

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: c1_cint()

*----------------------------------------------------------------------*/
void c1_call_mat(ELEMENT   *ele,
                 MATERIAL  *mat, 
                 int ip,       
                 double *stress,
                 double *strain,
                 double **d,
                 double *disd,
                 double g[6][6], 
                 double gi[6][6],
                 int istore,
                 int newval)
{
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("c1_call_mat");
#endif
/*------------------------------------------------ call material law ---*/
  switch(mat->mattyp)
  {
  case m_stvenant:/*------------------------------- linear elastic ---*/
    c1_mat_linel(mat->m.stvenant->youngs,
                 mat->m.stvenant->possionratio,
                 d);
    c1gld (d,g);/* transform local to global material matrix */
    c1mefm(strain, d, stress);
  break;
  case m_stvenpor:/*------------------------ porous linear elastic ---*/
    c1_mat_stvpor(mat, ele->e.c1->elewa->matdata, d);
    c1gld (d,g);/* transform local to global material matrix */
    c1mefm(strain, d, stress);
  break;
  case m_mfoc:/*-------------- open cell metal foam linear elastic ---*/
    c1_mat_mfoc(mat, ele->e.c1->elewa->matdata, d);
    c1gld (d,g);/* transform local to global material matrix */
    c1mefm(strain, d, stress);
  break;
  case m_mfcc:/*------------ closed cell metal foam linear elastic ---*/
    c1_mat_mfcc(mat, ele->e.c1->elewa->matdata, d);
    c1gld (d,g);/* transform local to global material matrix */
    c1mefm(strain, d, stress);
  break;
  case m_pl_mises:/*----------------------- von mises material law ---*/
    c1_mat_plast_mises(mat->m.pl_mises->youngs,
                       mat->m.pl_mises->possionratio,
                       mat->m.pl_mises->Sigy,
                       mat->m.pl_mises->Hard,
                       ele,
                        ip,
                       stress,
                       d,
                       disd,
                       g,
                       gi,
                       istore,
                       newval);
  break;
  case m_el_orth:/*-------------- elastic orthotropic material law ---*/
    c1_mat_elorth(mat->m.el_orth->emod1,
                  mat->m.el_orth->emod2,
                  mat->m.el_orth->emod3,
                  mat->m.el_orth->xnue23,
                  mat->m.el_orth->xnue13,
                  mat->m.el_orth->xnue12,
                  mat->m.el_orth->gmod12,
                  mat->m.el_orth->gmod23,
                  mat->m.el_orth->gmod13,
                  d);
    c1mefm(strain, d, stress);
  break;
  case m_pl_mises_ls:/*----- von mises material law - large strains---*/
    c1_mat_plast_mises_ls(
                       mat->m.pl_mises_ls->youngs,
                       mat->m.pl_mises_ls->possionratio,
                       mat->m.pl_mises_ls->Sigy,
                       mat->m.pl_mises_ls->Hard,
                       ele,
                       ip,
                       stress,
                       d,
                       disd,
                       istore,
                       newval);
  break;
  case m_pl_dp:/*------------------------ drucker prager material law ---*/
    dserror(" drucker prager material law not implemented for brick");
  break;
  case m_pl_epc:/*---------- elastoplastic concrete material law ---*/
    dserror(" elastoplastic concrete material law not implemented for brick");
  break;                                         
  default:
    dserror(" unknown type of material law");
  break;
  }
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of c1_call_mat */
/*----------------------------------------------------------------------*/

/*!----------------------------------------------------------------------
\brief calculates derivatives of material law

<pre>                                                              al 06/02
This routine calculates derivatives of material law.

</pre>
\param      ele ELEMENT * (i)   element data
\param      mat MATERIAL* (i)   material data
\param   stress double*   (o)   stress vector
\param   strain double*   (o)   strain vector
\param        d double**  (o)   DERIVATIVE of constitutive matrix
\param  g[6][6] double    (i)   transformation matrix s(glob)=g*s(loc)

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: opt_c1_cint()

*----------------------------------------------------------------------*/
void c1_call_matd(ELEMENT   *ele,
                 MATERIAL  *mat, 
                 double *stress,
                 double *strain,
                 double **d,
                 double g[6][6])
{
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("c1_call_matd");
#endif
/*------------------------------------------------ call material law ---*/
  switch(mat->mattyp)
  {
  case m_stvenant:/*------------------------------- linear elastic ---*/
    dserror(" material law not implemented for optimization");
  break;
  case m_stvenpor:/*------------------------ porous linear elastic ---*/
    c1_matd_stvpor(mat, ele->e.c1->elewa->matdata, d);
    c1gld (d,g);/* transform local to global material matrix */
    c1mefm(strain, d, stress);
  break;
  case m_mfoc:/*-------------- open cell metal foam linear elastic ---*/
    c1_matd_mfoc(mat, ele->e.c1->elewa->matdata, d);
    c1gld (d,g);/* transform local to global material matrix */
    c1mefm(strain, d, stress);
  break;
  case m_mfcc:/*------------ closed cell metal foam linear elastic ---*/
    c1_matd_mfcc(mat, ele->e.c1->elewa->matdata, d);
    c1gld (d,g);/* transform local to global material matrix */
    c1mefm(strain, d, stress);
  break;
  default:
    dserror(" material law not implemented for optimization");
  break;
  }
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of c1_call_matd */
/*!----------------------------------------------------------------------
\brief get density out of material law

<pre>                                                              al 06/02
This routine gives density out of material law.

</pre>
\param      mat MATERIAL*   (i)   material data
\param    density double*   (i)   density value 

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: opt_c1_cint()

*----------------------------------------------------------------------*/
void c1_getdensity(MATERIAL   *mat, double *density)
{
#ifdef DEBUG 
dstrc_enter("c1_getdensity");
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
case m_stvenpor:/*------------------------ porous linear elastic ---*/
   *density = mat->m.stvenpor->density;
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
} /* end of c1_getdensity */
/*----------------------------------------------------------------------*/
#endif
/*! @} (documentation module close)*/

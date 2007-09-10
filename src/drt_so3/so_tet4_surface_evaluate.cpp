/*!----------------------------------------------------------------------*###
\file so_tet4_surface.cpp
\brief

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
written by : Alexander Volf
			alexander.volf@mytum.de  
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOTET4
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "so_tet4.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_timecurve.H"

extern "C"
{
#include "../headers/standardtypes.h"
}
#include "../drt_lib/dstrc.H"

/*----------------------------------------------------------------------*
 * Integrate a Surface Neumann boundary condition (public)     vlf 08/07*
 * ---------------------------------------------------------------------*/
int DRT::Elements::Sotet4Surface::EvaluateNeumann(ParameterList&           params,
                                                DRT::Discretization&     discretization,
                                                DRT::Condition&          condition,
                                                vector<int>&             lm,
                                                Epetra_SerialDenseVector& elevec1)
{
  DSTraceHelper dst("Sotet4Surface::EvaluateNeumann");

  // get values and switches from the condition
  Epetra_SerialDenseMatrix* shapefct;
  Epetra_SerialDenseVector* weights;  //[NUMGPT_SOTET4_FACE]
  sotet4_surface_shapefunc(&shapefct,&weights);
  const vector<int>*    onoff = condition.Get<vector<int> >   ("onoff");
  const vector<double>* val   = condition.Get<vector<double> >("val"  );

  /*
  **    TIME CURVE BUSINESS
  */
  // find out whether we will use a time curve
  bool usetime = true;
  const double time = params.get("total time",-1.0);
  if (time<0.0) usetime = false;

  // find out whether we will use a time curve and get the factor
  const vector<int>* curve  = condition.Get<vector<int> >("curve");
  int curvenum = -1;
  if (curve) curvenum = (*curve)[0];
  double curvefac = 1.0;
  if (curvenum>=0 && usetime)
    curvefac = DRT::Utils::TimeCurveManager::Instance().Curve(curvenum).f(time);

  // element geometry
  const int numnod = 6;
  Epetra_SerialDenseMatrix xsrefe(numnod,NUMDIM_SOTET4);  // material coord. of element
  for (int i=0; i<numnod; i++){
    xsrefe(i,0) = Nodes()[i]->X()[0];
    xsrefe(i,1) = Nodes()[i]->X()[1];
    xsrefe(i,2) = Nodes()[i]->X()[2];
  }

  Epetra_SerialDenseVector A(NUMDIM_SOTET4);
  Epetra_SerialDenseVector B(NUMDIM_SOTET4);
  Epetra_SerialDenseVector C(NUMDIM_SOTET4);
  
  /*
   * to compute the Jacobian i compute the area of the triangle 
   * for that i get the boundary vectors A,B of the triangle 
   * the area is then |A x B| 
   */
  
  A(0)=xsrefe(1,0)-xsrefe(0,0);
  A(1)=xsrefe(1,1)-xsrefe(0,1);
  A(2)=xsrefe(1,2)-xsrefe(0,2);

  
  B(0)=xsrefe(2,0)-xsrefe(0,0);
  B(1)=xsrefe(2,1)-xsrefe(0,1);
  B(2)=xsrefe(2,2)-xsrefe(0,2);
  
  /* C = A x B  */
  C(0)=A(0)*B(1) - A(1)*B(0);
  C(1)=A(1)*B(2) - A(2)*B(1);
  C(2)=A(2)*B(0) - A(0)*B(2);
  
  /* detJ = |A x B|*/
  double detJ= C.Norm2()/2;

  /*
  ** Here, we integrate a 6-node surface with 3 Gauss Points
  */
  double fac = (*weights)(0) * detJ * curvefac;   // integration factor

  // gauss parameters
  for (int gpid = 0; gpid < NUMGPT_SOTET4_FACE; gpid++) {    // loop over intergration points
    // get shape functions and derivatives of element surface
    // distribute over element load vector
    for (int nodid=0; nodid < NUMNOD_SOTET4_FACE; nodid++) {
      for(int dim=0; dim < NUMDIM_SOTET4; dim++) {
        elevec1[nodid*NUMDIM_SOTET4 + dim] +=\
        	(*shapefct)(nodid,gpid) * (*onoff)[dim] * (*val)[dim] * fac;
      }
    }
  }
 
  return 0;
} //Sotet4Surface::EvaluateNeumann(..)


/*----------------------------------------------------------------------*
 * Get shape functions for a tet4 face					       vlf 08/07*
 * ---------------------------------------------------------------------*/
void DRT::Elements::Sotet4Surface::sotet4_surface_shapefunc(
      Epetra_SerialDenseMatrix** shapefct,  // pointer to pointer of shapefct
      Epetra_SerialDenseVector** weights)   // pointer to pointer of weights
{
  DSTraceHelper dst("Sotet4Surface::sotet4_surface_shapefunc");

  static Epetra_SerialDenseMatrix  f(NUMNOD_SOTET4_FACE,NUMGPT_SOTET4_FACE);  // shape functions
  static Epetra_SerialDenseVector weightfactors(NUMGPT_SOTET4_FACE);   // weights for each gp

 //Quadrature rule from Carlos A. Felippa: Adv. FEM  §17 
  const double gploc_alpha    = (double)0.5;    // gp sampling point value for quadr. fct
  const double w			  = (double)1;

  const double ksi1[NUMGPT_SOTET4_FACE] = {gploc_alpha };
  const double ksi2[NUMGPT_SOTET4_FACE] = {gploc_alpha };
  const double ksi3[NUMGPT_SOTET4_FACE] = {gploc_alpha };
  
  for (int i=0; i<NUMGPT_SOTET4_FACE; i++) {
      f(0,i) = ksi1[i];
      f(1,i) = ksi2[i];
      f(2,i) = ksi3[i];
      weightfactors[i] = w; // just for clarity how to get weight factors
   } 
   
   *weights  = &weightfactors;
   *shapefct = &f;

   return;

}


#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif // #ifdef D_SOTET4

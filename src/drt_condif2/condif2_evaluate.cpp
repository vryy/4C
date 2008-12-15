/*!----------------------------------------------------------------------
\file condif2_evaluate.cpp
\brief

<pre>
Maintainer: Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15245
</pre>

*----------------------------------------------------------------------*/
#ifdef D_FLUID2
#ifdef CCADISCRET

// This is just here to get the c++ mpi header, otherwise it would
// use the c version included inside standardtypes.h
#ifdef PARALLEL
#include "mpi.h"
#endif

#include "condif2.H"
#include "../drt_condif3/condif3_impl.H"
#include "../drt_mat/convecdiffus.H"
#include "../drt_mat/matlist.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H" // for CalculateFlux()
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_timecurve.H"
#include <Epetra_SerialDenseSolver.h>

#include "../drt_lib/drt_globalproblem.H"


using namespace DRT::UTILS;

/*---------------------------------------------------------------------*
|  converts a string into an action for this element(private) gjb 06/08|
 *---------------------------------------------------------------------*/
DRT::ELEMENTS::Condif2::ActionType DRT::ELEMENTS::Condif2::convertStringToActionType(
  const string& action) const
{
  dsassert(action != "none", "No action supplied");

  DRT::ELEMENTS::Condif2::ActionType act = Condif2::none;
  if (action == "calc_condif_systemmat_and_residual")
    act = Condif2::calc_condif_systemmat_and_residual;
  else if (action == "calc_initial_time_deriv")
    act = Condif2::calc_initial_time_deriv;
  else if (action == "calc_subgrid_diffusivity_matrix")
    act = Condif2::calc_subgrid_diffusivity_matrix;
  else if (action == "calc_condif_flux")
    act = Condif2::calc_condif_flux;
  else if (action == "calc_temp_and_dens")
    act = Condif2::calc_temp_and_dens;
  else
    dserror("Unknown type of action for Condif2");
  return act;
}


/*----------------------------------------------------------------------*
 |  evaluate the element (public)                               vg 05/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Condif2::Evaluate(
    ParameterList&            params,
    DRT::Discretization&      discretization,
    vector<int>&              lm,
    Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseMatrix& elemat2,
    Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2,
    Epetra_SerialDenseVector& elevec3)
{

  // all physics-related stuff is included in the implementation class that can
  // be used in principle inside any element (at the moment: condif3 and condif2)
  // If this element has special features/ methods that do not fit in the
  // generalized implementation class, you have to do a switch here in order to
  // call element-specific routines
  return DRT::ELEMENTS::Condif3ImplInterface::Impl(this)->Evaluate(
      this,
      params,
      discretization,
      lm,
      elemat1,
      elemat2,
      elevec1,
      elevec2,
      elevec3
  );

} // end of DRT::ELEMENTS::Condif2::Evaluate


/*----------------------------------------------------------------------*
 |  do nothing (public)                                         vg 08/07|
 |                                                                      |
 |  The function is just a dummy. For the condif2 elements, the         |
 |  integration of the surface neumann loads takes place in the element.|
 |  We need it there for the stabilization terms!                       |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Condif2::EvaluateNeumann(ParameterList&            params,
                                            DRT::Discretization&      discretization,
                                            DRT::Condition&           condition,
                                            vector<int>&              lm,
                                            Epetra_SerialDenseVector& elevec1)
{
  return 0;
}


/*----------------------------------------------------------------------*
 |  calculate mass flux                               (private) vg 09/08|
 *----------------------------------------------------------------------*/
Epetra_SerialDenseMatrix DRT::ELEMENTS::Condif2::CalculateFlux(
    vector<double>&           ephinp,
    struct _MATERIAL*         material,
    bool                      temperature,
    Epetra_SerialDenseVector& evel,
    Condif2::FluxType         fluxtype,
    const int&                dofindex
)
{
  /*------------------------------------------------- set element data */
  const int iel = NumNode();
  const DiscretizationType distype = Shape();
  const int nsd = 2;

  Epetra_SerialDenseMatrix xyze(nsd,iel);
  Epetra_SerialDenseMatrix flux(nsd,iel);

  // get node coordinates
  for (int i=0;i<iel;i++)
  {
    xyze(0,i)=Nodes()[i]->X()[0];
    xyze(1,i)=Nodes()[i]->X()[1];
  }

  // get diffusivity
  double diffus = 0;

  if (material->mattyp == m_matlist)
  {
    const int matid = material->m.matlist->matids[dofindex];
    const _MATERIAL& singlemat =  DRT::Problem::Instance()->Material(matid-1);

    if (singlemat.mattyp == m_condif)
      diffus = singlemat.m.condif->diffusivity;
    else if (singlemat.mattyp == m_ion)
      diffus = singlemat.m.ion->diffusivity;
    else
      dserror("type of material found in material list is not supported.");
  }
  else if (material->mattyp == m_condif)
  {
    dsassert(numdofpernode_==1,"more than 1 dof per node for condif material"); // paranoia?
    if (temperature)
      diffus = material->m.condif->diffusivity/material->m.condif->shc;
    else
      diffus = material->m.condif->diffusivity;
  }
  else
    dserror("Material type is not supported");

  /*----------------------------------------- declaration of variables ---*/
  Epetra_SerialDenseVector        funct(iel);
  Epetra_SerialDenseMatrix        deriv(nsd,iel);
  static Epetra_SerialDenseMatrix xjm(nsd,nsd);
  Epetra_SerialDenseMatrix        derxy(nsd,iel);

  vector< LINALG::Matrix<3,1> > nodecoords;
  nodecoords = DRT::UTILS::getEleNodeNumbering_nodes_reference(distype);

  if ((int) nodecoords.size() != iel) dserror("number of nodes does not match");

  // loop over all nodes
  for (int iquad=0; iquad<iel; ++iquad)
  {
    // coordiantes of the current integration point
    const double e1 = nodecoords[iquad](0);
    const double e2 = nodecoords[iquad](1);

    // shape functions and their derivatives
    DRT::UTILS::shape_function_2D(funct,e1,e2,distype);
    DRT::UTILS::shape_function_2D_deriv1(deriv,e1,e2,distype);

    /*----------------------------------------- compute Jacobian matrix */
    double dum;
    /*-------------------------------- determine jacobian at point r,s ---*/
    for (int i=0; i<nsd; i++)
    {
      for (int j=0; j<nsd; j++)
      {
        dum=0.0;
        for (int l=0; l<iel; l++)
        {
          dum += deriv(i,l)*xyze(j,l);
        }
        xjm(i,j)=dum;
      } /* end of loop j */
    } /* end of loop i */

    // The determinant is computed using Sarrus's rule:
    const double det = xjm(0,0)*xjm(1,1)-xjm(0,1)*xjm(1,0);

    if (det < 0.0)
      dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f", Id(), det);
    if (abs(det) < 1E-16)
      dserror("GLOBAL ELEMENT NO.%i\nZERO JACOBIAN DETERMINANT: %f", Id(), det);

    /*------------------------------------------------------------------*/
    /*                                         compute global derivates */
    /*------------------------------------------------------------------*/

    /*------------------------------------------------- initialization */
    for(int k=0;k<iel;k++)
    {
      derxy(0,k)=0.0;
      derxy(1,k)=0.0;
    } /* end of loop over k */

    // ---------------------------------------inverse of transposed jacobian
    static Epetra_SerialDenseMatrix       xij(nsd,nsd);
    xij(0,0) =  xjm(1,1)/det;
    xij(1,0) = -xjm(1,0)/det;
    xij(0,1) = -xjm(0,1)/det;
    xij(1,1) =  xjm(0,0)/det;

    // ---------------------------------------- calculate global derivatives
    for (int k=0;k<iel;k++)
    {
      derxy(0,k) +=  xij(0,0) * deriv(0,k) + xij(0,1) * deriv(1,k) ;
      derxy(1,k) +=  xij(1,0) * deriv(0,k) + xij(1,1) * deriv(1,k) ;
    }

    // add different flux contributions as specified by user input
    switch (fluxtype)
    {
      case Condif2::totalflux:
        //convective flux terms
        flux(0,iquad)+=evel[iquad*nsd]*ephinp[iquad*numdofpernode_+dofindex];
        flux(1,iquad)+=evel[1+iquad*nsd]*ephinp[iquad*numdofpernode_+dofindex];
        // no break statement here!
      case Condif2::diffusiveflux:
        //diffusive flux terms
        for (int k=0;k<iel;k++)
        {
          flux(0,iquad)+=-diffus*derxy(0,k)*ephinp[k*numdofpernode_+dofindex];
          flux(1,iquad)+=-diffus*derxy(1,k)*ephinp[k*numdofpernode_+dofindex];
        }
        break;
      case Condif2::noflux:
        dserror("received noflux flag inside CONDIF2 flux evaluation");
    };

  } // loop over nodes

  return flux;
} // Condif2::CalculateFlux


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID2

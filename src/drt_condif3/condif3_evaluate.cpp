/*!
\file condif3_evaluate.cpp
\brief

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>

*/
#ifdef D_FLUID3
#ifdef CCADISCRET

// This is just here to get the c++ mpi header, otherwise it would
// use the c version included inside standardtypes.h
#ifdef PARALLEL
#include "mpi.h"
#endif

#include "condif3.H"
#include "condif3_utils.H"
#include "condif3_impl.H"
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


/*----------------------------------------------------------------------*
 |  evaluate the element (public)                              gjb 06/08|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Condif3::Evaluate(ParameterList& params,
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

} //DRT::ELEMENTS::Condif3::Evaluate


/*----------------------------------------------------------------------*
 |  do nothing (public)                                        gjb 06/08|
 |                                                                      |
 |  The function is just a dummy. For the condif elements, the          |
 |  integration of the volume neumann (body forces) loads takes place   |
 |  in the element. We need it there for the stabilisation terms!       |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Condif3::EvaluateNeumann(ParameterList& params,
    DRT::Discretization&      discretization,
    DRT::Condition&           condition,
    vector<int>&              lm,
    Epetra_SerialDenseVector& elevec1)
{
  return 0;
}


/*----------------------------------------------------------------------*
 |  calculate mass flux                              (private) gjb 06/08|
 *----------------------------------------------------------------------*/
Epetra_SerialDenseMatrix DRT::ELEMENTS::Condif3::CalculateFlux(
    vector<double>&           ephinp,
    struct _MATERIAL*         material,
    const bool                temperature,
    const double              frt,
    Epetra_SerialDenseVector& evel,
    Condif3::FluxType         fluxtype,
    const int&                dofindex
)
{
  /*------------------------------------------------- set element data */
  const int iel = NumNode();
  const DiscretizationType distype = Shape();
  const int nsd = 3;

  Epetra_SerialDenseMatrix xyze(nsd,iel);
  Epetra_SerialDenseMatrix flux(nsd,iel);

  // get node coordinates
  for (int i=0;i<iel;i++)
  {
    xyze(0,i)=Nodes()[i]->X()[0];
    xyze(1,i)=Nodes()[i]->X()[1];
    xyze(2,i)=Nodes()[i]->X()[2];
  }

  // get diffusivity
  double diffus(0.0);
  double valence(0.0);
  double diffus_valence_frt(0.0);

  if (material->mattyp == m_matlist)
  {
    const int matid = material->m.matlist->matids[dofindex];
    const _MATERIAL& singlemat =  DRT::Problem::Instance()->Material(matid-1);

    if (singlemat.mattyp == m_condif)
      diffus = singlemat.m.condif->diffusivity;
    else if (singlemat.mattyp == m_ion)
    {
      diffus = singlemat.m.ion->diffusivity;
      valence = singlemat.m.ion->valence;
      diffus_valence_frt = diffus*valence*frt;
    }
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
  Epetra_SerialDenseMatrix        xjm(nsd,nsd);
  Epetra_SerialDenseMatrix        derxy(nsd,iel);

  LINALG::SerialDenseMatrix nodecoords;
  nodecoords = DRT::UTILS::getEleNodeNumbering_nodes_paramspace(distype);

  if ((int) nodecoords.N() != iel) dserror("number of nodes does not match");

  // loop over all nodes
  for (int iquad=0; iquad<iel; ++iquad)
  {
    // coordinates of the current integration point
    const double e1 = nodecoords(0, iquad);
    const double e2 = nodecoords(1, iquad);
    const double e3 = nodecoords(2, iquad);

    // shape functions and their derivatives
    DRT::UTILS::shape_function_3D(funct,e1,e2,e3,distype);
    DRT::UTILS::shape_function_3D_deriv1(deriv,e1,e2,e3,distype);

    /*----------------------------------------- compute Jacobian matrix */

    // get Jacobian matrix and determinant
    // actually compute its transpose....
    /*
    +-            -+ T      +-            -+
    | dx   dx   dx |        | dx   dy   dz |
    | --   --   -- |        | --   --   -- |
    | dr   ds   dt |        | dr   dr   dr |
    |              |        |              |
    | dy   dy   dy |        | dx   dy   dz |
    | --   --   -- |   =    | --   --   -- |
    | dr   ds   dt |        | ds   ds   ds |
    |              |        |              |
    | dz   dz   dz |        | dx   dy   dz |
    | --   --   -- |        | --   --   -- |
    | dr   ds   dt |        | dt   dt   dt |
    +-            -+        +-            -+
     */
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
    // ---------------------------------------- calculate determinant
    const double det = xjm(0,0)*xjm(1,1)*xjm(2,2)+
    xjm(0,1)*xjm(1,2)*xjm(2,0)+
    xjm(0,2)*xjm(1,0)*xjm(2,1)-
    xjm(0,2)*xjm(1,1)*xjm(2,0)-
    xjm(0,0)*xjm(1,2)*xjm(2,1)-
    xjm(0,1)*xjm(1,0)*xjm(2,2);

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
      derxy(2,k)=0.0;
    } /* end of loop over k */

    // ---------------------------------------inverse of transposed jacobian
    Epetra_SerialDenseMatrix       xij(nsd,nsd);
    double idet = 1./det;
    xij(0,0) = (  xjm(1,1)*xjm(2,2) - xjm(2,1)*xjm(1,2))*idet;
    xij(1,0) = (- xjm(1,0)*xjm(2,2) + xjm(2,0)*xjm(1,2))*idet;
    xij(2,0) = (  xjm(1,0)*xjm(2,1) - xjm(2,0)*xjm(1,1))*idet;
    xij(0,1) = (- xjm(0,1)*xjm(2,2) + xjm(2,1)*xjm(0,2))*idet;
    xij(1,1) = (  xjm(0,0)*xjm(2,2) - xjm(2,0)*xjm(0,2))*idet;
    xij(2,1) = (- xjm(0,0)*xjm(2,1) + xjm(2,0)*xjm(0,1))*idet;
    xij(0,2) = (  xjm(0,1)*xjm(1,2) - xjm(1,1)*xjm(0,2))*idet;
    xij(1,2) = (- xjm(0,0)*xjm(1,2) + xjm(1,0)*xjm(0,2))*idet;
    xij(2,2) = (  xjm(0,0)*xjm(1,1) - xjm(1,0)*xjm(0,1))*idet;
    /*---------------------------------------- calculate global derivatives */
    for (int k=0;k<iel;k++)
    {
      derxy(0,k) +=   xij(0,0) * deriv(0,k) + xij(0,1) * deriv(1,k) + xij(0,2) * deriv(2,k);
      derxy(1,k) +=   xij(1,0) * deriv(0,k) + xij(1,1) * deriv(1,k) + xij(1,2) * deriv(2,k);
      derxy(2,k) +=   xij(2,0) * deriv(0,k) + xij(2,1) * deriv(1,k) + xij(2,2) * deriv(2,k);
    } /* end of loop over k */

    // gradient of electric potential
    Epetra_SerialDenseVector gradpot(3);
    if (frt > 0.0) // ELCH
    {
      for (int k=0;k<iel;k++)
      {
        gradpot[0] += derxy(0,k)*ephinp[k*numdofpernode_+(numdofpernode_-1)];
        gradpot[1] += derxy(1,k)*ephinp[k*numdofpernode_+(numdofpernode_-1)];
        gradpot[2] += derxy(2,k)*ephinp[k*numdofpernode_+(numdofpernode_-1)];
      } /* end of loop over k */
    }

    const double ephinpatnode = ephinp[iquad*numdofpernode_+dofindex];
    // add different flux contributions as specified by user input
    switch (fluxtype)
    {
      case Condif3::totalflux:
        if (frt > 0.0) // ELCH
        {
          //migration flux terms
          flux(0,iquad)-=diffus_valence_frt*gradpot[0]*ephinpatnode;
          flux(1,iquad)-=diffus_valence_frt*gradpot[1]*ephinpatnode;
          flux(2,iquad)-=diffus_valence_frt*gradpot[2]*ephinpatnode;
        }
        //convective flux terms
        flux(0,iquad)+=evel[iquad*nsd]*ephinpatnode;
        flux(1,iquad)+=evel[1+iquad*nsd]*ephinpatnode;
        flux(2,iquad)+=evel[2+iquad*nsd]*ephinpatnode;
        // no break statement here!
      case Condif3::diffusiveflux:
        //diffusive flux terms
        for (int k=0;k<iel;k++)
        {
          flux(0,iquad)+=-diffus*derxy(0,k)*ephinp[k*numdofpernode_+dofindex];
          flux(1,iquad)+=-diffus*derxy(1,k)*ephinp[k*numdofpernode_+dofindex];
          flux(2,iquad)+=-diffus*derxy(2,k)*ephinp[k*numdofpernode_+dofindex];
        }
        break;
      case Condif3::noflux:
        dserror("received noflux flag inside CONDIF3 flux evaluation");
    };

  } // loop over nodes

  return flux;
} // Condif3::CalculateFlux



#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID3

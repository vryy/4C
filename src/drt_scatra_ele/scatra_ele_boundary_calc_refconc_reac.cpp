/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_boundary_calc_refconc_reac.cpp

\brief main file containing routines for calculation of scatra element formulated in reference concentrations
  and with advanced reaction terms

\level 3

 <pre>
   \maintainer Moritz Thon
               thon@mhpc.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-10364
 </pre>
 *----------------------------------------------------------------------*/

#include "scatra_ele_boundary_calc_refconc_reac.H"
#include "scatra_ele_parameter_std.H"
#include "../drt_geometry/position_array.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_fem_general/drt_utils_boundary_integration.H"


/*----------------------------------------------------------------------*
 |  Singleton access method                                  thon 02/16 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleBoundaryCalcRefConcReac<distype> * DRT::ELEMENTS::ScaTraEleBoundaryCalcRefConcReac<distype>::Instance(
    const int numdofpernode,
    const int numscal,
    const std::string& disname,
    const ScaTraEleBoundaryCalcRefConcReac* delete_me
    )
{
  static std::map<std::string,ScaTraEleBoundaryCalcRefConcReac<distype>* >  instances;

  if(delete_me == NULL)
  {
    if(instances.find(disname) == instances.end())
      instances[disname] = new ScaTraEleBoundaryCalcRefConcReac<distype>(numdofpernode,numscal,disname);
  }

  else
  {
    for( typename std::map<std::string,ScaTraEleBoundaryCalcRefConcReac<distype>* >::iterator i=instances.begin(); i!=instances.end(); ++i )
      if ( i->second == delete_me )
      {
        delete i->second;
        instances.erase(i);
        return NULL;
      }
  }

  return instances[disname];
}


/*----------------------------------------------------------------------*
 |  Clean up                                                 thon 02/16 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcRefConcReac<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
    Instance( 0, 0, "", this );
}


/*----------------------------------------------------------------------*
 |  Private constructor                                      thon 02/16 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleBoundaryCalcRefConcReac<distype>::ScaTraEleBoundaryCalcRefConcReac(const int numdofpernode, const int numscal,const std::string& disname)
  : DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::ScaTraEleBoundaryCalc(numdofpernode,numscal,disname)
{
  return;
}


/*---------------------------------------------------------------------------*
 | Factor needed for the calculation of reference concentrations  thon 02/16 |
 *---------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::ScaTraEleBoundaryCalcRefConcReac<distype>::FacForRefConc(
    const int                      iquad,             ///< current boundary integration point
    const DRT::FaceElement*        bele,              ///< current boundary element
    Teuchos::ParameterList&        params,            ///< parameter list
    DRT::Discretization&           discretization     ///< discretization
    )
{
  const DRT::Element* pele = bele->ParentElement();

  double J=1.0;
  // only 3D cases:
  if (bele->Shape()==DRT::Element::tri3)
  {
    if (pele->Shape()==DRT::Element::tet4)
      J = CalcJatIntPoint<DRT::Element::tri3,DRT::Element::tet4>(iquad,bele,pele,params,discretization);
    else if (pele->Shape()==DRT::Element::pyramid5)
      J = CalcJatIntPoint<DRT::Element::tri3,DRT::Element::pyramid5>(iquad,bele,pele,params,discretization);
    else
      dserror("Parent element not supported here!");

  }
  else if (bele->Shape()==DRT::Element::quad4)
  {
    if (pele->Shape()==DRT::Element::hex8)
      J = CalcJatIntPoint<DRT::Element::quad4,DRT::Element::hex8>(iquad,bele,pele,params,discretization);
    else if (pele->Shape()==DRT::Element::pyramid5)
      J = CalcJatIntPoint<DRT::Element::quad4,DRT::Element::pyramid5>(iquad,bele,pele,params,discretization);
    else
      dserror("Parent element not supported here!");
  }
  else
    dserror("Boundary element not supported here!");

  return 1.0/J;
}


/*---------------------------------------------------------------------------*
 | Factor needed for the calculation of reference concentrations  thon 02/16 |
 *---------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
template <DRT::Element::DiscretizationType bdistype, DRT::Element::DiscretizationType pdistype>
double DRT::ELEMENTS::ScaTraEleBoundaryCalcRefConcReac<distype>::CalcJatIntPoint(
    const int                      iquad,             ///< current boundary integration point
    const DRT::FaceElement*        bele,              ///< current boundary element
    const DRT::Element*            pele,              ///< current parent element
    Teuchos::ParameterList&        params,            ///< parameter list
    DRT::Discretization&           discretization     ///< discretization
    )
{
  // NOTE: we want to evaluate J=det(F) on the current gauß point of the current boundary element.
  // Since this does depend on ALL values of the involved element this is quite a hassle :(

  // number of parent spatial dimensions
  const int pnsd = DRT::UTILS::DisTypeToDim<pdistype>::dim;
  // number of boundary spatial dimensions
  const int bnsd = DRT::UTILS::DisTypeToDim<bdistype>::dim;

  if (pnsd != (my::nsd_+1))
    dserror("dimension do not match!");
  if (bnsd != my::nsd_)
    dserror("dimension do not match!");

  // number of parent element nodes
  const int pnen = DRT::UTILS::DisTypeToNumNodePerEle<pdistype>::numNodePerElement;
  // number of (boundary) element nodes
  static const int bnen = DRT::UTILS::DisTypeToNumNodePerEle<bdistype>::numNodePerElement;

  if (bnen != my::nen_)
    dserror("Number of element nodes do not match!");

  // get local node coordinates
  LINALG::Matrix<pnsd,pnen> pxyze(true);
  LINALG::Matrix<pnsd,pnen> pxyze0(true);
    GEO::fillInitialPositionArray<pdistype,pnsd,LINALG::Matrix<pnsd,pnen> >(pele,pxyze0);
  pxyze=pxyze0;

  if (my::scatraparams_->IsAle())
  {
    // get number of dof-set associated with displacement related dofs
    const int ndsdisp = params.get<int>("ndsdisp");

    Teuchos::RCP<const Epetra_Vector> dispnp = discretization.GetState(ndsdisp, "dispnp");
    if (dispnp==Teuchos::null)
      dserror("Cannot get state vector 'dispnp'");

    // parent element location array
    DRT::Element::LocationArray pla(discretization.NumDofSets());
    pele->LocationVector(discretization,pla,false);

    // determine number of velocity related dofs per node
    const int numdispdofpernode = pla[ndsdisp].lm_.size()/pnen;

    // construct location vector for velocity related dofs
    std::vector<int> plmdisp(pnsd*pnen,-1);
    for (int inode=0; inode<pnen; ++inode)
      for (int idim=0; idim<pnsd; ++idim)
        plmdisp[inode*pnsd+idim] = pla[ndsdisp].lm_[inode*numdispdofpernode+idim];

    // we deal with a (nsd_+1)-dimensional flow field
    LINALG::Matrix<pnsd,pnen>  pedispnp(true);

    // extract local values of convective velocity field from global state vector
    DRT::UTILS::ExtractMyValues<LINALG::Matrix<pnsd,pnen> >(*dispnp,pedispnp,plmdisp);

    // rotate the vector field in the case of rotationally symmetric boundary conditions
    //my::rotsymmpbc_->template RotateMyValuesIfNecessary<pnsd,pnen>(pedispnp);

    pxyze += pedispnp;
  }

  // get Gaussian integration points
  const DRT::UTILS::IntPointsAndWeights<pnsd>
      pintpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<pdistype>::rule);

  // get Gaussian integration points
  const DRT::UTILS::IntPointsAndWeights<bnsd>
      bintpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<bdistype>::rule);

  Epetra_SerialDenseMatrix gps(bintpoints.IP().nquad,bnsd);
  for (int biquad=0; biquad<bintpoints.IP().nquad; ++biquad)
  {
    const double* gpcoord = (bintpoints.IP().qxg)[biquad];
    for (int idim=0;idim<bnsd ;idim++)
    {
      gps(biquad,idim) = gpcoord[idim];
    }
  }

  // distinguish 2- and 3-D case
  Epetra_SerialDenseMatrix pqxg(pintpoints.IP().nquad,pnsd);
  if (pnsd==2)
     DRT::UTILS::BoundaryGPToParentGP2(pqxg,gps,pdistype,bdistype,bele->FaceMasterNumber());
  else if(pnsd==3)
    DRT::UTILS::BoundaryGPToParentGP3(pqxg,gps,pdistype,bdistype,bele->FaceMasterNumber());


  LINALG::Matrix<pnsd,1>  pxsi(true);
  LINALG::Matrix<pnsd,pnen>  pderiv(true);

  // reference coordinates of integration point from parent element
  for (int idim=0;idim<pnsd;idim++)
  {
    pxsi(idim) = pqxg(iquad,idim);
  }

  // parent element shape functions and local derivatives
  DRT::UTILS::shape_function_deriv1<pdistype>(pxsi,pderiv);

  // Jacobian matrix and determinant of parent element (including check)
  LINALG::Matrix<pnsd,pnsd>  dxds(true);
  dxds.MultiplyNT(pderiv,pxyze);
  const double detdxds = dxds.Determinant();

  // Jacobian matrix and determinant of parent element (including check)
  LINALG::Matrix<pnsd,pnsd>  dXds(true);
  dXds.MultiplyNT(pderiv,pxyze0);
  const double detdXds = dXds.Determinant();

  // deformation gradtient dx/dX = dx/ds * ds/dX = dx/ds * (dX/ds)^(-1)
  const double J=detdxds/detdXds;

  return J;
}


// template classes
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcRefConcReac<DRT::Element::quad4>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcRefConcReac<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcRefConcReac<DRT::Element::quad9>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcRefConcReac<DRT::Element::tri3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcRefConcReac<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcRefConcReac<DRT::Element::line2>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcRefConcReac<DRT::Element::line3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcRefConcReac<DRT::Element::nurbs3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcRefConcReac<DRT::Element::nurbs9>;

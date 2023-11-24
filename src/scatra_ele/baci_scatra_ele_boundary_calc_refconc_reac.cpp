/*----------------------------------------------------------------------*/
/*! \file
\brief main file containing routines for calculation of scatra element formulated in reference
concentrations and with advanced reaction terms

\level 3

 *----------------------------------------------------------------------*/

#include "baci_scatra_ele_boundary_calc_refconc_reac.H"

#include "baci_discretization_fem_general_utils_boundary_integration.H"
#include "baci_discretization_geometry_position_array.H"
#include "baci_lib_utils.H"
#include "baci_scatra_ele_parameter_std.H"
#include "baci_utils_singleton_owner.H"


/*----------------------------------------------------------------------*
 |  Singleton access method                                  thon 02/16 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
DRT::ELEMENTS::ScaTraEleBoundaryCalcRefConcReac<distype, probdim>*
DRT::ELEMENTS::ScaTraEleBoundaryCalcRefConcReac<distype, probdim>::Instance(
    const int numdofpernode, const int numscal, const std::string& disname)
{
  static auto singleton_map = CORE::UTILS::MakeSingletonMap<std::string>(
      [](const int numdofpernode, const int numscal, const std::string& disname)
      {
        return std::unique_ptr<ScaTraEleBoundaryCalcRefConcReac<distype, probdim>>(
            new ScaTraEleBoundaryCalcRefConcReac<distype, probdim>(
                numdofpernode, numscal, disname));
      });

  return singleton_map[disname].Instance(
      CORE::UTILS::SingletonAction::create, numdofpernode, numscal, disname);
}


/*----------------------------------------------------------------------*
 |  Private constructor                                      thon 02/16 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
DRT::ELEMENTS::ScaTraEleBoundaryCalcRefConcReac<distype, probdim>::ScaTraEleBoundaryCalcRefConcReac(
    const int numdofpernode, const int numscal, const std::string& disname)
    : DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype, probdim>::ScaTraEleBoundaryCalc(
          numdofpernode, numscal, disname)
{
  return;
}


/*---------------------------------------------------------------------------*
 | Factor needed for the calculation of reference concentrations  thon 02/16 |
 *---------------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
double DRT::ELEMENTS::ScaTraEleBoundaryCalcRefConcReac<distype, probdim>::FacForRefConc(
    const int iquad,                     ///< current boundary integration point
    const DRT::FaceElement* bele,        ///< current boundary element
    Teuchos::ParameterList& params,      ///< parameter list
    DRT::Discretization& discretization  ///< discretization
)
{
  const DRT::Element* pele = bele->ParentElement();

  double J = 1.0;
  // only 3D cases:
  if (bele->Shape() == CORE::FE::CellType::tri3)
  {
    if (pele->Shape() == CORE::FE::CellType::tet4)
      J = CalcJatIntPoint<CORE::FE::CellType::tri3, CORE::FE::CellType::tet4>(
          iquad, bele, pele, params, discretization);
    else if (pele->Shape() == CORE::FE::CellType::pyramid5)
      J = CalcJatIntPoint<CORE::FE::CellType::tri3, CORE::FE::CellType::pyramid5>(
          iquad, bele, pele, params, discretization);
    else
      dserror("Parent element not supported here!");
  }
  else if (bele->Shape() == CORE::FE::CellType::quad4)
  {
    if (pele->Shape() == CORE::FE::CellType::hex8)
      J = CalcJatIntPoint<CORE::FE::CellType::quad4, CORE::FE::CellType::hex8>(
          iquad, bele, pele, params, discretization);
    else if (pele->Shape() == CORE::FE::CellType::pyramid5)
      J = CalcJatIntPoint<CORE::FE::CellType::quad4, CORE::FE::CellType::pyramid5>(
          iquad, bele, pele, params, discretization);
    else
      dserror("Parent element not supported here!");
  }
  else
    dserror("Boundary element not supported here!");

  return 1.0 / J;
}


/*---------------------------------------------------------------------------*
 | Factor needed for the calculation of reference concentrations  thon 02/16 |
 *---------------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
template <CORE::FE::CellType bdistype, CORE::FE::CellType pdistype>
double DRT::ELEMENTS::ScaTraEleBoundaryCalcRefConcReac<distype, probdim>::CalcJatIntPoint(
    const int iquad,                     ///< current boundary integration point
    const DRT::FaceElement* bele,        ///< current boundary element
    const DRT::Element* pele,            ///< current parent element
    Teuchos::ParameterList& params,      ///< parameter list
    DRT::Discretization& discretization  ///< discretization
)
{
  // NOTE: we want to evaluate J=det(F) on the current gauss point of the current boundary element.
  // Since this does depend on ALL values of the involved element this is quite a hassle :(

  // number of parent spatial dimensions
  const int pnsd = CORE::FE::dim<pdistype>;
  // number of boundary spatial dimensions
  const int bnsd = CORE::FE::dim<bdistype>;

  if (pnsd != nsd_) dserror("dimension do not match!");
  if (bnsd != nsd_ele_) dserror("dimension do not match!");

  // number of parent element nodes
  const int pnen = CORE::FE::num_nodes<pdistype>;
  // number of (boundary) element nodes
  static const int bnen = CORE::FE::num_nodes<bdistype>;

  if (bnen != nen_) dserror("Number of element nodes do not match!");

  // get local node coordinates
  CORE::LINALG::Matrix<pnsd, pnen> pxyze(true);
  CORE::LINALG::Matrix<pnsd, pnen> pxyze0(true);
  CORE::GEO::fillInitialPositionArray<pdistype, pnsd, CORE::LINALG::Matrix<pnsd, pnen>>(
      pele, pxyze0);
  pxyze = pxyze0;

  if (my::scatraparams_->IsAle())
  {
    // get number of dof-set associated with displacement related dofs
    const int ndsdisp = my::scatraparams_->NdsDisp();

    Teuchos::RCP<const Epetra_Vector> dispnp = discretization.GetState(ndsdisp, "dispnp");
    if (dispnp == Teuchos::null) dserror("Cannot get state vector 'dispnp'");

    // parent element location array
    DRT::Element::LocationArray pla(discretization.NumDofSets());
    pele->LocationVector(discretization, pla, false);

    // determine number of velocity related dofs per node
    const int numdispdofpernode = pla[ndsdisp].lm_.size() / pnen;

    // construct location vector for velocity related dofs
    std::vector<int> plmdisp(pnsd * pnen, -1);
    for (int inode = 0; inode < pnen; ++inode)
      for (int idim = 0; idim < pnsd; ++idim)
        plmdisp[inode * pnsd + idim] = pla[ndsdisp].lm_[inode * numdispdofpernode + idim];

    // we deal with a nsd_-dimensional flow field
    CORE::LINALG::Matrix<pnsd, pnen> pedispnp(true);

    // extract local values of convective velocity field from global state vector
    DRT::UTILS::ExtractMyValues<CORE::LINALG::Matrix<pnsd, pnen>>(*dispnp, pedispnp, plmdisp);

    // rotate the vector field in the case of rotationally symmetric boundary conditions
    // my::rotsymmpbc_->template RotateMyValuesIfNecessary<pnsd,pnen>(pedispnp);

    pxyze += pedispnp;
  }

  // get Gaussian integration points
  const CORE::DRT::UTILS::IntPointsAndWeights<pnsd> pintpoints(
      DRT::ELEMENTS::DisTypeToOptGaussRule<pdistype>::rule);

  // get Gaussian integration points
  const CORE::DRT::UTILS::IntPointsAndWeights<bnsd> bintpoints(
      DRT::ELEMENTS::DisTypeToOptGaussRule<bdistype>::rule);

  CORE::LINALG::SerialDenseMatrix gps(bintpoints.IP().nquad, bnsd);
  for (int biquad = 0; biquad < bintpoints.IP().nquad; ++biquad)
  {
    const double* gpcoord = (bintpoints.IP().qxg)[biquad];
    for (int idim = 0; idim < bnsd; idim++)
    {
      gps(biquad, idim) = gpcoord[idim];
    }
  }

  // distinguish 2- and 3-D case
  CORE::LINALG::SerialDenseMatrix pqxg(pintpoints.IP().nquad, pnsd);
  if (pnsd == 2)
    CORE::DRT::UTILS::BoundaryGPToParentGP2(
        pqxg, gps, pdistype, bdistype, bele->FaceMasterNumber());
  else if (pnsd == 3)
    CORE::DRT::UTILS::BoundaryGPToParentGP3(
        pqxg, gps, pdistype, bdistype, bele->FaceMasterNumber());


  CORE::LINALG::Matrix<pnsd, 1> pxsi(true);
  CORE::LINALG::Matrix<pnsd, pnen> pderiv(true);

  // reference coordinates of integration point from parent element
  for (int idim = 0; idim < pnsd; idim++)
  {
    pxsi(idim) = pqxg(iquad, idim);
  }

  // parent element shape functions and local derivatives
  CORE::DRT::UTILS::shape_function_deriv1<pdistype>(pxsi, pderiv);

  // Jacobian matrix and determinant of parent element (including check)
  CORE::LINALG::Matrix<pnsd, pnsd> dxds(true);
  dxds.MultiplyNT(pderiv, pxyze);
  const double detdxds = dxds.Determinant();

  // Jacobian matrix and determinant of parent element (including check)
  CORE::LINALG::Matrix<pnsd, pnsd> dXds(true);
  dXds.MultiplyNT(pderiv, pxyze0);
  const double detdXds = dXds.Determinant();

  // deformation gradtient dx/dX = dx/ds * ds/dX = dx/ds * (dX/ds)^(-1)
  const double J = detdxds / detdXds;

  return J;
}


// template classes
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcRefConcReac<CORE::FE::CellType::quad4, 3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcRefConcReac<CORE::FE::CellType::quad8, 3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcRefConcReac<CORE::FE::CellType::quad9, 3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcRefConcReac<CORE::FE::CellType::tri3, 3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcRefConcReac<CORE::FE::CellType::tri6, 3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcRefConcReac<CORE::FE::CellType::line2, 2>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcRefConcReac<CORE::FE::CellType::line2, 3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcRefConcReac<CORE::FE::CellType::line3, 2>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcRefConcReac<CORE::FE::CellType::nurbs3, 2>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcRefConcReac<CORE::FE::CellType::nurbs9, 3>;

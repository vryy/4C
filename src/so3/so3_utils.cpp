/*----------------------------------------------------------------------*/
/*! \file

\brief A collection of helper methods for solid elements

\level 1
 *-----------------------------------------------------------------------*/

#include "so3_utils.H"
#include <algorithm>
#include "linalg_utils_densematrix_svd.H"
#include "lib_element.H"
#include "discretization_fem_general_utils_fem_shapefunctions.H"
#include "fiber_node.H"
#include "so3_prestress.H"

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::UTILS::CalcR(const DRT::Element* ele, const std::vector<double>& disp,
    CORE::LINALG::Matrix<CORE::DRT::UTILS::DisTypeToDim<distype>::dim,
        CORE::DRT::UTILS::DisTypeToDim<distype>::dim>& R)
{
  // number of nodes per element
  const int nen = CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

  // spatial dimension
  const int nsd = CORE::DRT::UTILS::DisTypeToDim<distype>::dim;

  if (disp.size() != nsd * nen) dserror("mismatch in dimensions");

  CORE::LINALG::Matrix<nsd, 1> xi_ele_center =
      CORE::DRT::UTILS::getLocalCenterPosition<nsd>(distype);  // depending on distype

  CORE::LINALG::Matrix<nen, nsd> xrefe;  // X, material coord. of element
  CORE::LINALG::Matrix<nen, nsd> xcurr;  // x, current  coord. of element
  for (int i = 0; i < nen; ++i)
  {
    for (int d = 0; d < nsd; ++d)
    {
      xrefe(i, d) = ele->Nodes()[i]->X()[d];
      xcurr(i, d) = ele->Nodes()[i]->X()[d] + disp[i * nsd + d];
    }
  }
  CORE::LINALG::Matrix<nsd, nen> deriv;
  CORE::DRT::UTILS::shape_function_deriv1<distype>(xi_ele_center, deriv);

  CORE::LINALG::Matrix<nsd, nsd> jac;
  CORE::LINALG::Matrix<nsd, nsd> defgrd;
  CORE::LINALG::Matrix<nsd, nen> deriv_xyz;
  jac.Multiply(deriv, xrefe);
  jac.Invert();
  deriv_xyz.Multiply(jac, deriv);
  defgrd.MultiplyTT(xcurr, deriv_xyz);

  // Calculate rotcurr from defgrd
  CORE::LINALG::Matrix<nsd, nsd> Q(true);
  CORE::LINALG::Matrix<nsd, nsd> S(true);
  CORE::LINALG::Matrix<nsd, nsd> VT(true);
  CORE::LINALG::SVD<nsd, nsd>(defgrd, Q, S, VT);
  R.MultiplyNN(Q, VT);
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::UTILS::GetTemperatureForStructuralMaterial(
    const CORE::LINALG::Matrix<CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement,
        1>& shapefctsGP,            // shape function of current Gauss-point
    Teuchos::ParameterList& params  // special material parameter e.g. scalartemp
)
{
  // initialise the temperature
  Teuchos::RCP<std::vector<double>> temperature_vector =
      params.get<Teuchos::RCP<std::vector<double>>>("nodal_tempnp", Teuchos::null);

  // current temperature vector is available
  if (temperature_vector != Teuchos::null)
  {
    double scalartemp = 0.0;
    for (int i = 0; i < CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement; ++i)
    {
      scalartemp += shapefctsGP(i) * (*temperature_vector)[i];
    }

    // insert current element temperature T_{n+1} into parameter list
    params.set<double>("scalartemp", scalartemp);
  }
}

template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::UTILS::ComputeDeformationGradient(
    CORE::LINALG::Matrix<probdim, probdim>& defgrd, DRT::Node** nodes,
    const CORE::LINALG::Matrix<probdim, 1>& xsi,
    const CORE::LINALG::Matrix<CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement,
        probdim>& xdisp)
{
  static CORE::LINALG::Matrix<CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement,
      probdim>
      xrefe, xcurr;

  EvaluateNodalCoordinates<distype, probdim>(nodes, xrefe);
  EvaluateCurrentNodalCoordinates<distype, probdim>(xrefe, xdisp, xcurr);

  CORE::LINALG::Matrix<probdim,
      CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement>
      N_rst(true);
  CORE::DRT::UTILS::shape_function_deriv1<distype>(xsi, N_rst);

  static CORE::LINALG::Matrix<probdim, probdim> inv_detFJ;
  inv_detFJ.Multiply(N_rst, xrefe);
  inv_detFJ.Invert();

  ComputeDeformationGradientStandard<distype, probdim>(defgrd, xcurr, N_rst, inv_detFJ);
}

template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::UTILS::ComputeDeformationGradient(
    CORE::LINALG::Matrix<probdim, probdim>& defgrd, DRT::Node** nodes,
    const CORE::LINALG::Matrix<probdim, 1>& xsi, const std::vector<double>& displacement)
{
  static CORE::LINALG::Matrix<CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement,
      probdim>
      xdisp;
  EvaluateNodalDisplacements<distype, probdim>(displacement, xdisp);

  ComputeDeformationGradient<distype, probdim>(defgrd, nodes, xsi, xdisp);
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::UTILS::ComputeDeformationGradient(
    CORE::LINALG::Matrix<CORE::DRT::UTILS::DisTypeToDim<distype>::dim,
        CORE::DRT::UTILS::DisTypeToDim<distype>::dim>& defgrd,
    const INPAR::STR::KinemType kinemType,
    const CORE::LINALG::Matrix<CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement,
        CORE::DRT::UTILS::DisTypeToDim<distype>::dim>& xdisp,
    const CORE::LINALG::Matrix<CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement,
        CORE::DRT::UTILS::DisTypeToDim<distype>::dim>& xcurr,
    const CORE::LINALG::Matrix<CORE::DRT::UTILS::DisTypeToDim<distype>::dim,
        CORE::DRT::UTILS::DisTypeToDim<distype>::dim>& inverseJacobian,
    const CORE::LINALG::Matrix<CORE::DRT::UTILS::DisTypeToDim<distype>::dim,
        CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement>& derivs,
    const INPAR::STR::PreStress prestressType,
    const Teuchos::RCP<DRT::ELEMENTS::PreStress> mulfHistory, const int gp)
{
  if (kinemType == INPAR::STR::kinem_linear)
  {
    defgrd.Clear();
    for (auto i = 0; i < CORE::DRT::UTILS::DisTypeToDim<distype>::dim; ++i)
    {
      defgrd(i, i) = 1.0;
    }
    return;
  }

  if (prestressType == INPAR::STR::PreStress::mulf)
  {
    ComputeDeformationGradientMulf<distype>(defgrd, xdisp, derivs, mulfHistory, gp);
    return;
  }

  ComputeDeformationGradientStandard<distype, CORE::DRT::UTILS::DisTypeToDim<distype>::dim>(
      defgrd, xcurr, derivs, inverseJacobian);
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::UTILS::ComputeDeformationGradientMulf(
    CORE::LINALG::Matrix<CORE::DRT::UTILS::DisTypeToDim<distype>::dim,
        CORE::DRT::UTILS::DisTypeToDim<distype>::dim>& defgrd,
    const CORE::LINALG::Matrix<CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement,
        CORE::DRT::UTILS::DisTypeToDim<distype>::dim>& xdisp,
    const CORE::LINALG::Matrix<CORE::DRT::UTILS::DisTypeToDim<distype>::dim,
        CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement>& derivs,
    const Teuchos::RCP<DRT::ELEMENTS::PreStress> mulfHistory, const int gp)
{
  // get Jacobian mapping wrt to the stored configuration
  CORE::LINALG::Matrix<CORE::DRT::UTILS::DisTypeToDim<distype>::dim,
      CORE::DRT::UTILS::DisTypeToDim<distype>::dim>
      invJdef;
  mulfHistory->StoragetoMatrix(gp, invJdef, mulfHistory->JHistory());

  // get derivatives wrt to last spatial configuration
  CORE::LINALG::Matrix<CORE::DRT::UTILS::DisTypeToDim<distype>::dim,
      CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement>
      N_xyz;
  N_xyz.Multiply(invJdef, derivs);

  // build multiplicative incremental defgrd
  CORE::LINALG::Matrix<CORE::DRT::UTILS::DisTypeToDim<distype>::dim,
      CORE::DRT::UTILS::DisTypeToDim<distype>::dim>
      Finc;
  Finc.MultiplyTT(xdisp, N_xyz);
  for (auto i = 0; i < CORE::DRT::UTILS::DisTypeToDim<distype>::dim; ++i)
  {
    defgrd(i, i) += 1.0;
  }

  // get stored old incremental F
  CORE::LINALG::Matrix<CORE::DRT::UTILS::DisTypeToDim<distype>::dim,
      CORE::DRT::UTILS::DisTypeToDim<distype>::dim>
      Fhist;
  mulfHistory->StoragetoMatrix(gp, Fhist, mulfHistory->FHistory());

  // build total defgrd = delta F * F_old
  defgrd.Multiply(Finc, Fhist);
}

template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::UTILS::ComputeDeformationGradientStandard(
    CORE::LINALG::Matrix<probdim, probdim>& defgrd,
    const CORE::LINALG::Matrix<CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement,
        probdim>& xcurr,
    const CORE::LINALG::Matrix<probdim,
        CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement>& derivs,
    const CORE::LINALG::Matrix<probdim, probdim>& inverseJacobian)
{
  CORE::LINALG::Matrix<probdim,
      CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement>
      N_XYZ(false);
  N_XYZ.Multiply(inverseJacobian, derivs);

  defgrd.MultiplyTT(xcurr, N_XYZ);
}

template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::UTILS::EvaluateNodalCoordinates(DRT::Node** nodes,
    CORE::LINALG::Matrix<CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement,
        probdim>& xrefe)
{
  for (auto i = 0; i < CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement; ++i)
  {
    const double* x = nodes[i]->X();
    for (auto dim = 0; dim < probdim; ++dim) xrefe(i, dim) = x[dim];
  }
}

template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::UTILS::EvaluateNodalDisplacements(const std::vector<double>& disp,
    CORE::LINALG::Matrix<CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement,
        probdim>& xdisp)
{
  for (auto i = 0; i < CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement; ++i)
  {
    for (auto dim = 0; dim < probdim; ++dim) xdisp(i, dim) = disp[i * probdim + dim];
  }
}

template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::UTILS::EvaluateCurrentNodalCoordinates(
    const CORE::LINALG::Matrix<CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement,
        probdim>& xrefe,
    const CORE::LINALG::Matrix<CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement,
        probdim>& xdisp,
    CORE::LINALG::Matrix<CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement,
        probdim>& xcurr)
{
  xcurr.Update(1.0, xrefe, 1.0, xdisp);
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::UTILS::EvaluateInverseJacobian(
    const CORE::LINALG::Matrix<CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement,
        CORE::DRT::UTILS::DisTypeToDim<distype>::dim>& xrefe,
    const CORE::LINALG::Matrix<CORE::DRT::UTILS::DisTypeToDim<distype>::dim,
        CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement>& derivs,
    CORE::LINALG::Matrix<CORE::DRT::UTILS::DisTypeToDim<distype>::dim,
        CORE::DRT::UTILS::DisTypeToDim<distype>::dim>& inverseJacobian)
{
  inverseJacobian.Multiply(1.0, derivs, xrefe, 0.0);
  inverseJacobian.Invert();
}

template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::UTILS::EvaluateShapeFunctionSpatialDerivative(
    CORE::LINALG::Matrix<probdim,
        CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement>& deriv_xyz,
    const CORE::LINALG::Matrix<CORE::DRT::UTILS::DisTypeToDim<distype>::dim,
        CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement>& deriv,
    const CORE::LINALG::Matrix<CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement,
        probdim>& xyze,
    const CORE::LINALG::Matrix<probdim, 1>& normal)
{
  const int nen = CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;
  const int nsd_ele = CORE::DRT::UTILS::DisTypeToDim<distype>::dim;

  static CORE::LINALG::Matrix<probdim, nen> deriv_full;

  // transform the derivatives and Jacobians to the higher dimensional coordinates (problem
  // dimension)
  static CORE::LINALG::Matrix<nsd_ele, probdim> dx_dr_red;
  static CORE::LINALG::Matrix<probdim, probdim> dx_dr, dr_dx;
  dx_dr_red.MultiplyNN(deriv, xyze);

  for (unsigned i = 0; i < probdim; ++i)
  {
    for (unsigned j = 0; j < nsd_ele; ++j) dx_dr(j, i) = dx_dr_red(j, i);
    dx_dr(nsd_ele, i) = normal(i);
  }

  for (unsigned i = 0; i < nen; ++i)
  {
    for (unsigned j = 0; j < nsd_ele; ++j) deriv_full(j, i) = deriv(j, i);
    deriv_full(nsd_ele, i) = 0.0;
  }

  // special case: 1D element embedded in 3D problem
  if (nsd_ele == 1 and probdim == 3)
  {
    // compute second unit normal
    const double normalvec2_0 = dx_dr_red(0, 1) * normal(2) - normal(1) * dx_dr_red(0, 2);
    const double normalvec2_1 = dx_dr_red(0, 2) * normal(0) - normal(2) * dx_dr_red(0, 0);
    const double normalvec2_2 = dx_dr_red(0, 0) * normal(1) - normal(0) * dx_dr_red(0, 1);

    // norm
    const double norm2 = std::sqrt(
        normalvec2_0 * normalvec2_0 + normalvec2_1 * normalvec2_1 + normalvec2_2 * normalvec2_2);

    dx_dr(2, 0) = normalvec2_0 / norm2;
    dx_dr(2, 1) = normalvec2_1 / norm2;
    dx_dr(2, 2) = normalvec2_2 / norm2;

    for (unsigned i = 0; i < nen; i++) deriv_full(2, i) = 0.0;
  }

  dr_dx.Invert(dx_dr);

  // compute global spatial derivatives
  deriv_xyz.Multiply(dr_dx, deriv_full);
}

void DRT::ELEMENTS::UTILS::ThrowErrorFDMaterialTangent(
    const Teuchos::ParameterList& sdyn, const std::string& eletype)
{
  bool doFDCheck = static_cast<bool>(DRT::INPUT::IntegralValue<int>(sdyn, "MATERIALTANGENT"));
  if (doFDCheck)
  {
    dserror(
        "Approximation of material tangent by finite differences not implemented by %s elements. "
        "Set parameter MATERIALTANGENT to analytical.",
        eletype.c_str());
  }
}

template void DRT::ELEMENTS::UTILS::CalcR<DRT::Element::tet10>(
    const DRT::Element*, const std::vector<double>&, CORE::LINALG::Matrix<3, 3>&);

template void DRT::ELEMENTS::UTILS::GetTemperatureForStructuralMaterial<DRT::Element::tet4>(
    const CORE::LINALG::Matrix<
        CORE::DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tet4>::numNodePerElement, 1>&
        shapefctsGP,
    Teuchos::ParameterList& params);

template void DRT::ELEMENTS::UTILS::GetTemperatureForStructuralMaterial<DRT::Element::hex27>(
    const CORE::LINALG::Matrix<
        CORE::DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex27>::numNodePerElement, 1>&
        shapefctsGP,
    Teuchos::ParameterList& params);

template void DRT::ELEMENTS::UTILS::GetTemperatureForStructuralMaterial<DRT::Element::hex8>(
    const CORE::LINALG::Matrix<
        CORE::DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement, 1>&
        shapefctsGP,
    Teuchos::ParameterList& params);

template void DRT::ELEMENTS::UTILS::GetTemperatureForStructuralMaterial<DRT::Element::nurbs27>(
    const CORE::LINALG::Matrix<
        CORE::DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::nurbs27>::numNodePerElement, 1>&
        shapefctsGP,
    Teuchos::ParameterList& params);

template void DRT::ELEMENTS::UTILS::GetTemperatureForStructuralMaterial<DRT::Element::tet10>(
    const CORE::LINALG::Matrix<
        CORE::DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tet10>::numNodePerElement, 1>&
        shapefctsGP,
    Teuchos::ParameterList& params);

template void DRT::ELEMENTS::UTILS::GetTemperatureForStructuralMaterial<DRT::Element::hex20>(
    const CORE::LINALG::Matrix<
        CORE::DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex20>::numNodePerElement, 1>&
        shapefctsGP,
    Teuchos::ParameterList& params);

template void DRT::ELEMENTS::UTILS::ComputeDeformationGradient<DRT::Element::hex8, 3>(
    CORE::LINALG::Matrix<3, 3>& defgrd, DRT::Node** nodes, const CORE::LINALG::Matrix<3, 1>& xsi,
    const CORE::LINALG::Matrix<8, 3>& xdisp);
template void DRT::ELEMENTS::UTILS::ComputeDeformationGradient<DRT::Element::tet4, 3>(
    CORE::LINALG::Matrix<3, 3>& defgrd, DRT::Node** nodes, const CORE::LINALG::Matrix<3, 1>& xsi,
    const CORE::LINALG::Matrix<4, 3>& xdisp);

template void DRT::ELEMENTS::UTILS::ComputeDeformationGradient<DRT::Element::hex8, 3>(
    CORE::LINALG::Matrix<3, 3>& defgrd, DRT::Node** nodes, const CORE::LINALG::Matrix<3, 1>& xsi,
    const std::vector<double>& displacement);
template void DRT::ELEMENTS::UTILS::ComputeDeformationGradient<DRT::Element::tet4, 3>(
    CORE::LINALG::Matrix<3, 3>& defgrd, DRT::Node** nodes, const CORE::LINALG::Matrix<3, 1>& xsi,
    const std::vector<double>& displacement);

template void DRT::ELEMENTS::UTILS::ComputeDeformationGradient<DRT::Element::hex8>(
    CORE::LINALG::Matrix<3, 3>& defgrd, const INPAR::STR::KinemType kinemType,
    const CORE::LINALG::Matrix<8, 3>& xdisp, const CORE::LINALG::Matrix<8, 3>& xcurr,
    const CORE::LINALG::Matrix<3, 3>& inverseJacobian, const CORE::LINALG::Matrix<3, 8>& derivs,
    const INPAR::STR::PreStress prestressType,
    const Teuchos::RCP<DRT::ELEMENTS::PreStress> mulfHistory, const int gp);
template void DRT::ELEMENTS::UTILS::ComputeDeformationGradient<DRT::Element::tet4>(
    CORE::LINALG::Matrix<3, 3>& defgrd, const INPAR::STR::KinemType kinemType,
    const CORE::LINALG::Matrix<4, 3>& xdisp, const CORE::LINALG::Matrix<4, 3>& xcurr,
    const CORE::LINALG::Matrix<3, 3>& inverseJacobian, const CORE::LINALG::Matrix<3, 4>& derivs,
    const INPAR::STR::PreStress prestressType,
    const Teuchos::RCP<DRT::ELEMENTS::PreStress> mulfHistory, const int gp);
template void DRT::ELEMENTS::UTILS::ComputeDeformationGradient<DRT::Element::tet10>(
    CORE::LINALG::Matrix<3, 3>& defgrd, const INPAR::STR::KinemType kinemType,
    const CORE::LINALG::Matrix<10, 3>& xdisp, const CORE::LINALG::Matrix<10, 3>& xcurr,
    const CORE::LINALG::Matrix<3, 3>& inverseJacobian, const CORE::LINALG::Matrix<3, 10>& derivs,
    const INPAR::STR::PreStress prestressType,
    const Teuchos::RCP<DRT::ELEMENTS::PreStress> mulfHistory, const int gp);

template void DRT::ELEMENTS::UTILS::ComputeDeformationGradientMulf<DRT::Element::hex8>(
    CORE::LINALG::Matrix<3, 3>& defgrd, const CORE::LINALG::Matrix<8, 3>& xdisp,
    const CORE::LINALG::Matrix<3, 8>& derivs,
    const Teuchos::RCP<DRT::ELEMENTS::PreStress> mulfHistory, const int gp);
template void DRT::ELEMENTS::UTILS::ComputeDeformationGradientMulf<DRT::Element::tet4>(
    CORE::LINALG::Matrix<3, 3>& defgrd, const CORE::LINALG::Matrix<4, 3>& xdisp,
    const CORE::LINALG::Matrix<3, 4>& derivs,
    const Teuchos::RCP<DRT::ELEMENTS::PreStress> mulfHistory, const int gp);
template void DRT::ELEMENTS::UTILS::ComputeDeformationGradientMulf<DRT::Element::tet10>(
    CORE::LINALG::Matrix<3, 3>& defgrd, const CORE::LINALG::Matrix<10, 3>& xdisp,
    const CORE::LINALG::Matrix<3, 10>& derivs,
    const Teuchos::RCP<DRT::ELEMENTS::PreStress> mulfHistory, const int gp);

template void DRT::ELEMENTS::UTILS::ComputeDeformationGradientStandard<DRT::Element::hex8, 3>(
    CORE::LINALG::Matrix<3, 3>& defgrd, const CORE::LINALG::Matrix<8, 3>& xcurr,
    const CORE::LINALG::Matrix<3, 8>& derivs, const CORE::LINALG::Matrix<3, 3>& inverseJacobian);
template void DRT::ELEMENTS::UTILS::ComputeDeformationGradientStandard<DRT::Element::tet4, 3>(
    CORE::LINALG::Matrix<3, 3>& defgrd, const CORE::LINALG::Matrix<4, 3>& xcurr,
    const CORE::LINALG::Matrix<3, 4>& derivs, const CORE::LINALG::Matrix<3, 3>& inverseJacobian);
template void DRT::ELEMENTS::UTILS::ComputeDeformationGradientStandard<DRT::Element::tet10, 3>(
    CORE::LINALG::Matrix<3, 3>& defgrd, const CORE::LINALG::Matrix<10, 3>& xcurr,
    const CORE::LINALG::Matrix<3, 10>& derivs, const CORE::LINALG::Matrix<3, 3>& inverseJacobian);

template void DRT::ELEMENTS::UTILS::EvaluateNodalCoordinates<DRT::Element::hex8, 3>(
    DRT::Node** nodes, CORE::LINALG::Matrix<8, 3>& xrefe);
template void DRT::ELEMENTS::UTILS::EvaluateNodalCoordinates<DRT::Element::tet4, 3>(
    DRT::Node** nodes, CORE::LINALG::Matrix<4, 3>& xrefe);
template void DRT::ELEMENTS::UTILS::EvaluateNodalCoordinates<DRT::Element::tet10, 3>(
    DRT::Node** nodes, CORE::LINALG::Matrix<10, 3>& xrefe);
template void DRT::ELEMENTS::UTILS::EvaluateNodalCoordinates<DRT::Element::quad4, 3>(
    DRT::Node** nodes, CORE::LINALG::Matrix<4, 3>& xrefe);
template void DRT::ELEMENTS::UTILS::EvaluateNodalCoordinates<DRT::Element::tri3, 3>(
    DRT::Node** nodes, CORE::LINALG::Matrix<3, 3>& xrefe);

template void DRT::ELEMENTS::UTILS::EvaluateNodalDisplacements<DRT::Element::hex8, 3>(
    const std::vector<double>&, CORE::LINALG::Matrix<8, 3>& xrefe);
template void DRT::ELEMENTS::UTILS::EvaluateNodalDisplacements<DRT::Element::tet4, 3>(
    const std::vector<double>&, CORE::LINALG::Matrix<4, 3>& xrefe);
template void DRT::ELEMENTS::UTILS::EvaluateNodalDisplacements<DRT::Element::tet10, 3>(
    const std::vector<double>&, CORE::LINALG::Matrix<10, 3>& xrefe);

template void DRT::ELEMENTS::UTILS::EvaluateCurrentNodalCoordinates<DRT::Element::hex8, 3>(
    const CORE::LINALG::Matrix<8, 3>& xrefe, const CORE::LINALG::Matrix<8, 3>& xdisp,
    CORE::LINALG::Matrix<8, 3>& xcurr);
template void DRT::ELEMENTS::UTILS::EvaluateCurrentNodalCoordinates<DRT::Element::tet4, 3>(
    const CORE::LINALG::Matrix<4, 3>& xrefe, const CORE::LINALG::Matrix<4, 3>& xdisp,
    CORE::LINALG::Matrix<4, 3>& xcurr);
template void DRT::ELEMENTS::UTILS::EvaluateCurrentNodalCoordinates<DRT::Element::tet10, 3>(
    const CORE::LINALG::Matrix<10, 3>& xrefe, const CORE::LINALG::Matrix<10, 3>& xdisp,
    CORE::LINALG::Matrix<10, 3>& xcurr);

template void DRT::ELEMENTS::UTILS::EvaluateInverseJacobian<DRT::Element::tet4>(
    const CORE::LINALG::Matrix<4, 3>& xrefe, const CORE::LINALG::Matrix<3, 4>& derivs,
    CORE::LINALG::Matrix<3, 3>& inverseJacobian);

template void DRT::ELEMENTS::UTILS::EvaluateShapeFunctionSpatialDerivative<DRT::Element::line2, 2>(
    CORE::LINALG::Matrix<2, 2>& deriv_xyz, const CORE::LINALG::Matrix<1, 2>& deriv,
    const CORE::LINALG::Matrix<2, 2>& xyze, const CORE::LINALG::Matrix<2, 1>& normal);
template void DRT::ELEMENTS::UTILS::EvaluateShapeFunctionSpatialDerivative<DRT::Element::line2, 3>(
    CORE::LINALG::Matrix<3, 2>& deriv_xyz, const CORE::LINALG::Matrix<1, 2>& deriv,
    const CORE::LINALG::Matrix<2, 3>& xyze, const CORE::LINALG::Matrix<3, 1>& normal);
template void DRT::ELEMENTS::UTILS::EvaluateShapeFunctionSpatialDerivative<DRT::Element::line3, 2>(
    CORE::LINALG::Matrix<2, 3>& deriv_xyz, const CORE::LINALG::Matrix<1, 3>& deriv,
    const CORE::LINALG::Matrix<3, 2>& xyze, const CORE::LINALG::Matrix<2, 1>& normal);
template void DRT::ELEMENTS::UTILS::EvaluateShapeFunctionSpatialDerivative<DRT::Element::nurbs3, 2>(
    CORE::LINALG::Matrix<2, 3>& deriv_xyz, const CORE::LINALG::Matrix<1, 3>& deriv,
    const CORE::LINALG::Matrix<3, 2>& xyze, const CORE::LINALG::Matrix<2, 1>& normal);
template void DRT::ELEMENTS::UTILS::EvaluateShapeFunctionSpatialDerivative<DRT::Element::nurbs9, 3>(
    CORE::LINALG::Matrix<3, 9>& deriv_xyz, const CORE::LINALG::Matrix<2, 9>& deriv,
    const CORE::LINALG::Matrix<9, 3>& xyze, const CORE::LINALG::Matrix<3, 1>& normal);
template void DRT::ELEMENTS::UTILS::EvaluateShapeFunctionSpatialDerivative<DRT::Element::quad4, 3>(
    CORE::LINALG::Matrix<3, 4>& deriv_xyz, const CORE::LINALG::Matrix<2, 4>& deriv,
    const CORE::LINALG::Matrix<4, 3>& xyze, const CORE::LINALG::Matrix<3, 1>& normal);
template void DRT::ELEMENTS::UTILS::EvaluateShapeFunctionSpatialDerivative<DRT::Element::quad8, 3>(
    CORE::LINALG::Matrix<3, 8>& deriv_xyz, const CORE::LINALG::Matrix<2, 8>& deriv,
    const CORE::LINALG::Matrix<8, 3>& xyze, const CORE::LINALG::Matrix<3, 1>& normal);
template void DRT::ELEMENTS::UTILS::EvaluateShapeFunctionSpatialDerivative<DRT::Element::quad9, 3>(
    CORE::LINALG::Matrix<3, 9>& deriv_xyz, const CORE::LINALG::Matrix<2, 9>& deriv,
    const CORE::LINALG::Matrix<9, 3>& xyze, const CORE::LINALG::Matrix<3, 1>& normal);
template void DRT::ELEMENTS::UTILS::EvaluateShapeFunctionSpatialDerivative<DRT::Element::tri3, 3>(
    CORE::LINALG::Matrix<3, 3>& deriv_xyz, const CORE::LINALG::Matrix<2, 3>& deriv,
    const CORE::LINALG::Matrix<3, 3>& xyze, const CORE::LINALG::Matrix<3, 1>& normal);
template void DRT::ELEMENTS::UTILS::EvaluateShapeFunctionSpatialDerivative<DRT::Element::tri6, 3>(
    CORE::LINALG::Matrix<3, 6>& deriv_xyz, const CORE::LINALG::Matrix<2, 6>& deriv,
    const CORE::LINALG::Matrix<6, 3>& xyze, const CORE::LINALG::Matrix<3, 1>& normal);

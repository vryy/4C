/*----------------------------------------------------------------------*/
/*! \file

\brief A collection of helper methods for solid elements

\level 1
 *-----------------------------------------------------------------------*/

#include "4C_so3_utils.hpp"

#include "4C_fem_general_element.hpp"
#include "4C_fem_general_fiber_node.hpp"
#include "4C_fem_general_utils_fem_shapefunctions.hpp"
#include "4C_linalg_utils_densematrix_svd.hpp"
#include "4C_so3_prestress.hpp"

#include <algorithm>

FOUR_C_NAMESPACE_OPEN

template <Core::FE::CellType distype>
void Discret::ELEMENTS::UTILS::CalcR(const Core::Elements::Element* ele,
    const std::vector<double>& disp,
    Core::LinAlg::Matrix<Core::FE::dim<distype>, Core::FE::dim<distype>>& R)
{
  // number of nodes per element
  const int nen = Core::FE::num_nodes<distype>;

  // spatial dimension
  const int nsd = Core::FE::dim<distype>;

  if (disp.size() != nsd * nen) FOUR_C_THROW("mismatch in dimensions");

  Core::LinAlg::Matrix<nsd, 1> xi_ele_center =
      Core::FE::getLocalCenterPosition<nsd>(distype);  // depending on distype

  Core::LinAlg::Matrix<nen, nsd> xrefe;  // X, material coord. of element
  Core::LinAlg::Matrix<nen, nsd> xcurr;  // x, current  coord. of element
  for (int i = 0; i < nen; ++i)
  {
    for (int d = 0; d < nsd; ++d)
    {
      xrefe(i, d) = ele->Nodes()[i]->X()[d];
      xcurr(i, d) = ele->Nodes()[i]->X()[d] + disp[i * nsd + d];
    }
  }
  Core::LinAlg::Matrix<nsd, nen> deriv;
  Core::FE::shape_function_deriv1<distype>(xi_ele_center, deriv);

  Core::LinAlg::Matrix<nsd, nsd> jac;
  Core::LinAlg::Matrix<nsd, nsd> defgrd;
  Core::LinAlg::Matrix<nsd, nen> deriv_xyz;
  jac.Multiply(deriv, xrefe);
  jac.Invert();
  deriv_xyz.Multiply(jac, deriv);
  defgrd.MultiplyTT(xcurr, deriv_xyz);

  // Calculate rotcurr from defgrd
  Core::LinAlg::Matrix<nsd, nsd> Q(true);
  Core::LinAlg::Matrix<nsd, nsd> S(true);
  Core::LinAlg::Matrix<nsd, nsd> VT(true);
  Core::LinAlg::SVD<nsd, nsd>(defgrd, Q, S, VT);
  R.MultiplyNN(Q, VT);
}

template <Core::FE::CellType distype>
void Discret::ELEMENTS::UTILS::get_temperature_for_structural_material(
    const Core::LinAlg::Matrix<Core::FE::num_nodes<distype>,
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
    for (int i = 0; i < Core::FE::num_nodes<distype>; ++i)
    {
      scalartemp += shapefctsGP(i) * (*temperature_vector)[i];
    }

    // insert current element temperature T_{n+1} into parameter list
    params.set<double>("scalartemp", scalartemp);
  }
}

template <Core::FE::CellType distype, int probdim>
void Discret::ELEMENTS::UTILS::compute_deformation_gradient(
    Core::LinAlg::Matrix<probdim, probdim>& defgrd, Core::Nodes::Node** nodes,
    const Core::LinAlg::Matrix<probdim, 1>& xsi,
    const Core::LinAlg::Matrix<Core::FE::num_nodes<distype>, probdim>& xdisp)
{
  static Core::LinAlg::Matrix<Core::FE::num_nodes<distype>, probdim> xrefe, xcurr;

  EvaluateNodalCoordinates<distype, probdim>(nodes, xrefe);
  EvaluateCurrentNodalCoordinates<distype, probdim>(xrefe, xdisp, xcurr);

  Core::LinAlg::Matrix<probdim, Core::FE::num_nodes<distype>> N_rst(true);
  Core::FE::shape_function_deriv1<distype>(xsi, N_rst);

  static Core::LinAlg::Matrix<probdim, probdim> inv_detFJ;
  inv_detFJ.Multiply(N_rst, xrefe);
  inv_detFJ.Invert();

  ComputeDeformationGradientStandard<distype, probdim>(defgrd, xcurr, N_rst, inv_detFJ);
}

template <Core::FE::CellType distype, int probdim>
void Discret::ELEMENTS::UTILS::compute_deformation_gradient(
    Core::LinAlg::Matrix<probdim, probdim>& defgrd, Core::Nodes::Node** nodes,
    const Core::LinAlg::Matrix<probdim, 1>& xsi, const std::vector<double>& displacement)
{
  static Core::LinAlg::Matrix<Core::FE::num_nodes<distype>, probdim> xdisp;
  EvaluateNodalDisplacements<distype, probdim>(displacement, xdisp);

  compute_deformation_gradient<distype, probdim>(defgrd, nodes, xsi, xdisp);
}

template <Core::FE::CellType distype>
void Discret::ELEMENTS::UTILS::compute_deformation_gradient(
    Core::LinAlg::Matrix<Core::FE::dim<distype>, Core::FE::dim<distype>>& defgrd,
    const Inpar::STR::KinemType kinemType,
    const Core::LinAlg::Matrix<Core::FE::num_nodes<distype>, Core::FE::dim<distype>>& xdisp,
    const Core::LinAlg::Matrix<Core::FE::num_nodes<distype>, Core::FE::dim<distype>>& xcurr,
    const Core::LinAlg::Matrix<Core::FE::dim<distype>, Core::FE::dim<distype>>& inverseJacobian,
    const Core::LinAlg::Matrix<Core::FE::dim<distype>, Core::FE::num_nodes<distype>>& derivs,
    const Inpar::STR::PreStress prestressType,
    const Teuchos::RCP<Discret::ELEMENTS::PreStress> mulfHistory, const int gp)
{
  if (kinemType == Inpar::STR::KinemType::linear)
  {
    defgrd.Clear();
    for (auto i = 0; i < Core::FE::dim<distype>; ++i)
    {
      defgrd(i, i) = 1.0;
    }
    return;
  }

  if (prestressType == Inpar::STR::PreStress::mulf)
  {
    ComputeDeformationGradientMulf<distype>(defgrd, xdisp, derivs, mulfHistory, gp);
    return;
  }

  ComputeDeformationGradientStandard<distype, Core::FE::dim<distype>>(
      defgrd, xcurr, derivs, inverseJacobian);
}

template <Core::FE::CellType distype>
void Discret::ELEMENTS::UTILS::ComputeDeformationGradientMulf(
    Core::LinAlg::Matrix<Core::FE::dim<distype>, Core::FE::dim<distype>>& defgrd,
    const Core::LinAlg::Matrix<Core::FE::num_nodes<distype>, Core::FE::dim<distype>>& xdisp,
    const Core::LinAlg::Matrix<Core::FE::dim<distype>, Core::FE::num_nodes<distype>>& derivs,
    const Teuchos::RCP<Discret::ELEMENTS::PreStress> mulfHistory, const int gp)
{
  // get Jacobian mapping wrt to the stored configuration
  Core::LinAlg::Matrix<Core::FE::dim<distype>, Core::FE::dim<distype>> invJdef;
  mulfHistory->StoragetoMatrix(gp, invJdef, mulfHistory->JHistory());

  // get derivatives wrt to last spatial configuration
  Core::LinAlg::Matrix<Core::FE::dim<distype>, Core::FE::num_nodes<distype>> N_xyz;
  N_xyz.Multiply(invJdef, derivs);

  // build multiplicative incremental defgrd
  Core::LinAlg::Matrix<Core::FE::dim<distype>, Core::FE::dim<distype>> Finc;
  Finc.MultiplyTT(xdisp, N_xyz);
  for (auto i = 0; i < Core::FE::dim<distype>; ++i)
  {
    defgrd(i, i) += 1.0;
  }

  // get stored old incremental F
  Core::LinAlg::Matrix<Core::FE::dim<distype>, Core::FE::dim<distype>> Fhist;
  mulfHistory->StoragetoMatrix(gp, Fhist, mulfHistory->FHistory());

  // build total defgrd = delta F * F_old
  defgrd.Multiply(Finc, Fhist);
}

template <Core::FE::CellType distype, int probdim>
void Discret::ELEMENTS::UTILS::ComputeDeformationGradientStandard(
    Core::LinAlg::Matrix<probdim, probdim>& defgrd,
    const Core::LinAlg::Matrix<Core::FE::num_nodes<distype>, probdim>& xcurr,
    const Core::LinAlg::Matrix<probdim, Core::FE::num_nodes<distype>>& derivs,
    const Core::LinAlg::Matrix<probdim, probdim>& inverseJacobian)
{
  Core::LinAlg::Matrix<probdim, Core::FE::num_nodes<distype>> N_XYZ(false);
  N_XYZ.Multiply(inverseJacobian, derivs);

  defgrd.MultiplyTT(xcurr, N_XYZ);
}

template <Core::FE::CellType distype, int probdim>
void Discret::ELEMENTS::UTILS::EvaluateNodalCoordinates(
    Core::Nodes::Node** nodes, Core::LinAlg::Matrix<Core::FE::num_nodes<distype>, probdim>& xrefe)
{
  for (auto i = 0; i < Core::FE::num_nodes<distype>; ++i)
  {
    const auto& x = nodes[i]->X();
    for (auto dim = 0; dim < probdim; ++dim) xrefe(i, dim) = x[dim];
  }
}

template <Core::FE::CellType distype, int probdim>
void Discret::ELEMENTS::UTILS::EvaluateNodalDisplacements(const std::vector<double>& disp,
    Core::LinAlg::Matrix<Core::FE::num_nodes<distype>, probdim>& xdisp)
{
  for (auto i = 0; i < Core::FE::num_nodes<distype>; ++i)
  {
    for (auto dim = 0; dim < probdim; ++dim) xdisp(i, dim) = disp[i * probdim + dim];
  }
}

template <Core::FE::CellType distype, int probdim>
void Discret::ELEMENTS::UTILS::EvaluateCurrentNodalCoordinates(
    const Core::LinAlg::Matrix<Core::FE::num_nodes<distype>, probdim>& xrefe,
    const Core::LinAlg::Matrix<Core::FE::num_nodes<distype>, probdim>& xdisp,
    Core::LinAlg::Matrix<Core::FE::num_nodes<distype>, probdim>& xcurr)
{
  xcurr.Update(1.0, xrefe, 1.0, xdisp);
}

template <Core::FE::CellType distype>
void Discret::ELEMENTS::UTILS::EvaluateInverseJacobian(
    const Core::LinAlg::Matrix<Core::FE::num_nodes<distype>, Core::FE::dim<distype>>& xrefe,
    const Core::LinAlg::Matrix<Core::FE::dim<distype>, Core::FE::num_nodes<distype>>& derivs,
    Core::LinAlg::Matrix<Core::FE::dim<distype>, Core::FE::dim<distype>>& inverseJacobian)
{
  inverseJacobian.Multiply(1.0, derivs, xrefe, 0.0);
  inverseJacobian.Invert();
}

void Discret::ELEMENTS::UTILS::ThrowErrorFDMaterialTangent(
    const Teuchos::ParameterList& sdyn, const std::string& eletype)
{
  bool doFDCheck = static_cast<bool>(Core::UTILS::IntegralValue<int>(sdyn, "MATERIALTANGENT"));
  if (doFDCheck)
  {
    FOUR_C_THROW(
        "Approximation of material tangent by finite differences not implemented by %s elements. "
        "Set parameter MATERIALTANGENT to analytical.",
        eletype.c_str());
  }
}

template void Discret::ELEMENTS::UTILS::CalcR<Core::FE::CellType::tet10>(
    const Core::Elements::Element*, const std::vector<double>&, Core::LinAlg::Matrix<3, 3>&);

template void
Discret::ELEMENTS::UTILS::get_temperature_for_structural_material<Core::FE::CellType::tet4>(
    const Core::LinAlg::Matrix<Core::FE::num_nodes<Core::FE::CellType::tet4>, 1>& shapefctsGP,
    Teuchos::ParameterList& params);

template void
Discret::ELEMENTS::UTILS::get_temperature_for_structural_material<Core::FE::CellType::hex27>(
    const Core::LinAlg::Matrix<Core::FE::num_nodes<Core::FE::CellType::hex27>, 1>& shapefctsGP,
    Teuchos::ParameterList& params);

template void
Discret::ELEMENTS::UTILS::get_temperature_for_structural_material<Core::FE::CellType::hex8>(
    const Core::LinAlg::Matrix<Core::FE::num_nodes<Core::FE::CellType::hex8>, 1>& shapefctsGP,
    Teuchos::ParameterList& params);

template void
Discret::ELEMENTS::UTILS::get_temperature_for_structural_material<Core::FE::CellType::nurbs27>(
    const Core::LinAlg::Matrix<Core::FE::num_nodes<Core::FE::CellType::nurbs27>, 1>& shapefctsGP,
    Teuchos::ParameterList& params);

template void
Discret::ELEMENTS::UTILS::get_temperature_for_structural_material<Core::FE::CellType::tet10>(
    const Core::LinAlg::Matrix<Core::FE::num_nodes<Core::FE::CellType::tet10>, 1>& shapefctsGP,
    Teuchos::ParameterList& params);

template void
Discret::ELEMENTS::UTILS::get_temperature_for_structural_material<Core::FE::CellType::hex20>(
    const Core::LinAlg::Matrix<Core::FE::num_nodes<Core::FE::CellType::hex20>, 1>& shapefctsGP,
    Teuchos::ParameterList& params);

template void Discret::ELEMENTS::UTILS::compute_deformation_gradient<Core::FE::CellType::hex8, 3>(
    Core::LinAlg::Matrix<3, 3>& defgrd, Core::Nodes::Node** nodes,
    const Core::LinAlg::Matrix<3, 1>& xsi, const Core::LinAlg::Matrix<8, 3>& xdisp);
template void Discret::ELEMENTS::UTILS::compute_deformation_gradient<Core::FE::CellType::tet4, 3>(
    Core::LinAlg::Matrix<3, 3>& defgrd, Core::Nodes::Node** nodes,
    const Core::LinAlg::Matrix<3, 1>& xsi, const Core::LinAlg::Matrix<4, 3>& xdisp);

template void Discret::ELEMENTS::UTILS::compute_deformation_gradient<Core::FE::CellType::hex8, 3>(
    Core::LinAlg::Matrix<3, 3>& defgrd, Core::Nodes::Node** nodes,
    const Core::LinAlg::Matrix<3, 1>& xsi, const std::vector<double>& displacement);
template void Discret::ELEMENTS::UTILS::compute_deformation_gradient<Core::FE::CellType::tet4, 3>(
    Core::LinAlg::Matrix<3, 3>& defgrd, Core::Nodes::Node** nodes,
    const Core::LinAlg::Matrix<3, 1>& xsi, const std::vector<double>& displacement);

template void Discret::ELEMENTS::UTILS::compute_deformation_gradient<Core::FE::CellType::hex8>(
    Core::LinAlg::Matrix<3, 3>& defgrd, const Inpar::STR::KinemType kinemType,
    const Core::LinAlg::Matrix<8, 3>& xdisp, const Core::LinAlg::Matrix<8, 3>& xcurr,
    const Core::LinAlg::Matrix<3, 3>& inverseJacobian, const Core::LinAlg::Matrix<3, 8>& derivs,
    const Inpar::STR::PreStress prestressType,
    const Teuchos::RCP<Discret::ELEMENTS::PreStress> mulfHistory, const int gp);
template void Discret::ELEMENTS::UTILS::compute_deformation_gradient<Core::FE::CellType::tet4>(
    Core::LinAlg::Matrix<3, 3>& defgrd, const Inpar::STR::KinemType kinemType,
    const Core::LinAlg::Matrix<4, 3>& xdisp, const Core::LinAlg::Matrix<4, 3>& xcurr,
    const Core::LinAlg::Matrix<3, 3>& inverseJacobian, const Core::LinAlg::Matrix<3, 4>& derivs,
    const Inpar::STR::PreStress prestressType,
    const Teuchos::RCP<Discret::ELEMENTS::PreStress> mulfHistory, const int gp);
template void Discret::ELEMENTS::UTILS::compute_deformation_gradient<Core::FE::CellType::tet10>(
    Core::LinAlg::Matrix<3, 3>& defgrd, const Inpar::STR::KinemType kinemType,
    const Core::LinAlg::Matrix<10, 3>& xdisp, const Core::LinAlg::Matrix<10, 3>& xcurr,
    const Core::LinAlg::Matrix<3, 3>& inverseJacobian, const Core::LinAlg::Matrix<3, 10>& derivs,
    const Inpar::STR::PreStress prestressType,
    const Teuchos::RCP<Discret::ELEMENTS::PreStress> mulfHistory, const int gp);

template void Discret::ELEMENTS::UTILS::ComputeDeformationGradientMulf<Core::FE::CellType::hex8>(
    Core::LinAlg::Matrix<3, 3>& defgrd, const Core::LinAlg::Matrix<8, 3>& xdisp,
    const Core::LinAlg::Matrix<3, 8>& derivs,
    const Teuchos::RCP<Discret::ELEMENTS::PreStress> mulfHistory, const int gp);
template void Discret::ELEMENTS::UTILS::ComputeDeformationGradientMulf<Core::FE::CellType::tet4>(
    Core::LinAlg::Matrix<3, 3>& defgrd, const Core::LinAlg::Matrix<4, 3>& xdisp,
    const Core::LinAlg::Matrix<3, 4>& derivs,
    const Teuchos::RCP<Discret::ELEMENTS::PreStress> mulfHistory, const int gp);
template void Discret::ELEMENTS::UTILS::ComputeDeformationGradientMulf<Core::FE::CellType::tet10>(
    Core::LinAlg::Matrix<3, 3>& defgrd, const Core::LinAlg::Matrix<10, 3>& xdisp,
    const Core::LinAlg::Matrix<3, 10>& derivs,
    const Teuchos::RCP<Discret::ELEMENTS::PreStress> mulfHistory, const int gp);

template void Discret::ELEMENTS::UTILS::ComputeDeformationGradientStandard<Core::FE::CellType::hex8,
    3>(Core::LinAlg::Matrix<3, 3>& defgrd, const Core::LinAlg::Matrix<8, 3>& xcurr,
    const Core::LinAlg::Matrix<3, 8>& derivs, const Core::LinAlg::Matrix<3, 3>& inverseJacobian);
template void Discret::ELEMENTS::UTILS::ComputeDeformationGradientStandard<Core::FE::CellType::tet4,
    3>(Core::LinAlg::Matrix<3, 3>& defgrd, const Core::LinAlg::Matrix<4, 3>& xcurr,
    const Core::LinAlg::Matrix<3, 4>& derivs, const Core::LinAlg::Matrix<3, 3>& inverseJacobian);
template void
Discret::ELEMENTS::UTILS::ComputeDeformationGradientStandard<Core::FE::CellType::tet10, 3>(
    Core::LinAlg::Matrix<3, 3>& defgrd, const Core::LinAlg::Matrix<10, 3>& xcurr,
    const Core::LinAlg::Matrix<3, 10>& derivs, const Core::LinAlg::Matrix<3, 3>& inverseJacobian);

template void Discret::ELEMENTS::UTILS::EvaluateNodalCoordinates<Core::FE::CellType::hex8, 3>(
    Core::Nodes::Node** nodes, Core::LinAlg::Matrix<8, 3>& xrefe);
template void Discret::ELEMENTS::UTILS::EvaluateNodalCoordinates<Core::FE::CellType::tet4, 3>(
    Core::Nodes::Node** nodes, Core::LinAlg::Matrix<4, 3>& xrefe);
template void Discret::ELEMENTS::UTILS::EvaluateNodalCoordinates<Core::FE::CellType::tet10, 3>(
    Core::Nodes::Node** nodes, Core::LinAlg::Matrix<10, 3>& xrefe);
template void Discret::ELEMENTS::UTILS::EvaluateNodalCoordinates<Core::FE::CellType::quad4, 3>(
    Core::Nodes::Node** nodes, Core::LinAlg::Matrix<4, 3>& xrefe);
template void Discret::ELEMENTS::UTILS::EvaluateNodalCoordinates<Core::FE::CellType::tri3, 3>(
    Core::Nodes::Node** nodes, Core::LinAlg::Matrix<3, 3>& xrefe);

template void Discret::ELEMENTS::UTILS::EvaluateNodalDisplacements<Core::FE::CellType::hex8, 3>(
    const std::vector<double>&, Core::LinAlg::Matrix<8, 3>& xrefe);
template void Discret::ELEMENTS::UTILS::EvaluateNodalDisplacements<Core::FE::CellType::tet4, 3>(
    const std::vector<double>&, Core::LinAlg::Matrix<4, 3>& xrefe);
template void Discret::ELEMENTS::UTILS::EvaluateNodalDisplacements<Core::FE::CellType::tet10, 3>(
    const std::vector<double>&, Core::LinAlg::Matrix<10, 3>& xrefe);

template void Discret::ELEMENTS::UTILS::EvaluateCurrentNodalCoordinates<Core::FE::CellType::hex8,
    3>(const Core::LinAlg::Matrix<8, 3>& xrefe, const Core::LinAlg::Matrix<8, 3>& xdisp,
    Core::LinAlg::Matrix<8, 3>& xcurr);
template void Discret::ELEMENTS::UTILS::EvaluateCurrentNodalCoordinates<Core::FE::CellType::tet4,
    3>(const Core::LinAlg::Matrix<4, 3>& xrefe, const Core::LinAlg::Matrix<4, 3>& xdisp,
    Core::LinAlg::Matrix<4, 3>& xcurr);
template void Discret::ELEMENTS::UTILS::EvaluateCurrentNodalCoordinates<Core::FE::CellType::tet10,
    3>(const Core::LinAlg::Matrix<10, 3>& xrefe, const Core::LinAlg::Matrix<10, 3>& xdisp,
    Core::LinAlg::Matrix<10, 3>& xcurr);

template void Discret::ELEMENTS::UTILS::EvaluateInverseJacobian<Core::FE::CellType::tet4>(
    const Core::LinAlg::Matrix<4, 3>& xrefe, const Core::LinAlg::Matrix<3, 4>& derivs,
    Core::LinAlg::Matrix<3, 3>& inverseJacobian);

FOUR_C_NAMESPACE_CLOSE

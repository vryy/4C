/*----------------------------------------------------------------------*/
/*! \file

\brief A collection of helper methods for solid elements

\level 1
 *-----------------------------------------------------------------------*/

#include "so_utils.H"
#include <algorithm>
#include <utility>
#include "linalg_utils_densematrix_svd.H"
#include "element.H"
#include "standardtypes_cpp.H"
#include "utils_fem_shapefunctions.H"
#include "fiber_node.H"
#include "nodal_fiber_holder.H"
#include "prestress.H"

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::UTILS::CalcR(const DRT::Element* ele, const std::vector<double>& disp,
    LINALG::Matrix<DRT::UTILS::DisTypeToDim<distype>::dim, DRT::UTILS::DisTypeToDim<distype>::dim>&
        R)
{
  // number of nodes per element
  const int nen = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

  // spatial dimension
  const int nsd = DRT::UTILS::DisTypeToDim<distype>::dim;

  if (disp.size() != nsd * nen) dserror("mismatch in dimensions");

  LINALG::Matrix<nsd, 1> xi_ele_center =
      DRT::UTILS::getLocalCenterPosition<nsd>(distype);  // depending on distype

  LINALG::Matrix<nen, nsd> xrefe;  // X, material coord. of element
  LINALG::Matrix<nen, nsd> xcurr;  // x, current  coord. of element
  for (int i = 0; i < nen; ++i)
  {
    for (int d = 0; d < nsd; ++d)
    {
      xrefe(i, d) = ele->Nodes()[i]->X()[d];
      xcurr(i, d) = ele->Nodes()[i]->X()[d] + disp[i * nsd + d];
    }
  }
  LINALG::Matrix<nsd, nen> deriv;
  DRT::UTILS::shape_function_deriv1<distype>(xi_ele_center, deriv);

  LINALG::Matrix<nsd, nsd> jac;
  LINALG::Matrix<nsd, nsd> defgrd;
  LINALG::Matrix<nsd, nen> deriv_xyz;
  jac.Multiply(deriv, xrefe);
  jac.Invert();
  deriv_xyz.Multiply(jac, deriv);
  defgrd.MultiplyTT(xcurr, deriv_xyz);

  // Calculate rotcurr from defgrd
  LINALG::Matrix<nsd, nsd> Q(true);
  LINALG::Matrix<nsd, nsd> S(true);
  LINALG::Matrix<nsd, nsd> VT(true);
  LINALG::SVD<nsd, nsd>(defgrd, Q, S, VT);
  R.MultiplyNN(Q, VT);
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::UTILS::GetTemperatureForStructuralMaterial(
    const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement, 1>&
        shapefctsGP,                // shape function of current Gauss-point
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
    for (int i = 0; i < DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement; ++i)
    {
      scalartemp += shapefctsGP(i) * (*temperature_vector)[i];
    }

    // insert current element temperature T_{n+1} into parameter list
    params.set<double>("scalartemp", scalartemp);
  }
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::UTILS::ComputeDeformationGradient(
    LINALG::Matrix<DRT::UTILS::DisTypeToDim<distype>::dim, DRT::UTILS::DisTypeToDim<distype>::dim>&
        defgrd,
    const INPAR::STR::KinemType kinemType,
    const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement,
        DRT::UTILS::DisTypeToDim<distype>::dim>& xdisp,
    const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement,
        DRT::UTILS::DisTypeToDim<distype>::dim>& xcurr,
    const LINALG::Matrix<DRT::UTILS::DisTypeToDim<distype>::dim,
        DRT::UTILS::DisTypeToDim<distype>::dim>& inverseJacobian,
    const LINALG::Matrix<DRT::UTILS::DisTypeToDim<distype>::dim,
        DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement>& derivs,
    const INPAR::STR::PreStress prestressType,
    const Teuchos::RCP<DRT::ELEMENTS::PreStress> mulfHistory, const int gp)
{
  if (kinemType == INPAR::STR::kinem_linear)
  {
    defgrd.Clear();
    for (auto i = 0; i < DRT::UTILS::DisTypeToDim<distype>::dim; ++i)
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

  ComputeDeformationGradientStandard<distype>(defgrd, xcurr, derivs, inverseJacobian);
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::UTILS::ComputeDeformationGradientMulf(
    LINALG::Matrix<DRT::UTILS::DisTypeToDim<distype>::dim, DRT::UTILS::DisTypeToDim<distype>::dim>&
        defgrd,
    const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement,
        DRT::UTILS::DisTypeToDim<distype>::dim>& xdisp,
    const LINALG::Matrix<DRT::UTILS::DisTypeToDim<distype>::dim,
        DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement>& derivs,
    const Teuchos::RCP<DRT::ELEMENTS::PreStress> mulfHistory, const int gp)
{
  // get Jacobian mapping wrt to the stored configuration
  LINALG::Matrix<DRT::UTILS::DisTypeToDim<distype>::dim, DRT::UTILS::DisTypeToDim<distype>::dim>
      invJdef;
  mulfHistory->StoragetoMatrix(gp, invJdef, mulfHistory->JHistory());

  // get derivatives wrt to last spatial configuration
  LINALG::Matrix<DRT::UTILS::DisTypeToDim<distype>::dim,
      DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement>
      N_xyz;
  N_xyz.Multiply(invJdef, derivs);

  // build multiplicative incremental defgrd
  LINALG::Matrix<DRT::UTILS::DisTypeToDim<distype>::dim, DRT::UTILS::DisTypeToDim<distype>::dim>
      Finc;
  Finc.MultiplyTT(xdisp, N_xyz);
  for (auto i = 0; i < DRT::UTILS::DisTypeToDim<distype>::dim; ++i)
  {
    defgrd(i, i) += 1.0;
  }

  // get stored old incremental F
  LINALG::Matrix<DRT::UTILS::DisTypeToDim<distype>::dim, DRT::UTILS::DisTypeToDim<distype>::dim>
      Fhist;
  mulfHistory->StoragetoMatrix(gp, Fhist, mulfHistory->FHistory());

  // build total defgrd = delta F * F_old
  defgrd.Multiply(Finc, Fhist);
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::UTILS::ComputeDeformationGradientStandard(
    LINALG::Matrix<DRT::UTILS::DisTypeToDim<distype>::dim, DRT::UTILS::DisTypeToDim<distype>::dim>&
        defgrd,
    const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement,
        DRT::UTILS::DisTypeToDim<distype>::dim>& xcurr,
    const LINALG::Matrix<DRT::UTILS::DisTypeToDim<distype>::dim,
        DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement>& derivs,
    const LINALG::Matrix<DRT::UTILS::DisTypeToDim<distype>::dim,
        DRT::UTILS::DisTypeToDim<distype>::dim>& inverseJacobian)
{
  LINALG::Matrix<DRT::UTILS::DisTypeToDim<distype>::dim,
      DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement>
      N_XYZ(false);
  N_XYZ.Multiply(inverseJacobian, derivs);

  defgrd.MultiplyTT(xcurr, N_XYZ);
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::UTILS::EvaluateNodalCoordinates(DRT::Node** nodes,
    LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement,
        DRT::UTILS::DisTypeToDim<distype>::dim>& xrefe)
{
  for (auto i = 0; i < DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement; ++i)
  {
    const double* x = nodes[i]->X();
    xrefe(i, 0) = x[0];
    xrefe(i, 1) = x[1];
    xrefe(i, 2) = x[2];
  }
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::UTILS::EvaluateNodalDisplacements(const std::vector<double>& disp,
    LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement,
        DRT::UTILS::DisTypeToDim<distype>::dim>& xdisp)
{
  for (auto i = 0; i < DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement; ++i)
  {
    xdisp(i, 0) = disp[i * DRT::UTILS::DisTypeToDim<distype>::dim + 0];
    xdisp(i, 1) = disp[i * DRT::UTILS::DisTypeToDim<distype>::dim + 1];
    xdisp(i, 2) = disp[i * DRT::UTILS::DisTypeToDim<distype>::dim + 2];
  }
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::UTILS::EvaluateCurrentNodalCoordinates(
    const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement,
        DRT::UTILS::DisTypeToDim<distype>::dim>& xrefe,
    const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement,
        DRT::UTILS::DisTypeToDim<distype>::dim>& xdisp,
    LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement,
        DRT::UTILS::DisTypeToDim<distype>::dim>& xcurr)
{
  xcurr.Update(1.0, xrefe, 1.0, xdisp);
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::UTILS::EvaluateInverseJacobian(
    const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement,
        DRT::UTILS::DisTypeToDim<distype>::dim>& xrefe,
    const LINALG::Matrix<DRT::UTILS::DisTypeToDim<distype>::dim,
        DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement>& derivs,
    LINALG::Matrix<DRT::UTILS::DisTypeToDim<distype>::dim, DRT::UTILS::DisTypeToDim<distype>::dim>&
        inverseJacobian)
{
  inverseJacobian.Multiply(1.0, derivs, xrefe, 0.0);
  inverseJacobian.Invert();
}

void DRT::ELEMENTS::UTILS::ThrowErrorFDMaterialTangent(
    const Teuchos::ParameterList& sdyn, const std::string& eletype)
{
  bool doFDCheck = static_cast<bool>(DRT::INPUT::IntegralValue<int>(sdyn, "MATERIALTANGENT"));
  if (doFDCheck)
    dserror(
        "Approximation of material tangent by finite differences not implemented by %s elements. "
        "Set parameter MATERIALTANGENT to analytical.",
        eletype.c_str());
}

template void DRT::ELEMENTS::UTILS::CalcR<DRT::Element::tet10>(
    const DRT::Element*, const std::vector<double>&, LINALG::Matrix<3, 3>&);

template void DRT::ELEMENTS::UTILS::GetTemperatureForStructuralMaterial<DRT::Element::tet4>(
    const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tet4>::numNodePerElement,
        1>& shapefctsGP,
    Teuchos::ParameterList& params);

template void DRT::ELEMENTS::UTILS::GetTemperatureForStructuralMaterial<DRT::Element::hex27>(
    const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex27>::numNodePerElement,
        1>& shapefctsGP,
    Teuchos::ParameterList& params);

template void DRT::ELEMENTS::UTILS::GetTemperatureForStructuralMaterial<DRT::Element::hex8>(
    const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement,
        1>& shapefctsGP,
    Teuchos::ParameterList& params);

template void DRT::ELEMENTS::UTILS::GetTemperatureForStructuralMaterial<DRT::Element::nurbs27>(
    const LINALG::Matrix<
        DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::nurbs27>::numNodePerElement, 1>&
        shapefctsGP,
    Teuchos::ParameterList& params);

template void DRT::ELEMENTS::UTILS::GetTemperatureForStructuralMaterial<DRT::Element::tet10>(
    const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tet10>::numNodePerElement,
        1>& shapefctsGP,
    Teuchos::ParameterList& params);

template void DRT::ELEMENTS::UTILS::GetTemperatureForStructuralMaterial<DRT::Element::hex20>(
    const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex20>::numNodePerElement,
        1>& shapefctsGP,
    Teuchos::ParameterList& params);


template void DRT::ELEMENTS::UTILS::ComputeDeformationGradient<DRT::Element::hex8>(
    LINALG::Matrix<3, 3>& defgrd, const INPAR::STR::KinemType kinemType,
    const LINALG::Matrix<8, 3>& xdisp, const LINALG::Matrix<8, 3>& xcurr,
    const LINALG::Matrix<3, 3>& inverseJacobian, const LINALG::Matrix<3, 8>& derivs,
    const INPAR::STR::PreStress prestressType,
    const Teuchos::RCP<DRT::ELEMENTS::PreStress> mulfHistory, const int gp);


template void DRT::ELEMENTS::UTILS::ComputeDeformationGradient<DRT::Element::tet4>(
    LINALG::Matrix<3, 3>& defgrd, const INPAR::STR::KinemType kinemType,
    const LINALG::Matrix<4, 3>& xdisp, const LINALG::Matrix<4, 3>& xcurr,
    const LINALG::Matrix<3, 3>& inverseJacobian, const LINALG::Matrix<3, 4>& derivs,
    const INPAR::STR::PreStress prestressType,
    const Teuchos::RCP<DRT::ELEMENTS::PreStress> mulfHistory, const int gp);


template void DRT::ELEMENTS::UTILS::ComputeDeformationGradient<DRT::Element::tet10>(
    LINALG::Matrix<3, 3>& defgrd, const INPAR::STR::KinemType kinemType,
    const LINALG::Matrix<10, 3>& xdisp, const LINALG::Matrix<10, 3>& xcurr,
    const LINALG::Matrix<3, 3>& inverseJacobian, const LINALG::Matrix<3, 10>& derivs,
    const INPAR::STR::PreStress prestressType,
    const Teuchos::RCP<DRT::ELEMENTS::PreStress> mulfHistory, const int gp);

template void DRT::ELEMENTS::UTILS::ComputeDeformationGradientMulf<DRT::Element::hex8>(
    LINALG::Matrix<3, 3>& defgrd, const LINALG::Matrix<8, 3>& xdisp,
    const LINALG::Matrix<3, 8>& derivs, const Teuchos::RCP<DRT::ELEMENTS::PreStress> mulfHistory,
    const int gp);
template void DRT::ELEMENTS::UTILS::ComputeDeformationGradientMulf<DRT::Element::tet4>(
    LINALG::Matrix<3, 3>& defgrd, const LINALG::Matrix<4, 3>& xdisp,
    const LINALG::Matrix<3, 4>& derivs, const Teuchos::RCP<DRT::ELEMENTS::PreStress> mulfHistory,
    const int gp);
template void DRT::ELEMENTS::UTILS::ComputeDeformationGradientMulf<DRT::Element::tet10>(
    LINALG::Matrix<3, 3>& defgrd, const LINALG::Matrix<10, 3>& xdisp,
    const LINALG::Matrix<3, 10>& derivs, const Teuchos::RCP<DRT::ELEMENTS::PreStress> mulfHistory,
    const int gp);

template void DRT::ELEMENTS::UTILS::ComputeDeformationGradientStandard<DRT::Element::hex8>(
    LINALG::Matrix<3, 3>& defgrd, const LINALG::Matrix<8, 3>& xcurr,
    const LINALG::Matrix<3, 8>& derivs, const LINALG::Matrix<3, 3>& inverseJacobian);
template void DRT::ELEMENTS::UTILS::ComputeDeformationGradientStandard<DRT::Element::tet4>(
    LINALG::Matrix<3, 3>& defgrd, const LINALG::Matrix<4, 3>& xcurr,
    const LINALG::Matrix<3, 4>& derivs, const LINALG::Matrix<3, 3>& inverseJacobian);
template void DRT::ELEMENTS::UTILS::ComputeDeformationGradientStandard<DRT::Element::tet10>(
    LINALG::Matrix<3, 3>& defgrd, const LINALG::Matrix<10, 3>& xcurr,
    const LINALG::Matrix<3, 10>& derivs, const LINALG::Matrix<3, 3>& inverseJacobian);

template void DRT::ELEMENTS::UTILS::EvaluateNodalCoordinates<DRT::Element::hex8>(
    DRT::Node** nodes, LINALG::Matrix<8, 3>& xrefe);
template void DRT::ELEMENTS::UTILS::EvaluateNodalCoordinates<DRT::Element::tet4>(
    DRT::Node** nodes, LINALG::Matrix<4, 3>& xrefe);
template void DRT::ELEMENTS::UTILS::EvaluateNodalCoordinates<DRT::Element::tet10>(
    DRT::Node** nodes, LINALG::Matrix<10, 3>& xrefe);

template void DRT::ELEMENTS::UTILS::EvaluateNodalDisplacements<DRT::Element::hex8>(
    const std::vector<double>&, LINALG::Matrix<8, 3>& xrefe);
template void DRT::ELEMENTS::UTILS::EvaluateNodalDisplacements<DRT::Element::tet4>(
    const std::vector<double>&, LINALG::Matrix<4, 3>& xrefe);
template void DRT::ELEMENTS::UTILS::EvaluateNodalDisplacements<DRT::Element::tet10>(
    const std::vector<double>&, LINALG::Matrix<10, 3>& xrefe);

template void DRT::ELEMENTS::UTILS::EvaluateCurrentNodalCoordinates<DRT::Element::hex8>(
    const LINALG::Matrix<8, 3>& xrefe, const LINALG::Matrix<8, 3>& xdisp,
    LINALG::Matrix<8, 3>& xcurr);
template void DRT::ELEMENTS::UTILS::EvaluateCurrentNodalCoordinates<DRT::Element::tet4>(
    const LINALG::Matrix<4, 3>& xrefe, const LINALG::Matrix<4, 3>& xdisp,
    LINALG::Matrix<4, 3>& xcurr);
template void DRT::ELEMENTS::UTILS::EvaluateCurrentNodalCoordinates<DRT::Element::tet10>(
    const LINALG::Matrix<10, 3>& xrefe, const LINALG::Matrix<10, 3>& xdisp,
    LINALG::Matrix<10, 3>& xcurr);

template void DRT::ELEMENTS::UTILS::EvaluateInverseJacobian<DRT::Element::tet4>(
    const LINALG::Matrix<4, 3>& xrefe, const LINALG::Matrix<3, 4>& derivs,
    LINALG::Matrix<3, 3>& inverseJacobian);
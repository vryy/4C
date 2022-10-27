/*----------------------------------------------------------------------*/
/*! \file

\brief main file containing routines for calculation of solid element
       simple displacement based
\level 1
*/
/*----------------------------------------------------------------------*/

#include "solid_ele_calc.H"
#include <Teuchos_ParameterList.hpp>
#include <memory>
#include <optional>
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/voigt_notation.H"
#include "solid_ele.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_mat/so3_material.H"
#include "solid_utils.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "../drt_fiber/drt_fiber_node.H"
#include "../drt_fiber/drt_fiber_utils.H"
#include "../drt_fiber/nodal_fiber_holder.H"

#include "../drt_structure_new/gauss_point_data_output_manager.H"
#include "../drt_so3/so_element_service.H"

namespace
{
  template <DRT::Element::DiscretizationType distype>
  std::unique_ptr<DRT::UTILS::GaussIntegration> CreateDefaultGaussIntegration()
  {
    return CreateDefaultGaussIntegration<distype>(
        DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);
  }

  template <DRT::Element::DiscretizationType distype>
  std::unique_ptr<DRT::UTILS::GaussIntegration> CreateGaussIntegration(
      const DRT::UTILS::IntPointsAndWeights<DRT::UTILS::DisTypeToDim<distype>::dim>& intpoints)
  {
    // format as DRT::UTILS::GaussIntegration
    Teuchos::RCP<DRT::UTILS::CollectedGaussPoints> gp =
        Teuchos::rcp(new DRT::UTILS::CollectedGaussPoints);

    std::array<double, 3> xi = {0., 0., 0.};
    for (int i = 0; i < intpoints.IP().nquad; ++i)
    {
      for (int d = 0; d < DRT::UTILS::DisTypeToDim<distype>::dim; ++d)
        xi[d] = intpoints.IP().qxg[i][d];
      gp->Append(xi[0], xi[1], xi[2], intpoints.IP().qwgt[i]);
    }

    return std::make_unique<DRT::UTILS::GaussIntegration>(gp);
  }
}  // namespace

template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::SolidEleCalc<distype>* DRT::ELEMENTS::SolidEleCalc<distype>::Instance(bool create)
{
  static SolidEleCalc<distype>* instance;
  if (create)
  {
    if (!instance) instance = new SolidEleCalc<distype>();
  }
  else
  {
    if (instance) delete instance;
    instance = nullptr;
  }
  return instance;
}

template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::SolidEleCalc<distype>::SolidEleCalc()
    : DRT::ELEMENTS::SolidEleInterface::SolidEleInterface()
{
  InitializeDefaultQuadrature();
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::SolidEleCalc<distype>::InitializeDefaultQuadrature()
{
  DRT::UTILS::GaussRule3D rule = DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule;
  if (distype == DRT::Element::DiscretizationType::tet10)
  {
    rule = DRT::UTILS::GaussRule3D::tet_4point;
  }
  // setup default integration
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(rule);

  // format as DRT::UTILS::GaussIntegration
  Teuchos::RCP<DRT::UTILS::CollectedGaussPoints> gp =
      Teuchos::rcp(new DRT::UTILS::CollectedGaussPoints);

  std::array<double, 3> xi = {0., 0., 0.};
  for (int i = 0; i < intpoints.IP().nquad; ++i)
  {
    for (int d = 0; d < nsd_; ++d) xi[d] = intpoints.IP().qxg[i][d];
    gp->Append(xi[0], xi[1], xi[2], intpoints.IP().qwgt[i]);
  }

  // save default integration rule
  default_integration_ = Teuchos::rcp(new DRT::UTILS::GaussIntegration(gp));
}


template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::SolidEleCalc<distype>::Evaluate(DRT::ELEMENTS::Solid* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization, std::vector<int>& lm,
    Epetra_SerialDenseMatrix& elemat1, Epetra_SerialDenseMatrix& elemat2,
    Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseVector& elevec2,
    Epetra_SerialDenseVector& elevec3)
{
  // get the action
  // start with "none"
  // todo: remove the string-comparison version of old structural time integration
  ELEMENTS::ActionType act = ELEMENTS::none;

  if (ele->IsParamsInterface())
  {
    act = ele->ParamsInterface().GetActionType();
  }
  else
  {
    // get the required action
    std::string action = params.get<std::string>("action", "none");
    if (action == "none")
      dserror("No action supplied");
    else if (action == "calc_struct_linstiff")
      act = ELEMENTS::struct_calc_linstiff;  // unused
    else if (action == "calc_struct_nlnstiff")
      act = ELEMENTS::struct_calc_nlnstiff;
    else if (action == "calc_struct_internalforce")
      act = ELEMENTS::struct_calc_internalforce;
    else if (action == "calc_struct_linstiffmass")
      act = ELEMENTS::struct_calc_linstiffmass;  // unused
    else if (action == "calc_struct_nlnstiffmass")
      act = ELEMENTS::struct_calc_nlnstiffmass;
    else if (action == "calc_struct_nlnstifflmass")
      act = ELEMENTS::struct_calc_nlnstifflmass;
    else if (action == "calc_struct_nlnstiff_gemm")
      act = ELEMENTS::struct_calc_nlnstiff_gemm;
    else if (action == "calc_struct_stress")
      act = ELEMENTS::struct_calc_stress;
    else if (action == "calc_struct_eleload")
      act = ELEMENTS::struct_calc_eleload;
    else if (action == "calc_struct_fsiload")
      act = ELEMENTS::struct_calc_fsiload;  // unused
    else if (action == "calc_struct_update_istep")
      act = ELEMENTS::struct_calc_update_istep;
    else if (action == "calc_struct_reset_istep")
      act = ELEMENTS::struct_calc_reset_istep;
    else if (action == "calc_struct_store_istep")
      act = ELEMENTS::struct_calc_store_istep;
    else if (action == "calc_struct_recover_istep")
      act = ELEMENTS::struct_calc_recover_istep;
    else if (action == "calc_struct_reset_all")
      act = ELEMENTS::struct_calc_reset_all;  // unused
    else if (action == "calc_struct_energy")
      act = ELEMENTS::struct_calc_energy;
    else if (action == "calc_struct_errornorms")
      act = ELEMENTS::struct_calc_errornorms;
    else if (action == "multi_eas_init")
      act = ELEMENTS::multi_init_eas;  // unused
    else if (action == "multi_eas_set")
      act = ELEMENTS::multi_set_eas;
    else if (action == "multi_readrestart")
      act = ELEMENTS::multi_readrestart;
    else if (action == "multi_calc_dens")
      act = ELEMENTS::multi_calc_dens;
    else if (action == "postprocess_stress")
      act = ELEMENTS::struct_postprocess_stress;
    else if (action == "calc_struct_prestress_update")
      act = ELEMENTS::struct_update_prestress;
    else if (action == "calc_struct_inversedesign_update")
      act = ELEMENTS::inversedesign_update;
    else if (action == "calc_struct_inversedesign_switch")
      act = ELEMENTS::inversedesign_switch;
    else if (action == "calc_global_gpstresses_map")
      act = ELEMENTS::struct_calc_global_gpstresses_map;
    else if (action == "interpolate_velocity_to_given_point")
      act = ELEMENTS::struct_interpolate_velocity_to_point;
    else if (action == "calc_struct_mass_volume")
      act = ELEMENTS::struct_calc_mass_volume;
    else if (action == "struct_calc_predict")
      act = ELEMENTS::struct_calc_predict;
    else
      dserror("Unknown type of action for So_hex8: %s", action.c_str());
  }

  // go through all the potential actions
  switch (act)
  {
    case DRT::ELEMENTS::struct_calc_nlnstiff:
      return nln_force_stiff_mass(
          ele, discretization, lm, params, &elemat1, nullptr, &elevec1, nullptr, nullptr);
      break;
    case struct_calc_internalforce:
      return nln_force_stiff_mass(
          ele, discretization, lm, params, nullptr, nullptr, &elevec1, nullptr, nullptr);
      break;
    case struct_calc_nlnstiffmass:
      return nln_force_stiff_mass(
          ele, discretization, lm, params, &elemat1, &elemat2, &elevec1, nullptr, nullptr);
      break;
    case struct_calc_nlnstifflmass:
    {
      int r = nln_force_stiff_mass(
          ele, discretization, lm, params, &elemat1, &elemat2, &elevec1, nullptr, nullptr);
      LumpMass(&elemat2);
      return r;
    }
    break;
    case struct_calc_nlnstiff_gemm:
      return nln_force_stiff_mass_gemm(
          ele, discretization, lm, params, &elemat1, nullptr, &elevec1, nullptr, nullptr);
      break;
    case struct_calc_stress:
      return CalcStress(
          ele, discretization, lm, params, nullptr, nullptr, nullptr, nullptr, nullptr);
      break;
    case struct_calc_update_istep:
      return UpdateElement(
          ele, discretization, lm, params, nullptr, nullptr, nullptr, nullptr, nullptr);
      break;
    case struct_calc_reset_istep:
      return ResetStep(
          ele, discretization, lm, params, nullptr, nullptr, nullptr, nullptr, nullptr);
      break;
    case struct_calc_store_istep:
      return StoreStep(
          ele, discretization, lm, params, nullptr, nullptr, nullptr, nullptr, nullptr);
      break;
    case struct_calc_recover_istep:
      return RecoverStep(
          ele, discretization, lm, params, nullptr, nullptr, nullptr, nullptr, nullptr);
      break;
    case struct_calc_energy:
      return CalcEnergy(
          ele, discretization, lm, params, nullptr, nullptr, nullptr, nullptr, nullptr);
      break;
    case struct_calc_errornorms:
      return CalcErrorNorms(
          ele, discretization, lm, params, nullptr, nullptr, nullptr, nullptr, nullptr);
      break;
    case multi_set_eas:
      return SetMultiEAS(
          ele, discretization, lm, params, nullptr, nullptr, nullptr, nullptr, nullptr);
      break;
    case multi_readrestart:
      return ReadMultiRestart(
          ele, discretization, lm, params, nullptr, nullptr, nullptr, nullptr, nullptr);
      break;
    case multi_calc_dens:
      return CalcMultiDens(
          ele, discretization, lm, params, nullptr, nullptr, nullptr, nullptr, nullptr);
      break;
    case struct_postprocess_stress:
      return PostProcessStress(
          ele, discretization, lm, params, nullptr, nullptr, nullptr, nullptr, nullptr);
      break;
    case struct_update_prestress:
      return UpdatePrestress(
          ele, discretization, lm, params, nullptr, nullptr, nullptr, nullptr, nullptr);
      break;
    case inversedesign_update:
      return UpdateInverseDesign(
          ele, discretization, lm, params, nullptr, nullptr, nullptr, nullptr, nullptr);
      break;
    case inversedesign_switch:
      return SwitchInverseDesign(
          ele, discretization, lm, params, nullptr, nullptr, nullptr, nullptr, nullptr);
      break;
    case struct_calc_global_gpstresses_map:
      return GlobalGPstressMap(
          ele, discretization, lm, params, nullptr, nullptr, nullptr, nullptr, nullptr);
      break;
    case struct_interpolate_velocity_to_point:
      return InterpolatePoroVelocityToPoint(
          ele, discretization, lm, params, nullptr, nullptr, nullptr, nullptr, nullptr);
      break;
    case struct_calc_mass_volume:
      return CalcMassVolume(
          ele, discretization, lm, params, nullptr, nullptr, nullptr, nullptr, nullptr);
      break;
    case struct_calc_recover:
      return RecoverCondensed(
          ele, discretization, lm, params, nullptr, nullptr, nullptr, nullptr, nullptr);
      break;
    case DRT::ELEMENTS::struct_calc_predict:
      // nothing to do here
      break;
    case DRT::ELEMENTS::struct_init_gauss_point_data_output:
      InitializeGaussPointDataOutput(*ele);
      return 0;
    case DRT::ELEMENTS::struct_gauss_point_data_output:
      EvaluateGaussPointDataOutput(*ele, std::nullopt);
      return 0;
    default:
      dserror("unknown action %s", DRT::ELEMENTS::ActionType2String(act).c_str());
      break;
  }

  return 0;
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::SolidEleCalc<distype>::EvaluateShapeDeriv(
    const DRT::UTILS::GaussIntegration& intpoints, const int gp)
{
  for (int d = 0; d < nsd_; ++d) xi_(d) = intpoints.Point(gp)[d];
  DRT::UTILS::shape_function<distype>(xi_, shapefunction_);
  DRT::UTILS::shape_function_deriv1<distype>(xi_, deriv_);
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::SolidEleCalc<distype>::JacobianMapping(
    const DRT::UTILS::GaussIntegration& intpoints, const int gp)
{
  invJ_.Multiply(deriv_, xrefe_);
  detJ_ = invJ_.Invert();
  fac_ = detJ_ * intpoints.Weight(gp);
  n_xyz_.Multiply(invJ_, deriv_);
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::SolidEleCalc<distype>::Strains()
{
  defgrd_.MultiplyTT(xcurr_, n_xyz_);
  rcg_.MultiplyTN(defgrd_, defgrd_);

  // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
  // todo: make this work for 2D
  gl_strain_(0) = 0.5 * (rcg_(0, 0) - 1.0);
  gl_strain_(1) = 0.5 * (rcg_(1, 1) - 1.0);
  gl_strain_(2) = 0.5 * (rcg_(2, 2) - 1.0);
  gl_strain_(3) = rcg_(0, 1);
  gl_strain_(4) = rcg_(1, 2);
  gl_strain_(5) = rcg_(2, 0);
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::SolidEleCalc<distype>::StrainGradient()
{
  // B-operator
  for (int i = 0; i < nen_; ++i)
  {
    for (int d = 0; d < nsd_; ++d)
      for (int e = 0; e < nsd_; ++e) bop_(d, nsd_ * i + e) = defgrd_(e, d) * n_xyz_(d, i);

    /* ~~~ */
    // todo: make this work for 2D
    bop_(3, nsd_ * i + 0) = defgrd_(0, 0) * n_xyz_(1, i) + defgrd_(0, 1) * n_xyz_(0, i);
    bop_(3, nsd_ * i + 1) = defgrd_(1, 0) * n_xyz_(1, i) + defgrd_(1, 1) * n_xyz_(0, i);
    bop_(3, nsd_ * i + 2) = defgrd_(2, 0) * n_xyz_(1, i) + defgrd_(2, 1) * n_xyz_(0, i);
    bop_(4, nsd_ * i + 0) = defgrd_(0, 1) * n_xyz_(2, i) + defgrd_(0, 2) * n_xyz_(1, i);
    bop_(4, nsd_ * i + 1) = defgrd_(1, 1) * n_xyz_(2, i) + defgrd_(1, 2) * n_xyz_(1, i);
    bop_(4, nsd_ * i + 2) = defgrd_(2, 1) * n_xyz_(2, i) + defgrd_(2, 2) * n_xyz_(1, i);
    bop_(5, nsd_ * i + 0) = defgrd_(0, 2) * n_xyz_(0, i) + defgrd_(0, 0) * n_xyz_(2, i);
    bop_(5, nsd_ * i + 1) = defgrd_(1, 2) * n_xyz_(0, i) + defgrd_(1, 0) * n_xyz_(2, i);
    bop_(5, nsd_ * i + 2) = defgrd_(2, 2) * n_xyz_(0, i) + defgrd_(2, 0) * n_xyz_(2, i);
  }
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::SolidEleCalc<distype>::AddInternalForce(
    LINALG::Matrix<nsd_ * nen_, 1>& force, const double& fac)
{
  force.MultiplyTN(fac, bop_, pk2_, 1.);
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::SolidEleCalc<distype>::AddElasticStiff(
    LINALG::Matrix<nsd_ * nen_, nsd_ * nen_>& stiff, const double& fac)
{
  // integrate `elastic' and `initial-displacement' stiffness matrix
  // keu = keu + (B^T . C . B) * detJ * w(gp)
  cb_.Multiply(cmat_, bop_);
  stiff.MultiplyTN(fac, bop_, cb_, 1.0);
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::SolidEleCalc<distype>::AddGeometricStiff(
    LINALG::Matrix<nsd_ * nen_, nsd_ * nen_>& stiff, const double& fac)
{
  // integrate `geometric' stiffness matrix and add to keu *****************
  std::array<double, 3> SmB_L;  // intermediate Sm.B_L
  // kgeo += (B_L^T . sigma . B_L) * detJ * w(gp)  with B_L = Ni,Xj see NiliFEM-Skript
  for (int inod = 0; inod < nen_; ++inod)
  {
    SmB_L[0] = pk2_(0) * n_xyz_(0, inod) + pk2_(3) * n_xyz_(1, inod) + pk2_(5) * n_xyz_(2, inod);
    SmB_L[1] = pk2_(3) * n_xyz_(0, inod) + pk2_(1) * n_xyz_(1, inod) + pk2_(4) * n_xyz_(2, inod);
    SmB_L[2] = pk2_(5) * n_xyz_(0, inod) + pk2_(4) * n_xyz_(1, inod) + pk2_(2) * n_xyz_(2, inod);
    for (int jnod = 0; jnod < nen_; ++jnod)
    {
      double bopstrbop = 0.0;  // intermediate value
      for (int idim = 0; idim < nsd_; ++idim) bopstrbop += n_xyz_(idim, jnod) * SmB_L[idim];
      for (int d = 0; d < nsd_; ++d) stiff(nsd_ * inod + d, nsd_ * jnod + d) += fac_ * bopstrbop;
    }
  }  // end of integrate `geometric' stiffness******************************
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::SolidEleCalc<distype>::AddMass(
    LINALG::Matrix<nsd_ * nen_, nsd_ * nen_>& mass, const double& fac)
{
  // integrate consistent mass matrix
  double ifactor, massfactor;
  for (int inod = 0; inod < nen_; ++inod)
  {
    ifactor = shapefunction_(inod) * fac;
    for (int jnod = 0; jnod < nen_; ++jnod)
    {
      massfactor = shapefunction_(jnod) * ifactor;  // intermediate factor
      for (int d = 0; d < nsd_; ++d) mass(nsd_ * inod + d, nsd_ * jnod + d) += massfactor;
    }
  }
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::SolidEleCalc<distype>::FillPositionArray(
    DRT::ELEMENTS::Solid* ele, DRT::Discretization& discretization, const std::vector<int>& lm)
{
  // get displacements
  Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
  std::vector<double> mydisp(lm.size());
  DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);

  for (int i = 0; i < nen_; ++i)
  {
    for (int d = 0; d < nsd_; ++d)
    {
      xrefe_(i, d) = ele->Nodes()[i]->X()[d];
      xcurr_(i, d) = xrefe_(i, d) + mydisp[i * nsd_ + d];
    }
  }
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::SolidEleCalc<distype>::Material(LINALG::Matrix<nsd_, nsd_>* defgrd,
    LINALG::Matrix<numstr_, 1>* glstrain, DRT::ELEMENTS::Solid* ele, Teuchos::ParameterList& params,
    int gp)
{
  // evaluate material
  params.set<int>("gp", gp);  // todo: don't use the params here but some other interface

  if (ele->IsParamsInterface())
  {
    ParamsInterfaceToList(ele->ParamsInterface(), params);
  }
  pk2_.Clear();
  cmat_.Clear();

  ele->SolidMaterial()->Evaluate(defgrd, glstrain, params, &pk2_, &cmat_, gp, ele->Id());
}

template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::SolidEleCalc<distype>::nln_force_stiff_mass(DRT::ELEMENTS::Solid* ele,
    DRT::Discretization& discretization, const std::vector<int>& lm, Teuchos::ParameterList& params,
    Epetra_SerialDenseMatrix* stiff_ep, Epetra_SerialDenseMatrix* mass_ep,
    Epetra_SerialDenseVector* force_ep, Epetra_SerialDenseVector* elevec2_epetra,
    Epetra_SerialDenseVector* elevec3_epetra, DRT::UTILS::GaussIntegration* intpoints)
{
  DRT::UTILS::GaussIntegration* my_intrule = default_integration_.getRawPtr();
  if (intpoints) my_intrule = intpoints;

  // get view
  Teuchos::RCP<LINALG::Matrix<nsd_ * nen_, nsd_* nen_>> stiff = Teuchos::null;
  Teuchos::RCP<LINALG::Matrix<nsd_ * nen_, nsd_* nen_>> mass = Teuchos::null;
  Teuchos::RCP<LINALG::Matrix<nsd_ * nen_, 1>> force;
  if (stiff_ep) stiff = Teuchos::rcp(new LINALG::Matrix<nsd_ * nen_, nsd_ * nen_>(*stiff_ep, true));
  if (mass_ep) mass = Teuchos::rcp(new LINALG::Matrix<nsd_ * nen_, nsd_ * nen_>(*mass_ep, true));
  if (force_ep) force = Teuchos::rcp(new LINALG::Matrix<nsd_ * nen_, 1>(*force_ep, true));

  FillPositionArray(ele, discretization, lm);

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  for (int gp = 0; gp < my_intrule->NumPoints(); ++gp)
  {
    EvaluateShapeDeriv(*my_intrule, gp);

    JacobianMapping(*my_intrule, gp);

    Strains();

    StrainGradient();

    Material(&defgrd_, &gl_strain_, ele, params, gp);

    if (force != Teuchos::null) AddInternalForce(*force, fac_);

    if (stiff != Teuchos::null)
    {
      AddElasticStiff(*stiff, fac_);
      AddGeometricStiff(*stiff, fac_);
    }

    if (mass != Teuchos::null) AddMass(*mass, fac_ * ele->Material()->Density());

  }  // gp loop

  return 0;
}

template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::SolidEleCalc<distype>::EvaluateNeumann(DRT::ELEMENTS::Solid* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization, DRT::Condition& condition,
    std::vector<int>& lm, const DRT::UTILS::GaussIntegration* intpoints,
    Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseMatrix* elemat1)
{
  dserror("not implemented");
  return 0;
}

template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::SolidEleCalc<distype>::UpdateElement(DRT::ELEMENTS::Solid* ele,
    DRT::Discretization& discretization, const std::vector<int>& lm, Teuchos::ParameterList& params,
    Epetra_SerialDenseMatrix* elemat1_epetra, Epetra_SerialDenseMatrix* elemat2_epetra,
    Epetra_SerialDenseVector* elevec1_epetra, Epetra_SerialDenseVector* elevec2_epetra,
    Epetra_SerialDenseVector* elevec3_epetra)
{
  FillPositionArray(ele, discretization, lm);

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  for (int gp = 0; gp < default_integration_->NumPoints(); ++gp)
  {
    EvaluateShapeDeriv(*default_integration_, gp);

    JacobianMapping(*default_integration_, gp);

    Strains();

    if (ele->IsParamsInterface())
    {
      ParamsInterfaceToList(ele->ParamsInterface(), params);
    }

    ele->SolidMaterial()->Update(defgrd_, gp, params, ele->Id());

  }  // gp loop

  return 0;
}

template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::SolidEleCalc<distype>::PostProcessStress(DRT::ELEMENTS::Solid* ele,
    DRT::Discretization& discretization, const std::vector<int>& lm, Teuchos::ParameterList& params,
    Epetra_SerialDenseMatrix* elemat1_epetra, Epetra_SerialDenseMatrix* elemat2_epetra,
    Epetra_SerialDenseVector* elevec1_epetra, Epetra_SerialDenseVector* elevec2_epetra,
    Epetra_SerialDenseVector* elevec3_epetra)
{
  const Teuchos::RCP<std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix>>> gpstressmap =
      params.get<Teuchos::RCP<std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix>>>>(
          "gpstressmap", Teuchos::null);
  if (gpstressmap == Teuchos::null) dserror("no gp stress/strain map available for postprocessing");
  std::string stresstype = params.get<std::string>("stresstype", "ndxyz");
  int gid = ele->Id();
  LINALG::Matrix<
      STR::UTILS::IntRuleToNquad<DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule>::ngp, numstr_>
      gpstress(((*gpstressmap)[gid])->A(), true);
  Teuchos::RCP<Epetra_MultiVector> poststress =
      params.get<Teuchos::RCP<Epetra_MultiVector>>("poststress", Teuchos::null);
  if (poststress == Teuchos::null) dserror("No element stress/strain vector available");
  if (stresstype == "ndxyz")
  {
    // TODO: IMPLEMENT
    // dserror("no stress extrapolation in new solid elements yet");
  }
  else if (stresstype == "cxyz")
  {
    const Epetra_BlockMap& elemap = poststress->Map();
    int lid = elemap.LID(ele->Id());
    if (lid != -1)
    {
      for (int i = 0; i < numstr_; ++i)
      {
        double& s = (*((*poststress)(i)))[lid];  // resolve pointer for faster access
        s = 0.;
        for (int j = 0;
             j <
             STR::UTILS::IntRuleToNquad<DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule>::ngp;
             ++j)
        {
          s += gpstress(j, i);
        }
        s *= 1.0 /
             STR::UTILS::IntRuleToNquad<DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule>::ngp;
      }
    }
  }
  else
    dserror("unknown type of stress/strain output on element level");

  return 0;
}

template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::SolidEleCalc<distype>::CalcStress(DRT::ELEMENTS::Solid* ele,
    DRT::Discretization& discretization, const std::vector<int>& lm, Teuchos::ParameterList& params,
    Epetra_SerialDenseMatrix* elemat1_epetra, Epetra_SerialDenseMatrix* elemat2_epetra,
    Epetra_SerialDenseVector* elevec1_epetra, Epetra_SerialDenseVector* elevec2_epetra,
    Epetra_SerialDenseVector* elevec3_epetra)
{
  if (discretization.Comm().MyPID() != ele->Owner()) return 0;

  FillPositionArray(ele, discretization, lm);

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  for (int gp = 0; gp < default_integration_->NumPoints(); ++gp)
  {
    EvaluateShapeDeriv(*default_integration_, gp);

    JacobianMapping(*default_integration_, gp);

    Strains();

    Material(&defgrd_, &gl_strain_, ele, params, gp);

    CalcStressStrainOutput(ele, gp);

  }  // gp loop

  WriteStressStrainOutput(ele);

  return 0;
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::SolidEleCalc<distype>::LumpMass(Epetra_SerialDenseMatrix* mass)
{
  if (!mass) return;

  // we assume mass is a square matrix
  for (int c = 0; c < (*mass).N(); ++c)  // parse columns
  {
    double d = 0.0;
    for (int r = 0; r < (*mass).M(); ++r)  // parse rows
    {
      d += (*mass)(r, c);  // accumulate row entries
      (*mass)(r, c) = 0.0;
    }
    (*mass)(c, c) = d;  // apply sum of row entries on diagonal
  }
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::SolidEleCalc<distype>::Setup(
    DRT::ELEMENTS::Solid* ele, DRT::INPUT::LineDefinition* linedef)
{
  ele->SolidMaterial()->Setup(default_integration_->NumPoints(), linedef);
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::SolidEleCalc<distype>::CalcStressStrainOutput(
    DRT::ELEMENTS::Solid* ele, const int gp)
{
  INPAR::STR::StressType iostress = ele->ParamsInterface().GetStressOutputType();
  INPAR::STR::StrainType iostrain = ele->ParamsInterface().GetStrainOutputType();

  switch (iostress)
  {
    case INPAR::STR::stress_2pk:
      for (int i = 0; i < numstr_; ++i) stress_output_(gp, i) = pk2_(i);
      break;
    case INPAR::STR::stress_cauchy:
    {
      static LINALG::Matrix<numstr_, 1> cauchy;
      STR::UTILS::Pk2ToCauchy(pk2_, defgrd_, cauchy);
      for (int i = 0; i < numstr_; ++i) stress_output_(gp, i) = cauchy(i);
      break;
    }
    case INPAR::STR::stress_none:
      break;
  }

  switch (iostrain)
  {
    case INPAR::STR::strain_gl:
      for (int i = 0; i < numstr_; ++i) strain_output_(gp, i) = gl_strain_(i);
      break;
    case INPAR::STR::strain_ea:
      static LINALG::Matrix<numstr_, 1> ea;
      static LINALG::Matrix<nsd_, nsd_> ea_m;
      ea_m.MultiplyNT(defgrd_, defgrd_);
      ::UTILS::VOIGT::Strains::MatrixToVector(ea_m, ea);
      for (int i = 0; i < numstr_; ++i) strain_output_(gp, i) = ea(i);
      break;
    case INPAR::STR::strain_none:
      break;
    default:
      dserror("strain type not supported");
      break;
  }
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::SolidEleCalc<distype>::WriteStressStrainOutput(DRT::ELEMENTS::Solid* ele)
{
  Teuchos::RCP<std::vector<char>> stressdata = ele->ParamsInterface().MutableStressDataPtr();
  Teuchos::RCP<std::vector<char>> straindata = ele->ParamsInterface().MutableStrainDataPtr();

  {
    DRT::PackBuffer data;
    ele->AddtoPack(data, stress_output_);
    data.StartPacking();
    ele->AddtoPack(data, stress_output_);
    std::copy(data().begin(), data().end(), std::back_inserter(*stressdata));
  }

  {
    DRT::PackBuffer data;
    ele->AddtoPack(data, strain_output_);
    data.StartPacking();
    ele->AddtoPack(data, strain_output_);
    std::copy(data().begin(), data().end(), std::back_inserter(*straindata));
  }
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::SolidEleCalc<distype>::MaterialPostSetup(const DRT::ELEMENTS::Solid& ele)
{
  Teuchos::ParameterList params{};
  if (DRT::FIBER::UTILS::HaveNodalFibers<distype>(ele.Nodes()))
  {
    // This element has fiber nodes.
    // Interpolate fibers to the Gauss points and pass them to the material

    // Get shape functions
    const static std::vector<LINALG::Matrix<nen_, 1>> shapefcts = std::invoke(
        [&]
        {
          const DRT::UTILS::GaussIntegration& integration_rule = *default_integration_;

          std::vector<LINALG::Matrix<nen_, 1>> shapefcns(integration_rule.NumPoints());
          for (int gp = 0; gp < integration_rule.NumPoints(); ++gp)
          {
            LINALG::Matrix<nsd_, 1> xi(integration_rule.Point(gp), true);
            DRT::UTILS::shape_function<distype>(xi, shapefcns[gp]);
          }
          return shapefcns;
        });

    // add fibers to the ParameterList
    DRT::FIBER::NodalFiberHolder fiberHolder;

    // Do the interpolation
    DRT::FIBER::UTILS::ProjectFibersToGaussPoints<distype>(ele.Nodes(), shapefcts, fiberHolder);

    params.set("fiberholder", fiberHolder);
  }

  // Call PostSetup of material
  ele.SolidMaterial()->PostSetup(params, ele.Id());
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::SolidEleCalc<distype>::InitializeGaussPointDataOutput(DRT::ELEMENTS::Solid& ele,
    const std::optional<DRT::UTILS::GaussIntegration>& param_intpoints) const
{
  const DRT::UTILS::GaussIntegration& integration = param_intpoints.value_or(*default_integration_);

  dsassert(ele.IsParamsInterface(),
      "This action type should only be called from the new time integration framework!");

  // Save number of Gauss of the element for gauss point data output
  ele.ParamsInterface().MutableGaussPointDataOutputManagerPtr()->AddElementNumberOfGaussPoints(
      integration.NumPoints());

  // holder for output quantity names and their size
  std::unordered_map<std::string, int> quantities_map{};

  // Ask material for the output quantity names and sizes
  ele.SolidMaterial()->RegisterVtkOutputDataNames(quantities_map);

  // Add quantities to the Gauss point output data manager (if they do not already exist)
  ele.ParamsInterface().MutableGaussPointDataOutputManagerPtr()->MergeQuantities(quantities_map);
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::SolidEleCalc<distype>::EvaluateGaussPointDataOutput(DRT::ELEMENTS::Solid& ele,
    const std::optional<DRT::UTILS::GaussIntegration>& param_intpoints) const
{
  dsassert(ele.IsParamsInterface(),
      "This action type should only be called from the new time integration framework!");

  const DRT::UTILS::GaussIntegration& integration = param_intpoints.value_or(*default_integration_);
  // Collection and assembly of gauss point data
  for (const auto& quantity :
      ele.ParamsInterface().MutableGaussPointDataOutputManagerPtr()->GetQuantities())
  {
    const std::string& quantity_name = quantity.first;
    const int quantity_size = quantity.second;

    // Step 1: Collect the data for each Gauss point for the material
    LINALG::SerialDenseMatrix gp_data(integration.NumPoints(), quantity_size, true);
    bool data_available = ele.SolidMaterial()->EvaluateVtkOutputData(quantity_name, gp_data);

    // Step 3: Assemble data based on output type (elecenter, postprocessed to nodes, Gauss
    // point)
    if (data_available)
    {
      switch (ele.ParamsInterface().MutableGaussPointDataOutputManagerPtr()->GetOutputType())
      {
        case INPAR::STR::GaussPointDataOutputType::element_center:
        {
          // compute average of the quantities
          Teuchos::RCP<Epetra_MultiVector> global_data =
              ele.ParamsInterface()
                  .MutableGaussPointDataOutputManagerPtr()
                  ->GetMutableElementCenterData()
                  .at(quantity_name);
          DRT::ELEMENTS::AssembleAveragedElementValues(*global_data, gp_data, ele);
          break;
        }
        case INPAR::STR::GaussPointDataOutputType::nodes:
        {
          Teuchos::RCP<Epetra_MultiVector> global_data =
              ele.ParamsInterface()
                  .MutableGaussPointDataOutputManagerPtr()
                  ->GetMutableNodalData()
                  .at(quantity_name);

          Epetra_IntVector& global_nodal_element_count =
              *ele.ParamsInterface()
                   .MutableGaussPointDataOutputManagerPtr()
                   ->GetMutableNodalDataCount()
                   .at(quantity_name);

          ExtrapolateGPQuantityToNodesAndAssemble<distype>(
              ele, gp_data, *global_data, false, integration);
          DRT::ELEMENTS::AssembleNodalElementCount(global_nodal_element_count, ele);
          break;
        }
        case INPAR::STR::GaussPointDataOutputType::gauss_points:
        {
          std::vector<Teuchos::RCP<Epetra_MultiVector>>& global_data =
              ele.ParamsInterface()
                  .MutableGaussPointDataOutputManagerPtr()
                  ->GetMutableGaussPointData()
                  .at(quantity_name);
          DRT::ELEMENTS::AssembleGaussPointValues(global_data, gp_data, ele);
          break;
        }
        case INPAR::STR::GaussPointDataOutputType::none:
          dserror(
              "You specified a Gauss point data output type of none, so you should not end up "
              "here.");
        default:
          dserror("Unknown Gauss point data output type.");
      }
    }
  }
}

// template classes
template class DRT::ELEMENTS::SolidEleCalc<DRT::Element::hex8>;
template class DRT::ELEMENTS::SolidEleCalc<DRT::Element::hex18>;
template class DRT::ELEMENTS::SolidEleCalc<DRT::Element::hex20>;
template class DRT::ELEMENTS::SolidEleCalc<DRT::Element::hex27>;
template class DRT::ELEMENTS::SolidEleCalc<DRT::Element::tet4>;
template class DRT::ELEMENTS::SolidEleCalc<DRT::Element::tet10>;
template class DRT::ELEMENTS::SolidEleCalc<DRT::Element::pyramid5>;

/*----------------------------------------------------------------------*/
/*! \file

\brief scatra_ele_calc_cardiac_monodomain.cpp

\level 2

 *----------------------------------------------------------------------*/


#include "4C_scatra_ele_calc_cardiac_monodomain.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_fem_general_element.hpp"
#include "4C_fem_general_extract_values.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_cardiac_monodomain.hpp"
#include "4C_mat_list.hpp"
#include "4C_mat_myocard.hpp"
#include "4C_scatra_ele_parameter_std.hpp"
#include "4C_scatra_ele_parameter_timint.hpp"
#include "4C_utils_singleton_owner.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
Discret::ELEMENTS::ScaTraEleCalcCardiacMonodomain<distype, probdim>::ScaTraEleCalcCardiacMonodomain(
    const int numdofpernode, const int numscal, const std::string& disname)
    : Discret::ELEMENTS::ScaTraEleCalc<distype, probdim>::ScaTraEleCalc(
          numdofpernode, numscal, disname),
      Discret::ELEMENTS::ScaTraEleCalcAniso<distype, probdim>::ScaTraEleCalcAniso(
          numdofpernode, numscal, disname),
      Discret::ELEMENTS::ScaTraEleCalcAdvReac<distype, probdim>::ScaTraEleCalcAdvReac(
          numdofpernode, numscal, disname)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
Discret::ELEMENTS::ScaTraEleCalcCardiacMonodomain<distype, probdim>*
Discret::ELEMENTS::ScaTraEleCalcCardiacMonodomain<distype, probdim>::Instance(
    const int numdofpernode, const int numscal, const std::string& disname)
{
  static auto singleton_map = Core::UTILS::MakeSingletonMap<std::string>(
      [](const int numdofpernode, const int numscal, const std::string& disname)
      {
        return std::unique_ptr<ScaTraEleCalcCardiacMonodomain<distype, probdim>>(
            new ScaTraEleCalcCardiacMonodomain<distype, probdim>(numdofpernode, numscal, disname));
      });

  return singleton_map[disname].Instance(
      Core::UTILS::SingletonAction::create, numdofpernode, numscal, disname);
}


/*----------------------------------------------------------------------*
 |  evaluate single material  (protected)                    ljag 06/14 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
void Discret::ELEMENTS::ScaTraEleCalcCardiacMonodomain<distype, probdim>::materials(
    const Teuchos::RCP<const Core::Mat::Material> material,  //!< pointer to current material
    const int k,                                             //!< id of current scalar
    double& densn,                                           //!< density at t_(n)
    double& densnp,  //!< density at t_(n+1) or t_(n+alpha_F)
    double& densam,  //!< density at t_(n+alpha_M)
    double& visc,    //!< fluid viscosity
    const int iquad  //!< id of current gauss point

)
{
  // safety check
  if (material->MaterialType() != Core::Materials::m_myocard)
    FOUR_C_THROW("Material type is not supported");

  // safety check
  Teuchos::RCP<Mat::Myocard> actmat = Teuchos::rcp_dynamic_cast<Mat::Myocard>(
      Teuchos::rcp_const_cast<Core::Mat::Material>(material));
  if (actmat->GetNumberOfGP() != 1 and not my::scatrapara_->MatGP())
  {
    actmat->SetGP(1);
    actmat->resize_internal_state_variables();
  }
  mat_myocard(material, k, densn, densnp, densam, visc, iquad);

  return;
}


/*----------------------------------------------------------------------*
 |  Material ScaTra                                          ljag 06/14 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
void Discret::ELEMENTS::ScaTraEleCalcCardiacMonodomain<distype, probdim>::mat_myocard(
    const Teuchos::RCP<const Core::Mat::Material> material,  //!< pointer to current material
    const int k,                                             //!< id of current scalar
    double& densn,                                           //!< density at t_(n)
    double& densnp,  //!< density at t_(n+1) or t_(n+alpha_F)
    double& densam,  //!< density at t_(n+alpha_M)
    double& visc,    //!< fluid viscosity
    const int iquad  //!< id of current gauss point (default = -1)
)
{
  const Teuchos::RCP<const Mat::Myocard>& actmat =
      Teuchos::rcp_dynamic_cast<const Mat::Myocard>(material);

  // dynamic cast to Advanced_Reaction-specific reaction manager
  Teuchos::RCP<ScaTraEleReaManagerAdvReac> advreamanager =
      Teuchos::rcp_dynamic_cast<ScaTraEleReaManagerAdvReac>(my::reamanager_);

  // dynamic cast to anisotropic diffusion manager
  Teuchos::RCP<ScaTraEleDiffManagerAniso<nsd_>> diffmanageraniso =
      Teuchos::rcp_dynamic_cast<ScaTraEleDiffManagerAniso<nsd_>>(my::diffmanager_);

  // get constant diffusivity
  Core::LinAlg::Matrix<nsd_, nsd_> difftensor(true);
  actmat->Diffusivity(difftensor);

  diffmanageraniso->SetAnisotropicDiff(difftensor, k);

  // clear
  advreamanager->Clear(my::numscal_);

  if (my::scatrapara_->SemiImplicit())
  {
    // get membrane potential at n at integration point
    const double phin = my::scatravarmanager_->Phin(k);
    const double phinp = my::scatravarmanager_->Phinp(k);
    // get reaction coefficient
    double react = -actmat->ReaCoeffN(phin, my::scatraparatimint_->Dt(), iquad);
    if (my::scatraparatimint_->IsGenAlpha())
      react *= my::scatraparatimint_->Dt() / my::scatraparatimint_->TimeFac();
    advreamanager->AddToReaBodyForce(react, k);
    advreamanager->add_to_rea_body_force_deriv_matrix(0.0, k, k);
    actmat->ReaCoeff(phinp, my::scatraparatimint_->Dt(), iquad);
  }
  else
  {
    // get membrane potential at n+1 or n+alpha_F at integration point
    const double phinp = my::scatravarmanager_->Phinp(k);
    // get reaction coefficient
    advreamanager->AddToReaBodyForce(
        -actmat->ReaCoeff(phinp, my::scatraparatimint_->Dt(), iquad), k);
    advreamanager->add_to_rea_body_force_deriv_matrix(
        -actmat->ReaCoeffDeriv(phinp, my::scatraparatimint_->Dt(), iquad), k, k);
  }

  return;
}  // ScaTraEleCalcCardiacMonodomain<distype>::MatMyocard


/*----------------------------------------------------------------------*
|  calculate system matrix and rhs for ep                 hoermann 06/16|
*----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
void Discret::ELEMENTS::ScaTraEleCalcCardiacMonodomain<distype, probdim>::sysmat(
    Core::Elements::Element* ele,               ///< the element whose matrix is calculated
    Core::LinAlg::SerialDenseMatrix& emat,      ///< element matrix to calculate
    Core::LinAlg::SerialDenseVector& erhs,      ///< element rhs to calculate
    Core::LinAlg::SerialDenseVector& subgrdiff  ///< subgrid-diff.-scaling vector
)
{
  // density at t_(n) (one per transported scalar)
  std::vector<double> densn(my::numscal_, 1.0);
  // density at t_(n+1) or t_(n+alpha_F) (one per transported scalar)
  std::vector<double> densnp(my::numscal_, 1.0);
  // density at t_(n+alpha_M) (one per transported scalar)
  std::vector<double> densam(my::numscal_, 1.0);

  // fluid viscosity
  double visc(0.0);

  // calculation of material parameter at element center
  if (not my::scatrapara_->MatGP())
  {
    advreac::eval_shape_func_and_derivs_at_ele_center();
    // set Gauss point variables needed for evaluation of mat and rhs
    my::set_internal_variables_for_mat_and_rhs();
    advreac::get_material_params(ele, densn, densnp, densam, visc);
  }
  // calculation of material at integration points (different number of integration points possible)
  else
  {
    int deg = 0;
    if (ele->Degree() == 1)
      deg = 4 * ele->Degree();
    else
      deg = 3 * ele->Degree();
    const Core::FE::IntPointsAndWeights<nsd_ele_> intpoints(
        ScaTra::DisTypeToMatGaussRule<distype>::get_gauss_rule(deg));

    // loop over integration points
    for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
    {
      const double fac = my::eval_shape_func_and_derivs_at_int_point(intpoints, iquad);

      // set gauss point variables needed for evaluation of mat and rhs
      my::set_internal_variables_for_mat_and_rhs();

      // get material parameters (evaluation at integration point)
      advreac::get_material_params(ele, densn, densnp, densam, visc, iquad);

      // loop all scalars
      for (int k = 0; k < my::numscal_; ++k)  // deal with a system of transported scalars
      {
        double rhsint(0.0);
        advreac::get_rhs_int(rhsint, densnp[k], k);

        Core::LinAlg::Matrix<nen_, 1> dummy(true);
        const double timefacfac = my::scatraparatimint_->TimeFac() * fac;

        // reactive terms on integration point on rhs
        my::compute_rhs_int(rhsint, densam[k], densnp[k], my::scatravarmanager_->Hist(k));

        // standard Galerkin transient, old part of rhs and bodyforce term
        my::calc_rhs_hist_and_source(erhs, k, fac, rhsint);

        // element matrix: reactive term
        advreac::calc_mat_react(emat, k, timefacfac, 0., 0., densnp[k], dummy, dummy);
      }
    }
  }

  //----------------------------------------------------------------------
  // integration loop for one element
  //----------------------------------------------------------------------
  // integration points and weights
  const Core::FE::IntPointsAndWeights<nsd_ele_> intpoints(
      ScaTra::DisTypeToOptGaussRule<distype>::rule);

  for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
  {
    const double fac = my::eval_shape_func_and_derivs_at_int_point(intpoints, iquad);

    // set gauss point variables needed for evaluation of mat and rhs
    my::set_internal_variables_for_mat_and_rhs();

    // loop all scalars
    for (int k = 0; k < my::numscal_; ++k)  // deal with a system of transported scalars
    {
      // compute rhs containing bodyforce
      double rhsint(0.0);
      advreac::get_rhs_int(rhsint, densnp[k], k);

      // integration factors
      const double timefacfac = my::scatraparatimint_->TimeFac() * fac;

      //----------------------------------------------------------------
      // 1) element matrix: stationary terms
      //----------------------------------------------------------------

      // calculation of diffusive element matrix
      aniso::calc_mat_diff(emat, k, timefacfac);

      //----------------------------------------------------------------
      // 2) element matrix: instationary terms
      //----------------------------------------------------------------

      if (not my::scatraparatimint_->IsStationary()) my::calc_mat_mass(emat, k, fac, densam[k]);

      //----------------------------------------------------------------
      // 3) element matrix: reactive term
      //----------------------------------------------------------------

      Core::LinAlg::Matrix<nen_, 1> dummy(true);
      if (not my::scatrapara_->MatGP())
        advreac::calc_mat_react(emat, k, timefacfac, 0., 0., densnp[k], dummy, dummy);

      //----------------------------------------------------------------
      // 5) element right hand side
      //----------------------------------------------------------------
      //----------------------------------------------------------------
      // computation of bodyforce (and potentially history) term,
      // residual, integration factors and standard Galerkin transient
      // term (if required) on right hand side depending on respective
      // (non-)incremental stationary or time-integration scheme
      //----------------------------------------------------------------
      double rhsfac = my::scatraparatimint_->TimeFacRhs() * fac;

      if (my::scatraparatimint_->IsIncremental() and not my::scatraparatimint_->IsStationary())
        my::calc_rhs_lin_mass(erhs, k, rhsfac, fac, densam[k], densnp[k]);


      if (not my::scatrapara_->MatGP())
      {
        my::compute_rhs_int(rhsint, densam[k], densnp[k], my::scatravarmanager_->Hist(k));
        // standard Galerkin transient, old part of rhs and bodyforce term
        my::calc_rhs_hist_and_source(erhs, k, fac, rhsint);
      }

      //----------------------------------------------------------------
      // standard Galerkin terms on right hand side
      //----------------------------------------------------------------

      // diffusive term
      aniso::calc_rhs_diff(erhs, k, rhsfac);

    }  // end loop all scalars
  }    // end loop Gauss points

  return;
}


/*----------------------------------------------------------------------*
 | extract element based or nodal values                 hoermann 06/16 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
void Discret::ELEMENTS::ScaTraEleCalcCardiacMonodomain<distype,
    probdim>::extract_element_and_node_values(Core::Elements::Element* ele,
    Teuchos::ParameterList& params, Discret::Discretization& discretization,
    Core::Elements::Element::LocationArray& la)
{
  my::extract_element_and_node_values(ele, params, discretization, la);

  // extract additional local values from global vector
  Teuchos::RCP<const Epetra_Vector> phin = discretization.GetState("phin");
  if (phin == Teuchos::null) FOUR_C_THROW("Cannot get state vector 'phin'");
  Core::FE::ExtractMyValues<Core::LinAlg::Matrix<nen_, 1>>(*phin, my::ephin_, la[0].lm_);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
// template classes
// 1D elements
template class Discret::ELEMENTS::ScaTraEleCalcCardiacMonodomain<Core::FE::CellType::line2, 1>;
template class Discret::ELEMENTS::ScaTraEleCalcCardiacMonodomain<Core::FE::CellType::line2, 2>;
template class Discret::ELEMENTS::ScaTraEleCalcCardiacMonodomain<Core::FE::CellType::line2, 3>;
template class Discret::ELEMENTS::ScaTraEleCalcCardiacMonodomain<Core::FE::CellType::line3, 1>;

// 2D elements
template class Discret::ELEMENTS::ScaTraEleCalcCardiacMonodomain<Core::FE::CellType::tri3, 2>;
template class Discret::ELEMENTS::ScaTraEleCalcCardiacMonodomain<Core::FE::CellType::tri3, 3>;
template class Discret::ELEMENTS::ScaTraEleCalcCardiacMonodomain<Core::FE::CellType::tri6, 2>;
template class Discret::ELEMENTS::ScaTraEleCalcCardiacMonodomain<Core::FE::CellType::quad4, 2>;
template class Discret::ELEMENTS::ScaTraEleCalcCardiacMonodomain<Core::FE::CellType::quad4, 3>;
// template class
// Discret::ELEMENTS::ScaTraEleCalcCardiacMonodomain<Core::FE::CellType::quad8>;
template class Discret::ELEMENTS::ScaTraEleCalcCardiacMonodomain<Core::FE::CellType::quad9, 2>;
template class Discret::ELEMENTS::ScaTraEleCalcCardiacMonodomain<Core::FE::CellType::nurbs9, 2>;

// 3D elements
template class Discret::ELEMENTS::ScaTraEleCalcCardiacMonodomain<Core::FE::CellType::hex8, 3>;
// template class
// Discret::ELEMENTS::ScaTraEleCalcCardiacMonodomain<Core::FE::CellType::hex20>;
template class Discret::ELEMENTS::ScaTraEleCalcCardiacMonodomain<Core::FE::CellType::hex27, 3>;
template class Discret::ELEMENTS::ScaTraEleCalcCardiacMonodomain<Core::FE::CellType::tet4, 3>;
template class Discret::ELEMENTS::ScaTraEleCalcCardiacMonodomain<Core::FE::CellType::tet10, 3>;
// template class
// Discret::ELEMENTS::ScaTraEleCalcCardiacMonodomain<Core::FE::CellType::wedge6>;
template class Discret::ELEMENTS::ScaTraEleCalcCardiacMonodomain<Core::FE::CellType::pyramid5, 3>;
// template class
// Discret::ELEMENTS::ScaTraEleCalcCardiacMonodomain<Core::FE::CellType::nurbs27>;

FOUR_C_NAMESPACE_CLOSE

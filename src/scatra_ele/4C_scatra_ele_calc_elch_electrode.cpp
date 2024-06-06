/*--------------------------------------------------------------------------*/
/*! \file

\brief evaluation of scatra elements for conservation of mass concentration and electronic charge
within isothermal electrodes

\level 2

*/
/*--------------------------------------------------------------------------*/
#include "4C_scatra_ele_calc_elch_electrode.hpp"

#include "4C_mat_material_factory.hpp"
#include "4C_material_base.hpp"
#include "4C_scatra_ele_parameter_std.hpp"
#include "4C_scatra_ele_parameter_timint.hpp"
#include "4C_scatra_ele_utils_elch_electrode.hpp"
#include "4C_utils_singleton_owner.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
Discret::ELEMENTS::ScaTraEleCalcElchElectrode<distype, probdim>*
Discret::ELEMENTS::ScaTraEleCalcElchElectrode<distype, probdim>::Instance(
    const int numdofpernode, const int numscal, const std::string& disname)
{
  static auto singleton_map = Core::UTILS::MakeSingletonMap<std::string>(
      [](const int numdofpernode, const int numscal, const std::string& disname)
      {
        return std::unique_ptr<ScaTraEleCalcElchElectrode<distype, probdim>>(
            new ScaTraEleCalcElchElectrode<distype, probdim>(numdofpernode, numscal, disname));
      });

  return singleton_map[disname].Instance(
      Core::UTILS::SingletonAction::create, numdofpernode, numscal, disname);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
Discret::ELEMENTS::ScaTraEleCalcElchElectrode<distype, probdim>::ScaTraEleCalcElchElectrode(
    const int numdofpernode, const int numscal, const std::string& disname)
    : myelch::ScaTraEleCalcElch(numdofpernode, numscal, disname)
{
  // replace elch diffusion manager by diffusion manager for electrodes
  my::diffmanager_ = Teuchos::rcp(new ScaTraEleDiffManagerElchElectrode(my::numscal_));

  // replace elch internal variable manager by internal variable manager for electrodes
  my::scatravarmanager_ =
      Teuchos::rcp(new ScaTraEleInternalVariableManagerElchElectrode<nsd_, nen_>(
          my::numscal_, myelch::elchparams_));

  // replace elch utility class by utility class for electrodes
  myelch::utils_ = Discret::ELEMENTS::ScaTraEleUtilsElchElectrode<distype>::Instance(
      numdofpernode, numscal, disname);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
void Discret::ELEMENTS::ScaTraEleCalcElchElectrode<distype, probdim>::calc_mat_and_rhs(
    Core::LinAlg::SerialDenseMatrix& emat, Core::LinAlg::SerialDenseVector& erhs, const int k,
    const double fac, const double timefacfac, const double rhsfac, const double taufac,
    const double timetaufac, const double rhstaufac, Core::LinAlg::Matrix<nen_, 1>& tauderpot,
    double& rhsint)
{
  //----------------------------------------------------------------------
  // 1) element matrix: instationary terms arising from transport equation
  //----------------------------------------------------------------------

  if (not my::scatraparatimint_->IsStationary())
    // 1a) element matrix: standard Galerkin mass term
    my::calc_mat_mass(emat, k, fac, 1.);

  //--------------------------------------------------------------------
  // 2) element matrix: stationary terms arising from transport equation
  //--------------------------------------------------------------------

  // 2a) element matrix: standard Galerkin diffusive term
  my::calc_mat_diff(emat, k, timefacfac);

  // 2b) element matrix: additional term arising from concentration dependency of diffusion
  // coefficient
  calc_mat_diff_coeff_lin(emat, k, timefacfac, var_manager()->GradPhi(k), 1.);

  // 2c) element matrix: conservative part of convective term, needed for deforming electrodes,
  //                     i.e., for scalar-structure interaction
  double vdiv(0.);
  if (my::scatrapara_->IsConservative())
  {
    my::get_divergence(vdiv, my::evelnp_);
    my::calc_mat_conv_add_cons(emat, k, timefacfac, vdiv, 1.);
  }

  //----------------------------------------------------------------------------
  // 3) element right hand side vector (negative residual of nonlinear problem):
  //    terms arising from transport equation
  //----------------------------------------------------------------------------

  // 3a) element rhs: standard Galerkin contributions from non-history part of instationary term if
  // needed
  if (not my::scatraparatimint_->IsStationary())
    my::calc_rhs_lin_mass(erhs, k, rhsfac, fac, 1., 1.);

  // 3b) element rhs: standard Galerkin contributions from rhsint vector (contains body force vector
  // and history vector) need to adapt rhsint vector to time integration scheme first
  my::compute_rhs_int(rhsint, 1., 1., var_manager()->Hist(k));
  my::calc_rhs_hist_and_source(erhs, k, fac, rhsint);

  // 3c) element rhs: standard Galerkin diffusion term
  my::calc_rhs_diff(erhs, k, rhsfac);

  // 3d) element rhs: conservative part of convective term, needed for deforming electrodes,
  //                  i.e., for scalar-structure interaction
  if (my::scatrapara_->IsConservative())
    calc_rhs_conservative_part_of_convective_term(erhs, k, rhsfac, vdiv);

  //----------------------------------------------------------------------------
  // 4) element matrix: stationary terms arising from potential equation
  // 5) element right hand side vector (negative residual of nonlinear problem):
  //    terms arising from potential equation
  //----------------------------------------------------------------------------
  // see function calc_mat_and_rhs_outside_scalar_loop()
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
void Discret::ELEMENTS::ScaTraEleCalcElchElectrode<distype,
    probdim>::calc_mat_and_rhs_outside_scalar_loop(Core::LinAlg::SerialDenseMatrix& emat,
    Core::LinAlg::SerialDenseVector& erhs, const double fac, const double timefacfac,
    const double rhsfac)
{
  //--------------------------------------------------------------------
  // 4) element matrix: stationary terms arising from potential equation
  //--------------------------------------------------------------------

  // element matrix: standard Galerkin terms from potential equation
  calc_mat_pot_equ_divi_ohm(emat, timefacfac, var_manager()->InvF(), var_manager()->GradPot(), 1.);

  //----------------------------------------------------------------------------
  // 5) element right hand side vector (negative residual of nonlinear problem):
  //    terms arising from potential equation
  //----------------------------------------------------------------------------

  // element rhs: standard Galerkin terms from potential equation
  calc_rhs_pot_equ_divi_ohm(erhs, rhsfac, var_manager()->InvF(), var_manager()->GradPot(), 1.);

  // safety check
  if (my::bodyforce_[my::numscal_].Dot(my::funct_) != 0.0)
    FOUR_C_THROW("body force not implemented for potential equation");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
void Discret::ELEMENTS::ScaTraEleCalcElchElectrode<distype, probdim>::calc_diff_od_mesh(
    Core::LinAlg::SerialDenseMatrix& emat, const int k, const int ndofpernodemesh,
    const double diffcoeff, const double fac, const double rhsfac, const double J,
    const Core::LinAlg::Matrix<nsd_, 1>& gradphi, const Core::LinAlg::Matrix<nsd_, 1>& convelint,
    const Core::LinAlg::Matrix<1, nsd_ * nen_>& dJ_dmesh)
{
  // safety check
  if (k != 0) FOUR_C_THROW("Invalid species index!");

  // call base class routine to compute linearizations of diffusion term w.r.t. structural
  // displacements
  my::calc_diff_od_mesh(
      emat, 0, ndofpernodemesh, diffcoeff, fac, rhsfac, J, gradphi, convelint, dJ_dmesh);

  // call base class routine again to compute linearizations of Ohmic overpotential w.r.t.
  // structural displacements
  my::calc_diff_od_mesh(emat, 1, ndofpernodemesh, var_manager()->InvF() * diff_manager()->GetCond(),
      fac, rhsfac, J, var_manager()->GradPot(), convelint, dJ_dmesh);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
void Discret::ELEMENTS::ScaTraEleCalcElchElectrode<distype, probdim>::calc_mat_diff_coeff_lin(
    Core::LinAlg::SerialDenseMatrix& emat, const int k, const double timefacfac,
    const Core::LinAlg::Matrix<nsd_, 1>& gradphi, const double scalar)
{
  // linearization of diffusion coefficient in ionic diffusion term (transport equation):
  //
  // (nabla w, D(D(c)) nabla c)
  //
  for (unsigned vi = 0; vi < nen_; ++vi)
  {
    for (unsigned ui = 0; ui < nen_; ++ui)
    {
      double laplawfrhs_gradphi(0.);
      my::get_laplacian_weak_form_rhs(laplawfrhs_gradphi, gradphi, vi);

      emat(vi * my::numdofpernode_ + k, ui * my::numdofpernode_ + k) +=
          scalar * timefacfac * diff_manager()->get_conc_deriv_iso_diff_coef(k, k) *
          laplawfrhs_gradphi * my::funct_(ui);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
void Discret::ELEMENTS::ScaTraEleCalcElchElectrode<distype, probdim>::calc_mat_pot_equ_divi_ohm(
    Core::LinAlg::SerialDenseMatrix& emat, const double timefacfac, const double invf,
    const Core::LinAlg::Matrix<nsd_, 1>& gradpot, const double scalar)
{
  for (unsigned vi = 0; vi < nen_; ++vi)
  {
    for (unsigned ui = 0; ui < nen_; ++ui)
    {
      double laplawf(0.);
      my::get_laplacian_weak_form(laplawf, ui, vi);

      // linearization of the ohmic term
      //
      // (grad w, 1/F kappa D(grad pot))
      //
      emat(vi * my::numdofpernode_ + my::numscal_, ui * my::numdofpernode_ + my::numscal_) +=
          scalar * timefacfac * invf * diff_manager()->GetCond() * laplawf;

      for (int iscal = 0; iscal < my::numscal_; ++iscal)
      {
        double laplawfrhs_gradpot(0.);
        my::get_laplacian_weak_form_rhs(laplawfrhs_gradpot, gradpot, vi);

        // linearization of the ohmic term with respect to conductivity
        //
        // (grad w, 1/F kappa D(grad pot))
        //
        emat(vi * my::numdofpernode_ + my::numscal_, ui * my::numdofpernode_ + iscal) +=
            scalar * timefacfac * invf * diff_manager()->GetConcDerivCond(iscal) * my::funct_(ui) *
            laplawfrhs_gradpot;
      }
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
void Discret::ELEMENTS::ScaTraEleCalcElchElectrode<distype,
    probdim>::calc_rhs_conservative_part_of_convective_term(Core::LinAlg::SerialDenseVector& erhs,
    const int k, const double rhsfac, const double vdiv)
{
  double vrhs = rhsfac * my::scatravarmanager_->Phinp(k) * vdiv;
  for (unsigned vi = 0; vi < nen_; ++vi) erhs[vi * my::numdofpernode_ + k] -= vrhs * my::funct_(vi);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
void Discret::ELEMENTS::ScaTraEleCalcElchElectrode<distype, probdim>::calc_rhs_pot_equ_divi_ohm(
    Core::LinAlg::SerialDenseVector& erhs, const double rhsfac, const double invf,
    const Core::LinAlg::Matrix<nsd_, 1>& gradpot, const double scalar)
{
  for (unsigned vi = 0; vi < nen_; ++vi)
  {
    double laplawfrhs_gradpot(0.);
    my::get_laplacian_weak_form_rhs(laplawfrhs_gradpot, gradpot, vi);

    erhs[vi * my::numdofpernode_ + my::numscal_] -=
        scalar * rhsfac * invf * diff_manager()->GetCond() * laplawfrhs_gradpot;
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
void Discret::ELEMENTS::ScaTraEleCalcElchElectrode<distype, probdim>::get_material_params(
    const Core::Elements::Element* ele, std::vector<double>& densn, std::vector<double>& densnp,
    std::vector<double>& densam, double& visc, const int iquad)
{
  // get material
  Teuchos::RCP<const Core::Mat::Material> material = ele->Material();

  // evaluate electrode material
  if (material->MaterialType() == Core::Materials::m_electrode)
  {
    utils()->mat_electrode(
        material, var_manager()->Phinp(0), var_manager()->Temperature(), diff_manager());
  }
  else
    FOUR_C_THROW("Material type not supported!");
}


// template classes
// 1D elements
template class Discret::ELEMENTS::ScaTraEleCalcElchElectrode<Core::FE::CellType::line2, 1>;
template class Discret::ELEMENTS::ScaTraEleCalcElchElectrode<Core::FE::CellType::line2, 2>;
template class Discret::ELEMENTS::ScaTraEleCalcElchElectrode<Core::FE::CellType::line2, 3>;
template class Discret::ELEMENTS::ScaTraEleCalcElchElectrode<Core::FE::CellType::line3, 1>;

// 2D elements
template class Discret::ELEMENTS::ScaTraEleCalcElchElectrode<Core::FE::CellType::tri3, 2>;
template class Discret::ELEMENTS::ScaTraEleCalcElchElectrode<Core::FE::CellType::tri3, 3>;
template class Discret::ELEMENTS::ScaTraEleCalcElchElectrode<Core::FE::CellType::tri6, 2>;
template class Discret::ELEMENTS::ScaTraEleCalcElchElectrode<Core::FE::CellType::quad4, 2>;
template class Discret::ELEMENTS::ScaTraEleCalcElchElectrode<Core::FE::CellType::quad4, 3>;
template class Discret::ELEMENTS::ScaTraEleCalcElchElectrode<Core::FE::CellType::quad9, 2>;
template class Discret::ELEMENTS::ScaTraEleCalcElchElectrode<Core::FE::CellType::nurbs9, 2>;

// 3D elements
template class Discret::ELEMENTS::ScaTraEleCalcElchElectrode<Core::FE::CellType::hex8, 3>;
template class Discret::ELEMENTS::ScaTraEleCalcElchElectrode<Core::FE::CellType::hex27, 3>;
template class Discret::ELEMENTS::ScaTraEleCalcElchElectrode<Core::FE::CellType::tet4, 3>;
template class Discret::ELEMENTS::ScaTraEleCalcElchElectrode<Core::FE::CellType::tet10, 3>;
template class Discret::ELEMENTS::ScaTraEleCalcElchElectrode<Core::FE::CellType::pyramid5, 3>;

FOUR_C_NAMESPACE_CLOSE

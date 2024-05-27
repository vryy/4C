/*--------------------------------------------------------------------------*/
/*! \file

\brief evaluation of scatra elements for isothermal space charge layer formation

\level 2

*/
/*--------------------------------------------------------------------------*/
#include "4C_scatra_ele_calc_elch_scl.hpp"

#include "4C_lib_discret.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_material_base.hpp"
#include "4C_scatra_ele_parameter_std.hpp"
#include "4C_scatra_ele_parameter_timint.hpp"
#include "4C_scatra_ele_utils_elch_scl.hpp"
#include "4C_utils_singleton_owner.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
DRT::ELEMENTS::ScaTraEleCalcElchScl<distype, probdim>*
DRT::ELEMENTS::ScaTraEleCalcElchScl<distype, probdim>::Instance(
    const int numdofpernode, const int numscal, const std::string& disname)
{
  static auto singleton_map = CORE::UTILS::MakeSingletonMap<std::string>(
      [](const int numdofpernode, const int numscal, const std::string& disname)
      {
        return std::unique_ptr<ScaTraEleCalcElchScl<distype, probdim>>(
            new ScaTraEleCalcElchScl<distype, probdim>(numdofpernode, numscal, disname));
      });

  return singleton_map[disname].Instance(
      CORE::UTILS::SingletonAction::create, numdofpernode, numscal, disname);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
DRT::ELEMENTS::ScaTraEleCalcElchScl<distype, probdim>::ScaTraEleCalcElchScl(
    const int numdofpernode, const int numscal, const std::string& disname)
    : DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype, probdim>::ScaTraEleCalcElchDiffCond(
          numdofpernode, numscal, disname),
      diffcondmat_(INPAR::ELCH::diffcondmat_undefined),
      diffcondparams_(DRT::ELEMENTS::ScaTraEleParameterElchDiffCond::Instance(disname))
{
  // replace diffusion manager for diffusion-conduciton formulation by diffusion manager for SCLs
  my::diffmanager_ = Teuchos::rcp(new ScaTraEleDiffManagerElchScl(my::numscal_));

  // replace internal variable manager for diffusion-conduction by internal variable manager for
  // SCL formulation
  my::scatravarmanager_ =
      Teuchos::rcp(new ScaTraEleInternalVariableManagerElchScl<my::nsd_, my::nen_>(
          my::numscal_, myelch::elchparams_, diffcondparams_));

  // replace utility class for diffusion-conduction formulation by utility class for SCLs
  myelch::utils_ =
      DRT::ELEMENTS::ScaTraEleUtilsElchScl<distype>::Instance(numdofpernode, numscal, disname);

  // safety checks for stabilization settings
  if (my::scatrapara_->StabType() != INPAR::SCATRA::stabtype_no_stabilization or
      my::scatrapara_->TauDef() != INPAR::SCATRA::tau_zero)
  {
    FOUR_C_THROW(
        "No stabilization available for the diffusion-conduction formulation, since we had no "
        "problems so far.");
  }
  if (not my::scatrapara_->MatGP() or not my::scatrapara_->TauGP())
  {
    FOUR_C_THROW(
        "Since most of the materials of the Diffusion-conduction formulation depend on the "
        "concentration, an evaluation of the material and the stabilization parameter at the "
        "element center is disabled.");
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
double DRT::ELEMENTS::ScaTraEleCalcElchScl<distype, probdim>::calc_free_charge(
    const double concentration)
{
  return diff_manager()->GetValence(0) * myelch::elchparams_->Faraday() *
         (concentration - diff_manager()->GetBulkConc());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
double DRT::ELEMENTS::ScaTraEleCalcElchScl<distype, probdim>::calc_free_charge_der_conc()
{
  return diff_manager()->GetValence(0) * myelch::elchparams_->Faraday();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElchScl<distype, probdim>::calc_mat_pot_coulomb(
    CORE::LINALG::SerialDenseMatrix& emat, const double fac, const double invf,
    const double scalefac, const CORE::LINALG::Matrix<my::nsd_, 1>& gradpot, const double epsilon)
{
  for (unsigned vi = 0; vi < my::nen_; ++vi)
  {
    for (unsigned ui = 0; ui < my::nen_; ++ui)
    {
      double laplawf(0.);
      my::get_laplacian_weak_form(laplawf, ui, vi);

      // linearization of the ohmic term
      //
      // (grad w, -epsilon D(grad pot))
      emat(vi * my::numdofpernode_ + my::numscal_, ui * my::numdofpernode_ + my::numscal_) +=
          fac * invf * scalefac * epsilon * laplawf;
    }
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElchScl<distype, probdim>::calc_rhs_pot_coulomb(
    CORE::LINALG::SerialDenseVector& erhs, const double fac, const double invf,
    const double cond_invperm, const CORE::LINALG::Matrix<my::nsd_, 1>& gradpot,
    const double epsilon)
{
  for (unsigned vi = 0; vi < my::nen_; ++vi)
  {
    double laplawfrhs_gradpot(0.);
    my::get_laplacian_weak_form_rhs(laplawfrhs_gradpot, gradpot, vi);

    erhs[vi * my::numdofpernode_ + my::numscal_] -=
        fac * invf * cond_invperm * epsilon * laplawfrhs_gradpot;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElchScl<distype, probdim>::calc_mat_pot_src(
    CORE::LINALG::SerialDenseMatrix& emat, const int k, const double timefacfac, const double invf,
    const double cond_invperm, const double z_k_F)
{
  for (unsigned vi = 0; vi < my::nen_; ++vi)
  {
    for (unsigned ui = 0; ui < my::nen_; ++ui)
    {
      emat(vi * my::numdofpernode_ + my::numscal_, ui * my::numdofpernode_ + k) +=
          -z_k_F * timefacfac * invf * cond_invperm * my::funct_(vi) * my::funct_(ui);
    }
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElchScl<distype, probdim>::calc_rhs_pot_src(
    CORE::LINALG::SerialDenseVector& erhs, const int k, const double fac, const double invf,
    const double cond_invperm, const double q_F)
{
  for (unsigned vi = 0; vi < my::nen_; ++vi)
  {
    erhs[vi * my::numdofpernode_ + my::numscal_] -=
        -fac * invf * cond_invperm * my::funct_(vi) * q_F;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElchScl<distype, probdim>::calc_rhs_diff_cur(
    CORE::LINALG::SerialDenseVector& erhs, const double rhsfac, const std::vector<double>& invfval,
    const std::vector<CORE::LINALG::Matrix<my::nsd_, 1>>& gradphi)
{
  if (diffcondmat_ != INPAR::ELCH::diffcondmat_scl)
    FOUR_C_THROW("Diffusion-Conduction material has to be SCL material");

  for (unsigned vi = 0; vi < my::nen_; ++vi)
  {
    for (unsigned idim = 0; idim < my::nsd_; ++idim)
    {
      for (int k = 0; k < my::numscal_; ++k)
      {
        erhs[vi * my::numdofpernode_ + (my::numscal_ + 1) + idim] -=
            rhsfac * diff_manager()->GetPhasePoroTort(0) * my::funct_(vi) *
            diff_manager()->GetIsotropicDiff(k) * gradphi[k](idim);
      }
    }
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElchScl<distype, probdim>::calc_mat_diff_cur(
    CORE::LINALG::SerialDenseMatrix& emat, const double timefacfac,
    const std::vector<double>& invfval,
    const std::vector<CORE::LINALG::Matrix<my::nsd_, 1>>& gradphi)
{
  for (unsigned vi = 0; vi < my::nen_; ++vi)
  {
    for (unsigned ui = 0; ui < my::nen_; ++ui)
    {
      // diffusive term
      // (grad w, D grad c)
      for (unsigned idim = 0; idim < my::nsd_; ++idim)
      {
        for (int k = 0; k < my::numscal_; ++k)
        {
          //  - D nabla c
          emat(vi * my::numdofpernode_ + (my::numscal_ + 1) + idim, ui * my::numdofpernode_ + k) +=
              timefacfac * diff_manager()->GetPhasePoroTort(0) * my::funct_(vi) *
              diff_manager()->GetIsotropicDiff(k) * my::derxy_(idim, ui);

          // linearization wrt DiffCoeff
          emat(vi * my::numdofpernode_ + (my::numscal_ + 1) + idim, ui * my::numdofpernode_ + k) +=
              timefacfac * diff_manager()->GetPhasePoroTort(0) *
              diff_manager()->get_conc_deriv_iso_diff_coef(k, k) * my::funct_(vi) *
              (gradphi[k])(idim)*my::funct_(ui);
        }
      }
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElchScl<distype, probdim>::calc_mat_and_rhs(
    CORE::LINALG::SerialDenseMatrix& emat, CORE::LINALG::SerialDenseVector& erhs, const int k,
    const double fac, const double timefacfac, const double rhsfac, const double taufac,
    const double timetaufac, const double rhstaufac, CORE::LINALG::Matrix<my::nen_, 1>& tauderpot,
    double& rhsint)
{
  //----------------------------------------------------------------
  // 1) element matrix: instationary terms
  //----------------------------------------------------------------
  if (not my::scatraparatimint_->IsStationary())
    my::calc_mat_mass(emat, k, fac, diff_manager()->GetPhasePoro(0));
  //----------------------------------------------------------------
  // 2) element matrix: stationary terms of ion-transport equation
  //----------------------------------------------------------------
  // 2b)  element matrix: diffusion term

  // current is not a solution variable
  if (not diffcondparams_->CurSolVar())
  {
    // i)  constant diffusion coefficient
    my::calc_mat_diff(emat, k, timefacfac * diff_manager()->GetPhasePoroTort(0));

    // ii) concentration depending diffusion coefficient
    mydiffcond::calc_mat_diff_coeff_lin(
        emat, k, timefacfac, var_manager()->GradPhi(k), diff_manager()->GetPhasePoroTort(0));


    // 2d) electrical conduction term (transport equation)
    //     i)  conduction term + ohmic overpotential
    //         (w_k, - t_k kappa nabla phi /(z_k F)) , transference number: const., unity
    mydiffcond::calc_mat_cond_ohm(
        emat, k, timefacfac, diff_manager()->InvFVal(k), var_manager()->GradPot());
  }
  // equation for current is solved independently: our case!!!
  else if (diffcondparams_->CurSolVar())
  // dc/dt + nabla N = 0
  {
    // current term (with current as a solution variable)
    mydiffcond::calc_mat_cond(
        emat, k, timefacfac, diff_manager()->InvFVal(k), var_manager()->CurInt());
  }

  //---------------------------------------------------------------------
  // 3)   governing equation for the electric potential field and free current
  //---------------------------------------------------------------------
  // see function calc_mat_and_rhs_outside_scalar_loop()

  //-----------------------------------------------------------------------
  // 4) element right hand side vector (neg. residual of nonlinear problem)
  //-----------------------------------------------------------------------
  if (my::scatraparatimint_->IsIncremental() and not my::scatraparatimint_->IsStationary())
  {
    my::calc_rhs_lin_mass(
        erhs, k, rhsfac, fac, diff_manager()->GetPhasePoro(0), diff_manager()->GetPhasePoro(0));
  }

  // adaption of rhs with respect to time integration: no sources
  // Evaluation at Gauss Points before spatial integration
  my::compute_rhs_int(rhsint, mydiffcond ::diff_manager()->GetPhasePoro(0),
      diff_manager()->GetPhasePoro(0), var_manager()->Hist(k));

  // add RHS and history contribution
  // Integrate RHS (@n) over element volume ==> total impact of timestep n
  my::calc_rhs_hist_and_source(erhs, k, fac, rhsint);

  if (not diffcondparams_->CurSolVar())  // not utilized in previous investigations
  {
    // diffusion term
    my::calc_rhs_diff(erhs, k, rhsfac * diff_manager()->GetPhasePoroTort(0));

    // electrical conduction term (transport equation)
    // equation for current is inserted in the mass transport equation
    mydiffcond::calc_rhs_cond_ohm(
        erhs, k, rhsfac, diff_manager()->InvFVal(k), var_manager()->GradPot());
  }
  // equation for current is solved independently: free current density!
  // nabla dot (i/z_k F)
  else if (diffcondparams_->CurSolVar())
  {
    // curint: current density at GP, InvFVal(k): 1/(z_k F)
    mydiffcond::calc_rhs_cond(erhs, k, rhsfac, diff_manager()->InvFVal(k), var_manager()->CurInt());
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElchScl<distype, probdim>::calc_mat_and_rhs_outside_scalar_loop(
    CORE::LINALG::SerialDenseMatrix& emat, CORE::LINALG::SerialDenseVector& erhs, const double fac,
    const double timefacfac, const double rhsfac)
{
  //----------------------------------------------------------------
  // 3)   governing equation for the electric potential field
  //----------------------------------------------------------------
  if (not diffcondparams_->CurSolVar())
  {
    // 3c) Laplace equation: nabla^2 Phi + sum (F z_k c_k) = 0

    // i) eps nabla^2 Phi = 0: Matrix
    calc_mat_pot_coulomb(emat, timefacfac, var_manager()->InvF(),
        diff_manager()->GetCond() / diff_manager()->GetPermittivity(), var_manager()->GradPot(),
        diff_manager()->GetPermittivity());

    //  RHS
    calc_rhs_pot_coulomb(erhs, rhsfac, var_manager()->InvF(),
        diff_manager()->GetCond() / diff_manager()->GetPermittivity(), var_manager()->GradPot(),
        diff_manager()->GetPermittivity());

    // ii) -sum (F z_k c_k) = 0 (use of defined function);

    // set this to zero (only laplace equation zero charge) ==> linear function
    for (int k = 0; k < my::numscal_; ++k)
    {
      calc_mat_pot_src(emat, k, timefacfac, var_manager()->InvF(),
          diff_manager()->GetCond() / diff_manager()->GetPermittivity(),
          calc_free_charge_der_conc());

      calc_rhs_pot_src(erhs, k, rhsfac, var_manager()->InvF(),
          diff_manager()->GetCond() / diff_manager()->GetPermittivity(),
          calc_free_charge(var_manager()->Phinp(k)));
    }
  }

  // 3c) Laplace equation based on free charge
  // i_F/(z_k F) = N+, ion flux!
  // equation for current is solved independently
  else if (diffcondparams_->CurSolVar())
  {
    //-----------------------------------------------------------------------
    // 5) equation for the current incl. rhs-terms
    //-----------------------------------------------------------------------

    // matrix terms
    // (xsi_i,Di)
    mydiffcond::calc_mat_cur_equ_cur(emat, timefacfac, var_manager()->InvF());

    // (xsi, -D(kappa phi))
    mydiffcond::calc_mat_cur_equ_ohm(
        emat, timefacfac, var_manager()->InvF(), var_manager()->GradPot());

    // (xsi, -D(z_k F D (c) nabla c)
    calc_mat_diff_cur(emat, timefacfac, diff_manager()->InvFVal(), var_manager()->GradPhi());

    // (xsi_i,Di): stays the same
    mydiffcond::calc_rhs_cur_equ_cur(erhs, rhsfac, var_manager()->InvF(), var_manager()->CurInt());

    // (xsi, -D(kappa phi)): stays the same, but local version of conductivity
    mydiffcond::calc_rhs_cur_equ_ohm(erhs, rhsfac, var_manager()->InvF(), var_manager()->GradPot());

    // (xsi, - D(z_k F D(c) nabla c)
    calc_rhs_diff_cur(erhs, rhsfac, diff_manager()->InvFVal(), var_manager()->GradPhi());

    //------------------------------------------------------------------------------------------
    // 3)   governing equation for the electric potential field and current (incl. rhs-terms)
    //------------------------------------------------------------------------------------------

    // i) eps nabla^2 Phi = 0: Matrix
    calc_mat_pot_coulomb(emat, timefacfac, var_manager()->InvF(),
        diff_manager()->GetCond() / diff_manager()->GetPermittivity(), var_manager()->GradPot(),
        diff_manager()->GetPermittivity());

    //  RHS
    calc_rhs_pot_coulomb(erhs, rhsfac, var_manager()->InvF(),
        diff_manager()->GetCond() / diff_manager()->GetPermittivity(), var_manager()->GradPot(),
        diff_manager()->GetPermittivity());
    // ii) -sum (F z_k c_k) = 0

    // set this to zero (only laplace equation zero charge) ==> linear function
    for (int k = 0; k < my::numscal_; ++k)
    {
      calc_mat_pot_src(emat, k, timefacfac, var_manager()->InvF(),
          diff_manager()->GetCond() / diff_manager()->GetPermittivity(),
          calc_free_charge_der_conc());

      calc_rhs_pot_src(erhs, k, rhsfac, var_manager()->InvF(),
          diff_manager()->GetCond() / diff_manager()->GetPermittivity(),
          calc_free_charge(var_manager()->Phinp(k)));
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElchScl<distype, probdim>::get_material_params(
    const DRT::Element* ele, std::vector<double>& densn, std::vector<double>& densnp,
    std::vector<double>& densam, double& visc, const int iquad)
{
  // extract material from element
  Teuchos::RCP<CORE::MAT::Material> material = ele->Material();

  // evaluate electrolyte material
  if (material->MaterialType() == CORE::Materials::m_elchmat)
  {
    utils()->MatElchMat(material, var_manager()->Phinp(), var_manager()->Temperature(),
        diff_manager(), diffcondmat_);
  }
  else
    FOUR_C_THROW("Invalid material type!");
}

// template classes
// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalcElchScl<CORE::FE::CellType::line2, 1>;
template class DRT::ELEMENTS::ScaTraEleCalcElchScl<CORE::FE::CellType::line2, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcElchScl<CORE::FE::CellType::line2, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcElchScl<CORE::FE::CellType::line3, 1>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalcElchScl<CORE::FE::CellType::tri3, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcElchScl<CORE::FE::CellType::tri3, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcElchScl<CORE::FE::CellType::tri6, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcElchScl<CORE::FE::CellType::quad4, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcElchScl<CORE::FE::CellType::quad4, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcElchScl<CORE::FE::CellType::quad9, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcElchScl<CORE::FE::CellType::nurbs9, 2>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcElchScl<CORE::FE::CellType::hex8, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcElchScl<CORE::FE::CellType::hex27, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcElchScl<CORE::FE::CellType::tet4, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcElchScl<CORE::FE::CellType::tet10, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcElchScl<CORE::FE::CellType::pyramid5, 3>;

FOUR_C_NAMESPACE_CLOSE

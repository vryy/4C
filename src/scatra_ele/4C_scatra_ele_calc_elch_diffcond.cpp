/*--------------------------------------------------------------------------*/
/*! \file

\brief evaluation of scatra elements for isothermal diffusion-conduction ion-transport equations

\level 2

*/
/*--------------------------------------------------------------------------*/
#include "4C_scatra_ele_calc_elch_diffcond.hpp"

#include "4C_discretization_fem_general_extract_values.hpp"
#include "4C_lib_discret.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_material_base.hpp"
#include "4C_scatra_ele_parameter_std.hpp"
#include "4C_scatra_ele_parameter_timint.hpp"
#include "4C_scatra_ele_utils_elch_diffcond.hpp"
#include "4C_utils_singleton_owner.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype, probdim>*
DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype, probdim>::Instance(
    const int numdofpernode, const int numscal, const std::string& disname)
{
  static auto singleton_map = CORE::UTILS::MakeSingletonMap<std::string>(
      [](const int numdofpernode, const int numscal, const std::string& disname)
      {
        return std::unique_ptr<ScaTraEleCalcElchDiffCond<distype, probdim>>(
            new ScaTraEleCalcElchDiffCond<distype, probdim>(numdofpernode, numscal, disname));
      });

  return singleton_map[disname].Instance(
      CORE::UTILS::SingletonAction::create, numdofpernode, numscal, disname);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype, probdim>::ScaTraEleCalcElchDiffCond(
    const int numdofpernode, const int numscal, const std::string& disname)
    : DRT::ELEMENTS::ScaTraEleCalcElchElectrode<distype, probdim>::ScaTraEleCalcElchElectrode(
          numdofpernode, numscal, disname),
      diffcondmat_(INPAR::ELCH::diffcondmat_undefined),
      diffcondparams_(DRT::ELEMENTS::ScaTraEleParameterElchDiffCond::Instance(disname))
{
  // replace diffusion manager for electrodes by diffusion manager for diffusion-conduction
  // formulation
  my::diffmanager_ = Teuchos::rcp(new ScaTraEleDiffManagerElchDiffCond(my::numscal_));

  // replace internal variable manager for electrodes by internal variable manager for
  // diffusion-conduction formulation
  my::scatravarmanager_ = Teuchos::rcp(new ScaTraEleInternalVariableManagerElchDiffCond<nsd_, nen_>(
      my::numscal_, myelch::elchparams_, diffcondparams_));

  // replace utility class for electrodes by utility class for diffusion-conduction formulation
  myelch::utils_ =
      DRT::ELEMENTS::ScaTraEleUtilsElchDiffCond<distype>::Instance(numdofpernode, numscal, disname);

  // safety check for closing equation
  switch (myelch::elchparams_->EquPot())
  {
    case INPAR::ELCH::equpot_divi:
    case INPAR::ELCH::equpot_enc:
      // valid closing equations for electric potential
      break;
    default:
    {
      FOUR_C_THROW("Invalid closing equation for electric potential!");
      break;
    }
  }

  // safety checks for stabilization settings
  if (my::scatrapara_->StabType() != INPAR::SCATRA::stabtype_no_stabilization or
      my::scatrapara_->TauDef() != INPAR::SCATRA::tau_zero)
  {
    FOUR_C_THROW(
        "No stabilization available for the diffusion-conduction formulation, since we had no "
        "problems so far.");
  }
  if (my::scatrapara_->MatGP() == false or my::scatrapara_->TauGP() == false)
  {
    FOUR_C_THROW(
        "Since most of the materials of the Diffusion-conduction formulation depend on the "
        "concentration,\n"
        "an evaluation of the material and the stabilization parameter at the element center is "
        "disabled.");
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype, probdim>::calc_mat_and_rhs(
    CORE::LINALG::SerialDenseMatrix& emat, CORE::LINALG::SerialDenseVector& erhs, const int k,
    const double fac, const double timefacfac, const double rhsfac, const double taufac,
    const double timetaufac, const double rhstaufac, CORE::LINALG::Matrix<nen_, 1>& tauderpot,
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

  // 2a)  element matrix: convective term
  my::calc_mat_conv(emat, k, timefacfac, diff_manager()->GetPhasePoro(0), var_manager()->SGConv());

  // 2b)  element matrix: diffusion term
  //      i)  constant diffusion coefficient
  my::calc_mat_diff(emat, k, timefacfac * diff_manager()->GetPhasePoroTort(0));

  //      ii) concentration depending diffusion coefficient
  //          (additional term for Newman material)
  if (diffcondmat_ == INPAR::ELCH::diffcondmat_newman)
    myelectrode::calc_mat_diff_coeff_lin(
        emat, k, timefacfac, var_manager()->GradPhi(k), diff_manager()->GetPhasePoroTort(0));

  // 2c) element matrix: conservative part of convective term, needed for deforming bodies,
  //                     i.e., for scalar-structure interaction
  double velocity_divergence(0.0);
  if (my::scatrapara_->IsConservative())
  {
    my::get_divergence(velocity_divergence, my::evelnp_);
    my::calc_mat_conv_add_cons(
        emat, k, timefacfac * diff_manager()->GetPhasePoro(0), velocity_divergence, 1.0);
  }

  // 2d) electrical conduction term (transport equation)
  //
  //     mass transport equation:
  //
  //               |     diffusion term      | |     conduction term    |
  //               |                         | |                        |
  //      dc_k/dt - nabla dot (D_k nabla c_k) + nabla dot (t_k i/(z_k F))
  //
  //    equation for current (based on Newman):
  //
  //          | ohmic overpot.   |         concentration overpotential           |
  //          |                  |                                               |
  //      i = - kappa nabla phi  + RT/F kappa (thermfactor) f(t_k) nabla ln c_k

  // equation for current is inserted in the mass transport equation
  if (not diffcondparams_->CurSolVar())
  {
    //    i)  conduction term + ohmic overpotential
    //        (w_k, - t_k kappa nabla phi /(z_k F))
    calc_mat_cond_ohm(emat, k, timefacfac, diff_manager()->InvFVal(k), var_manager()->GradPot());

    //    ii) conduction term + concentration overpotential
    //        (w_k, - t_k RT/F kappa (thermfactor) f(t_k) nabla ln c_k /(z_k F))
    if (diffcondmat_ == INPAR::ELCH::diffcondmat_newman)
    {
      calc_mat_cond_conc(emat, k, timefacfac,
          var_manager()->RTFFC() / diff_manager()->GetValence(k), diffcondparams_->NewmanConstA(),
          diffcondparams_->NewmanConstB(), var_manager()->GradPhi(k), var_manager()->ConIntInv());
    }
  }
  // equation for current is solved independently
  else
  {
    // current term (with current as a solution variable)
    calc_mat_cond(emat, k, timefacfac, diff_manager()->InvFVal(k), var_manager()->CurInt());

    // this coupling term cancels out for a 2 equation system
    if (diffcondmat_ == INPAR::ELCH::diffcondmat_ion)
      calc_mat_cond_diff(emat, k, timefacfac, diff_manager()->InvFVal(k), var_manager()->GradPhi());
  }  // end if(not diffcondparams_->CurSolVar())

  //---------------------------------------------------------------------
  // 3)   governing equation for the electric potential field and current
  //---------------------------------------------------------------------
  // see function calc_mat_and_rhs_outside_scalar_loop()

  //-----------------------------------------------------------------------
  // 4) element right hand side vector (neg. residual of nonlinear problem)
  //-----------------------------------------------------------------------

  if (my::scatraparatimint_->IsIncremental() and not my::scatraparatimint_->IsStationary())
    my::calc_rhs_lin_mass(
        erhs, k, rhsfac, fac, diff_manager()->GetPhasePoro(0), diff_manager()->GetPhasePoro(0));

  // adaption of rhs with respect to time integration
  my::compute_rhs_int(rhsint, diff_manager()->GetPhasePoro(0), diff_manager()->GetPhasePoro(0),
      var_manager()->Hist(k));

  // add RHS and history contribution
  my::calc_rhs_hist_and_source(erhs, k, fac, rhsint);

  // convective term
  my::calc_rhs_conv(erhs, k, rhsfac * diff_manager()->GetPhasePoro(0));

  // diffusion term
  my::calc_rhs_diff(erhs, k, rhsfac * diff_manager()->GetPhasePoroTort(0));

  // conservative part of convective term, needed for deforming bodies, i.e., for scalar-structure
  // interaction
  if (my::scatrapara_->IsConservative())
    myelectrode::calc_rhs_conservative_part_of_convective_term(
        erhs, k, rhsfac * diff_manager()->GetPhasePoro(0), velocity_divergence);

  // electrical conduction term (transport equation)
  // equation for current is inserted in the mass transport equation
  //
  // mass transport equation:
  //
  //               |     diffusion term      | |     conduction term    |
  //               |                         | |                        |
  //      dc_k/dt - nabla dot (D_k nabla c_k) + nabla dot (t_k i/(z_k F))
  //
  // equation for current:
  //
  //          | ohmic overpot.   |         concentration overpotential           |
  //          |                  |                                               |
  //      i = - kappa nabla phi  + RT/F kappa (thermfactor) f(t_k) nabla ln c_k
  if (not diffcondparams_->CurSolVar())
  {
    calc_rhs_cond_ohm(erhs, k, rhsfac, diff_manager()->InvFVal(k), var_manager()->GradPot());

    // if(diffcondmat_==INPAR::ELCH::diffcondmat_ion): all terms cancel out
    if (diffcondmat_ == INPAR::ELCH::diffcondmat_newman)
    {
      calc_rhs_cond_conc(erhs, k, rhsfac, var_manager()->RTFFC() / diff_manager()->GetValence(k),
          diffcondparams_->NewmanConstA(), diffcondparams_->NewmanConstB(),
          var_manager()->GradPhi(k), var_manager()->ConIntInv());
    }
  }

  // equation for current is solved independently
  else
  {
    calc_rhs_cond(erhs, k, rhsfac, diff_manager()->InvFVal(k), var_manager()->CurInt());

    if (diffcondmat_ == INPAR::ELCH::diffcondmat_ion)
      calc_rhs_cond_diff(erhs, k, rhsfac, var_manager()->GradPhi());
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype,
    probdim>::calc_mat_and_rhs_outside_scalar_loop(CORE::LINALG::SerialDenseMatrix& emat,
    CORE::LINALG::SerialDenseVector& erhs, const double fac, const double timefacfac,
    const double rhsfac)
{
  //----------------------------------------------------------------
  // 3)   governing equation for the electric potential field
  //----------------------------------------------------------------
  //
  // mass transport equation:
  //            |     diffusion term      | |     conduction term    |
  //            |                         | |                        |
  //   dc_k/dt - nabla dot (D_k nabla c_k) + nabla dot (t_k i/(z_k F)
  //
  // equation for current:
  //   i = - kappa nabla phi + RT/F kappa (thermfactor) f(t_k) nabla ln c_k
  //
  // equation for potential:
  //
  // a) nabla cdot i = 0
  // b) ENC

  // // equation for current is NOT solved independently
  if (not diffcondparams_->CurSolVar())
  {
    // equation for current is inserted in the mass transport equation
    // 3a)  nabla cdot i = 0
    if (myelch::elchparams_->EquPot() == INPAR::ELCH::equpot_divi)
    {
      //  i)  ohmic overpotential (implemented after the scalar loop)
      //      (w_k, - kappa nabla phi)
      myelectrode::calc_mat_pot_equ_divi_ohm(emat, timefacfac, var_manager()->InvF(),
          var_manager()->GradPot(), diff_manager()->GetPhasePoroTort(0));

      //
      myelectrode::calc_rhs_pot_equ_divi_ohm(erhs, rhsfac, var_manager()->InvF(),
          var_manager()->GradPot(), diff_manager()->GetPhasePoroTort(0));

      //  ii)  concentration overpotential
      //      (w_k, RT/F kappa (thermfactor) f(t_k) nabla ln c_k)
      for (int k = 0; k < my::numscal_; ++k)
      {
        //
        calc_mat_pot_equ_divi_conc(emat, k, timefacfac, var_manager()->RTFFC(),
            var_manager()->RTF(), var_manager()->InvF(), diffcondparams_->NewmanConstA(),
            diffcondparams_->NewmanConstB(), var_manager()->GradPhi(k),
            var_manager()->ConIntInv(k));

        //
        calc_rhs_pot_equ_divi_conc(erhs, k, rhsfac, var_manager()->RTF(), diff_manager()->InvFVal(),
            var_manager()->RTFFC(), diffcondparams_->NewmanConstA(),
            diffcondparams_->NewmanConstB(), var_manager()->GradPhi(k),
            var_manager()->ConIntInv(k));
      }
    }
    // 3b)  ENC
    else if (myelch::elchparams_->EquPot() == INPAR::ELCH::equpot_enc)
    {
      for (int k = 0; k < my::numscal_; ++k)
      {
        //
        myelch::calc_mat_pot_equ_enc(emat, k, fac, my::scatraparatimint_->AlphaF());

        //
        myelch::calc_rhs_pot_equ_enc(erhs, k, fac, var_manager()->Phinp(k));
      }
    }
    else
      FOUR_C_THROW("(div i, ENC) are the options available in the Diffusion-Conduction framework");
  }
  // equation for current is solved independently
  else
  {
    //-----------------------------------------------------------------------
    // 5) equation for the current incl. rhs-terms
    //-----------------------------------------------------------------------

    //   | cur  | ohmic overpot.   |         concentration overpotential           |
    //          |                  |                                               |
    //      i = - kappa nabla phi  + RT/F kappa (thermfactor) f(t_k) nabla ln c_k

    // matrix terms
    // (xsi_i,Di)
    calc_mat_cur_equ_cur(emat, timefacfac, var_manager()->InvF());

    // (xsi, -D(kappa phi))
    calc_mat_cur_equ_ohm(emat, timefacfac, var_manager()->InvF(), var_manager()->GradPot());

    // (xsi, -D(RT/F kappa (thermfactor) f(t_k) nabla ln c_k))
    calc_mat_cur_equ_conc(emat, timefacfac, var_manager()->RTF(), var_manager()->RTFFC(),
        diff_manager()->InvFVal(), diffcondparams_->NewmanConstA(), diffcondparams_->NewmanConstB(),
        var_manager()->GradPhi(), var_manager()->ConIntInv());

    // (xsi_i,Di)
    calc_rhs_cur_equ_cur(erhs, rhsfac, var_manager()->InvF(), var_manager()->CurInt());

    // (xsi, -D(kappa phi))
    calc_rhs_cur_equ_ohm(erhs, rhsfac, var_manager()->InvF(), var_manager()->GradPot());

    // (xsi, -D(RT/F kappa (thermfactor) f(t_k) nabla ln c_k))
    calc_rhs_cur_equ_conc(erhs, rhsfac, var_manager()->RTF(), diff_manager()->InvFVal(),
        var_manager()->RTFFC(), diffcondparams_->NewmanConstA(), diffcondparams_->NewmanConstB(),
        var_manager()->GradPhi(), var_manager()->ConIntInv());

    //------------------------------------------------------------------------------------------
    // 3)   governing equation for the electric potential field and current (incl. rhs-terms)
    //------------------------------------------------------------------------------------------

    if (myelch::elchparams_->EquPot() == INPAR::ELCH::equpot_divi)
    {
      calc_mat_pot_equ_divi(emat, timefacfac, var_manager()->InvF());

      calc_rhs_pot_equ_divi(erhs, rhsfac, var_manager()->InvF(), var_manager()->CurInt());
    }
    else if (myelch::elchparams_->EquPot() == INPAR::ELCH::equpot_enc)
    {
      for (int k = 0; k < my::numscal_; ++k)
      {
        //
        myelch::calc_mat_pot_equ_enc(emat, k, fac, my::scatraparatimint_->AlphaF());

        //
        myelch::calc_rhs_pot_equ_enc(erhs, k, fac, var_manager()->Phinp(k));
      }
    }
    else
      FOUR_C_THROW("(div i, ENC) are the options available in the Diffusion-Conduction framework");
  }
}

/*----------------------------------------------------------------------------------*
 *----------------------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype, probdim>::calc_mat_cond_ohm(
    CORE::LINALG::SerialDenseMatrix& emat, const int k, const double timefacfac,
    const double invfval, const CORE::LINALG::Matrix<nsd_, 1>& gradpot)
{
  for (unsigned vi = 0; vi < nen_; ++vi)
  {
    for (unsigned ui = 0; ui < nen_; ++ui)
    {
      double laplawf(0.0);
      my::get_laplacian_weak_form(laplawf, ui, vi);  // compute once, reuse below!

      // linearization of conduction term depending on the potential
      //
      // (grad w, t_k kappa/(F z_k) D(grad phi))
      //
      emat(vi * my::numdofpernode_ + k, ui * my::numdofpernode_ + my::numscal_) +=
          timefacfac * diff_manager()->GetPhasePoroTort(0) * diff_manager()->GetTransNum(k) *
          diff_manager()->GetCond() * invfval * laplawf;

      double laplawfrhs_gradpot(0.0);
      my::get_laplacian_weak_form_rhs(laplawfrhs_gradpot, gradpot, vi);

      for (int iscal = 0; iscal < my::numscal_; ++iscal)
      {
        // linearization of the conductivity in the conduction term depending on the potential
        //
        // (grad w, t_k D(kappa(c))/(F z_k) grad phi)
        //
        emat(vi * my::numdofpernode_ + k, ui * my::numdofpernode_ + iscal) +=
            timefacfac * diff_manager()->GetPhasePoroTort(0) * diff_manager()->GetTransNum(k) *
            invfval * diff_manager()->GetConcDerivCond(iscal) * my::funct_(ui) * laplawfrhs_gradpot;

        // linearization of the transference number in the conduction term depending on the
        // potential
        //
        // (grad w, D(t_k(c)) kappa/(F z_k) grad phi)
        //
        emat(vi * my::numdofpernode_ + k, ui * my::numdofpernode_ + iscal) +=
            timefacfac * diff_manager()->GetPhasePoroTort(0) *
            (diff_manager()->GetDerivTransNum(k, iscal)) * my::funct_(ui) * invfval *
            diff_manager()->GetCond() * laplawfrhs_gradpot;
      }
    }
  }
}

/*----------------------------------------------------------------------------------*
 *----------------------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype, probdim>::calc_mat_cond_conc(
    CORE::LINALG::SerialDenseMatrix& emat, const int k, const double timefacfac,
    const double rtffcval, const double newman_const_a, const double newman_const_b,
    const CORE::LINALG::Matrix<nsd_, 1>& gradphi, const std::vector<double>& conintinv)
{
  // additional safety check in the beginning for Newman materials
  if (k != 0)
    FOUR_C_THROW(
        "Material Newman is only valid for one scalar (binary electrolyte utilizing the ENC)");

  for (unsigned vi = 0; vi < nen_; ++vi)
  {
    for (unsigned ui = 0; ui < nen_; ++ui)
    {
      double laplawf(0.0);
      my::get_laplacian_weak_form(laplawf, ui, vi);

      // Material Newman is only valid for binary electrolyte utilizing the ENC:
      // -> the equations are solved only for one species (k=0)
      // -> all transport parameter only depend on a single species
      // -> all derivations of the transport parameters wrt this species
      //
      // additional safety check in the beginning of this material: k != 0
      // original safety check by method check_elch_element_parameter() in
      // scatra_ele_calc_service_elch.cpp

      // linearization of conduction term depending on the concentration
      //
      // (grad w, RT/(z_k F^2) kappa thermfac f(t_+) D(grad ln c_k))
      //
      emat(vi * my::numdofpernode_ + k, ui * my::numdofpernode_ + k) +=
          timefacfac * diff_manager()->GetPhasePoroTort(0) * rtffcval *
          diff_manager()->GetTransNum(k) * diff_manager()->GetCond() *
          diff_manager()->GetThermFac() *
          (newman_const_a + (newman_const_b * diff_manager()->GetTransNum(k))) * conintinv[k] *
          laplawf;

      // linearization of conduction term depending on the concentration is implemented
      // only for one species
      // otherwise you would need a second loop over the all scalars
      double laplawfrhs_gradc(0.0);
      my::get_laplacian_weak_form_rhs(laplawfrhs_gradc, gradphi, vi);

      // Linearization wrt ln c
      emat(vi * my::numdofpernode_ + k, ui * my::numdofpernode_ + k) +=
          -timefacfac * diff_manager()->GetPhasePoroTort(0) * rtffcval *
          diff_manager()->GetTransNum(k) * diff_manager()->GetCond() *
          diff_manager()->GetThermFac() *
          (newman_const_a + (newman_const_b * diff_manager()->GetTransNum(k))) * conintinv[k] *
          conintinv[k] * laplawfrhs_gradc * my::funct_(ui);

      // Linearization wrt kappa
      emat(vi * my::numdofpernode_ + k, ui * my::numdofpernode_ + k) +=
          timefacfac * diff_manager()->GetPhasePoroTort(0) * rtffcval *
          diff_manager()->GetTransNum(k) * diff_manager()->GetThermFac() *
          (newman_const_a + (newman_const_b * diff_manager()->GetTransNum(k))) * conintinv[k] *
          laplawfrhs_gradc * diff_manager()->GetConcDerivCond(k) * my::funct_(ui);

      // Linearization wrt transference number 1
      emat(vi * my::numdofpernode_ + k, ui * my::numdofpernode_ + k) +=
          timefacfac * diff_manager()->GetPhasePoroTort(0) * rtffcval * diff_manager()->GetCond() *
          conintinv[k] * laplawfrhs_gradc * diff_manager()->GetThermFac() *
          (newman_const_a + newman_const_b * diff_manager()->GetTransNum(k)) *
          diff_manager()->GetDerivTransNum(k, k) * my::funct_(ui);

      // Linearization wrt transference number 2
      emat(vi * my::numdofpernode_ + k, ui * my::numdofpernode_ + k) +=
          timefacfac * diff_manager()->GetPhasePoroTort(0) * rtffcval * diff_manager()->GetCond() *
          diff_manager()->GetThermFac() * conintinv[k] * laplawfrhs_gradc *
          diff_manager()->GetTransNum(k) * newman_const_b * diff_manager()->GetDerivTransNum(k, k) *
          my::funct_(ui);

      // Linearization wrt thermodynamic factor
      emat(vi * my::numdofpernode_ + k, ui * my::numdofpernode_ + k) +=
          timefacfac * diff_manager()->GetPhasePoroTort(0) * rtffcval *
          diff_manager()->GetTransNum(k) * diff_manager()->GetCond() *
          (newman_const_a + (newman_const_b * diff_manager()->GetTransNum(k))) * conintinv[k] *
          laplawfrhs_gradc * diff_manager()->GetDerivThermFac(k) * my::funct_(ui);
    }  // for ui
  }    // for vi
}

/*----------------------------------------------------------------------------------*
 *----------------------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype, probdim>::calc_mat_cond(
    CORE::LINALG::SerialDenseMatrix& emat, const int k, const double timefacfac,
    const double invfval, const CORE::LINALG::Matrix<nsd_, 1>& curint)
{
  for (unsigned vi = 0; vi < nen_; ++vi)
  {
    const int fvi = vi * my::numdofpernode_ + k;

    for (unsigned ui = 0; ui < nen_; ++ui)
    {
      for (unsigned idim = 0; idim < nsd_; ++idim)
      {
        const int fui = ui * my::numdofpernode_ + (my::numscal_ + 1) + idim;

        // linearization of conduction term depending on current flow
        //
        // (grad w, t_k/(F z_k) Di)
        //
        emat(fvi, fui) += -timefacfac * my::derxy_(idim, vi) * diff_manager()->GetTransNum(k) *
                          invfval * my::funct_(ui);

        // linearization of transference number in conduction term depending on current flow
        //
        // (grad w, Dt_k(c)/(F z_k) i)
        //
        for (int iscal = 0; iscal < my::numscal_; ++iscal)
        {
          emat(fvi, ui * my::numdofpernode_ + iscal) +=
              -timefacfac * my::derxy_(idim, vi) * (diff_manager()->GetDerivTransNum(k, iscal)) *
              my::funct_(ui) * invfval * curint(idim);
        }
      }
    }
  }
}

/*----------------------------------------------------------------------------------*
 *----------------------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype, probdim>::calc_mat_cond_diff(
    CORE::LINALG::SerialDenseMatrix& emat, const int k, const double timefacfac,
    const double invfval, const std::vector<CORE::LINALG::Matrix<nsd_, 1>>& gradphi)
{
  for (unsigned vi = 0; vi < nen_; ++vi)
  {
    for (unsigned ui = 0; ui < nen_; ++ui)
    {
      // compute once, reuse below!
      double laplawf(0.0);
      my::get_laplacian_weak_form(laplawf, ui, vi);

      for (int iscal = 0; iscal < my::numscal_; ++iscal)
      {
        // formulation a): plain ionic diffusion coefficients without using ENC
        //
        // (grad w, t_k/z_k*sum_i(D_i grad Dc))
        //
        emat(vi * my::numdofpernode_ + k, ui * my::numdofpernode_ + iscal) +=
            -timefacfac * diff_manager()->GetPhasePoroTort(0) * diff_manager()->GetTransNum(k) /
            diff_manager()->GetValence(k) * diff_manager()->GetValence(iscal) *
            diff_manager()->GetIsotropicDiff(iscal) * laplawf;
      }

      // linearization of transference number in the coupling term (transport equation)
      //
      // (grad w, Dt_k(c)/z_k (Sum_i z_i D_i grad c_i))
      //
      for (int iscal = 0; iscal < my::numscal_; ++iscal)
      {
        double term_vi = 0.0;
        for (int iscal2 = 0; iscal2 < my::numscal_; ++iscal2)
        {
          double laplawfrhs_gradphi = 0.0;
          my::get_laplacian_weak_form_rhs(laplawfrhs_gradphi, gradphi[iscal2], vi);

          term_vi += diff_manager()->GetValence(iscal2) * diff_manager()->GetIsotropicDiff(iscal2) *
                     laplawfrhs_gradphi;
        }

        emat(vi * my::numdofpernode_ + k, ui * my::numdofpernode_ + iscal) +=
            -timefacfac * diff_manager()->GetPhasePoroTort(0) *
            (diff_manager()->GetDerivTransNum(k, iscal)) * my::funct_(ui) /
            diff_manager()->GetValence(k) * term_vi;
      }  // for(iscal)
    }    // end for ui
  }      // end for vi
}

/*---------------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype, probdim>::calc_mat_pot_equ_divi_conc(
    CORE::LINALG::SerialDenseMatrix& emat, const int k, const double timefacfac, const double rtffc,
    const double rtf, const double invf, const double newman_const_a, const double newman_const_b,
    const CORE::LINALG::Matrix<nsd_, 1>& gradphi, const double conintinv)
{
  for (unsigned vi = 0; vi < nen_; ++vi)
  {
    for (unsigned ui = 0; ui < nen_; ++ui)
    {
      double laplawf(0.0);
      my::get_laplacian_weak_form(laplawf, ui, vi);

      if (diffcondmat_ == INPAR::ELCH::diffcondmat_newman)
      {
        // linearization of the diffusion overpotential term
        //
        // (grad w, RT/F^2 kappa (thermfactor) f(t_k) 1/c_k D nabla c_k)
        //
        emat(vi * my::numdofpernode_ + my::numscal_, ui * my::numdofpernode_ + k) +=
            timefacfac * diff_manager()->GetPhasePoroTort(0) * rtffc * diff_manager()->GetCond() *
            diff_manager()->GetThermFac() *
            (newman_const_a + (newman_const_b * diff_manager()->GetTransNum(k))) * conintinv *
            laplawf;

        // linearization of conduction term depending on the concentration is implemented
        // only for one species
        // otherwise you would need a second loop over the all scalars
        double laplawfrhs_gradphi(0.0);
        my::get_laplacian_weak_form_rhs(laplawfrhs_gradphi, gradphi, vi);

        // Linearization wrt ln c
        emat(vi * my::numdofpernode_ + my::numscal_, ui * my::numdofpernode_ + k) +=
            -timefacfac * diff_manager()->GetPhasePoroTort(0) * rtffc * diff_manager()->GetCond() *
            diff_manager()->GetThermFac() *
            (newman_const_a + (newman_const_b * diff_manager()->GetTransNum(k))) * conintinv *
            conintinv * laplawfrhs_gradphi * my::funct_(ui);

        // Linearization wrt kappa
        emat(vi * my::numdofpernode_ + my::numscal_, ui * my::numdofpernode_ + k) +=
            timefacfac * diff_manager()->GetPhasePoroTort(0) * rtffc *
            diff_manager()->GetThermFac() *
            (newman_const_a + (newman_const_b * diff_manager()->GetTransNum(k))) * conintinv *
            laplawfrhs_gradphi * diff_manager()->GetConcDerivCond(k) * my::funct_(ui);

        // Linearization wrt transference number
        emat(vi * my::numdofpernode_ + my::numscal_, ui * my::numdofpernode_ + k) +=
            timefacfac * diff_manager()->GetPhasePoroTort(0) * rtffc * diff_manager()->GetCond() *
            diff_manager()->GetThermFac() * conintinv * laplawfrhs_gradphi * newman_const_b *
            diff_manager()->GetDerivTransNum(k, k) * my::funct_(ui);

        // Linearization wrt thermodynamic factor
        emat(vi * my::numdofpernode_ + my::numscal_, ui * my::numdofpernode_ + k) +=
            timefacfac * diff_manager()->GetPhasePoroTort(0) * rtffc * diff_manager()->GetCond() *
            (newman_const_a + (newman_const_b * diff_manager()->GetTransNum(k))) * conintinv *
            laplawfrhs_gradphi * diff_manager()->GetDerivThermFac(k) * my::funct_(ui);
      }
      else if (diffcondmat_ == INPAR::ELCH::diffcondmat_ion)
      {
        if (diffcondparams_->DiffusionCoeffBased())
        {
          emat(vi * my::numdofpernode_ + my::numscal_, ui * my::numdofpernode_ + k) +=
              timefacfac * diff_manager()->GetPhasePoroTort(0) * diff_manager()->GetValence(k) *
              diff_manager()->GetIsotropicDiff(k) * laplawf;
        }
        else
        {
          // Attention:
          // Full linearization of transference number, conductivity, ... is still missing
          emat(vi * my::numdofpernode_ + my::numscal_, ui * my::numdofpernode_ + k) +=
              timefacfac * diff_manager()->GetPhasePoroTort(0) * rtf * invf /
              diff_manager()->GetValence(k) * diff_manager()->GetCond() *
              diff_manager()->GetTransNum(k) * conintinv * laplawf;
        }
      }
      else
        FOUR_C_THROW("Diffusion-Conduction material is not specified");
    }
  }
}

/*----------------------------------------------------------------------------------*
 *----------------------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype, probdim>::calc_mat_pot_equ_divi(
    CORE::LINALG::SerialDenseMatrix& emat, const double timefacfac, const double invf)
{
  for (unsigned vi = 0; vi < nen_; ++vi)
  {
    for (unsigned ui = 0; ui < nen_; ++ui)
    {
      for (unsigned idim = 0; idim < nsd_; ++idim)
      {
        const int fvi = my::numdofpernode_ * vi + my::numscal_;
        const int fui = my::numdofpernode_ * ui + (my::numscal_ + 1) + idim;
        /* current continuity term */
        /*
             /               \
            |                 |
            | w, nabla o Di   |
            |                 |
             \               /
        */
        // emat(fvi,fui) += timefacfac*funct_(vi);*derxy_(idim,ui);

        /* current continuity term */
        /*
             /               \
            |                 |
            | grad phi,  Di   |
            |                 |
             \               /
        */
        // version a: (grad phi,  Di)
        emat(fvi, fui) -= timefacfac * invf * my::derxy_(idim, vi) * my::funct_(ui);
      }  // end for(idim)
    }    // end for(ui)
  }      // end for(vi)
}

/*----------------------------------------------------------------------------------*
 *----------------------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype, probdim>::calc_mat_cur_equ_cur(
    CORE::LINALG::SerialDenseMatrix& emat, const double timefacfac, const double invf)
{
  // (v, i)
  for (unsigned vi = 0; vi < nen_; ++vi)
  {
    for (unsigned ui = 0; ui < nen_; ++ui)
    {
      for (unsigned idim = 0; idim < nsd_; ++idim)
      {
        const int fvi = vi * my::numdofpernode_ + (my::numscal_ + 1) + idim;
        const int fui = ui * my::numdofpernode_ + (my::numscal_ + 1) + idim;

        emat(fvi, fui) += timefacfac * invf * my::funct_(vi) * my::funct_(ui);
      }
    }
  }
}

/*----------------------------------------------------------------------------------*
 *----------------------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype, probdim>::calc_mat_cur_equ_ohm(
    CORE::LINALG::SerialDenseMatrix& emat, const double timefacfac, const double invf,
    const CORE::LINALG::Matrix<nsd_, 1>& gradpot)
{
  // (v, kappa grad phi)
  for (unsigned vi = 0; vi < nen_; ++vi)
  {
    for (unsigned ui = 0; ui < nen_; ++ui)
    {
      for (unsigned idim = 0; idim < nsd_; ++idim)
      {
        const int fvi = vi * my::numdofpernode_ + (my::numscal_ + 1) + idim;
        const int fui = ui * my::numdofpernode_ + my::numscal_;

        emat(fvi, fui) += timefacfac * invf * diff_manager()->GetPhasePoroTort(0) * my::funct_(vi) *
                          diff_manager()->GetCond() * my::derxy_(idim, ui);

        // linearization of conductivity in the ohmic resistance term (current equation)
        //
        // (w, D(kappa(c)) grad phi)
        //
        for (int k = 0; k < my::numscal_; ++k)
        {
          emat(fvi, ui * my::numdofpernode_ + k) +=
              timefacfac * invf * diff_manager()->GetPhasePoroTort(0) * my::funct_(vi) *
              diff_manager()->GetConcDerivCond(k) * my::funct_(ui) * gradpot(idim);
        }
      }
    }
  }
}

/*----------------------------------------------------------------------------------*
 *----------------------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype, probdim>::calc_mat_cur_equ_conc(
    CORE::LINALG::SerialDenseMatrix& emat, const double timefacfac, const double rtf,
    const double rtffc, const std::vector<double>& invfval, const double newman_const_a,
    const double newman_const_b, const std::vector<CORE::LINALG::Matrix<nsd_, 1>>& gradphi,
    const std::vector<double>& conintinv)
{
  for (unsigned vi = 0; vi < nen_; ++vi)
  {
    for (unsigned ui = 0; ui < nen_; ++ui)
    {
      // diffusive term
      // (grad w, D grad c)
      for (unsigned idim = 0; idim < nsd_; ++idim)
      {
        for (int k = 0; k < my::numscal_; ++k)
        {
          if (diffcondmat_ == INPAR::ELCH::diffcondmat_newman)
          {
            emat(
                vi * my::numdofpernode_ + (my::numscal_ + 1) + idim, ui * my::numdofpernode_ + k) +=
                timefacfac * rtffc * diff_manager()->GetPhasePoroTort(0) * my::funct_(vi) *
                diff_manager()->GetCond() *
                (newman_const_a + (newman_const_b * diff_manager()->GetTransNum(k))) *
                conintinv[k] * my::derxy_(idim, ui);
          }
          else if (diffcondmat_ == INPAR::ELCH::diffcondmat_ion)
          {
            if (diffcondparams_->DiffusionCoeffBased())
            {
              emat(vi * my::numdofpernode_ + (my::numscal_ + 1) + idim,
                  ui * my::numdofpernode_ + k) += timefacfac * diff_manager()->GetPhasePoroTort(0) *
                                                  my::funct_(vi) * diff_manager()->GetValence(k) *
                                                  diff_manager()->GetIsotropicDiff(k) *
                                                  my::derxy_(idim, ui);
            }
            else
            {
              // linearization wrt nabla c_k
              emat(vi * my::numdofpernode_ + (my::numscal_ + 1) + idim,
                  ui * my::numdofpernode_ + k) +=
                  timefacfac * diff_manager()->GetPhasePoroTort(0) * invfval[k] * rtf *
                  my::funct_(vi) * diff_manager()->GetCond() * diff_manager()->GetTransNum(k) *
                  conintinv[k] * my::derxy_(idim, ui);

              // linearization wrt 1/c_k
              emat(vi * my::numdofpernode_ + (my::numscal_ + 1) + idim,
                  ui * my::numdofpernode_ + k) +=
                  -timefacfac * diff_manager()->GetPhasePoroTort(0) * rtf * invfval[k] *
                  diff_manager()->GetCond() * diff_manager()->GetTransNum(k) * conintinv[k] *
                  conintinv[k] * my::funct_(vi) * (gradphi[k])(idim)*my::funct_(ui);

              // linearization wrt kappa
              double term_vi = 0.0;
              for (int iscal = 0; iscal < my::numscal_; ++iscal)
              {
                term_vi += invfval[iscal] * diff_manager()->GetTransNum(iscal) * conintinv[iscal] *
                           my::funct_(vi) * (gradphi[iscal])(idim)*my::funct_(ui);
              }
              emat(vi * my::numdofpernode_ + (my::numscal_ + 1) + idim,
                  ui * my::numdofpernode_ + k) += timefacfac * diff_manager()->GetPhasePoroTort(0) *
                                                  diff_manager()->GetConcDerivCond(k) * rtf *
                                                  term_vi;

              // linearization wrt transference number
              for (int iscal = 0; iscal < my::numscal_; ++iscal)
              {
                emat(vi * my::numdofpernode_ + (my::numscal_ + 1) + idim,
                    ui * my::numdofpernode_ + iscal) +=
                    timefacfac * diff_manager()->GetPhasePoroTort(0) * rtf * invfval[k] *
                    diff_manager()->GetCond() * (diff_manager()->GetDerivTransNum(k, iscal)) *
                    conintinv[k] * my::funct_(vi) * (gradphi[k])(idim)*my::funct_(ui);
              }
            }
          }
          else
            FOUR_C_THROW("Diffusion-Conduction material is not specified");
        }
      }
    }  // for ui
  }    // for vi
}

/*-------------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype, probdim>::calc_rhs_cond_ohm(
    CORE::LINALG::SerialDenseVector& erhs, const int k, const double rhsfac, const double invfval,
    const CORE::LINALG::Matrix<nsd_, 1>& gradpot)
{
  for (unsigned vi = 0; vi < nen_; ++vi)
  {
    // diffusive term
    double laplawfrhs_gradpot = 0.0;
    my::get_laplacian_weak_form_rhs(laplawfrhs_gradpot, gradpot, vi);
    erhs[vi * my::numdofpernode_ + k] -= rhsfac * diff_manager()->GetPhasePoroTort(0) *
                                         diff_manager()->GetTransNum(k) *
                                         diff_manager()->GetCond() * invfval * laplawfrhs_gradpot;
  }
}

/*-------------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype, probdim>::calc_rhs_cond_conc(
    CORE::LINALG::SerialDenseVector& erhs, const int k, const double rhsfac, const double rtffcval,
    const double newman_const_a, const double newman_const_b,
    const CORE::LINALG::Matrix<nsd_, 1>& gradphi, const std::vector<double>& conintinv)
{
  for (unsigned vi = 0; vi < nen_; ++vi)
  {
    for (int iscal = 0; iscal < my::numscal_; ++iscal)
    {
      // diffusive term second
      double laplawfrhs_gradphi(0.0);
      my::get_laplacian_weak_form_rhs(
          laplawfrhs_gradphi, gradphi, vi);  // compute once, reuse below!

      // formulation a): plain ionic diffusion coefficients without using ENC
      //
      // (grad w, sum(D_i grad Dc))
      erhs[vi * my::numdofpernode_ + k] -=
          rhsfac * diff_manager()->GetPhasePoroTort(0) * rtffcval * diff_manager()->GetTransNum(k) *
          diff_manager()->GetCond() * diff_manager()->GetThermFac() *
          (newman_const_a + (newman_const_b * diff_manager()->GetTransNum(iscal))) *
          conintinv[iscal] * laplawfrhs_gradphi;
    }
  }
}

/*-------------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype, probdim>::calc_rhs_cond(
    CORE::LINALG::SerialDenseVector& erhs, const int k, const double rhsfac, const double invfval,
    const CORE::LINALG::Matrix<nsd_, 1>& curint)
{
  for (unsigned vi = 0; vi < nen_; ++vi)
  {
    double laplawfrhs_cur = 0.0;
    my::get_laplacian_weak_form_rhs(laplawfrhs_cur, curint, vi);

    erhs[vi * my::numdofpernode_ + k] -=
        -rhsfac * diff_manager()->GetTransNum(k) * invfval * laplawfrhs_cur;
  }
}

/*-------------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype, probdim>::calc_rhs_cond_diff(
    CORE::LINALG::SerialDenseVector& erhs, const int k, const double rhsfac,
    const std::vector<CORE::LINALG::Matrix<nsd_, 1>>& gradphi)
{
  for (unsigned vi = 0; vi < nen_; ++vi)
  {
    for (int iscal = 0; iscal < my::numscal_; ++iscal)
    {
      // diffusive term second
      double laplawfrhs_gradphi(0.0);
      my::get_laplacian_weak_form_rhs(
          laplawfrhs_gradphi, gradphi[iscal], vi);  // compute once, reuse below!

      // formulation a:  plain ionic diffusion coefficients: sum (z_i D_i nabla c_i)
      erhs[vi * my::numdofpernode_ + k] -=
          -rhsfac * diff_manager()->GetPhasePoroTort(0) * diff_manager()->GetTransNum(k) *
          diff_manager()->GetValence(iscal) / diff_manager()->GetValence(k) *
          diff_manager()->GetIsotropicDiff(iscal) * laplawfrhs_gradphi;
    }
  }
}

/*-------------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype, probdim>::calc_rhs_pot_equ_divi_conc(
    CORE::LINALG::SerialDenseVector& erhs, const int k, const double rhsfac, const double rtf,
    const std::vector<double>& invfval, const double rtffc, const double newman_const_a,
    const double newman_const_b, const CORE::LINALG::Matrix<nsd_, 1>& gradphi,
    const double conintinv)
{
  for (unsigned vi = 0; vi < nen_; ++vi)
  {
    if (diffcondmat_ == INPAR::ELCH::diffcondmat_ion)
    {
      // diffusive term second
      double laplawf2(0.0);
      my::get_laplacian_weak_form_rhs(laplawf2, gradphi, vi);  // compute once, reuse below!

      if (diffcondparams_->DiffusionCoeffBased())
      {
        erhs[vi * my::numdofpernode_ + my::numscal_] -=
            rhsfac * diff_manager()->GetPhasePoroTort(0) * diff_manager()->GetValence(k) *
            diff_manager()->GetIsotropicDiff(k) * laplawf2;
      }
      else
      {
        erhs[vi * my::numdofpernode_ + my::numscal_] -=
            rhsfac * diff_manager()->GetPhasePoroTort(0) * rtf * invfval[k] *
            diff_manager()->GetCond() * diff_manager()->GetTransNum(k) * conintinv * laplawf2;
      }
    }
    // thermodynamic factor only implemented for Newman
    else if (diffcondmat_ == INPAR::ELCH::diffcondmat_newman)
    {
      // diffusive term second
      double laplawf2(0.0);
      my::get_laplacian_weak_form_rhs(laplawf2, gradphi, vi);  // compute once, reuse below!

      erhs[vi * my::numdofpernode_ + my::numscal_] -=
          rhsfac * diff_manager()->GetPhasePoroTort(0) * rtffc * diff_manager()->GetCond() *
          diff_manager()->GetThermFac() *
          (newman_const_a + (newman_const_b * diff_manager()->GetTransNum(k))) * conintinv *
          laplawf2;
    }
    else
      FOUR_C_THROW("Diffusion-Conduction material is not specified");
  }
}

/*-------------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype, probdim>::calc_rhs_pot_equ_divi(
    CORE::LINALG::SerialDenseVector& erhs, const double rhsfac, const double invf,
    const CORE::LINALG::Matrix<nsd_, 1>& curint)
{
  for (unsigned vi = 0; vi < nen_; ++vi)
  {
    double laplawf = 0.0;
    // version a: (grad phi,  Di)
    my::get_laplacian_weak_form_rhs(laplawf, curint, vi);
    erhs[vi * my::numdofpernode_ + my::numscal_] -= -rhsfac * invf * laplawf;
  }
}

/*-------------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype, probdim>::calc_rhs_cur_equ_cur(
    CORE::LINALG::SerialDenseVector& erhs, const double rhsfac, const double invf,
    const CORE::LINALG::Matrix<nsd_, 1>& curint)
{
  for (unsigned vi = 0; vi < nen_; ++vi)
  {
    for (unsigned idim = 0; idim < nsd_; ++idim)
    {
      // (v, i)
      erhs[vi * my::numdofpernode_ + (my::numscal_ + 1) + idim] -=
          rhsfac * invf * my::funct_(vi) * curint(idim);
    }
  }
}

/*-------------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype, probdim>::calc_rhs_cur_equ_ohm(
    CORE::LINALG::SerialDenseVector& erhs, const double rhsfac, const double invf,
    const CORE::LINALG::Matrix<nsd_, 1>& gradpot)
{
  for (unsigned vi = 0; vi < nen_; ++vi)
  {
    for (unsigned idim = 0; idim < nsd_; ++idim)
    {
      // (v, kappa grad phi)
      erhs[vi * my::numdofpernode_ + (my::numscal_ + 1) + idim] -=
          rhsfac * invf * diff_manager()->GetPhasePoroTort(0) * my::funct_(vi) *
          diff_manager()->GetCond() * gradpot(idim);
    }
  }
}

/*-------------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype, probdim>::calc_rhs_cur_equ_conc(
    CORE::LINALG::SerialDenseVector& erhs, const double rhsfac, const double rtf,
    const std::vector<double>& invfval, const double rtffc, const double newman_const_a,
    const double newman_const_b, const std::vector<CORE::LINALG::Matrix<nsd_, 1>>& gradphi,
    const std::vector<double>& conintinv)
{
  for (unsigned vi = 0; vi < nen_; ++vi)
  {
    for (unsigned idim = 0; idim < nsd_; ++idim)
    {
      for (int k = 0; k < my::numscal_; ++k)
      {
        if (diffcondmat_ == INPAR::ELCH::diffcondmat_newman)
        {
          erhs[vi * my::numdofpernode_ + (my::numscal_ + 1) + idim] -=
              rhsfac * diff_manager()->GetPhasePoroTort(0) * my::funct_(vi) * rtffc *
              diff_manager()->GetCond() *
              (newman_const_a + (newman_const_b * diff_manager()->GetTransNum(k))) * conintinv[k] *
              gradphi[k](idim);
        }
        else if (diffcondmat_ == INPAR::ELCH::diffcondmat_ion)
        {
          if (diffcondparams_->DiffusionCoeffBased())
          {
            erhs[vi * my::numdofpernode_ + (my::numscal_ + 1) + idim] -=
                rhsfac * diff_manager()->GetPhasePoroTort(0) * my::funct_(vi) *
                diff_manager()->GetValence(k) * diff_manager()->GetIsotropicDiff(k) *
                gradphi[k](idim);
          }
          else
          {
            erhs[vi * my::numdofpernode_ + (my::numscal_ + 1) + idim] -=
                rhsfac * diff_manager()->GetPhasePoroTort(0) * my::funct_(vi) * rtf *
                diff_manager()->GetCond() * invfval[k] * diff_manager()->GetTransNum(k) *
                conintinv[k] * gradphi[k](idim);
          }
        }
        else
          FOUR_C_THROW("Diffusion-Conduction material is not specified");
      }
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype, probdim>::correction_for_flux_across_dc(
    DRT::Discretization& discretization, const std::vector<int>& lm,
    CORE::LINALG::SerialDenseMatrix& emat, CORE::LINALG::SerialDenseVector& erhs)
{
  // get dirichlet toggle from the discretization
  // we always get the dirichet toggle:
  // in this function we check if the actual nodes have a dirichlet value
  Teuchos::RCP<const Epetra_Vector> dctoggle = discretization.GetState("dctoggle");
  std::vector<double> mydctoggle(lm.size());
  CORE::FE::ExtractMyValues(*dctoggle, mydctoggle, lm);

  double val = 0.0;
  for (unsigned vi = 0; vi < nen_; ++vi)
  {
    for (int k = 0; k < my::numscal_; ++k)
    {
      // here we check if the actual nodes have a dirichlet value
      if (mydctoggle[vi * my::numdofpernode_ + k] == 1)
      {
        const int fvi = vi * my::numdofpernode_ + k;
        // We use the fact, that the rhs vector value for boundary nodes
        // is equivalent to the integrated negative normal flux
        // due to diffusion and migration

        // scaling of div i results in a matrix with better condition number
        val = erhs[fvi];
        erhs[vi * my::numdofpernode_ + my::numscal_] += diff_manager()->GetValence(k) * (-val);
        // corresponding linearization
        for (unsigned ui = 0; ui < nen_; ++ui)
        {
          val = emat(vi * my::numdofpernode_ + k, ui * my::numdofpernode_ + k);
          emat(vi * my::numdofpernode_ + my::numscal_, ui * my::numdofpernode_ + k) +=
              diff_manager()->GetValence(k) * (-val);
          val = emat(vi * my::numdofpernode_ + k, ui * my::numdofpernode_ + my::numscal_);
          emat(vi * my::numdofpernode_ + my::numscal_, ui * my::numdofpernode_ + my::numscal_) +=
              diff_manager()->GetValence(k) * (-val);
        }
      }

      // Dirichlet conditions on the potential are only allowed for the newman material
      // since additional information about the reacting species is required. This is fulfilled
      // naturally for the Newman material since only one species is allowed in this case. Newman
      // material models binary electrolytes where the second species is condensed via the ENC!
      if (diffcondmat_ == INPAR::ELCH::diffcondmat_newman)
      {
        if (mydctoggle[vi * my::numdofpernode_ + my::numscal_] == 1)
        {
          // reacting species 0:
          // Newman material: reacting species is always the first species since there is only one
          // species other materials: one have to find a way to define the reacting species
          int l = 0;

          const int fvi = vi * my::numdofpernode_ + my::numscal_;
          // We use the fact, that the rhs vector value for boundary nodes
          // is equivalent to the integrated negative normal flux
          // due to diffusion and migration

          // scaling of div i results in a matrix with better condition number
          val = erhs[fvi];
          erhs[vi * my::numdofpernode_ + l] += 1.0 / diff_manager()->GetValence(l) * (-val);
          // corresponding linearization
          for (unsigned ui = 0; ui < nen_; ++ui)
          {
            val = emat(vi * my::numdofpernode_ + my::numscal_, ui * my::numdofpernode_ + l);
            emat(vi * my::numdofpernode_ + l, ui * my::numdofpernode_ + l) +=
                1.0 / diff_manager()->GetValence(l) * (-val);
            val = emat(
                vi * my::numdofpernode_ + my::numscal_, ui * my::numdofpernode_ + my::numscal_);
            emat(vi * my::numdofpernode_ + l, ui * my::numdofpernode_ + my::numscal_) +=
                1.0 / diff_manager()->GetValence(l) * (-val);
          }
        }
      }
    }  // for k
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype, probdim>::get_material_params(
    const DRT::Element* ele, std::vector<double>& densn, std::vector<double>& densnp,
    std::vector<double>& densam, double& visc, const int iquad)
{
  // extract material from element
  Teuchos::RCP<CORE::MAT::Material> material = ele->Material();

  // evaluate electrolyte material
  if (material->MaterialType() == CORE::Materials::m_elchmat)
  {
    utils()->MatElchMat(material, var_manager()->Phinp(), var_manager()->Temperature(),
        myelch::elchparams_->EquPot(), myelch::elchparams_->Faraday() * var_manager()->FRT(),
        diff_manager(), diffcondmat_);
  }
  else
    FOUR_C_THROW("Invalid material type!");
}  // DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::get_material_params


// template classes
// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<CORE::FE::CellType::line2, 1>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<CORE::FE::CellType::line2, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<CORE::FE::CellType::line2, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<CORE::FE::CellType::line3, 1>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<CORE::FE::CellType::tri3, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<CORE::FE::CellType::tri3, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<CORE::FE::CellType::tri6, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<CORE::FE::CellType::quad4, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<CORE::FE::CellType::quad4, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<CORE::FE::CellType::quad9, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<CORE::FE::CellType::nurbs9, 2>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<CORE::FE::CellType::hex8, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<CORE::FE::CellType::hex27, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<CORE::FE::CellType::tet4, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<CORE::FE::CellType::tet10, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<CORE::FE::CellType::pyramid5, 3>;

FOUR_C_NAMESPACE_CLOSE

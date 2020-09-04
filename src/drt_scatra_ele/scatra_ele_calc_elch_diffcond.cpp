/*--------------------------------------------------------------------------*/
/*! \file

\brief evaluation of scatra elements for isothermal diffusion-conduction ion-transport equations

\level 2

*/
/*--------------------------------------------------------------------------*/
#include "scatra_ele_calc_elch_diffcond.H"
#include "scatra_ele_parameter_std.H"
#include "scatra_ele_parameter_timint.H"
#include "scatra_ele_utils_elch_diffcond.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"

#include "../drt_mat/material.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>*
DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::Instance(const int numdofpernode,
    const int numscal, const std::string& disname, const ScaTraEleCalcElchDiffCond* delete_me)
{
  static std::map<std::string, ScaTraEleCalcElchDiffCond<distype>*> instances;

  if (delete_me == NULL)
  {
    if (instances.find(disname) == instances.end())
      instances[disname] = new ScaTraEleCalcElchDiffCond<distype>(numdofpernode, numscal, disname);
  }

  else
  {
    for (typename std::map<std::string, ScaTraEleCalcElchDiffCond<distype>*>::iterator i =
             instances.begin();
         i != instances.end(); ++i)
      if (i->second == delete_me)
      {
        delete i->second;
        instances.erase(i);
        return NULL;
      }
    dserror("Could not locate the desired instance. Internal error.");
  }

  return instances[disname];
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance(0, 0, "", this);

  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::ScaTraEleCalcElchDiffCond(
    const int numdofpernode, const int numscal, const std::string& disname)
    : DRT::ELEMENTS::ScaTraEleCalcElchElectrode<distype>::ScaTraEleCalcElchElectrode(
          numdofpernode, numscal, disname),
      diffcondmat_(INPAR::ELCH::diffcondmat_undefined),
      // parameter class for diffusion-conduction formulation
      diffcondparams_(DRT::ELEMENTS::ScaTraEleParameterElchDiffCond::Instance(disname))
{
  // replace diffusion manager for electrodes by diffusion manager for diffusion-conduction
  // formulation
  my::diffmanager_ = Teuchos::rcp(new ScaTraEleDiffManagerElchDiffCond(my::numscal_));

  // replace internal variable manager for electrodes by internal variable manager for
  // diffusion-conduction formulation
  my::scatravarmanager_ =
      Teuchos::rcp(new ScaTraEleInternalVariableManagerElchDiffCond<my::nsd_, my::nen_>(
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
      dserror("Invalid closing equation for electric potential!");
      break;
    }
  }

  // safety checks for stabilization settings
  if (my::scatrapara_->StabType() != INPAR::SCATRA::stabtype_no_stabilization or
      my::scatrapara_->TauDef() != INPAR::SCATRA::tau_zero)
    dserror(
        "No stabilization available for the diffusion-conduction formulation, since we had no "
        "problems so far.");
  if (my::scatrapara_->MatGP() == false or my::scatrapara_->TauGP() == false)
    dserror(
        "Since most of the materials of the Diffusion-conduction formulation depend on the "
        "concentration,\n"
        "an evaluation of the material and the stabilization parameter at the element center is "
        "disabled.");

  return;
}


/*----------------------------------------------------------------------*
|  calculate system matrix and rhs                           ehrl  02/14|
*----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::CalcMatAndRhs(
    Epetra_SerialDenseMatrix& emat,  //!< element matrix to calculate
    Epetra_SerialDenseVector& erhs,  //!< element rhs to calculate
    const int k,                     //!< index of current scalar
    const double fac,                //!< domain-integration factor
    const double timefacfac,         //!< domain-integration factor times time-integration factor
    const double rhsfac,      //!< time-integration factor for rhs times domain-integration factor
    const double taufac,      //!< tau times domain-integration factor
    const double timetaufac,  //!< domain-integration factor times tau times time-integration factor
    const double
        rhstaufac,  //!< time-integration factor for rhs times tau times domain-integration factor
    LINALG::Matrix<my::nen_, 1>&
        tauderpot,  //!< derivatives of stabilization parameter w.r.t. electric potential
    double& rhsint  //!< rhs at Gauss point
)
{
  //----------------------------------------------------------------
  // 1) element matrix: instationary terms
  //----------------------------------------------------------------

  if (not my::scatraparatimint_->IsStationary())
    my::CalcMatMass(emat, k, fac, DiffManager()->GetPhasePoro(0));

  //----------------------------------------------------------------
  // 2) element matrix: stationary terms of ion-transport equation
  //----------------------------------------------------------------

  // 2a)  element matrix: convective term
  my::CalcMatConv(emat, k, timefacfac, DiffManager()->GetPhasePoro(0), VarManager()->SGConv());

  // 2b)  element matrix: diffusion term
  //      i)  constant diffusion coefficient
  my::CalcMatDiff(emat, k, timefacfac * DiffManager()->GetPhasePoroTort(0));

  //      ii) concentration depending diffusion coefficient
  //          (additional term for Newman material)
  if (diffcondmat_ == INPAR::ELCH::diffcondmat_newman)
    myelectrode::CalcMatDiffCoeffLin(
        emat, k, timefacfac, VarManager()->GradPhi(k), DiffManager()->GetPhasePoroTort(0));

  // 2c) electrical conduction term (transport equation)
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
    CalcMatCondOhm(emat, k, timefacfac, DiffManager()->InvFVal(k), VarManager()->GradPot());

    //    ii) conduction term + concentration overpotential
    //        (w_k, - t_k RT/F kappa (thermfactor) f(t_k) nabla ln c_k /(z_k F))
    if (diffcondmat_ == INPAR::ELCH::diffcondmat_newman)
      CalcMatCondConc(emat, k, timefacfac, VarManager()->RTFFC() / DiffManager()->GetValence(k),
          diffcondparams_->NewmanConstA(), diffcondparams_->NewmanConstB(),
          VarManager()->GradPhi(k), VarManager()->ConIntInv());
  }
  // equation for current is solved independently
  else
  {
    // current term (with current as a solution variable)
    CalcMatCond(emat, k, timefacfac, DiffManager()->InvFVal(k), VarManager()->CurInt());

    // this coupling term cancels out for a 2 equation system
    if (diffcondmat_ == INPAR::ELCH::diffcondmat_ion)
      CalcMatCondDiff(emat, k, timefacfac, DiffManager()->InvFVal(k), VarManager()->GradPhi());
  }  // end if(not diffcondparams_->CurSolVar())

  //---------------------------------------------------------------------
  // 3)   governing equation for the electric potential field and current
  //---------------------------------------------------------------------
  // see function CalcMatAndRhsOutsideScalarLoop()

  //-----------------------------------------------------------------------
  // 4) element right hand side vector (neg. residual of nonlinear problem)
  //-----------------------------------------------------------------------

  if (my::scatraparatimint_->IsIncremental() and not my::scatraparatimint_->IsStationary())
    my::CalcRHSLinMass(
        erhs, k, rhsfac, fac, DiffManager()->GetPhasePoro(0), DiffManager()->GetPhasePoro(0));

  // adaption of rhs with respect to time integration
  my::ComputeRhsInt(rhsint, DiffManager()->GetPhasePoro(0), DiffManager()->GetPhasePoro(0),
      VarManager()->Hist(k));

  // add RHS and history contribution
  my::CalcRHSHistAndSource(erhs, k, fac, rhsint);

  // convective term
  my::CalcRHSConv(erhs, k, rhsfac * DiffManager()->GetPhasePoro(0));

  // diffusion term
  my::CalcRHSDiff(erhs, k, rhsfac * DiffManager()->GetPhasePoroTort(0));

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
    CalcRhsCondOhm(erhs, k, rhsfac, DiffManager()->InvFVal(k), VarManager()->GradPot());

    // if(diffcondmat_==INPAR::ELCH::diffcondmat_ion): all terms cancel out
    if (diffcondmat_ == INPAR::ELCH::diffcondmat_newman)
      CalcRhsCondConc(erhs, k, rhsfac, VarManager()->RTFFC() / DiffManager()->GetValence(k),
          diffcondparams_->NewmanConstA(), diffcondparams_->NewmanConstB(),
          VarManager()->GradPhi(k), VarManager()->ConIntInv());
  }

  // equation for current is solved independently
  else
  {
    CalcRhsCond(erhs, k, rhsfac, DiffManager()->InvFVal(k), VarManager()->CurInt());

    if (diffcondmat_ == INPAR::ELCH::diffcondmat_ion)
      CalcRhsCondDiff(erhs, k, rhsfac, VarManager()->GradPhi());
  }

  return;
}


/*----------------------------------------------------------------------*
|  calculate system matrix and rhs                           ehrl  02/14|
*----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::CalcMatAndRhsOutsideScalarLoop(
    Epetra_SerialDenseMatrix& emat,  //!< element matrix to calculate
    Epetra_SerialDenseVector& erhs,  //!< element rhs to calculate
    const double fac,                //!< domain-integration factor
    const double timefacfac,         //!< domain-integration factor times time-integration factor
    const double rhsfac  //!< time-integration factor for rhs times domain-integration factor
)
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
      myelectrode::CalcMatPotEquDiviOhm(emat, timefacfac, VarManager()->InvF(),
          VarManager()->GradPot(), DiffManager()->GetPhasePoroTort(0));

      //
      myelectrode::CalcRhsPotEquDiviOhm(erhs, rhsfac, VarManager()->InvF(), VarManager()->GradPot(),
          DiffManager()->GetPhasePoroTort(0));

      //  ii)  concentration overpotential
      //      (w_k, RT/F kappa (thermfactor) f(t_k) nabla ln c_k)
      for (int k = 0; k < my::numscal_; ++k)
      {
        //
        CalcMatPotEquDiviConc(emat, k, timefacfac, VarManager()->RTFFC(), VarManager()->RTF(),
            VarManager()->InvF(), diffcondparams_->NewmanConstA(), diffcondparams_->NewmanConstB(),
            VarManager()->GradPhi(k), VarManager()->ConIntInv(k));

        //
        CalcRhsPotEquDiviConc(erhs, k, rhsfac, VarManager()->RTF(), DiffManager()->InvFVal(),
            VarManager()->RTFFC(), diffcondparams_->NewmanConstA(), diffcondparams_->NewmanConstB(),
            VarManager()->GradPhi(k), VarManager()->ConIntInv(k));
      }
    }
    // 3b)  ENC
    else if (myelch::elchparams_->EquPot() == INPAR::ELCH::equpot_enc)
    {
      for (int k = 0; k < my::numscal_; ++k)
      {
        //
        myelch::CalcMatPotEquENC(emat, k, fac, my::scatraparatimint_->AlphaF());

        //
        myelch::CalcRhsPotEquENC(erhs, k, fac, VarManager()->Phinp(k));
      }
    }
    else
      dserror("(div i, ENC) are the options available in the Diffusion-Conduction framework");
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
    CalcMatCurEquCur(emat, timefacfac, VarManager()->InvF());

    // (xsi, -D(kappa phi))
    CalcMatCurEquOhm(emat, timefacfac, VarManager()->InvF(), VarManager()->GradPot());

    // (xsi, -D(RT/F kappa (thermfactor) f(t_k) nabla ln c_k))
    CalcMatCurEquConc(emat, timefacfac, VarManager()->RTF(), VarManager()->RTFFC(),
        DiffManager()->InvFVal(), diffcondparams_->NewmanConstA(), diffcondparams_->NewmanConstB(),
        VarManager()->GradPhi(), VarManager()->ConIntInv());

    // (xsi_i,Di)
    CalcRhsCurEquCur(erhs, rhsfac, VarManager()->InvF(), VarManager()->CurInt());

    // (xsi, -D(kappa phi))
    CalcRhsCurEquOhm(erhs, rhsfac, VarManager()->InvF(), VarManager()->GradPot());

    // (xsi, -D(RT/F kappa (thermfactor) f(t_k) nabla ln c_k))
    CalcRhsCurEquConc(erhs, rhsfac, VarManager()->RTF(), DiffManager()->InvFVal(),
        VarManager()->RTFFC(), diffcondparams_->NewmanConstA(), diffcondparams_->NewmanConstB(),
        VarManager()->GradPhi(), VarManager()->ConIntInv());

    //------------------------------------------------------------------------------------------
    // 3)   governing equation for the electric potential field and current (incl. rhs-terms)
    //------------------------------------------------------------------------------------------

    if (myelch::elchparams_->EquPot() == INPAR::ELCH::equpot_divi)
    {
      CalcMatPotEquDivi(emat, timefacfac, VarManager()->InvF());

      CalcRhsPotEquDivi(erhs, rhsfac, VarManager()->InvF(), VarManager()->CurInt());
    }
    else if (myelch::elchparams_->EquPot() == INPAR::ELCH::equpot_enc)
    {
      for (int k = 0; k < my::numscal_; ++k)
      {
        //
        myelch::CalcMatPotEquENC(emat, k, fac, my::scatraparatimint_->AlphaF());

        //
        myelch::CalcRhsPotEquENC(erhs, k, fac, VarManager()->Phinp(k));
      }
    }
    else
      dserror("(div i, ENC) are the options available in the Diffusion-Conduction framework");
  }

  return;
}


/*----------------------------------------------------------------------------------*
|  CalcMat: Conduction term with inserted current - ohmic overpotential  ehrl  02/14|
*-----------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::CalcMatCondOhm(
    Epetra_SerialDenseMatrix& emat,  //!< element matrix to be filled
    const int k,                     //!< index of current scalar
    const double timefacfac,         //!< domain-integration factor times time-integration factor
    const double invfval,            //!< 1/(F z_k)
    const LINALG::Matrix<my::nsd_, 1>& gradpot  //!< gradient of potential at GP
)
{
  for (unsigned vi = 0; vi < my::nen_; ++vi)
  {
    for (unsigned ui = 0; ui < my::nen_; ++ui)
    {
      double laplawf(0.0);
      my::GetLaplacianWeakForm(laplawf, ui, vi);  // compute once, reuse below!

      // linearization of conduction term depending on the potential
      //
      // (grad w, t_k kappa/(F z_k) D(grad phi))
      //
      emat(vi * my::numdofpernode_ + k, ui * my::numdofpernode_ + my::numscal_) +=
          timefacfac * DiffManager()->GetPhasePoroTort(0) * DiffManager()->GetTransNum(k) *
          DiffManager()->GetCond() * invfval * laplawf;

      double laplawfrhs_gradpot(0.0);
      my::GetLaplacianWeakFormRHS(laplawfrhs_gradpot, gradpot, vi);

      for (int iscal = 0; iscal < my::numscal_; ++iscal)
      {
        // linearization of the conductivity in the conduction term depending on the potential
        //
        // (grad w, t_k D(kappa(c))/(F z_k) grad phi)
        //
        emat(vi * my::numdofpernode_ + k, ui * my::numdofpernode_ + iscal) +=
            timefacfac * DiffManager()->GetPhasePoroTort(0) * DiffManager()->GetTransNum(k) *
            invfval * DiffManager()->GetDerivCond(iscal) * my::funct_(ui) * laplawfrhs_gradpot;

        // linearization of the transference number in the conduction term depending on the
        // potential
        //
        // (grad w, D(t_k(c)) kappa/(F z_k) grad phi)
        //
        emat(vi * my::numdofpernode_ + k, ui * my::numdofpernode_ + iscal) +=
            timefacfac * DiffManager()->GetPhasePoroTort(0) *
            (DiffManager()->GetDerivTransNum(k, iscal)) * my::funct_(ui) * invfval *
            DiffManager()->GetCond() * laplawfrhs_gradpot;
      }
    }
  }

  return;
}


/*----------------------------------------------------------------------------------*
|  CalcMat: Conduction term with inserted current - conc. overpotential  ehrl  02/14|
*-----------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::CalcMatCondConc(
    Epetra_SerialDenseMatrix& emat,  //!< element matrix to be filled
    const int k,                     //!< index of current scalar
    const double timefacfac,         //!< domain-integration factor times time-integration factor
    const double rtffcval,           //!< RT/F^2/Newman_const_c/z_k
    const double newman_const_a,     //!< Newman constant a
    const double newman_const_b,     //!< Newman constant b
    const LINALG::Matrix<my::nsd_, 1>& gradphi,  //!< gradient of concentration at GP
    const std::vector<double>& conintinv         //!< inverted concentration at GP
)
{
  // additional safety check in the beginning for Newman materials
  if (k != 0)
    dserror("Material Newman is only valid for one scalar (binary electrolyte utilizing the ENC)");

  for (unsigned vi = 0; vi < my::nen_; ++vi)
  {
    for (unsigned ui = 0; ui < my::nen_; ++ui)
    {
      double laplawf(0.0);
      my::GetLaplacianWeakForm(laplawf, ui, vi);

      // Material Newman is only valid for binary electrolyte utilizing the ENC:
      // -> the equations are solved only for one species (k=0)
      // -> all transport parameter only depend on a single species
      // -> all derivations of the transport parameters wrt this species
      //
      // additional safety check in the beginning of this material: k != 0
      // original safety check by method CheckElchElementParameter() in
      // scatra_ele_calc_service_elch.cpp

      // linearization of conduction term depending on the concentration
      //
      // (grad w, RT/(z_k F^2) kappa thermfac f(t_+) D(grad ln c_k))
      //
      emat(vi * my::numdofpernode_ + k, ui * my::numdofpernode_ + k) +=
          timefacfac * DiffManager()->GetPhasePoroTort(0) * rtffcval *
          DiffManager()->GetTransNum(k) * DiffManager()->GetCond() * DiffManager()->GetThermFac() *
          (newman_const_a + (newman_const_b * DiffManager()->GetTransNum(k))) * conintinv[k] *
          laplawf;

      // linearization of conduction term depending on the concentration is implemented
      // only for one species
      // otherwise you would need a second loop over the all scalars
      double laplawfrhs_gradc(0.0);
      my::GetLaplacianWeakFormRHS(laplawfrhs_gradc, gradphi, vi);

      // Linearization wrt ln c
      emat(vi * my::numdofpernode_ + k, ui * my::numdofpernode_ + k) +=
          -timefacfac * DiffManager()->GetPhasePoroTort(0) * rtffcval *
          DiffManager()->GetTransNum(k) * DiffManager()->GetCond() * DiffManager()->GetThermFac() *
          (newman_const_a + (newman_const_b * DiffManager()->GetTransNum(k))) * conintinv[k] *
          conintinv[k] * laplawfrhs_gradc * my::funct_(ui);

      // Linearization wrt kappa
      emat(vi * my::numdofpernode_ + k, ui * my::numdofpernode_ + k) +=
          timefacfac * DiffManager()->GetPhasePoroTort(0) * rtffcval *
          DiffManager()->GetTransNum(k) * DiffManager()->GetThermFac() *
          (newman_const_a + (newman_const_b * DiffManager()->GetTransNum(k))) * conintinv[k] *
          laplawfrhs_gradc * DiffManager()->GetDerivCond(k) * my::funct_(ui);

      // Linearization wrt transference number 1
      emat(vi * my::numdofpernode_ + k, ui * my::numdofpernode_ + k) +=
          timefacfac * DiffManager()->GetPhasePoroTort(0) * rtffcval * DiffManager()->GetCond() *
          conintinv[k] * laplawfrhs_gradc * DiffManager()->GetThermFac() *
          (newman_const_a + newman_const_b * DiffManager()->GetTransNum(k)) *
          DiffManager()->GetDerivTransNum(k, k) * my::funct_(ui);

      // Linearization wrt transference number 2
      emat(vi * my::numdofpernode_ + k, ui * my::numdofpernode_ + k) +=
          timefacfac * DiffManager()->GetPhasePoroTort(0) * rtffcval * DiffManager()->GetCond() *
          DiffManager()->GetThermFac() * conintinv[k] * laplawfrhs_gradc *
          DiffManager()->GetTransNum(k) * newman_const_b * DiffManager()->GetDerivTransNum(k, k) *
          my::funct_(ui);

      // Linearization wrt thermodynamic factor
      emat(vi * my::numdofpernode_ + k, ui * my::numdofpernode_ + k) +=
          timefacfac * DiffManager()->GetPhasePoroTort(0) * rtffcval *
          DiffManager()->GetTransNum(k) * DiffManager()->GetCond() *
          (newman_const_a + (newman_const_b * DiffManager()->GetTransNum(k))) * conintinv[k] *
          laplawfrhs_gradc * DiffManager()->GetDerivThermFac(k) * my::funct_(ui);
    }  // for ui
  }    // for vi

  return;
}


/*----------------------------------------------------------------------------------*
|  CalcMat: Conduction term without inserted current                     ehrl  02/14|
*-----------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::CalcMatCond(
    Epetra_SerialDenseMatrix& emat,  //!< element matrix to be filled
    const int k,                     //!< index of current scalar
    const double timefacfac,         //!< domain-integration factor times time-integration factor
    const double invfval,            //!< 1/(F z_k)
    const LINALG::Matrix<my::nsd_, 1>& curint  //!< current at GP
)
{
  for (unsigned vi = 0; vi < my::nen_; ++vi)
  {
    const int fvi = vi * my::numdofpernode_ + k;

    for (unsigned ui = 0; ui < my::nen_; ++ui)
    {
      for (unsigned idim = 0; idim < my::nsd_; ++idim)
      {
        const int fui = ui * my::numdofpernode_ + (my::numscal_ + 1) + idim;

        // linearization of conduction term depending on current flow
        //
        // (grad w, t_k/(F z_k) Di)
        //
        emat(fvi, fui) += -timefacfac * my::derxy_(idim, vi) * DiffManager()->GetTransNum(k) *
                          invfval * my::funct_(ui);

        // linearization of transference number in conduction term depending on current flow
        //
        // (grad w, Dt_k(c)/(F z_k) i)
        //
        for (int iscal = 0; iscal < my::numscal_; ++iscal)
          emat(fvi, ui * my::numdofpernode_ + iscal) +=
              -timefacfac * my::derxy_(idim, vi) * (DiffManager()->GetDerivTransNum(k, iscal)) *
              my::funct_(ui) * invfval * curint(idim);
      }
    }
  }

  return;
}


/*----------------------------------------------------------------------------------*
|  CalcMat: Additional diffusion term without inserted current           ehrl  02/14|
*-----------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::CalcMatCondDiff(
    Epetra_SerialDenseMatrix& emat,  //!< element matrix to be filled
    const int k,                     //!< index of current scalar
    const double timefacfac,         //!< domain-integration factor times time-integration factor
    const double invfval,            //!< 1/(F z_k)
    const std::vector<LINALG::Matrix<my::nsd_, 1>>& gradphi  //!< gradient of concentration at GP
)
{
  for (unsigned vi = 0; vi < my::nen_; ++vi)
  {
    for (unsigned ui = 0; ui < my::nen_; ++ui)
    {
      // compute once, reuse below!
      double laplawf(0.0);
      my::GetLaplacianWeakForm(laplawf, ui, vi);

      for (int iscal = 0; iscal < my::numscal_; ++iscal)
      {
        // formulation a): plain ionic diffusion coefficients without using ENC
        //
        // (grad w, t_k/z_k*sum_i(D_i grad Dc))
        //
        emat(vi * my::numdofpernode_ + k, ui * my::numdofpernode_ + iscal) +=
            -timefacfac * DiffManager()->GetPhasePoroTort(0) * DiffManager()->GetTransNum(k) /
            DiffManager()->GetValence(k) * DiffManager()->GetValence(iscal) *
            DiffManager()->GetIsotropicDiff(iscal) * laplawf;
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
          my::GetLaplacianWeakFormRHS(laplawfrhs_gradphi, gradphi[iscal2], vi);

          term_vi += DiffManager()->GetValence(iscal2) * DiffManager()->GetIsotropicDiff(iscal2) *
                     laplawfrhs_gradphi;
        }

        emat(vi * my::numdofpernode_ + k, ui * my::numdofpernode_ + iscal) +=
            -timefacfac * DiffManager()->GetPhasePoroTort(0) *
            (DiffManager()->GetDerivTransNum(k, iscal)) * my::funct_(ui) /
            DiffManager()->GetValence(k) * term_vi;
      }  // for(iscal)

      // formulation b):  plain ionic diffusion coefficients: one species eleminated by ENC
      //                  -> not implemented yet

      // formulation c):  plain ionic diffusion coefficients: sum (z_i D_i nabla c_i) replaced by
      // sum (t_i/c_i nabla c_i)
      //                  -> not activated
      //                  -> linearization is missing
      // emat(vi*my::numdofpernode_+k,ui*my::numdofpernode_+iscal)
      //    +=
      //    timefacfac*DiffManager()->GetPhasePoroTort(0)/frt/faraday/DiffManager()->GetValence(k)*DiffManager()->GetTransNum(k)*DiffManager()->GetCond()*DiffManager()->GetTransNum(iscal)/DiffManager()->GetValence(iscal]/conint[iscal]*laplawf;
    }  // end for ui
  }    // end for vi

  return;
}


/*---------------------------------------------------------------------------------------*
|  CalcMat: Potential equation div i inserted current - conc. overpotential   ehrl  02/14|
*----------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::CalcMatPotEquDiviConc(
    Epetra_SerialDenseMatrix& emat,  //!< element matrix to be filled
    const int k,                     //!< index of current scalar
    const double timefacfac,         //!< domain-integration factor times time-integration factor
    const double rtffc,              //!< RT/(F^2 Newman_const_c)
    const double rtf,                //!< RT/F
    const double invf,               //!< 1/F
    const double newman_const_a,     //!< Newman constant a
    const double newman_const_b,     //!< Newman constant b
    const LINALG::Matrix<my::nsd_, 1>& gradphi,  //!< gradient of concentration at GP
    const double conintinv                       //!< inverted concentration at GP
)
{
  for (unsigned vi = 0; vi < my::nen_; ++vi)
  {
    for (unsigned ui = 0; ui < my::nen_; ++ui)
    {
      double laplawf(0.0);
      my::GetLaplacianWeakForm(laplawf, ui, vi);

      if (diffcondmat_ == INPAR::ELCH::diffcondmat_newman)
      {
        // linearization of the diffusion overpotential term
        //
        // (grad w, RT/F^2 kappa (thermfactor) f(t_k) 1/c_k D nabla c_k)
        //
        emat(vi * my::numdofpernode_ + my::numscal_, ui * my::numdofpernode_ + k) +=
            timefacfac * DiffManager()->GetPhasePoroTort(0) * rtffc * DiffManager()->GetCond() *
            DiffManager()->GetThermFac() *
            (newman_const_a + (newman_const_b * DiffManager()->GetTransNum(k))) * conintinv *
            laplawf;

        // linearization of conduction term depending on the concentration is implemented
        // only for one species
        // otherwise you would need a second loop over the all scalars
        double laplawfrhs_gradphi(0.0);
        my::GetLaplacianWeakFormRHS(laplawfrhs_gradphi, gradphi, vi);

        // Linearization wrt ln c
        emat(vi * my::numdofpernode_ + my::numscal_, ui * my::numdofpernode_ + k) +=
            -timefacfac * DiffManager()->GetPhasePoroTort(0) * rtffc * DiffManager()->GetCond() *
            DiffManager()->GetThermFac() *
            (newman_const_a + (newman_const_b * DiffManager()->GetTransNum(k))) * conintinv *
            conintinv * laplawfrhs_gradphi * my::funct_(ui);

        // Linearization wrt kappa
        emat(vi * my::numdofpernode_ + my::numscal_, ui * my::numdofpernode_ + k) +=
            timefacfac * DiffManager()->GetPhasePoroTort(0) * rtffc * DiffManager()->GetThermFac() *
            (newman_const_a + (newman_const_b * DiffManager()->GetTransNum(k))) * conintinv *
            laplawfrhs_gradphi * DiffManager()->GetDerivCond(k) * my::funct_(ui);

        // Linearization wrt transference number
        emat(vi * my::numdofpernode_ + my::numscal_, ui * my::numdofpernode_ + k) +=
            timefacfac * DiffManager()->GetPhasePoroTort(0) * rtffc * DiffManager()->GetCond() *
            DiffManager()->GetThermFac() * conintinv * laplawfrhs_gradphi * newman_const_b *
            DiffManager()->GetDerivTransNum(k, k) * my::funct_(ui);

        // Linearization wrt thermodynamic factor
        emat(vi * my::numdofpernode_ + my::numscal_, ui * my::numdofpernode_ + k) +=
            timefacfac * DiffManager()->GetPhasePoroTort(0) * rtffc * DiffManager()->GetCond() *
            (newman_const_a + (newman_const_b * DiffManager()->GetTransNum(k))) * conintinv *
            laplawfrhs_gradphi * DiffManager()->GetDerivThermFac(k) * my::funct_(ui);
      }
      else if (diffcondmat_ == INPAR::ELCH::diffcondmat_ion)
      {
        if (diffcondparams_->DiffusionCoeffBased() == true)
          emat(vi * my::numdofpernode_ + my::numscal_, ui * my::numdofpernode_ + k) +=
              timefacfac * DiffManager()->GetPhasePoroTort(0) * DiffManager()->GetValence(k) *
              DiffManager()->GetIsotropicDiff(k) * laplawf;
        else
        {
          // Attention:
          // Full linearization of transference number, conductivity, ... is still missing
          emat(vi * my::numdofpernode_ + my::numscal_, ui * my::numdofpernode_ + k) +=
              timefacfac * DiffManager()->GetPhasePoroTort(0) * rtf * invf /
              DiffManager()->GetValence(k) * DiffManager()->GetCond() *
              DiffManager()->GetTransNum(k) * conintinv * laplawf;
        }
      }
      else
        dserror("Diffusion-Conduction material is not specified");
    }
  }

  return;
}


/*----------------------------------------------------------------------------------*
|  CalcMat: Potential equation div i without inserted current            ehrl  02/14|
*-----------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::CalcMatPotEquDivi(
    Epetra_SerialDenseMatrix& emat,  //!< element matrix to be filled
    const double timefacfac,         //!< domain-integration factor times time-integration factor
    const double invf                //!< 1/F
)
{
  for (unsigned vi = 0; vi < my::nen_; ++vi)
  {
    for (unsigned ui = 0; ui < my::nen_; ++ui)
    {
      for (unsigned idim = 0; idim < my::nsd_; ++idim)
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
        // version b: (phi, div Di) -> not partially integrated
        // emat(fvi,fui) += timefacfac*funct_(vi)*derxy_(idim,ui);
      }  // end for(idim)
    }    // end for(ui)
  }      // end for(vi)

  return;
}


/*----------------------------------------------------------------------------------*
|  CalcMat: Current equation current                                     ehrl  02/14|
*-----------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::CalcMatCurEquCur(
    Epetra_SerialDenseMatrix& emat,  //!< element matrix to be filled
    const double timefacfac,         //!< domain-integration factor times time-integration factor
    const double invf                //!< 1/F
)
{
  // (v, i)
  for (unsigned vi = 0; vi < my::nen_; ++vi)
  {
    for (unsigned ui = 0; ui < my::nen_; ++ui)
    {
      for (unsigned idim = 0; idim < my::nsd_; ++idim)
      {
        const int fvi = vi * my::numdofpernode_ + (my::numscal_ + 1) + idim;
        const int fui = ui * my::numdofpernode_ + (my::numscal_ + 1) + idim;

        emat(fvi, fui) += timefacfac * invf * my::funct_(vi) * my::funct_(ui);
      }
    }
  }

  return;
}


/*----------------------------------------------------------------------------------*
|  CalcMat: Current equation ohmic overpotential                         ehrl  02/14|
*-----------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::CalcMatCurEquOhm(
    Epetra_SerialDenseMatrix& emat,  //!< element matrix to be filled
    const double timefacfac,         //!< domain-integration factor times time-integration factor
    const double invf,               //!< 1/F
    const LINALG::Matrix<my::nsd_, 1>& gradpot  //!< gradient of potenial at GP
)
{
  // (v, kappa grad phi)
  for (unsigned vi = 0; vi < my::nen_; ++vi)
  {
    for (unsigned ui = 0; ui < my::nen_; ++ui)
    {
      for (unsigned idim = 0; idim < my::nsd_; ++idim)
      {
        const int fvi = vi * my::numdofpernode_ + (my::numscal_ + 1) + idim;
        const int fui = ui * my::numdofpernode_ + my::numscal_;

        emat(fvi, fui) += timefacfac * invf * DiffManager()->GetPhasePoroTort(0) * my::funct_(vi) *
                          DiffManager()->GetCond() * my::derxy_(idim, ui);

        // linearization of conductivity in the ohmic resistance term (current equation)
        //
        // (w, D(kappa(c)) grad phi)
        //
        for (int k = 0; k < my::numscal_; ++k)
          emat(fvi, ui * my::numdofpernode_ + k) +=
              timefacfac * invf * DiffManager()->GetPhasePoroTort(0) * my::funct_(vi) *
              DiffManager()->GetDerivCond(k) * my::funct_(ui) * gradpot(idim);
      }
    }
  }

  return;
}


/*----------------------------------------------------------------------------------*
|  CalcMat: Current equation concentration overpotential                 ehrl  02/14|
*-----------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::CalcMatCurEquConc(
    Epetra_SerialDenseMatrix& emat,  //!< element matrix to be filled
    const double timefacfac,         //!< domain-integration factor times time-integration factor
    const double rtf,                //!< RT/F
    const double rtffc,              //!< RT/(F^2 Newman_const_c)
    const std::vector<double>& invfval,                       //!< 1/(F z_k)
    const double newman_const_a,                              //!< Newman constant a
    const double newman_const_b,                              //!< Newman constant b
    const std::vector<LINALG::Matrix<my::nsd_, 1>>& gradphi,  //!< gradient of concentration at GP
    const std::vector<double>& conintinv                      //!< inverted concentration at GP
)
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
          if (diffcondmat_ == INPAR::ELCH::diffcondmat_newman)
          {
            emat(
                vi * my::numdofpernode_ + (my::numscal_ + 1) + idim, ui * my::numdofpernode_ + k) +=
                timefacfac * rtffc * DiffManager()->GetPhasePoroTort(0) * my::funct_(vi) *
                DiffManager()->GetCond() *
                (newman_const_a + (newman_const_b * DiffManager()->GetTransNum(k))) * conintinv[k] *
                my::derxy_(idim, ui);

            // Attention: Newman with current as solution variable
            // full linearization of transference number, conductivity, ... is still missing

            //  //linearization of coupling term in the current equation is still missing
            //  emat(vi*numdofpernode_+(numscal_+1)+idim,ui*numdofpernode_+k)
            //    +=
            //    timefacfac*my::funct_(vi)/frt_*DiffManager()->GetCond()*DiffManager()->GetTransNum(k)/DiffManager()->GetValence(k]/((-1)*conint_[k]*conint_[k])*my::funct_(ui)*gradphicoupling_[k](idim);
            //
            //  for(int k2=0; k2<numscal_;++k2)
            //  {
            //    //Check if necessary??
            //    emat(vi*numdofpernode_+(numscal_+1)+idim,ui*numdofpernode_+k)
            //       +=
            //       timefacfac*my::funct_(vi)/frt_*DiffManager()->GetDerivCond(k)*my::funct_(ui)*trans_[k2]*my::funct_(ui)/DiffManager()->GetValence(k2]/conint_[k2]*gradphicoupling_[k2](idim);
            //
            //
            //    emat(vi*numdofpernode_+(numscal_+1)+idim,ui*numdofpernode_+k)
            //      +=
            //      timefacfac*my::funct_(vi)/frt_*DiffManager()->GetCond()*(transderiv_[k2])[k]*my::funct_(ui)/DiffManager()->GetValence(k2]/conint_[k2]*gradphicoupling_[k2](idim);
            //  }
          }
          else if (diffcondmat_ == INPAR::ELCH::diffcondmat_ion)
          {
            if (diffcondparams_->DiffusionCoeffBased() == true)
              emat(vi * my::numdofpernode_ + (my::numscal_ + 1) + idim,
                  ui * my::numdofpernode_ + k) += timefacfac * DiffManager()->GetPhasePoroTort(0) *
                                                  my::funct_(vi) * DiffManager()->GetValence(k) *
                                                  DiffManager()->GetIsotropicDiff(k) *
                                                  my::derxy_(idim, ui);
            else
            {
              // linearization wrt nabla c_k
              emat(vi * my::numdofpernode_ + (my::numscal_ + 1) + idim,
                  ui * my::numdofpernode_ + k) +=
                  timefacfac * DiffManager()->GetPhasePoroTort(0) * invfval[k] * rtf *
                  my::funct_(vi) * DiffManager()->GetCond() * DiffManager()->GetTransNum(k) *
                  conintinv[k] * my::derxy_(idim, ui);

              // linearization wrt 1/c_k
              emat(vi * my::numdofpernode_ + (my::numscal_ + 1) + idim,
                  ui * my::numdofpernode_ + k) +=
                  -timefacfac * DiffManager()->GetPhasePoroTort(0) * rtf * invfval[k] *
                  DiffManager()->GetCond() * DiffManager()->GetTransNum(k) * conintinv[k] *
                  conintinv[k] * my::funct_(vi) * (gradphi[k])(idim)*my::funct_(ui);

              // linearization wrt kappa
              double term_vi = 0.0;
              for (int iscal = 0; iscal < my::numscal_; ++iscal)
              {
                term_vi += invfval[iscal] * DiffManager()->GetTransNum(iscal) * conintinv[iscal] *
                           my::funct_(vi) * (gradphi[iscal])(idim)*my::funct_(ui);
              }
              emat(vi * my::numdofpernode_ + (my::numscal_ + 1) + idim,
                  ui * my::numdofpernode_ + k) += timefacfac * DiffManager()->GetPhasePoroTort(0) *
                                                  DiffManager()->GetDerivCond(k) * rtf * term_vi;

              // linearization wrt transference number
              for (int iscal = 0; iscal < my::numscal_; ++iscal)
              {
                emat(vi * my::numdofpernode_ + (my::numscal_ + 1) + idim,
                    ui * my::numdofpernode_ + iscal) +=
                    timefacfac * DiffManager()->GetPhasePoroTort(0) * rtf * invfval[k] *
                    DiffManager()->GetCond() * (DiffManager()->GetDerivTransNum(k, iscal)) *
                    conintinv[k] * my::funct_(vi) * (gradphi[k])(idim)*my::funct_(ui);
              }
            }
          }
          else
            dserror("Diffusion-Conduction material is not specified");
        }
      }
    }  // for ui
  }    // for vi

  return;
}


/*-------------------------------------------------------------------------------------*
 |  CalcRhs: Conduction term with inserted current - ohmic overpotential    ehrl 11/13 |
 *-------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::CalcRhsCondOhm(
    Epetra_SerialDenseVector& erhs,  //!< element vector to be filled
    const int k,                     //!< index of current scalar
    const double rhsfac,   //!< time-integration factor for rhs times domain-integration factor
    const double invfval,  //!< 1/(F z_k)
    const LINALG::Matrix<my::nsd_, 1>& gradpot  //!< gradient of potenial at GP
)
{
  for (unsigned vi = 0; vi < my::nen_; ++vi)
  {
    // diffusive term
    double laplawfrhs_gradpot = 0.0;
    my::GetLaplacianWeakFormRHS(laplawfrhs_gradpot, gradpot, vi);
    erhs[vi * my::numdofpernode_ + k] -= rhsfac * DiffManager()->GetPhasePoroTort(0) *
                                         DiffManager()->GetTransNum(k) * DiffManager()->GetCond() *
                                         invfval * laplawfrhs_gradpot;
  }

  return;
}


/*-------------------------------------------------------------------------------------*
 |  CalcRhs: Conduction term with inserted current - conc. overpotential    ehrl 11/13 |
 *-------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::CalcRhsCondConc(
    Epetra_SerialDenseVector& erhs,  //!< element vector to be filled
    const int k,                     //!< index of current scalar
    const double rhsfac,    //!< time-integration factor for rhs times domain-integration factor
    const double rtffcval,  //!< RT/(F^2 Newman_const_c z_k)
    const double newman_const_a,                 //!< Newman constant a
    const double newman_const_b,                 //!< Newman constant b
    const LINALG::Matrix<my::nsd_, 1>& gradphi,  //!< gradient of concentration at GP
    const std::vector<double>& conintinv         //!< inverted concentration at GP
)
{
  for (unsigned vi = 0; vi < my::nen_; ++vi)
  {
    for (int iscal = 0; iscal < my::numscal_; ++iscal)
    {
      // diffusive term second
      double laplawfrhs_gradphi(0.0);
      my::GetLaplacianWeakFormRHS(laplawfrhs_gradphi, gradphi, vi);  // compute once, reuse below!

      // formulation a): plain ionic diffusion coefficients without using ENC
      //
      // (grad w, sum(D_i grad Dc))
      erhs[vi * my::numdofpernode_ + k] -=
          rhsfac * DiffManager()->GetPhasePoroTort(0) * rtffcval * DiffManager()->GetTransNum(k) *
          DiffManager()->GetCond() * DiffManager()->GetThermFac() *
          (newman_const_a + (newman_const_b * DiffManager()->GetTransNum(iscal))) *
          conintinv[iscal] * laplawfrhs_gradphi;

      // formulation b):  plain ionic diffusion coefficients: one species eleminated by ENC
      //                  -> not implemented yet
    }
  }

  return;
}

/*-------------------------------------------------------------------------------------*
 |  CalcRhs: Conduction term without inserted current                       ehrl 11/13 |
 *-------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::CalcRhsCond(
    Epetra_SerialDenseVector& erhs,  //!< element vector to be filled
    const int k,                     //!< index of current scalar
    const double rhsfac,   //!< time-integration factor for rhs times domain-integration factor
    const double invfval,  //!< 1/(F z_k)
    const LINALG::Matrix<my::nsd_, 1>& curint  //!< current at GP
)
{
  for (unsigned vi = 0; vi < my::nen_; ++vi)
  {
    double laplawfrhs_cur = 0.0;
    my::GetLaplacianWeakFormRHS(laplawfrhs_cur, curint, vi);

    erhs[vi * my::numdofpernode_ + k] -=
        -rhsfac * DiffManager()->GetTransNum(k) * invfval * laplawfrhs_cur;
  }

  return;
}


/*-------------------------------------------------------------------------------------*
 |  CalcRhs: Additional diffusion term without inserted current             ehrl 11/13 |
 *-------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::CalcRhsCondDiff(
    Epetra_SerialDenseVector& erhs,  //!< element vector to be filled
    const int k,                     //!< index of current scalar
    const double rhsfac,  //!< time-integration factor for rhs times domain-integration factor
    const std::vector<LINALG::Matrix<my::nsd_, 1>>& gradphi  //!< gradient of concentration at GP
)
{
  for (unsigned vi = 0; vi < my::nen_; ++vi)
  {
    for (int iscal = 0; iscal < my::numscal_; ++iscal)
    {
      // diffusive term second
      double laplawfrhs_gradphi(0.0);
      my::GetLaplacianWeakFormRHS(
          laplawfrhs_gradphi, gradphi[iscal], vi);  // compute once, reuse below!

      // formulation a:  plain ionic diffusion coefficients: sum (z_i D_i nabla c_i)
      erhs[vi * my::numdofpernode_ + k] -=
          -rhsfac * DiffManager()->GetPhasePoroTort(0) * DiffManager()->GetTransNum(k) *
          DiffManager()->GetValence(iscal) / DiffManager()->GetValence(k) *
          DiffManager()->GetIsotropicDiff(iscal) * laplawfrhs_gradphi;

      // formulation b):  plain ionic diffusion coefficients: one species eleminated by ENC
      //                  -> not implemented yet
      // formulation c:  plain ionic diffusion coefficients: sum (z_i D_i nabla c_i) replaced by sum
      // (t_i/c_i nabla c_i)
      //                  -> not activated
      // erhs[fvi]
      //   -=
      //   rhsfac*DiffManager()->GetPhasePoroTort(0)/frt/faraday/DiffManager()->GetValence(k)*DiffManager()->GetTransNum(k)*DiffManager()->GetCond()*DiffManager()->GetTransNum(iscal)/DiffManager()->GetValence(iscal]/conint[iscal]*laplawf2;
    }
  }

  return;
}


/*-------------------------------------------------------------------------------------*
 |CalcRhs: Potential equation div i inserted current - conc. overpotential  ehrl 11/13 |
 *-------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::CalcRhsPotEquDiviConc(
    Epetra_SerialDenseVector& erhs,  //!< element vector to be filled
    const int k,                     //!< index of current scalar
    const double rhsfac,  //!< time-integration factor for rhs times domain-integration factor
    const double rtf,     //!< RT/F
    const std::vector<double>& invfval,          //!< 1/(F z_k)
    const double rtffc,                          //!< RT/(F^2 Newman_const_c)
    const double newman_const_a,                 //!< Newman constant a
    const double newman_const_b,                 //!< Newman constant b
    const LINALG::Matrix<my::nsd_, 1>& gradphi,  //!< gradient of concentration at GP
    const double conintinv                       //!< inverted concentration at GP
)
{
  for (unsigned vi = 0; vi < my::nen_; ++vi)
  {
    if (diffcondmat_ == INPAR::ELCH::diffcondmat_ion)
    {
      // diffusive term second
      double laplawf2(0.0);
      my::GetLaplacianWeakFormRHS(laplawf2, gradphi, vi);  // compute once, reuse below!

      if (diffcondparams_->DiffusionCoeffBased() == true)
        erhs[vi * my::numdofpernode_ + my::numscal_] -=
            rhsfac * DiffManager()->GetPhasePoroTort(0) * DiffManager()->GetValence(k) *
            DiffManager()->GetIsotropicDiff(k) * laplawf2;
      else
        erhs[vi * my::numdofpernode_ + my::numscal_] -=
            rhsfac * DiffManager()->GetPhasePoroTort(0) * rtf * invfval[k] *
            DiffManager()->GetCond() * DiffManager()->GetTransNum(k) * conintinv * laplawf2;
    }
    // thermodynamic factor only implemented for Newman
    else if (diffcondmat_ == INPAR::ELCH::diffcondmat_newman)
    {
      // diffusive term second
      double laplawf2(0.0);
      my::GetLaplacianWeakFormRHS(laplawf2, gradphi, vi);  // compute once, reuse below!

      erhs[vi * my::numdofpernode_ + my::numscal_] -=
          rhsfac * DiffManager()->GetPhasePoroTort(0) * rtffc * DiffManager()->GetCond() *
          DiffManager()->GetThermFac() *
          (newman_const_a + (newman_const_b * DiffManager()->GetTransNum(k))) * conintinv *
          laplawf2;
    }
    else
      dserror("Diffusion-Conduction material is not specified");
  }

  return;
}


/*-------------------------------------------------------------------------------------*
 |  CalcRhs: Potential equation divi without inserted current               ehrl 11/13 |
 *-------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::CalcRhsPotEquDivi(
    Epetra_SerialDenseVector& erhs,  //!< element vector to be filled
    const double rhsfac,  //!< time-integration factor for rhs times domain-integration factor
    const double invf,    //!< 1/F
    const LINALG::Matrix<my::nsd_, 1>& curint  //!< current at GP
)
{
  for (unsigned vi = 0; vi < my::nen_; ++vi)
  {
    double laplawf = 0.0;
    // version a: (grad phi,  Di)
    my::GetLaplacianWeakFormRHS(laplawf, curint, vi);
    erhs[vi * my::numdofpernode_ + my::numscal_] -= -rhsfac * invf * laplawf;
    // version b: (phi, div Di) -> not partially integrated
    // erhs[vi*my::numdofpernode_+my::numscal_] -= rhsfac*my::funct_(vi)*divi;
  }

  return;
}


/*-------------------------------------------------------------------------------------*
 |  CalcRhs: Current equation - current                                     ehrl 11/13 |
 *-------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::CalcRhsCurEquCur(
    Epetra_SerialDenseVector& erhs,  //!< element vector to be filled
    const double rhsfac,  //!< time-integration factor for rhs times domain-integration factor
    const double invf,    //!< 1/F
    const LINALG::Matrix<my::nsd_, 1>& curint  //!< current at GP
)
{
  for (unsigned vi = 0; vi < my::nen_; ++vi)
  {
    for (unsigned idim = 0; idim < my::nsd_; ++idim)
    {
      // (v, i)
      erhs[vi * my::numdofpernode_ + (my::numscal_ + 1) + idim] -=
          rhsfac * invf * my::funct_(vi) * curint(idim);
    }
  }

  return;
}

/*-------------------------------------------------------------------------------------*
 |  CalcRhs: Current equation - ohmic overpotential                         ehrl 11/13 |
 *-------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::CalcRhsCurEquOhm(
    Epetra_SerialDenseVector& erhs,  //!< element vector to be filled
    const double rhsfac,  //!< time-integration factor for rhs times domain-integration factor
    const double invf,    //!< 1/F
    const LINALG::Matrix<my::nsd_, 1>& gradpot  //!< gradient of potenial at GP
)
{
  for (unsigned vi = 0; vi < my::nen_; ++vi)
  {
    for (unsigned idim = 0; idim < my::nsd_; ++idim)
    {
      // (v, kappa grad phi)
      erhs[vi * my::numdofpernode_ + (my::numscal_ + 1) + idim] -=
          rhsfac * invf * DiffManager()->GetPhasePoroTort(0) * my::funct_(vi) *
          DiffManager()->GetCond() * gradpot(idim);
    }
  }

  return;
}


/*-------------------------------------------------------------------------------------*
 |  CalcRhs: Current equation - concentration overpotential                 ehrl 11/13 |
 *-------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::CalcRhsCurEquConc(
    Epetra_SerialDenseVector& erhs,  //!< element vector to be filled
    const double rhsfac,  //!< time-integration factor for rhs times domain-integration factor
    const double rtf,     //!< RT/F
    const std::vector<double>& invfval,  //!< 1/(F z_k)
    const double rtffc,                  //!< RT/(F^2 Newman_const_c)
    const double newman_const_a,         //!< Newman constant a
    const double newman_const_b,         //!< Newman constant b
    const std::vector<LINALG::Matrix<my::nsd_, 1>>&
        gradphi,                          //!< vector of gradient of concentration at GP
    const std::vector<double>& conintinv  //!< inverted concentration at GP
)
{
  for (unsigned vi = 0; vi < my::nen_; ++vi)
  {
    for (unsigned idim = 0; idim < my::nsd_; ++idim)
    {
      for (int k = 0; k < my::numscal_; ++k)
      {
        if (diffcondmat_ == INPAR::ELCH::diffcondmat_newman)
        {
          erhs[vi * my::numdofpernode_ + (my::numscal_ + 1) + idim] -=
              rhsfac * DiffManager()->GetPhasePoroTort(0) * my::funct_(vi) * rtffc *
              DiffManager()->GetCond() *
              (newman_const_a + (newman_const_b * DiffManager()->GetTransNum(k))) * conintinv[k] *
              gradphi[k](idim);
        }
        else if (diffcondmat_ == INPAR::ELCH::diffcondmat_ion)
        {
          if (diffcondparams_->DiffusionCoeffBased() == true)
            erhs[vi * my::numdofpernode_ + (my::numscal_ + 1) + idim] -=
                rhsfac * DiffManager()->GetPhasePoroTort(0) * my::funct_(vi) *
                DiffManager()->GetValence(k) * DiffManager()->GetIsotropicDiff(k) *
                gradphi[k](idim);
          else
          {
            erhs[vi * my::numdofpernode_ + (my::numscal_ + 1) + idim] -=
                rhsfac * DiffManager()->GetPhasePoroTort(0) * my::funct_(vi) * rtf *
                DiffManager()->GetCond() * invfval[k] * DiffManager()->GetTransNum(k) *
                conintinv[k] * gradphi[k](idim);
          }
        }
        else
          dserror("Diffusion-Conduction material is not specified");
      }
    }
  }

  return;
}


/*----------------------------------------------------------------------*
|  Correct sysmat for fluxes across DC                       ehrl  02/14|
*----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::CorrectionForFluxAcrossDC(
    DRT::Discretization& discretization, const std::vector<int>& lm, Epetra_SerialDenseMatrix& emat,
    Epetra_SerialDenseVector& erhs)
{
  // get dirichlet toggle from the discretization
  // we always get the dirichet toggle:
  // in this function we check if the actual nodes have a dirichlet value
  Teuchos::RCP<const Epetra_Vector> dctoggle = discretization.GetState("dctoggle");
  std::vector<double> mydctoggle(lm.size());
  DRT::UTILS::ExtractMyValues(*dctoggle, mydctoggle, lm);

  double val = 0.0;
  for (unsigned vi = 0; vi < my::nen_; ++vi)
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
        erhs[vi * my::numdofpernode_ + my::numscal_] += DiffManager()->GetValence(k) * (-val);
        // corresponding linearization
        for (unsigned ui = 0; ui < my::nen_; ++ui)
        {
          val = emat(vi * my::numdofpernode_ + k, ui * my::numdofpernode_ + k);
          emat(vi * my::numdofpernode_ + my::numscal_, ui * my::numdofpernode_ + k) +=
              DiffManager()->GetValence(k) * (-val);
          val = emat(vi * my::numdofpernode_ + k, ui * my::numdofpernode_ + my::numscal_);
          emat(vi * my::numdofpernode_ + my::numscal_, ui * my::numdofpernode_ + my::numscal_) +=
              DiffManager()->GetValence(k) * (-val);
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
          int k = 0;

          const int fvi = vi * my::numdofpernode_ + my::numscal_;
          // We use the fact, that the rhs vector value for boundary nodes
          // is equivalent to the integrated negative normal flux
          // due to diffusion and migration

          // scaling of div i results in a matrix with better condition number
          val = erhs[fvi];
          erhs[vi * my::numdofpernode_ + k] += 1.0 / DiffManager()->GetValence(k) * (-val);
          // corresponding linearization
          for (unsigned ui = 0; ui < my::nen_; ++ui)
          {
            val = emat(vi * my::numdofpernode_ + my::numscal_, ui * my::numdofpernode_ + k);
            emat(vi * my::numdofpernode_ + k, ui * my::numdofpernode_ + k) +=
                1.0 / DiffManager()->GetValence(k) * (-val);
            val = emat(
                vi * my::numdofpernode_ + my::numscal_, ui * my::numdofpernode_ + my::numscal_);
            emat(vi * my::numdofpernode_ + k, ui * my::numdofpernode_ + my::numscal_) +=
                1.0 / DiffManager()->GetValence(k) * (-val);
          }
        }
      }
    }  // for k
  }

  return;
}


/*----------------------------------------------------------------------*
 | get material parameters                                   fang 07/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::GetMaterialParams(const DRT::Element* ele,
    std::vector<double>& densn, std::vector<double>& densnp, std::vector<double>& densam,
    double& visc, const int iquad)
{
  // extract material from element
  Teuchos::RCP<MAT::Material> material = ele->Material();

  // evaluate electrolyte material
  if (material->MaterialType() == INPAR::MAT::m_elchmat)
    Utils()->MatElchMat(material, VarManager()->Phinp(), VarManager()->Temperature(),
        myelch::elchparams_->EquPot(), myelch::elchparams_->Faraday() * VarManager()->FRT(),
        DiffManager(), diffcondmat_);
  else
    dserror("Invalid material type!");
}  // DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::GetMaterialParams


// template classes
// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<DRT::Element::line2>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<DRT::Element::line3>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<DRT::Element::tri3>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<DRT::Element::quad4>;
// template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<DRT::Element::quad9>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<DRT::Element::nurbs9>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<DRT::Element::hex8>;
// template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<DRT::Element::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<DRT::Element::hex27>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<DRT::Element::tet4>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<DRT::Element::tet10>;
// template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<DRT::Element::wedge6>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<DRT::Element::pyramid5>;
// template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<DRT::Element::nurbs27>;

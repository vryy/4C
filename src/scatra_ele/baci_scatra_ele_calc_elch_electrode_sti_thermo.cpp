/*--------------------------------------------------------------------------*/
/*! \file

\brief evaluation of scatra elements for conservation of mass concentration and electronic charge
within thermodynamic electrodes

\level 2

*/
/*--------------------------------------------------------------------------*/
#include "baci_scatra_ele_calc_elch_electrode_sti_thermo.H"

#include "baci_scatra_ele_parameter_timint.H"
#include "baci_utils_singleton_owner.H"

/*----------------------------------------------------------------------*
 | singleton access method                                   fang 11/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<distype>*
DRT::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<distype>::Instance(
    const int numdofpernode, const int numscal, const std::string& disname)
{
  static auto singleton_map = CORE::UTILS::MakeSingletonMap<std::string>(
      [](const int numdofpernode, const int numscal, const std::string& disname)
      {
        return std::unique_ptr<ScaTraEleCalcElchElectrodeSTIThermo<distype>>(
            new ScaTraEleCalcElchElectrodeSTIThermo<distype>(numdofpernode, numscal, disname));
      });

  return singleton_map[disname].Instance(
      CORE::UTILS::SingletonAction::create, numdofpernode, numscal, disname);
}


/*----------------------------------------------------------------------*
 | extract quantities for element evaluation                 fang 11/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<distype>::ExtractElementAndNodeValues(
    DRT::Element* ele, Teuchos::ParameterList& params, DRT::Discretization& discretization,
    DRT::Element::LocationArray& la)
{
  // call base class routine to extract scatra-related quantities
  myelch::ExtractElementAndNodeValues(ele, params, discretization, la);

  // call base class routine to extract thermo-related quantitites
  mythermo::ExtractElementAndNodeValues(ele, params, discretization, la);
}


/*----------------------------------------------------------------------*
 | get material parameters                                   fang 11/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<distype>::GetMaterialParams(
    const DRT::Element* ele, std::vector<double>& densn, std::vector<double>& densnp,
    std::vector<double>& densam, double& visc, const int iquad)
{
  // Set GP values to MatElectrode
  myelectrode::Utils()->MatElectrode(
      ele->Material(), VarManager()->Phinp(0), VarManager()->Temp(), myelectrode::DiffManager());

  // get parameters of secondary, thermodynamic electrolyte material
  Teuchos::RCP<const MAT::Material> material = ele->Material(1);
  materialtype_ = material->MaterialType();
  if (materialtype_ == INPAR::MAT::m_soret) mythermo::MatSoret(material);
}  // DRT::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<distype>::GetMaterialParams


/*--------------------------------------------------------------------------*
 | calculate element matrix and element right-hand side vector   fang 11/15 |
 *--------------------------------------------------------------------------*/

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<distype>::CalcMatAndRhs(
    CORE::LINALG::SerialDenseMatrix& emat, CORE::LINALG::SerialDenseVector& erhs, const int k,
    const double fac, const double timefacfac, const double rhsfac, const double taufac,
    const double timetaufac, const double rhstaufac, CORE::LINALG::Matrix<nen_, 1>& tauderpot,
    double& rhsint)
{
  // call base class routine for isothermal problems
  myelectrode::CalcMatAndRhs(
      emat, erhs, k, fac, timefacfac, rhsfac, taufac, timetaufac, rhstaufac, tauderpot, rhsint);

  if (materialtype_ == INPAR::MAT::m_soret)
  {
    // matrix and vector contributions arising from additional, thermodynamic term for Soret effect
    mythermo::CalcMatSoret(emat, timefacfac, VarManager()->Phinp(0),
        myelectrode::DiffManager()->GetIsotropicDiff(0),
        myelectrode::DiffManager()->GetConcDerivIsoDiffCoef(0, 0), VarManager()->Temp(),
        VarManager()->GradTemp(), my::funct_, my::derxy_);
    mythermo::CalcRHSSoret(erhs, VarManager()->Phinp(0),
        myelectrode::DiffManager()->GetIsotropicDiff(0), rhsfac, VarManager()->Temp(),
        VarManager()->GradTemp(), my::derxy_);
  }
}


/*----------------------------------------------------------------------*
 | evaluate action for off-diagonal system matrix block      fang 11/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<distype>::EvaluateActionOD(DRT::Element* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization,
    const SCATRA::Action& action, DRT::Element::LocationArray& la,
    CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
    CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
    CORE::LINALG::SerialDenseVector& elevec1_epetra,
    CORE::LINALG::SerialDenseVector& elevec2_epetra,
    CORE::LINALG::SerialDenseVector& elevec3_epetra)
{
  // determine and evaluate action
  switch (action)
  {
    case SCATRA::Action::calc_scatra_mono_odblock_scatrathermo:
    {
      SysmatODScatraThermo(ele, elemat1_epetra);

      break;
    }

    default:
    {
      // call base class routine
      my::EvaluateActionOD(ele, params, discretization, action, la, elemat1_epetra, elemat2_epetra,
          elevec1_epetra, elevec2_epetra, elevec3_epetra);

      break;
    }
  }  // switch(action)

  return 0;
}


/*------------------------------------------------------------------------------------------------------*
 | fill element matrix with linearizations of discrete scatra residuals w.r.t. thermo dofs   fang
 11/15 |
 *------------------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<distype>::SysmatODScatraThermo(
    DRT::Element* ele, CORE::LINALG::SerialDenseMatrix& emat)
{
  // integration points and weights
  CORE::DRT::UTILS::IntPointsAndWeights<nsd_ele_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // loop over integration points
  for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
  {
    // evaluate shape functions, their derivatives, and domain integration factor at current
    // integration point
    const double fac = my::EvalShapeFuncAndDerivsAtIntPoint(intpoints, iquad);

    // evaluate overall integration factor
    const double timefacfac = my::scatraparatimint_->TimeFac() * fac;

    // evaluate internal variables at current integration point
    SetInternalVariablesForMatAndRHS();

    // evaluate material parameters at current integration point
    double dummy(0.);
    std::vector<double> dummyvec(my::numscal_, 0.);
    GetMaterialParams(ele, dummyvec, dummyvec, dummyvec, dummy, iquad);

    // calculating the off diagonal for the temperature derivative of concentration and electric
    // potential
    mythermo::CalcMatDiffThermoOD(emat, my::numdofpernode_, timefacfac, VarManager()->InvF(),
        VarManager()->GradPhi(0), VarManager()->GradPot(),
        myelectrode::DiffManager()->GetTempDerivIsoDiffCoef(0, 0),
        myelectrode::DiffManager()->GetTempDerivCond(0), my::funct_, my::derxy_, 1.);

    if (materialtype_ == INPAR::MAT::m_soret)
    {
      // provide element matrix with linearizations of Soret term in discrete scatra residuals
      // w.r.t. thermo dofs
      mythermo::CalcMatSoretOD(emat, timefacfac, VarManager()->Phinp(0),
          myelectrode::DiffManager()->GetIsotropicDiff(0), VarManager()->Temp(),
          VarManager()->GradTemp(), my::funct_, my::derxy_);
    }
  }
}


/*------------------------------------------------------------------------------*
 | set internal variables for element evaluation                     fang 11/15 |
 *------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<distype>::SetInternalVariablesForMatAndRHS()
{
  // set internal variables for element evaluation
  VarManager()->SetInternalVariables(my::funct_, my::derxy_, my::ephinp_, my::ephin_,
      mythermo::etempnp_, my::econvelnp_, my::ehist_);
}

/*----------------------------------------------------------------------*
 | private constructor for singletons                        fang 11/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<distype>::ScaTraEleCalcElchElectrodeSTIThermo(
    const int numdofpernode, const int numscal, const std::string& disname)
    :  // constructors of base classes
      ScaTraEleCalcElchElectrode<distype>::ScaTraEleCalcElchElectrode(
          numdofpernode, numscal, disname),
      ScaTraEleSTIThermo<distype>::ScaTraEleSTIThermo(numscal)
{
  // safety check
  if (numscal != 1 or numdofpernode != 2)
    dserror("Invalid number of transported scalars or degrees of freedom per node!");

  // replace internal variable manager for isothermal electrodes by internal variable manager for
  // thermodynamic electrodes
  my::scatravarmanager_ =
      Teuchos::rcp(new ScaTraEleInternalVariableManagerElchElectrodeSTIThermo<nsd_, nen_>(
          my::numscal_, myelch::elchparams_));
}


// template classes
// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<
    DRT::Element::DiscretizationType::line2>;
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<
    DRT::Element::DiscretizationType::line3>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<
    DRT::Element::DiscretizationType::tri3>;
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<
    DRT::Element::DiscretizationType::tri6>;
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<
    DRT::Element::DiscretizationType::quad4>;
// template class
// DRT::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<DRT::Element::DiscretizationType::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<
    DRT::Element::DiscretizationType::quad9>;
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<
    DRT::Element::DiscretizationType::nurbs9>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<
    DRT::Element::DiscretizationType::hex8>;
// template class
// DRT::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<DRT::Element::DiscretizationType::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<
    DRT::Element::DiscretizationType::hex27>;
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<
    DRT::Element::DiscretizationType::tet4>;
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<
    DRT::Element::DiscretizationType::tet10>;
// template class
// DRT::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<DRT::Element::DiscretizationType::wedge6>;
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<
    DRT::Element::DiscretizationType::pyramid5>;
// template class
// DRT::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<DRT::Element::DiscretizationType::nurbs27>;

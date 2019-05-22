/*--------------------------------------------------------------------------*/
/*!

\brief evaluation of scatra elements for thermodynamic diffusion-conduction ion-transport equations

\level 2

\maintainer Christoph Schmidt
*/
/*--------------------------------------------------------------------------*/
#include "scatra_ele_calc_elch_diffcond_sti_thermo.H"

#include "scatra_ele_parameter_timint.H"
#include "scatra_ele_parameter_elch.H"

#include "../drt_mat/material.H"

/*----------------------------------------------------------------------*
 | singleton access method                                   fang 11/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcElchDiffCondSTIThermo<distype>*
DRT::ELEMENTS::ScaTraEleCalcElchDiffCondSTIThermo<distype>::Instance(const int numdofpernode,
    const int numscal, const std::string& disname,
    const ScaTraEleCalcElchDiffCondSTIThermo* delete_me)
{
  static std::map<std::string, ScaTraEleCalcElchDiffCondSTIThermo<distype>*> instances;

  if (delete_me == NULL)
  {
    if (instances.find(disname) == instances.end())
      instances[disname] =
          new ScaTraEleCalcElchDiffCondSTIThermo<distype>(numdofpernode, numscal, disname);
  }

  else
  {
    for (typename std::map<std::string, ScaTraEleCalcElchDiffCondSTIThermo<distype>*>::iterator i =
             instances.begin();
         i != instances.end(); ++i)
      if (i->second == delete_me)
      {
        delete i->second;
        instances.erase(i);
        return NULL;
      }
  }

  return instances[disname];
}


/*----------------------------------------------------------------------*
 | singleton destruction                                     fang 11/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCondSTIThermo<distype>::Done()
{
  // delete singleton
  Instance(0, 0, "", this);

  return;
}


/*----------------------------------------------------------------------*
 | extract quantities for element evaluation                 fang 11/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCondSTIThermo<distype>::ExtractElementAndNodeValues(
    DRT::Element* ele,                    //!< current element
    Teuchos::ParameterList& params,       //!< parameter list
    DRT::Discretization& discretization,  //!< discretization
    DRT::Element::LocationArray& la       //!< location array
)
{
  // call base class routine to extract scatra-related quantities
  mydiffcond::ExtractElementAndNodeValues(ele, params, discretization, la);

  // call base class routine to extract thermo-related quantitites
  mythermo::ExtractElementAndNodeValues(ele, params, discretization, la);

  return;
}


/*----------------------------------------------------------------------*
 | get material parameters                                   fang 11/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCondSTIThermo<distype>::GetMaterialParams(
    const DRT::Element* ele,      //!< current element
    std::vector<double>& densn,   //!< density at t_(n)
    std::vector<double>& densnp,  //!< density at t_(n+1) or t_(n+alpha_F)
    std::vector<double>& densam,  //!< density at t_(n+alpha_M)
    double& visc,                 //!< fluid viscosity
    const int iquad               //!< ID of current integration point
)
{
  // call base class routine to get parameters of primary, isothermal electrolyte material
  mydiffcond::GetMaterialParams(ele, densn, densnp, densam, visc, iquad);

  // get parameters of secondary, thermodynamic electrolyte material
  Teuchos::RCP<const MAT::Material> material = ele->Material(1);
  if (material->MaterialType() == INPAR::MAT::m_soret)
    mythermo::MatSoret(material);
  else
    dserror("Invalid electrolyte material!");

  return;
}  // DRT::ELEMENTS::ScaTraEleCalcElchDiffCondSTIThermo<distype>::GetMaterialParams


/*--------------------------------------------------------------------------*
 | calculate element matrix and element right-hand side vector   fang 11/15 |
 *--------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCondSTIThermo<distype>::CalcMatAndRhs(
    Epetra_SerialDenseMatrix& emat,  //!< element matrix
    Epetra_SerialDenseVector& erhs,  //!< element right-hand side vector
    const int k,                     //!< index of current scalar
    const double fac,                //!< domain integration factor
    const double timefacfac,         //!< domain integration factor times time integration factor
    const double rhsfac,      //!< domain integration factor times time integration factor for
                              //!< right-hand side vector
    const double taufac,      //!< domain integration factor times stabilization parameter
    const double timetaufac,  //!< domain integration factor times stabilization parameter times
                              //!< time integration factor
    const double rhstaufac,  //!< domain integration factor times stabilization parameter times time
                             //!< integration factor for right-hand side vector
    LINALG::Matrix<my::nen_, 1>&
        tauderpot,  //!< derivatives of stabilization parameter w.r.t. electric potential
    double& rhsint  //!< body force value
)
{
  // call base class routine for isothermal problems
  mydiffcond::CalcMatAndRhs(
      emat, erhs, k, fac, timefacfac, rhsfac, taufac, timetaufac, rhstaufac, tauderpot, rhsint);

  // extract variables and parameters
  const double& concentration = VarManager()->Phinp(0);
  const LINALG::Matrix<my::nsd_, 1>& gradtemp = VarManager()->GradTemp();
  const double& kappa = mydiffcond::DiffManager()->GetCond();
  const double& kappaderiv = mydiffcond::DiffManager()->GetDerivCond(0);
  const double faraday = DRT::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->Faraday();
  const double invffval = mydiffcond::DiffManager()->InvFVal(0) / faraday;
  const double& invfval = mydiffcond::DiffManager()->InvFVal(0);
  const double& R = DRT::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->GasConstant();
  const double& t = mydiffcond::DiffManager()->GetTransNum(0);
  const double& tderiv = mydiffcond::DiffManager()->GetDerivTransNum(0, 0);

  // matrix and vector contributions arising from additional, thermodynamic term in expression for
  // current density
  for (unsigned vi = 0; vi < my::nen_; ++vi)
  {
    // recurring indices
    const int rowconc(vi * 2);
    const int rowpot(vi * 2 + 1);

    // gradient of test function times gradient of temperature
    double laplawfrhs_temp(0.);
    my::GetLaplacianWeakFormRHS(laplawfrhs_temp, gradtemp, vi);

    for (unsigned ui = 0; ui < my::nen_; ++ui)
    {
      // recurring index
      const int colconc(ui * 2);

      // linearizations of contributions for concentration residuals w.r.t. concentration dofs
      emat(rowconc, colconc) -=
          timefacfac * my::funct_(ui) * laplawfrhs_temp * pow(invfval, 2.0) *
          (tderiv * kappa * R * log(concentration) + t * kappaderiv * R * log(concentration) +
              t * kappa * R / concentration);

      // linearizations of contributions for electric potential residuals w.r.t. concentration dofs
      emat(rowpot, colconc) -= timefacfac * my::funct_(ui) * laplawfrhs_temp * invffval *
                               (kappaderiv * R * log(concentration) + kappa * R / concentration);

      // linearizations w.r.t. electric potential dofs are zero
    }

    // contribution for concentration residual
    erhs[rowconc] +=
        rhsfac * laplawfrhs_temp * kappa * pow(invfval, 2.0) * t * R * log(concentration);

    // contribution for electric potential residual
    erhs[rowpot] += rhsfac * laplawfrhs_temp * kappa * invffval * R * log(concentration);
  }

  // matrix and vector contributions arising from additional, thermodynamic term for Soret effect
  mythermo::CalcMatSoret(emat, timefacfac, VarManager()->Phinp(0),
      mydiffcond::DiffManager()->GetIsotropicDiff(0),
      mydiffcond::DiffManager()->GetDerivIsoDiffCoef(0, 0), VarManager()->Temp(),
      VarManager()->GradTemp(), my::funct_, my::derxy_);
  mythermo::CalcRHSSoret(erhs, VarManager()->Phinp(0),
      mydiffcond::DiffManager()->GetIsotropicDiff(0), rhsfac, VarManager()->Temp(),
      VarManager()->GradTemp(), my::derxy_);

  return;
}


/*----------------------------------------------------------------------*
 | evaluate action for off-diagonal system matrix block      fang 11/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ScaTraEleCalcElchDiffCondSTIThermo<distype>::EvaluateActionOD(
    DRT::Element* ele,                         //!< current element
    Teuchos::ParameterList& params,            //!< parameter list
    DRT::Discretization& discretization,       //!< discretization
    const SCATRA::Action& action,              //!< action parameter
    DRT::Element::LocationArray& la,           //!< location array
    Epetra_SerialDenseMatrix& elemat1_epetra,  //!< element matrix 1
    Epetra_SerialDenseMatrix& elemat2_epetra,  //!< element matrix 2
    Epetra_SerialDenseVector& elevec1_epetra,  //!< element right-hand side vector 1
    Epetra_SerialDenseVector& elevec2_epetra,  //!< element right-hand side vector 2
    Epetra_SerialDenseVector& elevec3_epetra   //!< element right-hand side vector 3
)
{
  // determine and evaluate action
  switch (action)
  {
    case SCATRA::calc_scatra_mono_odblock_scatrathermo:
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
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCondSTIThermo<distype>::SysmatODScatraThermo(
    DRT::Element* ele,              //!< current element
    Epetra_SerialDenseMatrix& emat  //!< element matrix
)
{
  // integration points and weights
  DRT::UTILS::IntPointsAndWeights<my::nsd_ele_> intpoints(
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

    // extract variables and parameters
    const double& concentration = VarManager()->Phinp(0);
    const LINALG::Matrix<my::nsd_, 1>& gradconc = VarManager()->GradPhi(0);
    const double faraday = DRT::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->Faraday();
    const double& invffval = mydiffcond::DiffManager()->InvFVal(0) / faraday;
    const double& invfval = mydiffcond::DiffManager()->InvFVal(0);
    const double& kappa = mydiffcond::DiffManager()->GetCond();
    const double& R = DRT::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->GasConstant();
    const double& t = mydiffcond::DiffManager()->GetTransNum(0);

    // matrix contributions arising from additional, thermodynamic term in expression for current
    // density
    for (unsigned vi = 0; vi < my::nen_; ++vi)
    {
      // recurring indices
      const int rowconc(vi * 2);
      const int rowpot(vi * 2 + 1);

      // gradient of test function times gradient of concentration
      double laplawfrhs_conc(0.);
      my::GetLaplacianWeakFormRHS(laplawfrhs_conc, gradconc, vi);

      for (unsigned ui = 0; ui < my::nen_; ++ui)
      {
        // gradient of test function times gradient of shape function
        double laplawf(0.);
        my::GetLaplacianWeakForm(laplawf, vi, ui);

        // gradient of test function times derivative of current density w.r.t. temperature
        const double di_dT = 2 * kappa * (1 - t) * invfval * R / concentration * laplawfrhs_conc;

        // formal, symbolic derivative of current density w.r.t. temperature gradient
        const double di_dgradT = kappa * invfval * R * log(concentration);

        // linearizations of contributions for concentration residuals w.r.t. thermo dofs
        emat(rowconc, ui) -=
            timefacfac * (t * invfval * di_dT * my::funct_(ui) + t * invfval * di_dgradT * laplawf);

        // linearizations of contributions for electric potential residuals w.r.t. thermo dofs
        emat(rowpot, ui) -= timefacfac * (laplawf * kappa * invffval * R * log(concentration) +
                                             my::funct_(ui) * kappa * (1 - t) * invffval * R /
                                                 concentration * 2 * laplawfrhs_conc);
      }
    }

    // provide element matrix with linearizations of Soret term in discrete scatra residuals w.r.t.
    // thermo dofs
    mythermo::CalcMatSoretOD(emat, timefacfac, VarManager()->Phinp(0),
        mydiffcond::DiffManager()->GetIsotropicDiff(0), VarManager()->Temp(),
        VarManager()->GradTemp(), my::funct_, my::derxy_);
  }

  return;
}


/*------------------------------------------------------------------------------*
 | set internal variables for element evaluation                     fang 11/15 |
 *------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCondSTIThermo<distype>::SetInternalVariablesForMatAndRHS()
{
  // set internal variables for element evaluation
  VarManager()->SetInternalVariables(my::funct_, my::derxy_, mythermo::etempnp_, my::ephinp_,
      my::ephin_, my::econvelnp_, my::ehist_);

  return;
}


/*----------------------------------------------------------------------*
 | private constructor for singletons                        fang 11/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcElchDiffCondSTIThermo<distype>::ScaTraEleCalcElchDiffCondSTIThermo(
    const int numdofpernode, const int numscal, const std::string& disname)
    :  // constructors of base classes
      mydiffcond::ScaTraEleCalcElchDiffCond(numdofpernode, numscal, disname),
      mythermo::ScaTraEleSTIThermo(numscal)
{
  // safety check
  if (numscal != 1 or numdofpernode != 2)
    dserror("Invalid number of transported scalars or degrees of freedom per node!");

  // replace internal variable manager for isothermal diffusion-conduction formulation by internal
  // variable manager for thermodynamic diffusion-conduction formulation
  my::scatravarmanager_ =
      Teuchos::rcp(new ScaTraEleInternalVariableManagerElchDiffCondSTIThermo<my::nsd_, my::nen_>(
          my::numscal_, myelch::elchparams_, mydiffcond::diffcondparams_));

  return;
}


// template classes
// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCondSTIThermo<DRT::Element::line2>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCondSTIThermo<DRT::Element::line3>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCondSTIThermo<DRT::Element::tri3>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCondSTIThermo<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCondSTIThermo<DRT::Element::quad4>;
// template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCondSTIThermo<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCondSTIThermo<DRT::Element::quad9>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCondSTIThermo<DRT::Element::nurbs9>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCondSTIThermo<DRT::Element::hex8>;
// template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCondSTIThermo<DRT::Element::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCondSTIThermo<DRT::Element::hex27>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCondSTIThermo<DRT::Element::tet4>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCondSTIThermo<DRT::Element::tet10>;
// template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCondSTIThermo<DRT::Element::wedge6>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCondSTIThermo<DRT::Element::pyramid5>;
// template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCondSTIThermo<DRT::Element::nurbs27>;

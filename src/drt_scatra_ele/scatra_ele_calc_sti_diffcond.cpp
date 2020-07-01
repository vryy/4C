/*--------------------------------------------------------------------------*/
/*! \file

\brief evaluate heat transport within binary, concentrated electrolytes on element level

\level 2

*/
/*--------------------------------------------------------------------------*/
#include "scatra_ele_calc_sti_diffcond.H"

#include "scatra_ele_calc_elch_diffcond.H"
#include "scatra_ele_parameter_timint.H"
#include "scatra_ele_sti_thermo.H"
#include "scatra_ele_utils_elch_diffcond.H"

#include "../drt_mat/soret.H"

/*----------------------------------------------------------------------*
 | singleton access method                                   fang 11/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcSTIDiffCond<distype>*
DRT::ELEMENTS::ScaTraEleCalcSTIDiffCond<distype>::Instance(const int numdofpernode,
    const int numscal, const std::string& disname, const ScaTraEleCalcSTIDiffCond* delete_me)
{
  static std::map<std::string, ScaTraEleCalcSTIDiffCond<distype>*> instances;

  if (delete_me == NULL)
  {
    if (instances.find(disname) == instances.end())
      instances[disname] = new ScaTraEleCalcSTIDiffCond<distype>(numdofpernode, numscal, disname);
  }

  else
  {
    for (typename std::map<std::string, ScaTraEleCalcSTIDiffCond<distype>*>::iterator i =
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
void DRT::ELEMENTS::ScaTraEleCalcSTIDiffCond<distype>::Done()
{
  // delete singleton
  Instance(0, 0, "", this);

  return;
}


/*--------------------------------------------------------------------------*
 | calculate element matrix and element right-hand side vector   fang 11/15 |
 *--------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcSTIDiffCond<distype>::Sysmat(
    DRT::Element* ele,                   ///< current element
    Epetra_SerialDenseMatrix& emat,      ///< element matrix
    Epetra_SerialDenseVector& erhs,      ///< element right-hand side vector
    Epetra_SerialDenseVector& subgrdiff  ///< subgrid diffusivity scaling vector
)
{
  // density at t_(n)
  std::vector<double> densn(my::numscal_, 1.0);
  // density at t_(n+1) or t_(n+alpha_F)
  std::vector<double> densnp(my::numscal_, 1.0);
  // density at t_(n+alpha_M)
  std::vector<double> densam(my::numscal_, 1.0);

  // dummy variable
  double dummy(0.);

  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<my::nsd_ele_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // loop over integration points
  for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
  {
    // evaluate shape functions, their derivatives, and domain integration factor at current
    // integration point
    const double fac = my::EvalShapeFuncAndDerivsAtIntPoint(intpoints, iquad);

    // evaluate overall integration factors
    double timefacfac = my::scatraparatimint_->TimeFac() * fac;
    double rhsfac = my::scatraparatimint_->TimeFacRhs() * fac;

    // evaluate internal variables at current integration point
    SetInternalVariablesForMatAndRHS();

    // evaluate material parameters at current integration point
    GetMaterialParams(ele, densn, densnp, densam, dummy, iquad);

    // matrix and vector contributions arising from mass term
    if (not my::scatraparatimint_->IsStationary())
    {
      my::CalcMatMass(emat, 0, fac, densam[0]);
      my::CalcRHSLinMass(erhs, 0, rhsfac, fac, densam[0], densnp[0]);
    }

    // vector contributions arising from history value
    // need to adapt history value to time integration scheme first
    double rhsint(0.0);
    my::ComputeRhsInt(rhsint, densam[0], densnp[0], my::scatravarmanager_->Hist(0));
    my::CalcRHSHistAndSource(erhs, 0, fac, rhsint);

    // matrix and vector contributions arising from diffusion term
    my::CalcMatDiff(emat, 0, timefacfac);
    my::CalcRHSDiff(erhs, 0, rhsfac);

    // matrix and vector contributions arising from source terms
    mystielch::CalcMatAndRhsSource(emat, erhs, timefacfac, rhsfac);
  }  // loop over integration points

  return;
}


/*------------------------------------------------------------------------------------------------*
 | element matrix and right-hand side vector contributions arising from Joule's heat   fang 11/15 |
 *------------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcSTIDiffCond<distype>::CalcMatAndRhsJoule(
    Epetra_SerialDenseMatrix& emat,  //!< element matrix
    Epetra_SerialDenseVector& erhs,  //!< element right-hand side vector
    const double& timefacfac,        //!< domain integration factor times time integration factor
    const double& rhsfac  //!< domain integration factor times time integration factor for
                          //!< right-hand side vector
)
{
  // extract variables and parameters
  const double& concentration = VarManager()->Conc();
  const double& invfval =
      1. / (diffmanagerdiffcond_->GetValence(0) *
               DRT::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->Faraday());
  const double& kappa = diffmanagerdiffcond_->GetCond();
  const double& R = DRT::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->GasConstant();
  const double& t = diffmanagerdiffcond_->GetTransNum(0);

  // current density
  LINALG::Matrix<my::nsd_, 1> i = VarManager()->GradPot();
  i.Update((1 - t) * invfval * 2. * R * my::scatravarmanager_->Phinp(0) / concentration,
      VarManager()->GradConc(), invfval * R * log(concentration), my::scatravarmanager_->GradPhi(0),
      -1.);
  i.Scale(kappa);

  // derivative of current density w.r.t. temperature
  LINALG::Matrix<my::nsd_, 1> di_dT = VarManager()->GradConc();
  di_dT.Scale(kappa * (1 - t) * invfval * 2. * R / concentration);

  // formal, symbolic derivative of current density w.r.t. temperature gradient
  const double di_dgradT = kappa * invfval * R * log(concentration);

  // derivative of square of current density w.r.t. temperature gradient
  LINALG::Matrix<my::nsd_, 1> di2_dgradT = i;
  di2_dgradT.Scale(2. * di_dgradT);

  for (unsigned vi = 0; vi < my::nen_; ++vi)
  {
    for (unsigned ui = 0; ui < my::nen_; ++ui)
    {
      // gradient of shape function times derivative of square of current density w.r.t. temperature
      // gradient
      double di2_dgradT_gradN(0.);
      my::GetLaplacianWeakFormRHS(di2_dgradT_gradN, di2_dgradT, ui);

      // linearizations of Joule's heat term in thermo residuals w.r.t. thermo dofs
      emat(vi, ui) -= timefacfac * my::funct_(vi) / kappa *
                      (di2_dgradT_gradN + 2. * i.Dot(di_dT) * my::funct_(ui));
    }

    // contributions of Joule's heat term to thermo residuals
    erhs[vi] += rhsfac * my::funct_(vi) * i.Dot(i) / kappa;
  }

  return;
}


/*--------------------------------------------------------------------------------------------------*
 | element matrix and right-hand side vector contributions arising from heat of mixing   fang 11/15
 |
 *--------------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcSTIDiffCond<distype>::CalcMatAndRhsMixing(
    Epetra_SerialDenseMatrix& emat,  //!< element matrix
    Epetra_SerialDenseVector& erhs,  //!< element right-hand side vector
    const double& timefacfac,        //!< domain integration factor times time integration factor
    const double& rhsfac  //!< domain integration factor times time integration factor for
                          //!< right-hand side vector
)
{
  // extract variables and parameters
  const double& concentration = VarManager()->Conc();
  const double& diffcoeff = diffmanagerdiffcond_->GetIsotropicDiff(0);
  const LINALG::Matrix<my::nsd_, 1>& gradtemp = my::scatravarmanager_->GradPhi(0);
  const double& soret = DiffManager()->GetSoret();
  const double& temperature = my::scatravarmanager_->Phinp(0);
  const double gasconstant =
      DRT::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->GasConstant();

  // gradient of concentration plus scaled gradient of temperature
  LINALG::Matrix<my::nsd_, 1> a = VarManager()->GradConc();
  a.Update(concentration * soret / temperature, gradtemp, 1.);

  // square of abovementioned gradient
  const double a2 = a.Dot(a);

  for (unsigned vi = 0; vi < my::nen_; ++vi)
  {
    for (unsigned ui = 0; ui < my::nen_; ++ui)
    {
      // abovementioned gradient times gradient of shape function
      double laplawfrhs_a(0.);
      my::GetLaplacianWeakFormRHS(laplawfrhs_a, a, ui);

      // intermediate terms
      const double term1 = 1. / concentration * a2 * my::funct_(ui);
      const double term2 =
          -2. * temperature * a.Dot(gradtemp) * soret * pow(1 / temperature, 2) * my::funct_(ui);
      const double term3 = 2. * temperature * laplawfrhs_a * soret / temperature;

      // linearizations of heat of mixing term in thermo residuals w.r.t. thermo dofs
      emat(vi, ui) -= timefacfac * my::funct_(vi) * pow(diffcoeff, 2) * 2. * gasconstant *
                      (term1 + term2 + term3);
    }

    // contributions of heat of mixing term to thermo residuals
    erhs[vi] += rhsfac * my::funct_(vi) * pow(diffcoeff, 2) * gasconstant * 2. * temperature /
                concentration * a2;
  }

  return;
}


/*------------------------------------------------------------------------------------------------*
 | element matrix and right-hand side vector contributions arising from Soret effect   fang 11/15 |
 *------------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcSTIDiffCond<distype>::CalcMatAndRhsSoret(
    Epetra_SerialDenseMatrix& emat,  //!< element matrix
    Epetra_SerialDenseVector& erhs,  //!< element right-hand side vector
    const double& timefacfac,        //!< domain integration factor times time integration factor
    const double& rhsfac  //!< domain integration factor times time integration factor for
                          //!< right-hand side vector
)
{
  // extract variables and parameters
  const double& concentration = VarManager()->Conc();
  const double& diffcoeff = diffmanagerdiffcond_->GetIsotropicDiff(0);
  const LINALG::Matrix<my::nsd_, 1>& gradtemp = my::scatravarmanager_->GradPhi(0);
  const double& R = DRT::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->GasConstant();
  const double& soret = DiffManager()->GetSoret();
  const double& temperature = my::scatravarmanager_->Phinp(0);

  // gradient of concentration plus scaled gradient of temperature
  LINALG::Matrix<my::nsd_, 1> a = VarManager()->GradConc();
  a.Update(concentration * soret / temperature, gradtemp, 1.);

  for (unsigned vi = 0; vi < my::nen_; ++vi)
  {
    // abovementioned gradient times gradient of test function
    double laplawfrhs_a_vi(0.);
    my::GetLaplacianWeakFormRHS(laplawfrhs_a_vi, a, vi);

    // temperature gradient times gradient of test function
    double laplawfrhs_gradtemp_vi(0.);
    my::GetLaplacianWeakFormRHS(laplawfrhs_gradtemp_vi, gradtemp, vi);

    for (unsigned ui = 0; ui < my::nen_; ++ui)
    {
      // abovementioned gradient times gradient of shape function
      double laplawfrhs_a_ui(0.);
      my::GetLaplacianWeakFormRHS(laplawfrhs_a_ui, a, ui);

      // temperature gradient times gradient of shape function
      double laplawfrhs_gradtemp_ui(0.);
      my::GetLaplacianWeakFormRHS(laplawfrhs_gradtemp_ui, gradtemp, ui);

      // gradient of test function times gradient of shape function
      double laplawf(0.);
      my::GetLaplacianWeakForm(laplawf, ui, vi);

      // intermediate terms
      const double term1 = -gradtemp.Dot(gradtemp) * soret / pow(temperature, 2) * my::funct_(ui);
      const double term2 = laplawfrhs_a_ui / concentration;
      const double term3 = laplawfrhs_gradtemp_ui * soret / temperature;
      const double term4 = my::funct_(ui) * laplawfrhs_a_vi / concentration;
      const double term5 = -soret / temperature * laplawfrhs_gradtemp_vi * my::funct_(ui);
      const double term6 = soret * laplawf;

      // linearizations of Soret effect term in thermo residuals w.r.t. thermo dofs
      emat(vi, ui) += timefacfac * diffcoeff * concentration * 2. * R * soret *
                      ((term1 + term2 + term3) * my::funct_(vi) + term4 + term5 + term6);
    }

    // contributions of Soret effect term to thermo residuals
    erhs[vi] -= rhsfac * diffcoeff * 2. * R * soret *
                (a.Dot(gradtemp) * my::funct_(vi) + temperature * laplawfrhs_a_vi);
  }

  return;
}


/*----------------------------------------------------------------------*
 | evaluate action for off-diagonal system matrix block      fang 11/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ScaTraEleCalcSTIDiffCond<distype>::EvaluateActionOD(
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
    case SCATRA::calc_scatra_mono_odblock_thermoscatra:
    {
      SysmatODThermoScatra(ele, elemat1_epetra);

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
 | fill element matrix with linearizations of discrete thermo residuals w.r.t. scatra dofs   fang
 11/15 |
 *------------------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcSTIDiffCond<distype>::SysmatODThermoScatra(
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

    // evaluate internal variables at current integration point
    SetInternalVariablesForMatAndRHS();

    // evaluate material parameters at current integration point
    std::vector<double> dummy(my::numscal_, 0.);
    double dummy2(0.);
    GetMaterialParams(ele, dummy, dummy, dummy, dummy2, iquad);

    // provide element matrix with linearizations of source terms in discrete thermo residuals
    // w.r.t. scatra dofs
    mystielch::CalcMatSourceOD(emat, my::scatraparatimint_->TimeFac() * fac);
  }

  return;
}


/*------------------------------------------------------------------------------------------------------------------------------*
 | provide element matrix with linearizations of Joule's heat term in discrete thermo residuals
 w.r.t. scatra dofs   fang 11/15 |
 *------------------------------------------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcSTIDiffCond<distype>::CalcMatJouleOD(
    Epetra_SerialDenseMatrix& emat,  //!< element matrix
    const double& timefacfac         //!< domain integration factor times time integration factor
)
{
  // extract variables and parameters
  const double& concentration = VarManager()->Conc();
  const LINALG::Matrix<my::nsd_, 1>& gradconc = VarManager()->GradConc();
  const LINALG::Matrix<my::nsd_, 1>& gradpot = VarManager()->GradPot();
  const LINALG::Matrix<my::nsd_, 1>& gradtemp = my::scatravarmanager_->GradPhi(0);
  const double invfval =
      1. / (diffmanagerdiffcond_->GetValence(0) *
               DRT::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->Faraday());
  const double& kappa = diffmanagerdiffcond_->GetCond();
  const double& kappaderiv = diffmanagerdiffcond_->GetDerivCond(0);
  const double& R = DRT::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->GasConstant();
  const double& t = diffmanagerdiffcond_->GetTransNum(0);
  const double& temperature = my::scatravarmanager_->Phinp(0);

  // current density
  LINALG::Matrix<my::nsd_, 1> i = gradpot;
  i.Update((1 - t) * invfval * 2. * R * temperature / concentration, gradconc,
      invfval * R * log(concentration), gradtemp, -1.);
  i.Scale(kappa);

  // derivative of current density w.r.t. concentration
  LINALG::Matrix<my::nsd_, 1> di_dc = gradpot;
  di_dc.Update(kappaderiv * (1 - t) * invfval * 2. * R * temperature / concentration -
                   kappa * diffmanagerdiffcond_->GetDerivTransNum(0, 0) * invfval * 2. * R *
                       temperature / concentration -
                   kappa * (1 - t) * invfval * 2. * R * temperature / pow(concentration, 2),
      gradconc, kappaderiv * invfval * R * log(concentration) + kappa * invfval * R / concentration,
      gradtemp, -kappaderiv);

  // formal, symbolic derivative of current density w.r.t. concentration gradient
  const double di_dgradc = kappa * (1 - t) * invfval * 2. * R * temperature / concentration;

  // square of current density
  const double i2 = i.Dot(i);

  // derivative of square of current density w.r.t. concentration
  const double di2_dc = 2. * i.Dot(di_dc);

  // derivative of square of current density w.r.t. concentration gradient
  LINALG::Matrix<my::nsd_, 1> di2_dgradc = i;
  di2_dgradc.Scale(2. * di_dgradc);

  // derivative of square of current density w.r.t. gradient of electric potential
  LINALG::Matrix<my::nsd_, 1> di2_dgradpot = i;
  di2_dgradpot.Scale(-2. * kappa);

  for (unsigned vi = 0; vi < my::nen_; ++vi)
  {
    for (unsigned ui = 0; ui < my::nen_; ++ui)
    {
      // gradient of shape function times derivative of square of current density w.r.t.
      // concentration gradient
      double di2_dgradc_gradN(0.);
      my::GetLaplacianWeakFormRHS(di2_dgradc_gradN, di2_dgradc, ui);

      // gradient of shape function times derivative of square of current density w.r.t. gradient of
      // electric potential
      double di2_dgradpot_gradN(0.0);
      my::GetLaplacianWeakFormRHS(di2_dgradpot_gradN, di2_dgradpot, ui);

      // intermediate terms
      const double term1 = my::funct_(ui) * di2_dc;
      const double term2 = di2_dgradc_gradN;
      const double term3 = -my::funct_(ui) * kappaderiv * i2 / kappa;

      // linearizations of Joule's heat term in thermo residuals w.r.t. concentration dofs
      emat(vi, ui * 2) -= timefacfac * my::funct_(vi) * (term1 + term2 + term3) / kappa;

      // linearizations of Joule's heat term in thermo residuals w.r.t. electric potential dofs
      emat(vi, ui * 2 + 1) -= timefacfac * my::funct_(vi) * di2_dgradpot_gradN / kappa;
    }
  }

  return;
}


/*--------------------------------------------------------------------------------------------------------------------------------*
 | provide element matrix with linearizations of heat of mixing term in discrete thermo residuals
 w.r.t. scatra dofs   fang 11/15 |
 *--------------------------------------------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcSTIDiffCond<distype>::CalcMatMixingOD(
    Epetra_SerialDenseMatrix& emat,  //!< element matrix
    const double& timefacfac         //!< domain integration factor times time integration factor
)
{
  // extract variables and parameters
  const double& concentration = VarManager()->Conc();
  const double& diffcoeff = diffmanagerdiffcond_->GetIsotropicDiff(0);
  const LINALG::Matrix<my::nsd_, 1>& gradtemp = my::scatravarmanager_->GradPhi(0);
  const double& soret = DiffManager()->GetSoret();
  const double& temperature = my::scatravarmanager_->Phinp(0);
  const double gasconstant =
      DRT::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->GasConstant();

  // gradient of concentration plus scaled gradient of temperature
  LINALG::Matrix<my::nsd_, 1> a = VarManager()->GradConc();
  a.Update(concentration * soret / temperature, gradtemp, 1.);

  // square of abovementioned gradient
  const double a2 = a.Dot(a);

  for (unsigned vi = 0; vi < my::nen_; ++vi)
  {
    for (unsigned ui = 0; ui < my::nen_; ++ui)
    {
      // abovementioned gradient times gradient of shape function
      double laplawfrhs_a(0.);
      my::GetLaplacianWeakFormRHS(laplawfrhs_a, a, ui);

      // intermediate terms
      const double term1 = 2. * diffcoeff / concentration *
                           diffmanagerdiffcond_->GetDerivIsoDiffCoef(0, 0) * a2 * my::funct_(ui);
      const double term2 = -pow(diffcoeff, 2) / pow(concentration, 2) * a2 * my::funct_(ui);
      const double term3 = 2. * pow(diffcoeff, 2) / concentration * a.Dot(gradtemp) * soret /
                           temperature * my::funct_(ui);
      const double term4 = 2. * pow(diffcoeff, 2) / concentration * laplawfrhs_a;

      // linearizations of heat of mixing term in thermo residuals w.r.t. concentration dofs
      emat(vi, ui * 2) += -timefacfac * my::funct_(vi) * gasconstant * temperature * 2. *
                          (term1 + term2 + term3 + term4);

      // linearizations of heat of mixing term in thermo residuals w.r.t. electric potential dofs
      // are zero
    }
  }

  return;
}


/*------------------------------------------------------------------------------------------------------------------------------*
 | provide element matrix with linearizations of Soret effect term in discrete thermo residuals
 w.r.t. scatra dofs   fang 11/15 |
 *------------------------------------------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcSTIDiffCond<distype>::CalcMatSoretOD(
    Epetra_SerialDenseMatrix& emat,  //!< element matrix
    const double& timefacfac         //!< domain integration factor times time integration factor
)
{
  // extract variables and parameters
  const double& concentration = VarManager()->Conc();
  const double& diffcoeff = diffmanagerdiffcond_->GetIsotropicDiff(0);
  const double& diffcoeffderiv = diffmanagerdiffcond_->GetDerivIsoDiffCoef(0, 0);
  const LINALG::Matrix<my::nsd_, 1>& gradtemp = my::scatravarmanager_->GradPhi(0);
  const double& soret = DiffManager()->GetSoret();
  const double& temperature = my::scatravarmanager_->Phinp(0);
  const double gasconstant =
      DRT::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->GasConstant();

  // square of temperature gradient
  const double gradtemp2 = gradtemp.Dot(gradtemp);

  // gradient of concentration plus scaled gradient of temperature
  LINALG::Matrix<my::nsd_, 1> a = VarManager()->GradConc();
  a.Update(concentration * soret / temperature, gradtemp, 1.);

  // abovementioned gradient times temperature gradient
  const double gradtemp_a = gradtemp.Dot(a);

  for (unsigned vi = 0; vi < my::nen_; ++vi)
  {
    // gradient of test function times abovementioned gradient
    double laplawfrhs_a(0.);
    my::GetLaplacianWeakFormRHS(laplawfrhs_a, a, vi);

    // gradient of test function times temperature gradient
    double laplawfrhs_gradtemp_vi(0.);
    my::GetLaplacianWeakFormRHS(laplawfrhs_gradtemp_vi, gradtemp, vi);

    for (unsigned ui = 0; ui < my::nen_; ++ui)
    {
      // gradient of shape function times temperature gradient
      double laplawfrhs_gradtemp_ui(0.);
      my::GetLaplacianWeakFormRHS(laplawfrhs_gradtemp_ui, gradtemp, ui);

      // gradient of test function times gradient of shape function
      double laplawf(0.);
      my::GetLaplacianWeakForm(laplawf, vi, ui);

      // linearizations of Soret effect term in thermo residuals w.r.t. concentration dofs
      emat(vi, ui * 2) +=
          timefacfac * soret * 2. * gasconstant *
          (my::funct_(ui) *
                  (diffcoeffderiv * (gradtemp_a * my::funct_(vi) + temperature * laplawfrhs_a) +
                      diffcoeff * soret *
                          (gradtemp2 * my::funct_(vi) / temperature + laplawfrhs_gradtemp_vi)) +
              diffcoeff * (laplawfrhs_gradtemp_ui * my::funct_(vi) + temperature * laplawf));

      // linearizations of Soret effect term in thermo residuals w.r.t. electric potential dofs are
      // zero
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 | extract quantities for element evaluation                 fang 11/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcSTIDiffCond<distype>::ExtractElementAndNodeValues(
    DRT::Element* ele,                    //!< current element
    Teuchos::ParameterList& params,       //!< parameter list
    DRT::Discretization& discretization,  //!< discretization
    DRT::Element::LocationArray& la       //!< location array
)
{
  // call base class routine to extract thermo-related quantities
  my::ExtractElementAndNodeValues(ele, params, discretization, la);

  // call base class routine to extract scatra-related quantities
  mystielch::ExtractElementAndNodeValues(ele, params, discretization, la);

  return;
}


/*----------------------------------------------------------------------*
 | get material parameters                                   fang 11/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcSTIDiffCond<distype>::GetMaterialParams(
    const DRT::Element* ele,      //!< current element
    std::vector<double>& densn,   //!< density at t_(n)
    std::vector<double>& densnp,  //!< density at t_(n+1) or t_(n+alpha_F)
    std::vector<double>& densam,  //!< density at t_(n+alpha_M)
    double& visc,                 //!< fluid viscosity
    const int iquad               //!< ID of current integration point
)
{
  // get parameters of primary, thermal material
  Teuchos::RCP<const MAT::Material> material = ele->Material();
  if (material->MaterialType() == INPAR::MAT::m_soret)
    MatSoret(material, densn[0], densnp[0], densam[0]);
  else
    dserror("Invalid thermal material!");

  // get parameters of secondary, scatra material
  material = ele->Material(1);
  if (material->MaterialType() == INPAR::MAT::m_elchmat)
  {
    std::vector<double> concentrations(1, VarManager()->Conc());
    INPAR::ELCH::DiffCondMat dummy(INPAR::ELCH::diffcondmat_undefined);
    utils_->MatElchMat(material, concentrations, INPAR::ELCH::equpot_undefined,
        pow(DRT::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->Faraday(), 2) /
            (DRT::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->GasConstant() *
                VarManager()->Phinp(0)),
        diffmanagerdiffcond_, dummy);
  }
  else
    dserror("Invalid scalar transport material!");

  return;
}  // DRT::ELEMENTS::ScaTraEleCalcSTIDiffCond<distype>::GetMaterialParams


/*----------------------------------------------------------------------*
 | evaluate Soret material                                   fang 11/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcSTIDiffCond<distype>::MatSoret(
    const Teuchos::RCP<const MAT::Material> material,  //!< Soret material
    double& densn,                                     //!< density at time t_(n)
    double& densnp,                                    //!< density at time t_(n+1) or t_(n+alpha_F)
    double& densam                                     //!< density at time t_(n+alpha_M)
)
{
  // extract material parameters from Soret material
  const Teuchos::RCP<const MAT::Soret> matsoret =
      Teuchos::rcp_static_cast<const MAT::Soret>(material);
  densn = densnp = densam = matsoret->Capacity();
  DiffManager()->SetIsotropicDiff(matsoret->Conductivity(), 0);
  DiffManager()->SetSoret(matsoret->SoretCoefficient());

  return;
}  // DRT::ELEMENTS::ScaTraEleCalcSTIDiffCond<distype>::MatSoret


/*------------------------------------------------------------------------------*
 | set internal variables for element evaluation                     fang 11/15 |
 *------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcSTIDiffCond<distype>::SetInternalVariablesForMatAndRHS()
{
  // set internal variables for element evaluation
  VarManager()->SetInternalVariablesSTIElch(my::funct_, my::derxy_, my::ephinp_, my::ephin_,
      mystielch::econcnp_, mystielch::epotnp_, my::econvelnp_, my::ehist_);

  return;
}  // DRT::ELEMENTS::ScaTraEleCalcSTIDiffCond<distype>::SetInternalVariablesForMatAndRHS


/*----------------------------------------------------------------------*
 | private constructor for singletons                        fang 11/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcSTIDiffCond<distype>::ScaTraEleCalcSTIDiffCond(
    const int numdofpernode, const int numscal, const std::string& disname)
    :  // constructors of base classes
      ScaTraEleCalc<distype>::ScaTraEleCalc(numdofpernode, numscal, disname),
      ScaTraEleSTIElch<distype>::ScaTraEleSTIElch(numdofpernode, numscal, disname),

      // diffusion manager for diffusion-conduction formulation
      diffmanagerdiffcond_(Teuchos::rcp(new ScaTraEleDiffManagerElchDiffCond(my::numscal_))),

      // utility class supporting element evaluation for diffusion-conduction formulation
      utils_(DRT::ELEMENTS::ScaTraEleUtilsElchDiffCond<distype>::Instance(
          numdofpernode, numscal, disname))
{
  // safety check
  if (numscal != 1 or numdofpernode != 1)
    dserror("Invalid number of transported scalars or degrees of freedom per node!");

  // replace diffusion manager for standard scalar transport by thermo diffusion manager
  my::diffmanager_ = Teuchos::rcp(new ScaTraEleDiffManagerSTIThermo(my::numscal_));

  // replace internal variable manager for standard scalar transport by internal variable manager
  // for heat transport within electrochemical substances
  my::scatravarmanager_ =
      Teuchos::rcp(new ScaTraEleInternalVariableManagerSTIElch<my::nsd_, my::nen_>(my::numscal_));

  return;
}


// template classes
// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalcSTIDiffCond<DRT::Element::line2>;
template class DRT::ELEMENTS::ScaTraEleCalcSTIDiffCond<DRT::Element::line3>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalcSTIDiffCond<DRT::Element::tri3>;
template class DRT::ELEMENTS::ScaTraEleCalcSTIDiffCond<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleCalcSTIDiffCond<DRT::Element::quad4>;
// template class DRT::ELEMENTS::ScaTraEleCalcSTIDiffCond<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcSTIDiffCond<DRT::Element::quad9>;
template class DRT::ELEMENTS::ScaTraEleCalcSTIDiffCond<DRT::Element::nurbs9>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcSTIDiffCond<DRT::Element::hex8>;
// template class DRT::ELEMENTS::ScaTraEleCalcSTIDiffCond<DRT::Element::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalcSTIDiffCond<DRT::Element::hex27>;
template class DRT::ELEMENTS::ScaTraEleCalcSTIDiffCond<DRT::Element::tet4>;
template class DRT::ELEMENTS::ScaTraEleCalcSTIDiffCond<DRT::Element::tet10>;
// template class DRT::ELEMENTS::ScaTraEleCalcSTIDiffCond<DRT::Element::wedge6>;
template class DRT::ELEMENTS::ScaTraEleCalcSTIDiffCond<DRT::Element::pyramid5>;
// template class DRT::ELEMENTS::ScaTraEleCalcSTIDiffCond<DRT::Element::nurbs27>;

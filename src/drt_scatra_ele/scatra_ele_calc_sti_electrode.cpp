/*--------------------------------------------------------------------------*/
/*!
\file scatra_ele_calc_sti_electrode.cpp

\brief evaluate heat transport within electrodes on element level

\level 2

<pre>
\maintainer Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>
*/
/*--------------------------------------------------------------------------*/
#include "scatra_ele_calc_sti_electrode.H"

#include "scatra_ele_calc_elch_electrode.H"
#include "scatra_ele_parameter_std.H"
#include "scatra_ele_parameter_timint.H"
#include "scatra_ele_sti_thermo.H"

#include "../drt_mat/soret.H"

/*----------------------------------------------------------------------*
 | singleton access method                                   fang 11/15 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcSTIElectrode<distype>* DRT::ELEMENTS::ScaTraEleCalcSTIElectrode<distype>::Instance(
    const int numdofpernode,
    const int numscal,
    const std::string& disname,
    const ScaTraEleCalcSTIElectrode* delete_me
    )
{
  static std::map<std::string,ScaTraEleCalcSTIElectrode<distype>* > instances;

  if(delete_me == NULL)
  {
    if(instances.find(disname) == instances.end())
      instances[disname] = new ScaTraEleCalcSTIElectrode<distype>(numdofpernode,numscal,disname);
  }

  else
  {
    for(typename std::map<std::string,ScaTraEleCalcSTIElectrode<distype>* >::iterator i=instances.begin(); i!=instances.end(); ++i)
    if(i->second == delete_me)
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
void DRT::ELEMENTS::ScaTraEleCalcSTIElectrode<distype>::Done()
{
  // delete singleton
  Instance(0,0,"",this);

  return;
}


/*--------------------------------------------------------------------------*
 | calculate element matrix and element right-hand side vector   fang 11/15 |
 *--------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcSTIElectrode<distype>::Sysmat(
    DRT::Element*               ele,        ///< current element
    Epetra_SerialDenseMatrix&   emat,       ///< element matrix
    Epetra_SerialDenseVector&   erhs,       ///< element right-hand side vector
    Epetra_SerialDenseVector&   subgrdiff   ///< subgrid diffusivity scaling vector
    )
{
  // density at time t_(n+1) or t_(n+alpha_F)
  double densnp(0.);

  // density at time t_(n+alpha_M)
  double densam(0.);

  // dummy variable
  double dummy(0.);

  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<my::nsd_ele_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // loop over integration points
  for(int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
  {
    // evaluate shape functions, their derivatives, and domain integration factor at current integration point
    const double fac = my::EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad);

    // evaluate overall integration factors
    double timefacfac = my::scatraparatimint_->TimeFac()*fac;
    double rhsfac = my::scatraparatimint_->TimeFacRhs()*fac;

    // evaluate internal variables at current integration point
    SetInternalVariablesForMatAndRHS();

    // evaluate material parameters at current integration point
    GetMaterialParams(ele,dummy,densnp,densam,dummy,iquad);

    // matrix and vector contributions arising from mass term
    if(not my::scatraparatimint_->IsStationary())
    {
      my::CalcMatMass(emat,0,fac,densam);
      my::CalcRHSLinMass(erhs,0,rhsfac,fac,densam,densnp);
    }

    // vector contributions arising from history value
    // need to adapt history value to time integration scheme first
    double rhsint(0.0);
    my::ComputeRhsInt(rhsint,densam,densnp,my::scatravarmanager_->Hist(0));
    my::CalcRHSHistAndSource(erhs,0,fac,rhsint);

    // matrix and vector contributions arising from diffusion term
    my::CalcMatDiff(emat,0,timefacfac);
    my::CalcRHSDiff(erhs,0,rhsfac);

    // matrix and vector contributions arising from source terms
    mystielch::CalcMatAndRhsSource(emat,erhs,timefacfac,rhsfac);
  } // loop over integration points

  return;
}


/*------------------------------------------------------------------------------------------------*
 | element matrix and right-hand side vector contributions arising from Joule's heat   fang 11/15 |
 *------------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcSTIElectrode<distype>::CalcMatAndRhsJoule(
    Epetra_SerialDenseMatrix&   emat,         //!< element matrix
    Epetra_SerialDenseVector&   erhs,         //!< element right-hand side vector
    const double&               timefacfac,   //!< domain integration factor times time integration factor
    const double&               rhsfac        //!< domain integration factor times time integration factor for right-hand side vector
    )
{
  // square of gradient of electric potential
  const double gradpot2 = VarManager()->GradPot().Dot(VarManager()->GradPot());

  for(int vi=0; vi<my::nen_; ++vi)
    // linearizations of Joule's heat term in thermo residuals w.r.t. thermo dofs are zero
    // contributions of Joule's heat term to thermo residuals
    erhs[vi] += rhsfac*my::funct_(vi)*gradpot2*diffmanagerstielectrode_->GetCond();

  return;
}


/*--------------------------------------------------------------------------------------------------*
 | element matrix and right-hand side vector contributions arising from heat of mixing   fang 11/15 |
 *--------------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcSTIElectrode<distype>::CalcMatAndRhsMixing(
    Epetra_SerialDenseMatrix&   emat,         //!< element matrix
    Epetra_SerialDenseVector&   erhs,         //!< element right-hand side vector
    const double&               timefacfac,   //!< domain integration factor times time integration factor
    const double&               rhsfac        //!< domain integration factor times time integration factor for right-hand side vector
    )
{
  // extract variables and parameters
  const double& concentration = VarManager()->Conc();
  const double& diffcoeff = diffmanagerstielectrode_->GetIsotropicDiff(0);
  const double& F = INPAR::ELCH::faraday_const;
  const LINALG::Matrix<my::nsd_,1>& gradtemp = my::scatravarmanager_->GradPhi(0);
  const double& soret = DiffManager()->GetSoret();
  const double& temperature = my::scatravarmanager_->Phinp(0);

  // ionic flux density
  LINALG::Matrix<my::nsd_,1> n = VarManager()->GradConc();
  n.Update(-diffcoeff*concentration*soret/temperature,gradtemp,-diffcoeff);

  // derivative of square of ionic flux density w.r.t. temperature
  const double dn2_dT = 2.*n.Dot(gradtemp)*diffcoeff*concentration*soret/pow(temperature,2);

  // formal, symbolic derivative of square of ionic flux density w.r.t. temperature gradient
  LINALG::Matrix<my::nsd_,1> dn2_dgradT = n;
  dn2_dgradT.Scale(-2.*diffcoeff*concentration*soret/temperature);

  for(int vi=0; vi<my::nen_; ++vi)
  {
    for(int ui=0; ui<my::nen_; ++ui)
    {
      // gradient of shape function times derivative of square of ionic flux density w.r.t. temperature gradient
      double dn2_dgradT_ui(0.);
      my::GetLaplacianWeakFormRHS(dn2_dgradT_ui,dn2_dgradT,ui);

      // linearizations of heat of mixing term in thermo residuals w.r.t. thermo dofs
      emat(vi,ui) += timefacfac*my::funct_(vi)*F*diffmanagerstielectrode_->GetOCPDeriv()/diffcoeff*(my::funct_(ui)*dn2_dT+dn2_dgradT_ui);
    }

    // contributions of heat of mixing term to thermo residuals
    erhs[vi] -= rhsfac*diffmanagerstielectrode_->GetOCPDeriv()*F*n.Dot(n)/diffcoeff*my::funct_(vi);
  }

  return;
}


/*------------------------------------------------------------------------------------------------*
 | element matrix and right-hand side vector contributions arising from Soret effect   fang 11/15 |
 *------------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcSTIElectrode<distype>::CalcMatAndRhsSoret(
    Epetra_SerialDenseMatrix&   emat,         //!< element matrix
    Epetra_SerialDenseVector&   erhs,         //!< element right-hand side vector
    const double&               timefacfac,   //!< domain integration factor times time integration factor
    const double&               rhsfac        //!< domain integration factor times time integration factor for right-hand side vector
    )
{
  // extract variables and parameters
  const double& concentration = VarManager()->Conc();
  const double& diffcoeff = diffmanagerstielectrode_->GetIsotropicDiff(0);
  const double& F = INPAR::ELCH::faraday_const;
  const LINALG::Matrix<my::nsd_,1>& gradtemp = my::scatravarmanager_->GradPhi(0);
  const double& soret = DiffManager()->GetSoret();
  const double& temperature = my::scatravarmanager_->Phinp(0);

  // ionic flux density
  LINALG::Matrix<my::nsd_,1> n = VarManager()->GradConc();
  n.Update(-diffcoeff*concentration*soret/temperature,gradtemp,-diffcoeff);

  // derivative of ionic flux density w.r.t. temperature
  LINALG::Matrix<my::nsd_,1> dn_dT = gradtemp;
  dn_dT.Scale(diffcoeff*concentration*soret/pow(temperature,2));

  for(int vi=0; vi<my::nen_; ++vi)
  {
    // gradient of test function times ionic flux density
    double laplawfrhs_n_vi(0.);
    my::GetLaplacianWeakFormRHS(laplawfrhs_n_vi,n,vi);

    // gradient of test function times derivative of ionic flux density w.r.t. temperature
    double laplawfrhs_dndT(0.);
    my::GetLaplacianWeakFormRHS(laplawfrhs_dndT,dn_dT,vi);

    for(int ui=0; ui<my::nen_; ++ui)
    {
      // gradient of shape function times gradient of test function
      double laplawf(0.);
      my::GetLaplacianWeakForm(laplawf,vi,ui);

      // gradient of shape function times gradient of temperature
      double laplawfrhs_gradtemp(0.);
      my::GetLaplacianWeakFormRHS(laplawfrhs_gradtemp,gradtemp,ui);

      // gradient of shape function times gradient of ionic flux density
      double laplawfrhs_n_ui(0.);
      my::GetLaplacianWeakFormRHS(laplawfrhs_n_ui,n,ui);

      // linearizations of Soret effect term in thermo residuals w.r.t. thermo dofs
      emat(vi,ui) += timefacfac*my::funct_(vi)*F*concentration*soret*diffmanagerstielectrode_->GetOCPDeriv()*(laplawfrhs_n_ui/temperature
                     + laplawfrhs_gradtemp/temperature*(-diffcoeff)*concentration*soret/temperature
                     - gradtemp.Dot(n)*pow(1/temperature,2.0)*my::funct_(ui)
                     + gradtemp.Dot(dn_dT)*my::funct_(ui)/temperature)
                     + timefacfac*F*concentration*soret*diffmanagerstielectrode_->GetOCPDeriv()*(-laplawf*diffcoeff*concentration*soret/temperature+laplawfrhs_dndT*my::funct_(ui));
    }

    // contributions of Soret effect term to thermo residuals
    erhs[vi] -= rhsfac*concentration*diffmanagerstielectrode_->GetOCPDeriv()*F*soret*(my::funct_(vi)*n.Dot(gradtemp)/temperature+laplawfrhs_n_vi);
  }

  return;
}


/*----------------------------------------------------------------------*
 | evaluate action for off-diagonal system matrix block      fang 11/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ScaTraEleCalcSTIElectrode<distype>::EvaluateActionOD(
    DRT::Element*                  ele,              //!< current element
    Teuchos::ParameterList&        params,           //!< parameter list
    DRT::Discretization&           discretization,   //!< discretization
    const SCATRA::Action&          action,           //!< action parameter
    DRT::Element::LocationArray&   la,               //!< location array
    Epetra_SerialDenseMatrix&      elemat1_epetra,   //!< element matrix 1
    Epetra_SerialDenseMatrix&      elemat2_epetra,   //!< element matrix 2
    Epetra_SerialDenseVector&      elevec1_epetra,   //!< element right-hand side vector 1
    Epetra_SerialDenseVector&      elevec2_epetra,   //!< element right-hand side vector 2
    Epetra_SerialDenseVector&      elevec3_epetra    //!< element right-hand side vector 3
    )
{
  // determine and evaluate action
  switch(action)
  {
    case SCATRA::calc_scatra_mono_odblock_thermoscatra:
    {
      SysmatODThermoScatra(
          ele,
          elemat1_epetra
          );

      break;
    }

    default:
    {
      // call base class routine
      my::EvaluateActionOD(
          ele,
          params,
          discretization,
          action,
          la,
          elemat1_epetra,
          elemat2_epetra,
          elevec1_epetra,
          elevec2_epetra,
          elevec3_epetra
          );

      break;
    }
  } // switch(action)

  return 0;
}


/*------------------------------------------------------------------------------------------------------*
 | fill element matrix with linearizations of discrete thermo residuals w.r.t. scatra dofs   fang 11/15 |
 *------------------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcSTIElectrode<distype>::SysmatODThermoScatra(
    DRT::Element*               ele,   //!< current element
    Epetra_SerialDenseMatrix&   emat   //!< element matrix
    )
{
  // integration points and weights
  DRT::UTILS::IntPointsAndWeights<my::nsd_ele_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // loop over integration points
  for(int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
  {
    // evaluate shape functions, their derivatives, and domain integration factor at current integration point
    const double fac = my::EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad);

    // evaluate internal variables at current integration point
    SetInternalVariablesForMatAndRHS();

    // evaluate material parameters at current integration point
    double dummy(0.);
    GetMaterialParams(ele,dummy,dummy,dummy,dummy,iquad);

    // provide element matrix with linearizations of source terms in discrete thermo residuals w.r.t. scatra dofs
    mystielch::CalcMatSourceOD(emat,my::scatraparatimint_->TimeFac()*fac);
  }

  return;
}


/*------------------------------------------------------------------------------------------------------------------------------*
 | provide element matrix with linearizations of Joule's heat term in discrete thermo residuals w.r.t. scatra dofs   fang 11/15 |
 *------------------------------------------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcSTIElectrode<distype>::CalcMatJouleOD(
    Epetra_SerialDenseMatrix&   emat,        //!< element matrix
    const double&               timefacfac   //!< domain integration factor times time integration factor
    )
{
  // extract variables and parameters
  const LINALG::Matrix<my::nsd_,1>& gradpot = VarManager()->GradPot();
  const double gradpot2 = gradpot.Dot(gradpot);

  for(int vi=0; vi<my::nen_; ++vi)
  {
    for(int ui=0; ui<my::nen_; ++ui)
    {
      // gradient of shape function times gradient of electric potential
      double laplawfrhs_gradpot(0.0);
      my::GetLaplacianWeakFormRHS(laplawfrhs_gradpot,gradpot,ui);

      // linearizations of Joule's heat term in thermo residuals w.r.t. concentration dofs
      emat(vi,ui*2) -= timefacfac*my::funct_(vi)*diffmanagerstielectrode_->GetDerivCond(0)*gradpot2*my::funct_(ui);

      // linearizations of Joule's heat term in thermo residuals w.r.t. electric potential dofs
      emat(vi,ui*2+1) -= timefacfac*my::funct_(vi)*2.*diffmanagerstielectrode_->GetCond()*laplawfrhs_gradpot;
    }
  }

  return;
}


/*--------------------------------------------------------------------------------------------------------------------------------*
 | provide element matrix with linearizations of heat of mixing term in discrete thermo residuals w.r.t. scatra dofs   fang 11/15 |
 *--------------------------------------------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcSTIElectrode<distype>::CalcMatMixingOD(
    Epetra_SerialDenseMatrix&   emat,        //!< element matrix
    const double&               timefacfac   //!< domain integration factor times time integration factor
    )
{
  // extract variables and parameters
  const double& concentration = VarManager()->Conc();
  const double& diffcoeff = diffmanagerstielectrode_->GetIsotropicDiff(0);
  const double& diffcoeffderiv = diffmanagerstielectrode_->GetDerivIsoDiffCoef(0,0);
  const LINALG::Matrix<my::nsd_,1>& gradconc = VarManager()->GradConc();
  const LINALG::Matrix<my::nsd_,1>& gradtemp = my::scatravarmanager_->GradPhi(0);
  const double& soret = DiffManager()->GetSoret();
  const double& temperature = my::scatravarmanager_->Phinp(0);

  // ionic flux density
  LINALG::Matrix<my::nsd_,1> n = gradconc;
  n.Update(-diffcoeff*concentration*soret/temperature,gradtemp,-diffcoeff);

  // derivative of ionic flux density w.r.t. concentration
  LINALG::Matrix<my::nsd_,1> dn_dc = gradconc;
  dn_dc.Update(-diffcoeffderiv*concentration*soret/temperature-diffcoeff*soret/temperature,gradtemp,-diffcoeffderiv);

  // square of ionic flux density
  const double n2 = n.Dot(n);

  // derivative of square of ionic flux density w.r.t. concentration
  double dn2_dc = 2.*n.Dot(dn_dc);

  // derivative of square of ionic flux density w.r.t. concentration gradient
  LINALG::Matrix<my::nsd_,1> dn2_dgradc = n;
  dn2_dgradc.Scale(-2.*diffcoeff);

  for(int vi=0; vi<my::nen_; ++vi)
  {
    for(int ui=0; ui<my::nen_; ++ui)
    {
      // gradient of shape function times derivative of square of ionic flux density w.r.t. concentration gradient
      double dn2_dgradc_ui(0.);
      my::GetLaplacianWeakFormRHS(dn2_dgradc_ui,dn2_dgradc,ui);

      // intermediate terms
      const double term1 = diffmanagerstielectrode_->GetOCPDeriv2()*n2*my::funct_(ui);
      const double term2 = diffmanagerstielectrode_->GetOCPDeriv()*dn2_dc*my::funct_(ui);
      const double term3 = diffmanagerstielectrode_->GetOCPDeriv()*dn2_dgradc_ui;
      const double term4 = -diffmanagerstielectrode_->GetOCPDeriv()*n2/diffcoeff*diffcoeffderiv*my::funct_(ui);

      // linearizations of heat of mixing term in thermo residuals w.r.t. concentration dofs
      emat(vi,ui*2) += timefacfac*my::funct_(vi)*INPAR::ELCH::faraday_const/diffcoeff*(term1+term2+term3+term4);

      // linearizations of heat of mixing term in thermo residuals w.r.t. electric potential dofs are zero
    }
  }

  return;
}


/*------------------------------------------------------------------------------------------------------------------------------*
 | provide element matrix with linearizations of Soret effect term in discrete thermo residuals w.r.t. scatra dofs   fang 11/15 |
 *------------------------------------------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcSTIElectrode<distype>::CalcMatSoretOD(
    Epetra_SerialDenseMatrix&   emat,        //!< element matrix
    const double&               timefacfac   //!< domain integration factor times time integration factor
    )
{
  // extract variables and parameters
  const double& concentration = VarManager()->Conc();
  const double& diffcoeff = diffmanagerstielectrode_->GetIsotropicDiff(0);
  const double& diffcoeffderiv = diffmanagerstielectrode_->GetDerivIsoDiffCoef(0,0);
  const double& F = INPAR::ELCH::faraday_const;
  const LINALG::Matrix<my::nsd_,1>& gradconc = VarManager()->GradConc();
  const LINALG::Matrix<my::nsd_,1>& gradtemp = my::scatravarmanager_->GradPhi(0);
  const double& soret = DiffManager()->GetSoret();
  const double& temperature = my::scatravarmanager_->Phinp(0);

  // ionic flux density
  LINALG::Matrix<my::nsd_,1> n = gradconc;
  n.Update(-diffcoeff*concentration*soret/temperature,gradtemp,-diffcoeff);

  // derivative of ionic flux density w.r.t. concentration
  LINALG::Matrix<my::nsd_,1> dn_dc = gradconc;
  dn_dc.Update(-diffcoeffderiv*concentration*soret/temperature-diffcoeff*soret/temperature,gradtemp,-diffcoeffderiv);

  for(int vi=0; vi<my::nen_; ++vi)
  {
    // gradient of test function times ionic flux density
    double laplawfrhs_n(0.);
    my::GetLaplacianWeakFormRHS(laplawfrhs_n,n,vi);

    // gradient of test function times derivative of ionic flux density w.r.t. concentration
    double laplawfrhs_dndc(0.);
    my::GetLaplacianWeakFormRHS(laplawfrhs_dndc,dn_dc,vi);

    for(int ui=0; ui<my::nen_; ++ui)
    {
      // gradient of shape function times gradient of test function
      double laplawf(0.);
      my::GetLaplacianWeakForm(laplawf,vi,ui);

      // gradient of shape function times temperature gradient
      double laplawfrhs_gradtemp(0.);
      my::GetLaplacianWeakFormRHS(laplawfrhs_gradtemp,gradtemp,ui);

      // linearizations of Soret effect term in thermo residuals w.r.t. concentration dofs
      emat(vi,ui*2) += timefacfac*my::funct_(vi)*F*soret/temperature*(diffmanagerstielectrode_->GetOCPDeriv()*gradtemp.Dot(n)*my::funct_(ui)
                       + concentration*diffmanagerstielectrode_->GetOCPDeriv2()*gradtemp.Dot(n)*my::funct_(ui)
                       + concentration*diffmanagerstielectrode_->GetOCPDeriv()*gradtemp.Dot(dn_dc)*my::funct_(ui)
                       + concentration*diffmanagerstielectrode_->GetOCPDeriv()*(-diffcoeff)*laplawfrhs_gradtemp)
                       + timefacfac*F*soret*(my::funct_(ui)*diffmanagerstielectrode_->GetOCPDeriv()*laplawfrhs_n
                       + my::funct_(ui)*concentration*diffmanagerstielectrode_->GetOCPDeriv2()*laplawfrhs_n
                       + my::funct_(ui)*concentration*diffmanagerstielectrode_->GetOCPDeriv()*laplawfrhs_dndc
                       + laplawf*concentration*diffmanagerstielectrode_->GetOCPDeriv()*(-diffcoeff));

      // linearizations of Soret effect term in thermo residuals w.r.t. electric potential dofs are zero
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 | extract quantities for element evaluation                 fang 11/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcSTIElectrode<distype>::ExtractElementAndNodeValues(
    DRT::Element*                  ele,              //!< current element
    Teuchos::ParameterList&        params,           //!< parameter list
    DRT::Discretization&           discretization,   //!< discretization
    DRT::Element::LocationArray&   la                //!< location array
    )
{
  // call base class routine to extract thermo-related quantities
  my::ExtractElementAndNodeValues(ele,params,discretization,la);

  // call base class routine to extract scatra-related quantities
  mystielch::ExtractElementAndNodeValues(ele,params,discretization,la);

  return;
}


/*----------------------------------------------------------------------*
 | get material parameters                                   fang 11/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcSTIElectrode<distype>::GetMaterialParams(
    const DRT::Element*   ele,           //!< current element
    double&               densn,         //!< density at time t_(n)
    double&               densnp,        //!< density at time t_(n+1) or t_(n+alpha_F)
    double&               densam,        //!< density at time t_(n+alpha_M)
    double&               visc,          //!< fluid viscosity
    const int             iquad          //!< ID of current integration point
    )
{
  // get parameters of primary, thermal material
  Teuchos::RCP<const MAT::Material> material = ele->Material();
  if(material->MaterialType() == INPAR::MAT::m_soret)
    MatSoret(material,densn,densnp,densam);
  else
    dserror("Invalid thermal material!");

  // get parameters of secondary, scatra material
  material = ele->Material(1);
  if(material->MaterialType() == INPAR::MAT::m_electrode)
  {
    utils_->MatElectrode(material,VarManager()->Conc(),diffmanagerstielectrode_);
    diffmanagerstielectrode_->SetOCPAndDerivs(ele,VarManager()->Conc(),my::scatravarmanager_->Phinp(0));
  }
  else
    dserror("Invalid scalar transport material!");

  return;
} // DRT::ELEMENTS::ScaTraEleCalcSTIElectrode<distype>::GetMaterialParams


/*----------------------------------------------------------------------*
 | evaluate Soret material                                   fang 11/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcSTIElectrode<distype>::MatSoret(
    const Teuchos::RCP<const MAT::Material>   material,   //!< Soret material
    double&                                   densn,      //!< density at time t_(n)
    double&                                   densnp,     //!< density at time t_(n+1) or t_(n+alpha_F)
    double&                                   densam      //!< density at time t_(n+alpha_M)
)
{
  // extract material parameters from Soret material
  const Teuchos::RCP<const MAT::Soret> matsoret = Teuchos::rcp_static_cast<const MAT::Soret>(material);
  densn = densnp = densam = matsoret->Capacity();
  DiffManager()->SetIsotropicDiff(matsoret->Conductivity(),0);
  DiffManager()->SetSoret(matsoret->SoretCoefficient());

  return;
} // DRT::ELEMENTS::ScaTraEleCalcSTIElectrode<distype>::MatSoret


/*------------------------------------------------------------------------------*
 | set internal variables for element evaluation                     fang 11/15 |
 *------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcSTIElectrode<distype>::SetInternalVariablesForMatAndRHS()
{
  // set internal variables for element evaluation
  VarManager()->SetInternalVariablesSTIElch(my::funct_,my::derxy_,my::ephinp_,my::ephin_,mystielch::econcnp_,mystielch::epotnp_,my::econvelnp_,my::ehist_);

  return;
}


/*----------------------------------------------------------------------*
 | private constructor for singletons                        fang 11/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcSTIElectrode<distype>::ScaTraEleCalcSTIElectrode(
    const int numdofpernode,
    const int numscal,
    const std::string& disname
    ) :
    // constructors of base classes
    ScaTraEleCalc<distype>::ScaTraEleCalc(numdofpernode,numscal,disname),
    ScaTraEleSTIElch<distype>::ScaTraEleSTIElch(numdofpernode,numscal,disname),

    // diffusion manager for electrodes
    diffmanagerstielectrode_(Teuchos::rcp(new ScaTraEleDiffManagerSTIElchElectrode(my::numscal_))),

    // utility class supporting element evaluation for electrodes
    utils_(DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<distype>::Instance(numdofpernode,numscal,disname))
{
  // safety check
  if(numscal != 1 or numdofpernode != 1)
    dserror("Invalid number of transported scalars or degrees of freedom per node!");

  // replace diffusion manager for standard scalar transport by thermo diffusion manager
  my::diffmanager_ = Teuchos::rcp(new ScaTraEleDiffManagerSTIThermo(my::numscal_));

  // replace internal variable manager for standard scalar transport by internal variable manager for heat transport within electrochemical substances
  my::scatravarmanager_ = Teuchos::rcp(new ScaTraEleInternalVariableManagerSTIElch<my::nsd_, my::nen_>(my::numscal_));

  return;
}


// template classes
// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalcSTIElectrode<DRT::Element::line2>;
template class DRT::ELEMENTS::ScaTraEleCalcSTIElectrode<DRT::Element::line3>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalcSTIElectrode<DRT::Element::tri3>;
template class DRT::ELEMENTS::ScaTraEleCalcSTIElectrode<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleCalcSTIElectrode<DRT::Element::quad4>;
//template class DRT::ELEMENTS::ScaTraEleCalcSTIElectrode<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcSTIElectrode<DRT::Element::quad9>;
template class DRT::ELEMENTS::ScaTraEleCalcSTIElectrode<DRT::Element::nurbs9>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcSTIElectrode<DRT::Element::hex8>;
//template class DRT::ELEMENTS::ScaTraEleCalcSTIElectrode<DRT::Element::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalcSTIElectrode<DRT::Element::hex27>;
template class DRT::ELEMENTS::ScaTraEleCalcSTIElectrode<DRT::Element::tet4>;
template class DRT::ELEMENTS::ScaTraEleCalcSTIElectrode<DRT::Element::tet10>;
//template class DRT::ELEMENTS::ScaTraEleCalcSTIElectrode<DRT::Element::wedge6>;
template class DRT::ELEMENTS::ScaTraEleCalcSTIElectrode<DRT::Element::pyramid5>;
//template class DRT::ELEMENTS::ScaTraEleCalcSTIElectrode<DRT::Element::nurbs27>;

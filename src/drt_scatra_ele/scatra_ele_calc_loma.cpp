/*--------------------------------------------------------------------------*/
/*!
\file scatra_ele_calc_loma.cpp

\brief Element evaluations for loma problems

\level 2

<pre>
\maintainer Ursula Rasthofer/Volker Gravemeier
            {rasthofer,vgravem}@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15236/45
</pre>
*/
/*--------------------------------------------------------------------------*/

#include "scatra_ele_calc_loma.H"
#include "scatra_ele.H"
#include "scatra_ele_parameter_std.H"
#include "scatra_ele_parameter_timint.H"
#include "scatra_ele_parameter_turbulence.H"

#include "../drt_fem_general/drt_utils_nurbs_shapefunctions.H"

#include "../drt_geometry/position_array.H"

#include "../drt_inpar/inpar_fluid.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"

#include "../drt_mat/mixfrac.H"
#include "../drt_mat/sutherland.H"
#include "../drt_mat/arrhenius_spec.H"
#include "../drt_mat/arrhenius_temp.H"
#include "../drt_mat/arrhenius_pv.H"
#include "../drt_mat/ferech_pv.H"
#include "../drt_mat/yoghurt.H"

#include "../drt_nurbs_discret/drt_nurbs_utils.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcLoma<distype> * DRT::ELEMENTS::ScaTraEleCalcLoma<distype>::Instance(
  const int numdofpernode,
  const int numscal,
  const std::string& disname,
  const ScaTraEleCalcLoma* delete_me )
{
  static std::map<std::string,ScaTraEleCalcLoma<distype>* >  instances;

  if(delete_me == NULL)
  {
    if(instances.find(disname) == instances.end())
      instances[disname] = new ScaTraEleCalcLoma<distype>(numdofpernode,numscal,disname);
  }

  else
  {
    for( typename std::map<std::string,ScaTraEleCalcLoma<distype>* >::iterator i=instances.begin(); i!=instances.end(); ++i )
      if ( i->second == delete_me )
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
void DRT::ELEMENTS::ScaTraEleCalcLoma<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance( 0, 0, "", this );
}


/*----------------------------------------------------------------------*
 | private constructor for singletons                        fang 02/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcLoma<distype>::ScaTraEleCalcLoma(
    const int numdofpernode,
    const int numscal,
    const std::string& disname
    ) :
    DRT::ELEMENTS::ScaTraEleCalc<distype>::ScaTraEleCalc(numdofpernode,numscal,disname),
    ephiam_(my::numscal_),
    densgradfac_(my::numscal_,0.0),
    thermpressnp_(0.0),
    thermpressam_(0.0),
    thermpressdt_(0.0),
    shc_(1.0)
{
  // set appropriate reaction manager
  my::reamanager_ = Teuchos::rcp(new ScaTraEleReaManagerLoma(my::numscal_));

  // safety check
  if(my::turbparams_->MfsConservative())
    dserror("Conservative formulation not supported for loma!");

  return;
}


/*----------------------------------------------------------------------*
 |  evaluate single loma material  (protected)                vg 12/13  |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcLoma<distype>::Materials(
  const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
  const int                               k,        //!< id of current scalar
  double&                                 densn,    //!< density at t_(n)
  double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
  double&                                 densam,   //!< density at t_(n+alpha_M)
  double&                                 visc,          //!< fluid viscosity
  const int                               iquad         //!< id of current gauss point
  )
{
  if (material->MaterialType() == INPAR::MAT::m_mixfrac)
    MatMixFrac(material,k,densn,densnp,densam,visc);
  else if (material->MaterialType() == INPAR::MAT::m_sutherland)
   MatSutherland(material,k,densn,densnp,densam,visc);
  else if (material->MaterialType() == INPAR::MAT::m_arrhenius_pv)
    MatArrheniusPV(material,k,densn,densnp,densam,visc);
  else if (material->MaterialType() == INPAR::MAT::m_arrhenius_spec)
    MatArrheniusSpec(material,k,densn,densnp,densam,visc);
  else if (material->MaterialType() == INPAR::MAT::m_arrhenius_temp)
    MatArrheniusTemp(material,k,densn,densnp,densam,visc);
  else if (material->MaterialType() == INPAR::MAT::m_ferech_pv)
    MatArrheniusPV(material,k,densn,densnp,densam,visc);
  else if (material->MaterialType() == INPAR::MAT::m_yoghurt)
    MatYoghurt(material,k,densn,densnp,densam,visc);
  else dserror("Material type is not supported");

  return;
}


/*----------------------------------------------------------------------*
 | material mixfrac                                            vg 12/13 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcLoma<distype>::MatMixFrac(
  const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
  const int                               k,        //!< id of current scalar
  double&                                 densn,    //!< density at t_(n)
  double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
  double&                                 densam,   //!< density at t_(n+alpha_M)
  double&                                 visc      //!< fluid viscosity
  )
{
  const Teuchos::RCP<const MAT::MixFrac>& actmat
    = Teuchos::rcp_dynamic_cast<const MAT::MixFrac>(material);

  // compute mixture fraction at n+1 or n+alpha_F
  const double mixfracnp = my::scatravarmanager_->Phinp(0);

  // compute dynamic diffusivity at n+1 or n+alpha_F based on mixture fraction
  my::diffmanager_->SetIsotropicDiff(actmat->ComputeDiffusivity(mixfracnp),k);

  // compute density at n+1 or n+alpha_F based on mixture fraction
  densnp = actmat->ComputeDensity(mixfracnp);

  // set specific heat capacity at constant pressure to 1.0
  shc_ = 1.0;

  if (my::scatraparatimint_->IsGenAlpha())
  {
    // compute density at n+alpha_M
    const double mixfracam = my::funct_.Dot(ephiam_[0]);
    densam = actmat->ComputeDensity(mixfracam);

    if (not my::scatraparatimint_->IsIncremental())
    {
      // compute density at n
      const double mixfracn = my::scatravarmanager_->Phin(0);
      densn = actmat->ComputeDensity(mixfracn);
    }
    else densn = 1.0;
  }
  else densam = densnp;

  // factor for density gradient
  densgradfac_[0] = -densnp*densnp*actmat->EosFacA();

  // get also fluid viscosity if subgrid-scale velocity is to be included
  // or multifractal subgrid-scales are used
  if (my::scatrapara_->RBSubGrVel() or my::turbparams_->TurbModel() == INPAR::FLUID::multifractal_subgrid_scales)
    visc = actmat->ComputeViscosity(mixfracnp);

  return;
}


/*----------------------------------------------------------------------*
 | material Sutherland                                         vg 12/13 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcLoma<distype>::MatSutherland(
  const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
  const int                               k,        //!< id of current scalar
  double&                                 densn,    //!< density at t_(n)
  double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
  double&                                 densam,   //!< density at t_(n+alpha_M)
  double&                                 visc      //!< fluid viscosity
  )
{
  const Teuchos::RCP<const MAT::Sutherland>& actmat
    = Teuchos::rcp_dynamic_cast<const MAT::Sutherland>(material);

  // get specific heat capacity at constant pressure
  shc_ = actmat->Shc();

  // compute temperature at n+1 or n+alpha_F and check whether it is positive
  const double tempnp = my::scatravarmanager_->Phinp(0);
  if (tempnp < 0.0)
    dserror("Negative temperature in ScaTra Sutherland material evaluation!");

  // compute diffusivity according to Sutherland's law
  my::diffmanager_->SetIsotropicDiff(actmat->ComputeDiffusivity(tempnp),k);

  // compute density at n+1 or n+alpha_F based on temperature
  // and thermodynamic pressure
  densnp = actmat->ComputeDensity(tempnp,thermpressnp_);

  if (my::scatraparatimint_->IsGenAlpha())
  {
    // compute density at n+alpha_M
    const double tempam = my::funct_.Dot(ephiam_[0]);
    densam = actmat->ComputeDensity(tempam,thermpressam_);

    if (not my::scatraparatimint_->IsIncremental())
    {
      // compute density at n (thermodynamic pressure approximated at n+alpha_M)
      const double tempn = my::scatravarmanager_->Phin(0);
      densn = actmat->ComputeDensity(tempn,thermpressam_);
    }
    else densn = 1.0;
  }
  else densam = densnp;

  // factor for density gradient
  densgradfac_[0] = -densnp/tempnp;

  // get also fluid viscosity if subgrid-scale velocity is to be included
  // or multifractal subgrid-scales are used
  if (my::scatrapara_->RBSubGrVel() or my::turbparams_->TurbModel() == INPAR::FLUID::multifractal_subgrid_scales)
    visc = actmat->ComputeViscosity(tempnp);

  return;
}


/*----------------------------------------------------------------------*
 | material Arrhenius PV                                       vg 12/13 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcLoma<distype>::MatArrheniusPV(
  const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
  const int                               k,        //!< id of current scalar
  double&                                 densn,    //!< density at t_(n)
  double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
  double&                                 densam,   //!< density at t_(n+alpha_M)
  double&                                 visc      //!< fluid viscosity
  )
{
  const Teuchos::RCP<const MAT::ArrheniusPV>& actmat
    = Teuchos::rcp_dynamic_cast<const MAT::ArrheniusPV>(material);

  // get progress variable at n+1 or n+alpha_F
  const double provarnp = my::scatravarmanager_->Phinp(0);

  // get specific heat capacity at constant pressure and
  // compute temperature based on progress variable
  shc_ = actmat->ComputeShc(provarnp);

  // compute temperature at n+1 or n+alpha_F and check whether it is positive
  const double tempnp = actmat->ComputeTemperature(provarnp);
  if (tempnp < 0.0)
    dserror("Negative temperature in ScaTra Arrhenius progress-variable material evaluation!");

  // compute density at n+1 or n+alpha_F
  densnp = actmat->ComputeDensity(provarnp);

  if (my::scatraparatimint_->IsGenAlpha())
  {
    // compute density at n+alpha_M
    const double provaram = my::funct_.Dot(ephiam_[0]);
    densam = actmat->ComputeDensity(provaram);

    if (not my::scatraparatimint_->IsIncremental())
    {
      // compute density at n
      const double provarn = my::scatravarmanager_->Phin(0);
      densn = actmat->ComputeDensity(provarn);
    }
    else densn = 1.0;
  }
  else densam = densnp;

  // factor for density gradient
  densgradfac_[0] = -densnp*actmat->ComputeFactor(provarnp);

  // compute diffusivity according to
  my::diffmanager_->SetIsotropicDiff(actmat->ComputeDiffusivity(tempnp),0);

  // compute reaction coefficient for progress variable
  const double reacoef = actmat->ComputeReactionCoeff(tempnp);

  // set different reaction terms in the reaction manager
  my::reamanager_->SetReaCoeff(reacoef,0);

  // compute right-hand side contribution for progress variable
  // -> equal to reaction coefficient
  Teuchos::rcp_dynamic_cast<ScaTraEleReaManagerLoma>(my::reamanager_)->SetReaTempRhs(reacoef,0);

  // get also fluid viscosity if subgrid-scale velocity is to be included
  // or multifractal subgrid-scales are used
  if (my::scatrapara_->RBSubGrVel() or my::turbparams_->TurbModel() == INPAR::FLUID::multifractal_subgrid_scales)
    visc = actmat->ComputeViscosity(tempnp);

  return;
}


/*----------------------------------------------------------------------*
 | material Arrhenius Spec                                     vg 12/13 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcLoma<distype>::MatArrheniusSpec(
  const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
  const int                               k,        //!< id of current scalar
  double&                                 densn,    //!< density at t_(n)
  double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
  double&                                 densam,   //!< density at t_(n+alpha_M)
  double&                                 visc      //!< fluid viscosity
  )
{
  const Teuchos::RCP<const MAT::ArrheniusSpec>& actmat
    = Teuchos::rcp_dynamic_cast<const MAT::ArrheniusSpec>(material);

  // compute temperature at n+1 or n+alpha_F and check whether it is positive
  const double tempnp = my::scatravarmanager_->Phinp(my::numscal_-1);
  if (tempnp < 0.0)
    dserror("Negative temperature in ScaTra Arrhenius species material evaluation!");

  // compute diffusivity according to Sutherland's law
  my::diffmanager_->SetIsotropicDiff(actmat->ComputeDiffusivity(tempnp),k);

  // compute density at n+1 or n+alpha_F based on temperature
  // and thermodynamic pressure
  densnp = actmat->ComputeDensity(tempnp,thermpressnp_);

  if (my::scatraparatimint_->IsGenAlpha())
  {
    // compute density at n+alpha_M
    const double tempam = my::funct_.Dot(ephiam_[my::numscal_-1]);
    densam = actmat->ComputeDensity(tempam,thermpressam_);

    if (not my::scatraparatimint_->IsIncremental())
    {
      // compute density at n (thermodynamic pressure approximated at n+alpha_M)
      const double tempn = my::scatravarmanager_->Phin(my::numscal_-1);
      densn = actmat->ComputeDensity(tempn,thermpressam_);
    }
    else densn = 1.0;
  }
  else densam = densnp;

  // factor for density gradient
  densgradfac_[k] = -densnp/tempnp;

  // compute reaction coefficient for species equation and set in reaction manager
  const double reacoef = actmat->ComputeReactionCoeff(tempnp);
  my::reamanager_->SetReaCoeff(reacoef,k);

  // get also fluid viscosity if subgrid-scale velocity is to be included
  // or multifractal subgrid-scales are used
  if (my::scatrapara_->RBSubGrVel() or my::turbparams_->TurbModel() == INPAR::FLUID::multifractal_subgrid_scales)
    visc = actmat->ComputeViscosity(tempnp);

  return;
}

/*----------------------------------------------------------------------*
 | material Arrhenius Spec                                     vg 12/13 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcLoma<distype>::MatArrheniusTemp(
  const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
  const int                               k,        //!< id of current scalar
  double&                                 densn,    //!< density at t_(n)
  double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
  double&                                 densam,   //!< density at t_(n+alpha_M)
  double&                                 visc      //!< fluid viscosity
  )
{
  if (k != my::numscal_-1)
    dserror("Temperature always needs to be the last variable for this Arrhenius-type system!");

  const Teuchos::RCP<const MAT::ArrheniusTemp>& actmat
    = Teuchos::rcp_dynamic_cast<const MAT::ArrheniusTemp>(material);

  // get specific heat capacity at constant pressure
  shc_ = actmat->Shc();

  // compute species mass fraction and temperature at n+1 or n+alpha_F, including
  // check whether temperature is positive
  // (only two-equation systems, for the time being, such that only one species
  //  mass fraction possible)
  const double spmfnp = my::scatravarmanager_->Phinp(0);
  const double tempnp = my::scatravarmanager_->Phinp(k);
  if (tempnp < 0.0)
    dserror("Negative temperature in ScaTra Arrhenius temperature material evaluation!");

  // compute diffusivity according to Sutherland's law
  my::diffmanager_->SetIsotropicDiff(actmat->ComputeDiffusivity(tempnp),k);

  // compute density at n+1 or n+alpha_F based on temperature
  // and thermodynamic pressure
  densnp = actmat->ComputeDensity(tempnp,thermpressnp_);

  if (my::scatraparatimint_->IsGenAlpha())
  {
    // compute density at n+alpha_M
    const double tempam = my::funct_.Dot(ephiam_[k]);
    densam = actmat->ComputeDensity(tempam,thermpressam_);

    if (not my::scatraparatimint_->IsIncremental())
    {
      // compute density at n (thermodynamic pressure approximated at n+alpha_M)
      const double tempn = my::scatravarmanager_->Phin(k);
      densn = actmat->ComputeDensity(tempn,thermpressam_);
    }
    else densn = 1.0;
  }
  else densam = densnp;

  // factor for density gradient
  densgradfac_[k] = -densnp/tempnp;

  // compute sum of reaction rates for temperature equation divided by specific
  // heat capacity -> will be considered as a right-hand side contribution
  const double reatemprhs = actmat->ComputeReactionRHS(spmfnp,tempnp);
  Teuchos::rcp_dynamic_cast<ScaTraEleReaManagerLoma>(my::reamanager_)->SetReaTempRhs(reatemprhs,k);

  // get also fluid viscosity if subgrid-scale velocity is to be included
  // or multifractal subgrid-scales are used
  if (my::scatrapara_->RBSubGrVel() or my::turbparams_->TurbModel() == INPAR::FLUID::multifractal_subgrid_scales)
    visc = actmat->ComputeViscosity(tempnp);

  return;
}


/*----------------------------------------------------------------------*
 | material Ferech PV                                          vg 12/13 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcLoma<distype>::MatFerechPV(
  const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
  const int                               k,        //!< id of current scalar
  double&                                 densn,    //!< density at t_(n)
  double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
  double&                                 densam,   //!< density at t_(n+alpha_M)
  double&                                 visc      //!< fluid viscosity
  )
{
  const Teuchos::RCP<const MAT::FerEchPV>& actmat
    = Teuchos::rcp_dynamic_cast<const MAT::FerEchPV>(material);

  // get progress variable at n+1 or n+alpha_F
  const double provarnp = my::scatravarmanager_->Phinp(0);

  // get specific heat capacity at constant pressure
  shc_ = actmat->ComputeShc(provarnp);

  // compute temperature at n+1 or n+alpha_F and check whether it is positive
  const double tempnp = actmat->ComputeTemperature(provarnp);
  if (tempnp < 0.0)
    dserror("Negative temperature in ScaTra Ferziger and Echekki progress-variable material evaluation!");

  // compute density at n+1 or n+alpha_F
  densnp = actmat->ComputeDensity(provarnp);

  if (my::scatraparatimint_->IsGenAlpha())
  {
    // compute density at n+alpha_M
    const double provaram = my::funct_.Dot(ephiam_[0]);
    densam = actmat->ComputeDensity(provaram);

    if (not my::scatraparatimint_->IsIncremental())
    {
      // compute density at n
      const double provarn = my::scatravarmanager_->Phin(0);
      densn = actmat->ComputeDensity(provarn);
    }
    else densn = 1.0;
  }
  else densam = densnp;

  // factor for density gradient
  densgradfac_[0] = -densnp*actmat->ComputeFactor(provarnp);

  // compute diffusivity according to Ferech law
  my::diffmanager_->SetIsotropicDiff(actmat->ComputeDiffusivity(tempnp),0);

  // compute reaction coefficient for progress variable
  const double reacoef = actmat->ComputeReactionCoeff(provarnp);

  // set different reaction terms in the reaction manager
  my::reamanager_->SetReaCoeff(reacoef,0);

  // compute right-hand side contribution for progress variable
  // -> equal to reaction coefficient
  Teuchos::rcp_dynamic_cast<ScaTraEleReaManagerLoma>(my::reamanager_)->SetReaTempRhs(reacoef,0);

   // get also fluid viscosity if subgrid-scale velocity is to be included
   // or multifractal subgrid-scales are used
   if (my::scatrapara_->RBSubGrVel() or my::turbparams_->TurbModel() == INPAR::FLUID::multifractal_subgrid_scales)
     visc = actmat->ComputeViscosity(tempnp);

  return;
}


/*----------------------------------------------------------------------*
 | material Yoghurt                                            vg 12/13 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcLoma<distype>::MatYoghurt(
  const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
  const int                               k,        //!< id of current scalar
  double&                                 densn,    //!< density at t_(n)
  double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
  double&                                 densam,   //!< density at t_(n+alpha_M)
  double&                                 visc      //!< fluid viscosity
  )
{
  const Teuchos::RCP<const MAT::Yoghurt>& actmat
     = Teuchos::rcp_dynamic_cast<const MAT::Yoghurt>(material);

  // get specific heat capacity at constant pressure
  shc_ = actmat->Shc();

  // compute diffusivity
  my::diffmanager_->SetIsotropicDiff(actmat->ComputeDiffusivity(),0);
  //diffus_[0] = actmat->ComputeDiffusivity();

  // get constant density
  densnp = actmat->Density();
  densam = densnp;
  densn = densnp;

  // get also fluid viscosity if subgrid-scale velocity is to be included
  // or multifractal subgrid-scales are used
  if (my::scatrapara_->RBSubGrVel() or my::turbparams_->TurbModel() == INPAR::FLUID::multifractal_subgrid_scales)
  {
    // compute temperature at n+1 or n+alpha_F and check whether it is positive
    const double tempnp = my::scatravarmanager_->Phinp(0);
    if (tempnp < 0.0)
      dserror("Negative temperature in ScaTra yoghurt material evaluation!");

    // compute rate of strain
    double rateofstrain = -1.0e30;
    rateofstrain = this->GetStrainRate(my::evelnp_);

    // compute viscosity for Yoghurt-like flows according to Afonso et al. (2003)
    visc = actmat->ComputeViscosity(rateofstrain,tempnp);
  }

  return;
}


/*-----------------------------------------------------------------------------*
 | compute rhs containing bodyforce                                 ehrl 11/13 |
 *-----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcLoma<distype>::GetRhsInt(
  double&      rhsint, //!< rhs containing bodyforce at Gauss point
  const double densnp, //!< density at t_(n+1)
  const int    k      //!< index of current scalar
  )
{
  // get reatemprhs of species k from the reaction manager
  const double reatemprhs = Teuchos::rcp_dynamic_cast<ScaTraEleReaManagerLoma>(my::reamanager_)->GetReaTempRhs(k);

  // Three cases have to be distinguished for computing the rhs:
  // 1) reactive temperature equation: reaction-rate term
  //    (divided by specific heat capacity)
  // 2) non-reactive temperature equation: heat-source term and
  //    temporal derivative of thermodynamic pressure
  //    (both divided by specific heat capacity)
  // 3) species equation: only potential body force (usually zero)
  const double tol=1e-8;
  if ((reatemprhs < (0.0-tol)) or (reatemprhs > (0.0+tol)))
    rhsint = densnp*reatemprhs/shc_;
  else
  {
    if (k == my::numscal_-1)
    {
      rhsint  = my::bodyforce_[k].Dot(my::funct_)/shc_;
      rhsint += thermpressdt_/shc_;
    }
    else rhsint  = my::bodyforce_[k].Dot(my::funct_);
  }

  return;
} // GetRhsInt


/*------------------------------------------------------------------------------------------*
 |  re-implementatio: calculation of convective element matrix: add conservative contributions  ehrl 11/13 |
 *------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcLoma<distype>::CalcMatConvAddCons(
  Epetra_SerialDenseMatrix&     emat,
  const int                     k,
  const double                  timefacfac,
  const LINALG::Matrix<my::nsd_,1>& convelint,
  const LINALG::Matrix<my::nsd_,1>& gradphi,
  const double                  vdiv,
  const double                  densnp,
  const double                  visc
  )
{
  // convective term using current scalar value
  const double cons_conv_phi = convelint.Dot(gradphi);

  const double consfac = timefacfac*(densnp*vdiv+densgradfac_[k]*cons_conv_phi);
  for (int vi=0; vi<my::nen_; ++vi)
  {
    const double v = consfac*my::funct_(vi);
    const int fvi = vi*my::numdofpernode_+k;

    for (int ui=0; ui<my::nen_; ++ui)
    {
      const int fui = ui*my::numdofpernode_+k;

      emat(fvi,fui) += v*my::funct_(ui);
    }
  }
  return;
}


/*------------------------------------------------------------------- *
 | re-implementatio: adaption of convective term for rhs   ehrl 11/13 |
 *--------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcLoma<distype>::RecomputeConvPhiForRhs(
  double&                       conv_phi,
  const int                     k,
  const LINALG::Matrix<my::nsd_,1>& sgvelint,
  const LINALG::Matrix<my::nsd_,1>& gradphi,
  const double                  densnp,
  const double                  densn,
  const double                  phinp,
  const double                  phin,
  const double                  vdiv
  )
{
  if (my::scatraparatimint_->IsIncremental())
  {
    // addition to convective term due to subgrid-scale velocity
    // (not included in residual)
    double sgconv_phi = sgvelint.Dot(gradphi);
    conv_phi += sgconv_phi;

    // addition to convective term for conservative form
    // (not included in residual)
    if (my::scatrapara_->IsConservative())
    {
      // convective term in conservative form
      conv_phi += phinp*(vdiv+(densgradfac_[k]/densnp)*conv_phi);
    }

    // multiply convective term by density
    conv_phi *= densnp;
  }
  else if (not my::scatraparatimint_->IsIncremental() and my::scatraparatimint_->IsGenAlpha())
  {
    // addition to convective term due to subgrid-scale velocity
    // (not included in residual)
    double sgconv_phi = sgvelint.Dot(gradphi);
    conv_phi += sgconv_phi;

    // addition to convective term for conservative form
    // (not included in residual)
    if (my::scatrapara_->IsConservative())
    {
      // convective term in conservative form
      // caution: velocity divergence is for n+1 and not for n!
      // -> hopefully, this inconsistency is of small amount
      conv_phi += phin*(vdiv+(densgradfac_[k]/densnp)*conv_phi);
    }

    // multiply convective term by density
    conv_phi *= densn;
  }

  return;
}


// template classes

// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalcLoma<DRT::Element::line2>;
template class DRT::ELEMENTS::ScaTraEleCalcLoma<DRT::Element::line3>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalcLoma<DRT::Element::tri3>;
template class DRT::ELEMENTS::ScaTraEleCalcLoma<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleCalcLoma<DRT::Element::quad4>;
//template class DRT::ELEMENTS::ScaTraEleCalcLoma<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcLoma<DRT::Element::quad9>;
template class DRT::ELEMENTS::ScaTraEleCalcLoma<DRT::Element::nurbs9>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcLoma<DRT::Element::hex8>;
//template class DRT::ELEMENTS::ScaTraEleCalcLoma<DRT::Element::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalcLoma<DRT::Element::hex27>;
template class DRT::ELEMENTS::ScaTraEleCalcLoma<DRT::Element::tet4>;
template class DRT::ELEMENTS::ScaTraEleCalcLoma<DRT::Element::tet10>;
//template class DRT::ELEMENTS::ScaTraEleCalcLoma<DRT::Element::wedge6>;
template class DRT::ELEMENTS::ScaTraEleCalcLoma<DRT::Element::pyramid5>;
//template class DRT::ELEMENTS::ScaTraEleCalcLoma<DRT::Element::nurbs27>;

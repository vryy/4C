/*--------------------------------------------------------------------------*/
/*!
\file scatra_ele_calc_loma.H

\brief Element evaluations for loma problems

<pre>
Maintainer: Ursula Rasthofer/Volker Gravemeier
            {rasthofer,vgravem}@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15236/45
</pre>
*/
/*--------------------------------------------------------------------------*/

#include "scatra_ele_calc_loma.H"
#include "scatra_ele.H"
#include "scatra_ele_parameter.H"
#include "scatra_ele_parameter_timint.H"

#include "../drt_geometry/position_array.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_fem_general/drt_utils_nurbs_shapefunctions.H"
#include "../drt_nurbs_discret/drt_nurbs_utils.H"

#include "../drt_mat/mixfrac.H"
#include "../drt_mat/sutherland.H"
#include "../drt_mat/arrhenius_spec.H"
#include "../drt_mat/arrhenius_temp.H"
#include "../drt_mat/arrhenius_pv.H"
#include "../drt_mat/ferech_pv.H"
#include "../drt_mat/yoghurt.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcLoma<distype> * DRT::ELEMENTS::ScaTraEleCalcLoma<distype>::Instance(
  const int numdofpernode,
  const int numscal,
  bool create )
{
  static ScaTraEleCalcLoma<distype> * instance;
  if ( create )
  {
    if ( instance==NULL )
    {
      instance = new ScaTraEleCalcLoma<distype>(numdofpernode,numscal);
    }
  }
  else
  {
    if ( instance!=NULL )
      delete instance;
    instance = NULL;
  }
  return instance;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcLoma<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance( 0, 0, false );
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcLoma<distype>::ScaTraEleCalcLoma(const int numdofpernode,const int numscal)
  : DRT::ELEMENTS::ScaTraEleCalc<distype>::ScaTraEleCalc(numdofpernode,numscal),
    ephiam_(my::numscal_),
    densgradfac_(my::numscal_,0.0),
    thermpressnp_(0.0),
    thermpressam_(0.0),
    thermpressdt_(0.0),
    shc_(1.0)
{
  // set appropriate reaction manager
  my::reamanager_ = Teuchos::rcp(new ScaTraEleReaManagerLoma(my::numscal_));
}


/*----------------------------------------------------------------------*
 | Action type: Evaluate                                rasthofer 12/13 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ScaTraEleCalcLoma<distype>::Evaluate(
  DRT::ELEMENTS::Transport*  ele,
  Teuchos::ParameterList&    params,
  DRT::Discretization&       discretization,
  const std::vector<int>&    lm,
  Epetra_SerialDenseMatrix&  elemat1_epetra,
  Epetra_SerialDenseMatrix&  elemat2_epetra,
  Epetra_SerialDenseVector&  elevec1_epetra,
  Epetra_SerialDenseVector&  elevec2_epetra,
  Epetra_SerialDenseVector&  elevec3_epetra
  )
{
  //--------------------------------------------------------------------------------
  // preparations for element
  //--------------------------------------------------------------------------------

  //get element coordinates
  GEO::fillInitialPositionArray<distype,my::nsd_,LINALG::Matrix<my::nsd_,my::nen_> >(ele,my::xyze_);

  // Now do the nurbs specific stuff (for isogeometric elements)
  if(DRT::NURBS::IsNurbs(distype))
  {
    // access knots and weights for this element
    bool zero_size = DRT::NURBS::GetMyNurbsKnotsAndWeights(discretization,ele,my::myknots_,my::weights_);

    // if we have a zero sized element due to a interpolated point -> exit here
    if(zero_size)
      return(0);
  } // Nurbs specific stuff

  //--------------------------------------------------------------------------------
  // extract element based or nodal values
  //--------------------------------------------------------------------------------

  // extract standard quantities first
  my::ExtractElementAndNodeValues(ele,params,discretization,lm);

  // add further loma-specific values
  if (my::scatraparatimint_->IsGenAlpha())
  {
    // extract additional local values from global vector
    Teuchos::RCP<const Epetra_Vector> phiam = discretization.GetState("phiam");
    if (phiam==Teuchos::null) dserror("Cannot get state vector 'phiam'");
    std::vector<double> myphiam(lm.size());
    DRT::UTILS::ExtractMyValues(*phiam,myphiam,lm);

    // fill element array
    for (int i=0;i<my::nen_;++i)
    {
      for (int k = 0; k< my::numscal_; ++k)
      {
        // split for each transported scalar, insert into element arrays
        ephiam_[k](i,0) = myphiam[k+(i*my::numdofpernode_)];
      }
    } // for i
  }

  // get thermodynamic pressure
  thermpressnp_ = params.get<double>("thermodynamic pressure");
  thermpressdt_ = params.get<double>("time derivative of thermodynamic pressure");
  if (my::scatraparatimint_->IsGenAlpha())
    thermpressam_ = params.get<double>("thermodynamic pressure at n+alpha_M");

  //--------------------------------------------------------------------------------
  // prepare turbulence models
  //--------------------------------------------------------------------------------

  int nlayer = 0;
  my::ExtractTurbulenceApproach(ele,params,discretization,lm,nlayer);

  //--------------------------------------------------------------------------------
  // calculate element coefficient matrix and rhs
  //--------------------------------------------------------------------------------

  // calculate element coefficient matrix and rhs
  my::Sysmat(
    ele,
    elemat1_epetra,
    elevec1_epetra,
    elevec2_epetra);

  // ---------------------------------------------------------------------
  // output values of Prt, diffeff and Cs_delta_sq_Prt (channel flow only)
  // ---------------------------------------------------------------------

  if (my::scatrapara_->TurbModel() == INPAR::FLUID::dynamic_smagorinsky and my::scatrapara_->CsAv())
  {
    Teuchos::ParameterList& turbulencelist = params.sublist("TURBULENCE MODEL");
    my::StoreModelParametersForOutput(ele,ele->Owner() == discretization.Comm().MyPID(),turbulencelist,nlayer);
  }

  return 0;
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
  Teuchos::RCP<ScaTraEleDiffManager>      diffmanager,  //!< diffusion manager handling diffusivity / diffusivities (in case of systems) or (thermal conductivity/specific heat) in case of loma
  Teuchos::RCP<ScaTraEleReaManager>       reamanager,   //!< reaction manager
  double&                                 visc,          //!< fluid viscosity
  const int                               iquad         //!< id of current gauss point
  )
{
  if (material->MaterialType() == INPAR::MAT::m_mixfrac)
    MatMixFrac(material,k,densn,densnp,densam,diffmanager,reamanager,visc);
  else if (material->MaterialType() == INPAR::MAT::m_sutherland)
   MatSutherland(material,k,densn,densnp,densam,diffmanager,reamanager,visc);
  else if (material->MaterialType() == INPAR::MAT::m_arrhenius_pv)
    MatArrheniusPV(material,k,densn,densnp,densam,diffmanager,reamanager,visc);
  else if (material->MaterialType() == INPAR::MAT::m_arrhenius_spec)
    MatArrheniusSpec(material,k,densn,densnp,densam,diffmanager,reamanager,visc);
  else if (material->MaterialType() == INPAR::MAT::m_arrhenius_temp)
    MatArrheniusTemp(material,k,densn,densnp,densam,diffmanager,reamanager,visc);
  else if (material->MaterialType() == INPAR::MAT::m_ferech_pv)
    MatArrheniusPV(material,k,densn,densnp,densam,diffmanager,reamanager,visc);
  else if (material->MaterialType() == INPAR::MAT::m_yoghurt)
    MatYoghurt(material,k,densn,densnp,densam,diffmanager,reamanager,visc);
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
  Teuchos::RCP<ScaTraEleDiffManager>      diffmanager,  //!< diffusion manager handling diffusivity / diffusivities (in case of systems) or (thermal conductivity/specific heat) in case of loma
  Teuchos::RCP<ScaTraEleReaManager>       reamanager,   //!< reaction manager
  double&                                 visc      //!< fluid viscosity
  )
{
  const Teuchos::RCP<const MAT::MixFrac>& actmat
    = Teuchos::rcp_dynamic_cast<const MAT::MixFrac>(material);

  if(my::numdofpernode_!=1)
    dserror("more than 1 dof per node for progress-variable material!");

  // compute mixture fraction at n+1 or n+alpha_F
  const double mixfracnp = my::funct_.Dot(my::ephinp_[0]);

  // compute dynamic diffusivity at n+1 or n+alpha_F based on mixture fraction
  diffmanager->SetIsotropicDiff(actmat->ComputeDiffusivity(mixfracnp),k);

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
      const double mixfracn = my::funct_.Dot(my::ephin_[0]);
      densn = actmat->ComputeDensity(mixfracn);
    }
    else densn = 1.0;
  }
  else densam = densnp;

  // factor for density gradient
  densgradfac_[0] = -densnp*densnp*actmat->EosFacA();

  // get also fluid viscosity if subgrid-scale velocity is to be included
  // or multifractal subgrid-scales are used
  if (my::scatrapara_->RBSubGrVel() or my::scatrapara_->TurbModel() == INPAR::FLUID::multifractal_subgrid_scales)
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
  Teuchos::RCP<ScaTraEleDiffManager>      diffmanager,  //!< diffusion manager handling diffusivity / diffusivities (in case of systems) or (thermal conductivity/specific heat) in case of loma
  Teuchos::RCP<ScaTraEleReaManager>       reamanager,   //!< reaction manager
  double&                                 visc      //!< fluid viscosity
  )
{
  const Teuchos::RCP<const MAT::Sutherland>& actmat
    = Teuchos::rcp_dynamic_cast<const MAT::Sutherland>(material);

  if(my::numdofpernode_!=1)
    dserror("more than 1 dof per node for progress-variable material!");

  // get specific heat capacity at constant pressure
  shc_ = actmat->Shc();

  // compute temperature at n+1 or n+alpha_F
  const double tempnp = my::funct_.Dot(my::ephinp_[0]);
  if (tempnp < 0.0)
    dserror("Negative temperature occurred! Sutherland's law is defined for positive temperatures, only!");

  // compute diffusivity according to material sutherland
  diffmanager->SetIsotropicDiff(actmat->ComputeDiffusivity(tempnp),k);

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
      const double tempn = my::funct_.Dot(my::ephin_[0]);
      densn = actmat->ComputeDensity(tempn,thermpressam_);
    }
    else densn = 1.0;
  }
  else densam = densnp;

  // factor for density gradient
  densgradfac_[0] = -densnp/tempnp;

  // get also fluid viscosity if subgrid-scale velocity is to be included
  // or multifractal subgrid-scales are used
  if (my::scatrapara_->RBSubGrVel() or my::scatrapara_->TurbModel() == INPAR::FLUID::multifractal_subgrid_scales)
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
  Teuchos::RCP<ScaTraEleDiffManager>      diffmanager,  //!< diffusion manager handling diffusivity / diffusivities (in case of systems) or (thermal conductivity/specific heat) in case of loma
  Teuchos::RCP<ScaTraEleReaManager>       reamanager,   //!< reaction manager
  double&                                 visc      //!< fluid viscosity
  )
{
  const Teuchos::RCP<const MAT::ArrheniusPV>& actmat
    = Teuchos::rcp_dynamic_cast<const MAT::ArrheniusPV>(material);

  if(my::numdofpernode_!=1)
    dserror("more than 1 dof per node for progress-variable material!");

  // get progress variable at n+1 or n+alpha_F
  const double provarnp = my::funct_.Dot(my::ephinp_[0]);

  // get specific heat capacity at constant pressure and
  // compute temperature based on progress variable
  shc_ = actmat->ComputeShc(provarnp);

  // compute temperature at n+1 or n+alpha_F
  const double tempnp = actmat->ComputeTemperature(provarnp);

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
      const double provarn = my::funct_.Dot(my::ephin_[0]);
      densn = actmat->ComputeDensity(provarn);
    }
    else densn = 1.0;
  }
  else densam = densnp;

  // factor for density gradient
  densgradfac_[0] = -densnp*actmat->ComputeFactor(provarnp);

  // compute diffusivity according to
  diffmanager->SetIsotropicDiff(actmat->ComputeDiffusivity(tempnp),0);

  // compute reaction coefficient for progress variable
  const double reacoef = actmat->ComputeReactionCoeff(tempnp);

  // set different reaction terms in the reaction manager
  reamanager->SetReaCoeff(reacoef,0);

  // compute right-hand side contribution for progress variable
  // -> equal to reaction coefficient
  Teuchos::rcp_dynamic_cast<ScaTraEleReaManagerLoma>(my::reamanager_)->SetReaTempRhs(reacoef,0);

  // get also fluid viscosity if subgrid-scale velocity is to be included
  // or multifractal subgrid-scales are used
  if (my::scatrapara_->RBSubGrVel() or my::scatrapara_->TurbModel() == INPAR::FLUID::multifractal_subgrid_scales)
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
  Teuchos::RCP<ScaTraEleDiffManager>      diffmanager,  //!< diffusion manager handling diffusivity / diffusivities (in case of systems) or (thermal conductivity/specific heat) in case of loma
  Teuchos::RCP<ScaTraEleReaManager>       reamanager,   //!< reaction manager
  double&                                 visc      //!< fluid viscosity
  )
{
  if (k != my::numscal_-1)
    dserror("Temperature equation always needs to be the last variable for reactive equation system!");

  const Teuchos::RCP<const MAT::ArrheniusSpec>& actmat
    = Teuchos::rcp_dynamic_cast<const MAT::ArrheniusSpec>(material);

  // compute temperature
  const double tempnp = my::funct_.Dot(my::ephinp_[my::numscal_-1]);

  // compute diffusivity according to
  diffmanager->SetIsotropicDiff(actmat->ComputeDiffusivity(tempnp),k);

  // compute reaction coefficient for species equation
  const double reacoef = actmat->ComputeReactionCoeff(tempnp);

  // set different reaction terms in the reaction manager
  reamanager->SetReaCoeff(reacoef,0);

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
  Teuchos::RCP<ScaTraEleDiffManager>      diffmanager,  //!< diffusion manager handling diffusivity / diffusivities (in case of systems) or (thermal conductivity/specific heat) in case of loma
  Teuchos::RCP<ScaTraEleReaManager>       reamanager,   //!< reaction manager
  double&                                 visc      //!< fluid viscosity
  )
{
  if (k != my::numscal_-1)
    dserror("Temperature equation always needs to be the last variable for reactive equation system!");

  const Teuchos::RCP<const MAT::ArrheniusTemp>& actmat
    = Teuchos::rcp_dynamic_cast<const MAT::ArrheniusTemp>(material);

  // get specific heat capacity at constant pressure
  shc_ = actmat->Shc();

  // compute species mass fraction and temperature
  const double spmf   = my::funct_.Dot(my::ephinp_[0]);
  const double tempnp = my::funct_.Dot(my::ephinp_[k]);

  // compute diffusivity according to
  diffmanager->SetIsotropicDiff(actmat->ComputeDiffusivity(tempnp),k);
  //diffus_[k] = actsinglemat->ComputeDiffusivity(tempnp);

  // compute density based on temperature and thermodynamic pressure
  densnp = actmat->ComputeDensity(tempnp,thermpressnp_);

  if (my::scatraparatimint_->IsGenAlpha())
  {
    // compute density at n+alpha_M
    const double tempam = my::funct_.Dot(ephiam_[k]);
    densam = actmat->ComputeDensity(tempam,thermpressam_);

    if (not my::scatraparatimint_->IsIncremental())
    {
      // compute density at n (thermodynamic pressure approximated at n+alpha_M)
      const double tempn = my::funct_.Dot(my::ephin_[k]);
      densn = actmat->ComputeDensity(tempn,thermpressam_);
    }
    else densn = 1.0;
  }
  else densam = densnp;

  // factor for density gradient
  densgradfac_[k] = -densnp/tempnp;

  // compute sum of reaction rates for temperature equation divided by specific
  // heat capacity -> will be considered a right-hand side contribution
  const double reatemprhs = actmat->ComputeReactionRHS(spmf,tempnp)/shc_;

  Teuchos::rcp_dynamic_cast<ScaTraEleReaManagerLoma>(my::reamanager_)->SetReaTempRhs(reatemprhs,0);

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
  Teuchos::RCP<ScaTraEleDiffManager>      diffmanager,  //!< diffusion manager handling diffusivity / diffusivities (in case of systems) or (thermal conductivity/specific heat) in case of loma
  Teuchos::RCP<ScaTraEleReaManager>       reamanager,   //!< reaction manager
  double&                                 visc      //!< fluid viscosity
  )
{
  const Teuchos::RCP<const MAT::FerEchPV>& actmat
    = Teuchos::rcp_dynamic_cast<const MAT::FerEchPV>(material);

  dsassert(my::numdofpernode_==1,"more than 1 dof per node for progress-variable material");

  if(my::numdofpernode_!=1)
    dserror("more than 1 dof per node for progress-variable material!");

  // get progress variable at n+1 or n+alpha_F
  const double provarnp = my::funct_.Dot(my::ephinp_[0]);

  // get specific heat capacity at constant pressure and
  // compute temperature based on progress variable
  shc_ = actmat->ComputeShc(provarnp);
  const double tempnp = actmat->ComputeTemperature(provarnp);

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
      const double provarn = my::funct_.Dot(my::ephin_[0]);
      densn = actmat->ComputeDensity(provarn);
    }
    else densn = 1.0;
  }
  else densam = densnp;

  // factor for density gradient
  densgradfac_[0] = -densnp*actmat->ComputeFactor(provarnp);

  // compute diffusivity according to Ferech law
  diffmanager->SetIsotropicDiff(actmat->ComputeDiffusivity(tempnp),0);

  // compute reaction coefficient for progress variable
  const double reacoef = actmat->ComputeReactionCoeff(provarnp);

  // set different reaction terms in the reaction manager
  reamanager->SetReaCoeff(reacoef,0);

  // compute right-hand side contribution for progress variable
  // -> equal to reaction coefficient
  Teuchos::rcp_dynamic_cast<ScaTraEleReaManagerLoma>(my::reamanager_)->SetReaTempRhs(reacoef,0);

   // get also fluid viscosity if subgrid-scale velocity is to be included
   // or multifractal subgrid-scales are used
   if (my::scatrapara_->RBSubGrVel() or my::scatrapara_->TurbModel() == INPAR::FLUID::multifractal_subgrid_scales)
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
  Teuchos::RCP<ScaTraEleDiffManager>      diffmanager,  //!< diffusion manager handling diffusivity / diffusivities (in case of systems) or (thermal conductivity/specific heat) in case of loma
  Teuchos::RCP<ScaTraEleReaManager>       reamanager,   //!< reaction manager
  double&                                 visc      //!< fluid viscosity
  )
{
  const Teuchos::RCP<const MAT::Yoghurt>& actmat
     = Teuchos::rcp_dynamic_cast<const MAT::Yoghurt>(material);

  if(my::numdofpernode_!=1)
    dserror("more than 1 dof per node for progress-variable material!");

  // get specific heat capacity at constant pressure
  shc_ = actmat->Shc();

  // compute diffusivity
  diffmanager->SetIsotropicDiff(actmat->ComputeDiffusivity(),0);
  //diffus_[0] = actmat->ComputeDiffusivity();

  // get constant density
  densnp = actmat->Density();
  densam = densnp;
  densn = densnp;

  // get also fluid viscosity if subgrid-scale velocity is to be included
  // or multifractal subgrid-scales are used
  if (my::scatrapara_->RBSubGrVel() or my::scatrapara_->TurbModel() == INPAR::FLUID::multifractal_subgrid_scales)
  {
    // compute temperature at n+1 or n+alpha_F
    const double tempnp = my::funct_.Dot(my::ephinp_[0]);

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

  // compute rhs containing bodyforce (divided by specific heat capacity) and,
  // for temperature equation, the time derivative of thermodynamic pressure,
  // if not constant, and for temperature equation of a reactive
  // equation system, the reaction-rate term
  rhsint = my::bodyforce_[k].Dot(my::funct_)/shc_;
  rhsint += thermpressdt_/shc_;
  rhsint += densnp*reatemprhs;

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
//template class DRT::ELEMENTS::ScaTraEleCalcLoma<DRT::Element::tri3>;
//template class DRT::ELEMENTS::ScaTraEleCalcLoma<DRT::Element::tri6>;
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
//template class DRT::ELEMENTS::ScaTraEleCalcLoma<DRT::Element::pyramid5>;
//template class DRT::ELEMENTS::ScaTraEleCalcLoma<DRT::Element::nurbs27>;




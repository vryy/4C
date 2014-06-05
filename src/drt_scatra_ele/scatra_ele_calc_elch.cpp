/*--------------------------------------------------------------------------*/
/*!
\file scatra_ele_calc_elch.cpp

\brief evaluation of ScaTra elements for ion-transport equation

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15252
</pre>
*/
/*--------------------------------------------------------------------------*/

#include "scatra_ele_calc_elch.H"
#include "scatra_ele.H"
#include "scatra_ele_parameter_elch.H"

#include "../drt_geometry/position_array.H"
#include "../drt_fem_general/drt_utils_nurbs_shapefunctions.H"
#include "../drt_nurbs_discret/drt_nurbs_utils.H"
#include "../drt_lib/drt_utils.H"

#include "../drt_mat/ion.H"
#include "../drt_mat/matlist.H"
#include "../drt_inpar/inpar_elch.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcElch<distype>::ScaTraEleCalcElch(const int numdofpernode,const int numscal)
  : DRT::ELEMENTS::ScaTraEleCalc<distype>::ScaTraEleCalc(numdofpernode,numscal),
    epotnp_(my::numscal_)
{
  // initialization of diffusion manager
  my::diffmanager_ = Teuchos::rcp(new ScaTraEleDiffManagerElch(my::numscal_));

  // set appropriate parameter list
  my::scatrapara_ = DRT::ELEMENTS::ScaTraEleParameterElch::Instance();
  elchpara_ = dynamic_cast<DRT::ELEMENTS::ScaTraEleParameterElch*>(my::scatrapara_);

  // initialize internal variable manager
  varmanager_ = Teuchos::rcp(new ScaTraEleInternalVariableManagerElch<my::nsd_, my::nen_>(my::numscal_,my::nsd_,elchpara_));

  if(not my::scatraparatimint_->IsIncremental())
    dserror("Since the ion-transport equations are non-linear, it can be solved only incrementally!!");

  return;
}


/*----------------------------------------------------------------------*
 | Action type: Evaluate                                     ehrl 01/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ScaTraEleCalcElch<distype>::Evaluate(
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

  // get myphi vector for potential and current
  const std::vector<double> myphinp = my::ExtractElementAndNodeValues(ele,params,discretization,lm);

  // get additional values for el. potential at element nodes
  for (int ien=0;ien<my::nen_;++ien)
    epotnp_(ien) = myphinp[ien*my::numdofpernode_+my::numscal_];

  //--------------------------------------------------------------------------------
  // prepare turbulence models
  //--------------------------------------------------------------------------------

  int nlayer = 0;
  my::ExtractTurbulenceApproach(ele,params,discretization,lm,nlayer);

  //--------------------------------------------------------------------------------
  // calculate element coefficient matrix and rhs
  //--------------------------------------------------------------------------------

  Sysmat(
    ele,
    elemat1_epetra,
    elevec1_epetra,
    elevec2_epetra);

//  if(my::eid_==10)
//  {
//  std::cout << "matrix: " << std::endl;
//  std::cout <<  elemat1_epetra << std::endl;
//  std::cout << "vector: " << std::endl;
//  std::cout <<  elevec1_epetra << std::endl;
//  }

  // for certain ELCH problem formulations we have to provide
  // additional flux terms / currents across Dirichlet boundaries
  CorrectionForFluxAcrossDC(discretization,lm,elemat1_epetra,elevec1_epetra);

#if 0
  // for debugging of matrix entries
  if((ele->Id()==2) and (time < 1.3 and time > 1.1))
  {
    FDcheck(
      ele,
      elemat1_epetra,
      elevec1_epetra,
      elevec2_epetra,
      time,
      dt,
      timefac,
      alphaF,
      whichassgd,
      whichfssgd,
      assgd,
      fssgd,
      turbmodel_,
      Cs,
      tpn,
      frt,
      scatratype);
  }
#endif
  return 0;
}

/*----------------------------------------------------------------------*
|  calculate system matrix and rhs (public)                 ehrl  08/08|
*----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElch<distype>::Sysmat(
  DRT::Element*                         ele, ///< the element those matrix is calculated
  Epetra_SerialDenseMatrix&             emat,///< element matrix to calculate
  Epetra_SerialDenseVector&             erhs, ///< element rhs to calculate
  Epetra_SerialDenseVector&             subgrdiff ///< subgrid-diff.-scaling vector
  )
{
  // dynamic cast to elch-specific diffusion manager
  Teuchos::RCP<ScaTraEleDiffManagerElch> dme = Teuchos::rcp_dynamic_cast<ScaTraEleDiffManagerElch>(my::diffmanager_);

  //----------------------------------------------------------------------
  // calculation of element volume both for tau at ele. cent. and int. pt.
  //----------------------------------------------------------------------
  my::EvalShapeFuncAndDerivsAtEleCenter();

  //----------------------------------------------------------------------
  // get material and stabilization parameters (evaluation at element center)
  //----------------------------------------------------------------------
  // density at t_(n)
  double densn(1.0);
  // density at t_(n+1) or t_(n+alpha_F)
  double densnp(1.0);
  // density at t_(n+alpha_M)
  double densam(1.0);

  // fluid viscosity
  double visc(0.0);

  // material parameter at the element center
  if (not my::scatrapara_->MatGP())
    this->GetMaterialParams(ele,densn,densnp,densam,my::diffmanager_,my::reamanager_,visc);

  //----------------------------------------------------------------------
  // integration loop for one element
  //----------------------------------------------------------------------
  // integrations points and weights
  DRT::UTILS::IntPointsAndWeights<my::nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
  {
    const double fac = my::EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad);

    //----------------------------------------------------------------------
    // get material parameters (evaluation at integration point)
    //----------------------------------------------------------------------
    if (my::scatrapara_->MatGP())
      this->GetMaterialParams(ele,densn,densnp,densam,my::diffmanager_,my::reamanager_,visc,iquad);

    // set internal variables
    varmanager_->SetInternalVariablesElch(my::funct_,my::derxy_,my::ephinp_,epotnp_,my::econvelnp_);

    SetFormulationSpecificInternalVariables(dme,varmanager_);

    const double timefacfac = my::scatraparatimint_->TimeFac()*fac;
    const double rhsfac    = my::scatraparatimint_->TimeFacRhs()*fac;

    // loop all scalars
    // deal with a system of transported scalars
    for (int k=0;k<my::numscal_;++k)
    {
      // get history data (or acceleration)
      double hist(0.0);
      hist = my::funct_.Dot(my::ehist_[k]);

      // compute rhs containing bodyforce (divided by specific heat capacity) and,
      // for temperature equation, the time derivative of thermodynamic pressure,
      // if not constant, and for temperature equation of a reactive
      // equation system, the reaction-rate term
      double rhsint(0.0);
      my::GetRhsInt(rhsint,densnp,k);

      // Compute element matrix and rhs
      CalMatAndRhs(varmanager_,emat,erhs,k,fac,timefacfac,rhsfac,dme,rhsint,hist);
    }  // end loop over scalar

    // Compute element matrix and rhs
    CalMatAndRhsOutsideScalarLoop(varmanager_,emat,erhs,fac,timefacfac,rhsfac,dme);
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Material ION                                             ehrl 11/13 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElch<distype>::MatIon(
  const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
  const int                               k,        //!< id of current scalar
  double&                                 densn,    //!< density at t_(n)
  double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
  double&                                 densam,   //!< density at t_(n+alpha_M)
  Teuchos::RCP<ScaTraEleDiffManagerElch>  diffmanager,  //!< diffusion manager handling diffusivity / diffusivities (in case of systems) or (thermal conductivity/specific heat) in case of loma
  Teuchos::RCP<ScaTraEleReaManager>       reamanager,   //!< reaction manager
  double&                                 visc,     //!< fluid viscosity
  const int                               iquad     //!< id of current gauss point
  )
{
  const MAT::Ion* actmat = static_cast<const MAT::Ion*>(material.get());

  // valence of ionic species
  diffmanager->SetValence(actmat->Valence(),k);

  // concentration depending diffusion coefficient
  diffmanager->SetIsotropicDiff(actmat->Diffusivity(),k);

  // Loop over materials is finished - now all material parameter are set
  if(k==(my::numscal_-1))
  {
    // Material data of eliminated ion species is read from the LAST ion material
    // in the matlist!
    if(elchpara_->EquPot()==INPAR::ELCH::equpot_enc_pde_elim)
    {
      diffmanager->IncreaseLengthVector(k, my::numscal_);

      // valence of ionic species
      diffmanager->SetValence(actmat->ElimValence(),my::numscal_);

      // concentration depending diffusion coefficient
      diffmanager->SetIsotropicDiff(actmat->ElimDiffusivity(),my::numscal_);
    }
  }

  return;
}

/*----------------------------------------------------------------------------------*
|  CalcMat: Potential equation ENC                                       ehrl  02/14|
*-----------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElch<distype>::CalcMatPotEquENC(
  Epetra_SerialDenseMatrix&                 emat,     //!< element matrix to be filled
  const int                                 k,        //!< index of current scalar
  const double                              fac,      //!< domain-integration factor
  const double                              alphaf,   //!< time factor for ENC
  Teuchos::RCP<ScaTraEleDiffManagerElch>&   dme       //!< diffusion manager

)
{
  for (int vi=0; vi<my::nen_; ++vi)
  {
    for (int ui=0; ui<my::nen_; ++ui)
    {
      // linearization of the transference number in the conduction term (transport equation)
      //
      // (w, sum(z_k c_k))
      //
      // electroneutrality condition (only derivative w.r.t. concentration c_k)
      emat(vi*my::numdofpernode_+my::numscal_, ui*my::numdofpernode_+k) += alphaf*dme->GetValence(k)*fac*my::funct_(vi)*my::funct_(ui);
    }
  }

  return;
}

/*-------------------------------------------------------------------------------------*
 |  CalcRhs: Potential equation ENC                                         ehrl 11/13 |
 *-------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElch<distype>::CalcRhsPotEquENC(
  Epetra_SerialDenseVector&                 erhs,    //!< element vector to be filled
  const int                                 k,       //!< index of current scalar
  const double                              fac,     //!< time-integration factor for rhs times domain-integration factor
  Teuchos::RCP<ScaTraEleDiffManagerElch>&   dme,     //!< diffusion manager
  const double                              conint   //!< concentration at GP
  )
{
  for (int vi=0; vi<my::nen_; ++vi)
    // electroneutrality condition
    // for incremental formulation, there is the residuum on the rhs! : 0-sum(z_k c_k)
    erhs[vi*my::numdofpernode_+my::numscal_] -= dme->GetValence(k)*fac*my::funct_(vi)*conint;

  return;
}


// template classes

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::line2>;
//template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::tri3>;
//template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::quad4>;
//template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::quad9>;
template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::nurbs9>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::hex8>;
//template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::hex27>;
template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::tet4>;
template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::tet10>;
//template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::wedge6>;
//template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::pyramid5>;
//template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::nurbs27>;




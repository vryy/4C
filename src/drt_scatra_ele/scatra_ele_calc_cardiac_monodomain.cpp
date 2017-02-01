/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_calc_cardiac_monodomain.cpp

\brief scatra_ele_calc_cardiac_monodomain.cpp

\level 2

<pre>
  \maintainer Lasse Jagschies
              jagschies@mhpc.mw.tum.de
              http://www.lnm.mw.tum.de
              089 - 289-10365
</pre>
 *----------------------------------------------------------------------*/


#include "scatra_ele_calc_cardiac_monodomain.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_element.H"

#include "scatra_ele_parameter_timint.H"

#include "../drt_mat/myocard.H"
#include "../drt_mat/matlist.H"

#include "../drt_inpar/inpar_cardiac_monodomain.H"
#include "scatra_ele_parameter_std.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,int probdim>
DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<distype,probdim>::ScaTraEleCalcCardiacMonodomain(const int numdofpernode,const int numscal,const std::string& disname)
  : DRT::ELEMENTS::ScaTraEleCalc<distype,probdim>::ScaTraEleCalc(numdofpernode,numscal,disname),
    DRT::ELEMENTS::ScaTraEleCalcAniso<distype,probdim>::ScaTraEleCalcAniso(numdofpernode,numscal,disname),
    DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype,probdim>::ScaTraEleCalcAdvReac(numdofpernode,numscal,disname)
{

}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype,int probdim>
DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<distype,probdim> * DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<distype,probdim>::Instance(
  const int numdofpernode,
  const int numscal,
  const std::string& disname,
  const ScaTraEleCalcCardiacMonodomain* delete_me )
{
  static std::map<std::string,ScaTraEleCalcCardiacMonodomain<distype,probdim>* >  instances;

  if(delete_me == NULL)
  {
    if(instances.find(disname) == instances.end())
      instances[disname] = new ScaTraEleCalcCardiacMonodomain<distype,probdim>(numdofpernode,numscal,disname);
  }

  else
  {
    for( typename std::map<std::string,ScaTraEleCalcCardiacMonodomain<distype,probdim>* >::iterator i=instances.begin(); i!=instances.end(); ++i )
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
template <DRT::Element::DiscretizationType distype,int probdim>
void DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<distype,probdim>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance( 0, 0, "", this );
}


/*----------------------------------------------------------------------*
 |  evaluate single material  (protected)                    ljag 06/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,int probdim>
void DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<distype,probdim>::Materials(
  const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
  const int                               k,        //!< id of current scalar
  double&                                 densn,    //!< density at t_(n)
  double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
  double&                                 densam,   //!< density at t_(n+alpha_M)
  double&                                 visc,         //!< fluid viscosity
  const int                               iquad         //!< id of current gauss point

  )
{
  // safety check
  if (material->MaterialType() != INPAR::MAT::m_myocard)
    dserror("Material type is not supported");

  // safety check
  Teuchos::RCP<MAT::Myocard> actmat = Teuchos::rcp_dynamic_cast<MAT::Myocard>(Teuchos::rcp_const_cast<MAT::Material>(material));
  if(actmat->GetNumberOfGP() != 1 and not my::scatrapara_->MatGP())
  {
    actmat->SetGP(1);
    actmat->ResizeInternalStateVariables();
  }
  MatMyocard(material,k,densn,densnp,densam,visc,iquad);

  return;
}


/*----------------------------------------------------------------------*
 |  Material ScaTra                                          ljag 06/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,int probdim>
void DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<distype,probdim>::MatMyocard(
  const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
  const int                               k,        //!< id of current scalar
  double&                                 densn,    //!< density at t_(n)
  double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
  double&                                 densam,   //!< density at t_(n+alpha_M)
  double&                                 visc,     //!< fluid viscosity
  const int                               iquad   //!< id of current gauss point (default = -1)
  )
{
  const Teuchos::RCP<const MAT::Myocard>& actmat
    = Teuchos::rcp_dynamic_cast<const MAT::Myocard>(material);

  // dynamic cast to Advanced_Reaction-specific reaction manager
  Teuchos::RCP<ScaTraEleReaManagerAdvReac> advreamanager = Teuchos::rcp_dynamic_cast<ScaTraEleReaManagerAdvReac>(my::reamanager_);

  // dynamic cast to anisotropic diffusion manager
  Teuchos::RCP<ScaTraEleDiffManagerAniso<my::nsd_> > diffmanageraniso = Teuchos::rcp_dynamic_cast<ScaTraEleDiffManagerAniso<my::nsd_> >(my::diffmanager_);

  // get constant diffusivity
  LINALG::Matrix<my::nsd_,my::nsd_> difftensor(true);
  actmat->Diffusivity(difftensor);

  diffmanageraniso->SetAnisotropicDiff(difftensor,k);

    //clear
  advreamanager->Clear(my::numscal_);

  if(my::scatrapara_->SemiImplicit())
  {
    // get membrane potential at n at integration point
    const double phin = my::scatravarmanager_->Phin(k);
    const double phinp = my::scatravarmanager_->Phinp(k);
    // get reaction coefficient
    double react = -actmat->ReaCoeffN(phin, my::scatraparatimint_->Dt(),iquad);
    if (my::scatraparatimint_->IsGenAlpha())
      react *= my::scatraparatimint_->Dt()/my::scatraparatimint_->TimeFac();
    advreamanager->AddToReaBodyForce(react,k);
    advreamanager->AddToReaBodyForceDerivMatrix(0.0,k,k);
    actmat->ReaCoeff(phinp, my::scatraparatimint_->Dt(),iquad);
  }
  else
  {
    // get membrane potential at n+1 or n+alpha_F at integration point
    const double phinp = my::scatravarmanager_->Phinp(k);
    // get reaction coefficient
    advreamanager->AddToReaBodyForce(-actmat->ReaCoeff(phinp, my::scatraparatimint_->Dt(),iquad),k);
    advreamanager->AddToReaBodyForceDerivMatrix(-actmat->ReaCoeffDeriv(phinp, my::scatraparatimint_->Dt(),iquad),k,k);
  }

  return;
} // ScaTraEleCalcCardiacMonodomain<distype>::MatMyocard


/*----------------------------------------------------------------------*
|  calculate system matrix and rhs for ep                 hoermann 06/16|
*----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,int probdim>
void DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<distype,probdim>::Sysmat(
  DRT::Element*                         ele,        ///< the element whose matrix is calculated
  Epetra_SerialDenseMatrix&             emat,       ///< element matrix to calculate
  Epetra_SerialDenseVector&             erhs,       ///< element rhs to calculate
  Epetra_SerialDenseVector&             subgrdiff   ///< subgrid-diff.-scaling vector
  )
{
  // density at t_(n) (one per transported scalar)
  std::vector<double> densn(my::numscal_,1.0);
  // density at t_(n+1) or t_(n+alpha_F) (one per transported scalar)
  std::vector<double> densnp(my::numscal_,1.0);
  // density at t_(n+alpha_M) (one per transported scalar)
  std::vector<double> densam(my::numscal_,1.0);

  // fluid viscosity
  double visc(0.0);

  // calculation of material parameter at element center
  if(not my::scatrapara_->MatGP())
  {
    advreac::EvalShapeFuncAndDerivsAtEleCenter();
    //set Gauss point variables needed for evaluation of mat and rhs
    my::SetInternalVariablesForMatAndRHS();
    advreac::GetMaterialParams(ele,densn,densnp,densam,visc);
  }
  // calculation of material at integration points (different number of integration points possible)
  else
  {
    const DRT::UTILS::IntPointsAndWeights<my::nsd_ele_> intpoints(SCATRA::DisTypeToMatGaussRule<distype>::GetGaussRule(3*ele->Degree()));

    //loop over integration points
    for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
    {
      const double fac = my::EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad);

      //set gauss point variables needed for evaluation of mat and rhs
      my::SetInternalVariablesForMatAndRHS();

      // get material parameters (evaluation at integration point)
      advreac::GetMaterialParams(ele,densn,densnp,densam,visc,iquad);

      // loop all scalars
      for (int k=0;k<my::numscal_;++k) // deal with a system of transported scalars
      {
        double rhsint(0.0);
        advreac::GetRhsInt(rhsint,densnp[k],k);

        LINALG::Matrix<my::nen_,1> dummy(true);
        const double timefacfac = my::scatraparatimint_->TimeFac()*fac;

        // reactive terms on integration point on rhs
        my::ComputeRhsInt(rhsint,densam[k],densnp[k],my::scatravarmanager_->Hist(k));

        // standard Galerkin transient, old part of rhs and bodyforce term
        my::CalcRHSHistAndSource(erhs,k,fac,rhsint);

        // element matrix: reactive term
        advreac::CalcMatReact(emat,k,timefacfac,0.,0.,densnp[k],dummy,dummy);
      }
    }
  }

  //----------------------------------------------------------------------
  // integration loop for one element
  //----------------------------------------------------------------------
  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<my::nsd_ele_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
  {
    const double fac = my::EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad);

    //set gauss point variables needed for evaluation of mat and rhs
    my::SetInternalVariablesForMatAndRHS();

    // loop all scalars
    for (int k=0;k<my::numscal_;++k) // deal with a system of transported scalars
    {

      // compute rhs containing bodyforce
      double rhsint(0.0);
      advreac::GetRhsInt(rhsint,densnp[k],k);

      // integration factors
      const double timefacfac = my::scatraparatimint_->TimeFac()*fac;

      //----------------------------------------------------------------
      // 1) element matrix: stationary terms
      //----------------------------------------------------------------

      // calculation of diffusive element matrix
      aniso::CalcMatDiff(emat,k,timefacfac);

      //----------------------------------------------------------------
      // 2) element matrix: instationary terms
      //----------------------------------------------------------------

      if (not my::scatraparatimint_->IsStationary())
        my::CalcMatMass(emat,k,fac,densam[k]);

      //----------------------------------------------------------------
      // 3) element matrix: reactive term
      //----------------------------------------------------------------

      LINALG::Matrix<my::nen_,1> dummy(true);
      if(not my::scatrapara_->MatGP())
        advreac::CalcMatReact(emat,k,timefacfac,0.,0.,densnp[k],dummy,dummy);

      //----------------------------------------------------------------
      // 5) element right hand side
      //----------------------------------------------------------------
      //----------------------------------------------------------------
      // computation of bodyforce (and potentially history) term,
      // residual, integration factors and standard Galerkin transient
      // term (if required) on right hand side depending on respective
      // (non-)incremental stationary or time-integration scheme
      //----------------------------------------------------------------
      double rhsfac = my::scatraparatimint_->TimeFacRhs() * fac;

      if (my::scatraparatimint_->IsIncremental() and not my::scatraparatimint_->IsStationary())
        my::CalcRHSLinMass(erhs,k,rhsfac,fac,densam[k],densnp[k]);


      if(not my::scatrapara_->MatGP())
      {
        my::ComputeRhsInt(rhsint,densam[k],densnp[k],my::scatravarmanager_->Hist(k));
        // standard Galerkin transient, old part of rhs and bodyforce term
        my::CalcRHSHistAndSource(erhs,k,fac,rhsint);
      }

      //----------------------------------------------------------------
      // standard Galerkin terms on right hand side
      //----------------------------------------------------------------

      // diffusive term
      aniso::CalcRHSDiff(erhs,k,rhsfac);

    }// end loop all scalars
  }// end loop Gauss points

  return;
}


/*----------------------------------------------------------------------*
 | extract element based or nodal values                 hoermann 06/16 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,int probdim>
void DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<distype,probdim>::ExtractElementAndNodeValues(
    DRT::Element*                 ele,
    Teuchos::ParameterList&       params,
    DRT::Discretization&          discretization,
    DRT::Element::LocationArray&  la
)
{
  my::ExtractElementAndNodeValues(ele,params,discretization,la);

  // extract additional local values from global vector
  Teuchos::RCP<const Epetra_Vector> phin = discretization.GetState("phin");
  if (phin==Teuchos::null) dserror("Cannot get state vector 'phin'");
  DRT::UTILS::ExtractMyValues<LINALG::Matrix<my::nen_,1> >(*phin,my::ephin_,la[0].lm_);

}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
// template classes
// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<DRT::Element::line2,1>;
template class DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<DRT::Element::line2,2>;
template class DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<DRT::Element::line2,3>;
template class DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<DRT::Element::line3,1>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<DRT::Element::tri3,2>;
template class DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<DRT::Element::tri3,3>;
template class DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<DRT::Element::tri6,2>;
template class DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<DRT::Element::quad4,2>;
template class DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<DRT::Element::quad4,3>;
//template class DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<DRT::Element::quad9,2>;
template class DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<DRT::Element::nurbs9,2>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<DRT::Element::hex8,3>;
//template class DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<DRT::Element::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<DRT::Element::hex27,3>;
template class DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<DRT::Element::tet4,3>;
template class DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<DRT::Element::tet10,3>;
//template class DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<DRT::Element::wedge6>;
template class DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<DRT::Element::pyramid5,3>;
//template class DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<DRT::Element::nurbs27>;

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
  if (material->MaterialType() == INPAR::MAT::m_myocard)
    MatMyocard(material,k,densn,densnp,densam,visc,iquad);
  else dserror("Material type is not supported");

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

  // get membrane potential at n+1 or n+alpha_F at integration point
  const double phinp = my::scatravarmanager_->Phinp(k);

  //clear
  advreamanager->Clear(my::numscal_);
  // get reaction coefficient
  advreamanager->AddToReaBodyForce(-actmat->ReaCoeff(phinp, my::scatraparatimint_->Dt()),k);
  advreamanager->AddToReaBodyForceDerivMatrix(-actmat->ReaCoeffDeriv(phinp, my::scatraparatimint_->Dt()),k,k);
  advreamanager->SetReaCoeff(0.0, k);
  advreamanager->SetReaCoeffDerivMatrix(0.0,k,k);

  return;
} // ScaTraEleCalcCardiacMonodomain<distype>::MatMyocard



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

/*----------------------------------------------------------------------*/
/*!
 \file scatra_ele_calc_poro_reac.cpp

 \brief

 <pre>
   Maintainer: Moritz Thon
               thon@mhpc.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-10364
 </pre>
 *----------------------------------------------------------------------*/
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_element.H"

#include "../drt_mat/structporo.H"
#include "../drt_mat/structporo_reaction_ecm.H"
#include "../drt_mat/scatra_mat.H"

#include "scatra_ele_parameter.H"
#include "scatra_ele_calc_poro_reac.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcPoroReac<distype>::ScaTraEleCalcPoroReac(const int numdofpernode,const int numscal)
  : DRT::ELEMENTS::ScaTraEleCalc<distype>::ScaTraEleCalc(numdofpernode,numscal),
    DRT::ELEMENTS::ScaTraEleCalcPoro<distype>::ScaTraEleCalcPoro(numdofpernode,numscal),
    DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype>::ScaTraEleCalcAdvReac(numdofpernode,numscal)
{
  // safety check
  if(not my::scatrapara_->TauGP())
    dserror("For poro reactions, tau needs to be evaluated by integration-point evaluations!");

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcPoroReac<distype> * DRT::ELEMENTS::ScaTraEleCalcPoroReac<distype>::Instance(
  const int numdofpernode,
  const int numscal,
  bool create )
{
  static ScaTraEleCalcPoroReac<distype> * instance;
  if ( create )
  {
    if ( instance==NULL )
    {
      instance = new ScaTraEleCalcPoroReac<distype>(numdofpernode,numscal);
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
void DRT::ELEMENTS::ScaTraEleCalcPoroReac<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance( 0, 0, false );
}

/*----------------------------------------------------------------------*
 |  evaluate single material  (protected)                    vuong 10/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcPoroReac<distype>::GetMaterialParams(
  const DRT::Element* ele,       //!< the element we are dealing with
  double&             densn,     //!< density at t_(n)
  double&             densnp,    //!< density at t_(n+1) or t_(n+alpha_F)
  double&             densam,    //!< density at t_(n+alpha_M)
  double&             visc,      //!< fluid viscosity
  const int           iquad      //!< id of current gauss point
  )
{
  //call poro base class to compute porosity
  poro::ComputePorosity(ele);

  //call advreac base class to evaluate porosity
  advreac::GetMaterialParams(ele,densn,densnp,densam,visc,iquad);

  return;
}

/*----------------------------------------------------------------------*
 |  evaluate single material  (protected)                    thon 02/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcPoroReac<distype>::Materials(
  const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
  const int                               k,        //!< id of current scalar
  double&                                 densn,    //!< density at t_(n)
  double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
  double&                                 densam,   //!< density at t_(n+alpha_M)
  double&                                 visc,         //!< fluid viscosity
  const int                               iquad         //!< id of current gauss point
  )
{
  switch(material->MaterialType())
  {
  case INPAR::MAT::m_scatra:
    MatScaTra(material,k,densn,densnp,densam,visc,iquad);
    break;
  case INPAR::MAT::m_scatra_poroECM:
    MatPoroECM(material,k,densn,densnp,densam,visc,iquad);
    break;
  default:
    dserror("Material type %i is not supported",material->MaterialType());
   break;
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Material ScaTra                                          thon 02/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcPoroReac<distype>::MatScaTra(
  const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
  const int                               k,        //!< id of current scalar
  double&                                 densn,    //!< density at t_(n)
  double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
  double&                                 densam,   //!< density at t_(n+alpha_M)
  double&                                 visc,     //!< fluid viscosity
  const int                               iquad     //!< id of current gauss point
  )
{
  poro::MatScaTra(material,k,densn,densnp,densam,visc,iquad);

  return;
} // ScaTraEleCalcPoroReac<distype>::MatScaTra

/*----------------------------------------------------------------------*
 |  Material ScaTra                                          thon 02/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcPoroReac<distype>::MatPoroECM(
  const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
  const int                               k,        //!< id of current scalar
  double&                                 densn,    //!< density at t_(n)
  double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
  double&                                 densam,   //!< density at t_(n+alpha_M)
  double&                                 visc,     //!< fluid viscosity
  const int                               iquad     //!< id of current gauss point
  )
{
  MatScaTra(material,k,densn,densnp,densam,visc,iquad);

  //Todo: clean up
  //access structure discretization
  Teuchos::RCP<DRT::Discretization> structdis = Teuchos::null;
  structdis = DRT::Problem::Instance()->GetDis("structure");
  //get corresponding fluid element (it has the same global ID as the scatra element)
  DRT::Element* structele = structdis->gElement(my::eid_);
  if (structele == NULL)
    dserror("Structure element %i not on local processor", my::eid_);

  const Teuchos::RCP<const MAT::StructPoroReactionECM>& structmat
            = Teuchos::rcp_dynamic_cast<const MAT::StructPoroReactionECM>(structele->Material());
  if(structmat == Teuchos::null)
    dserror("invalid structure material for reactive poroelasticity model for ECM");

  // dynamic cast to Advanced_Reaction-specific reaction manager
  Teuchos::RCP<ScaTraEleReaManagerAdvReac> reamanageradvreac = advreac::ReaManager();

  const double porosity = poro::DiffManager()->GetPorosity();
  const double bodyforce = structmat->BodyForceTerm(porosity);
  const double bodyforce_old = reamanageradvreac->GetReaBodyForce(k);
  reamanageradvreac->AddToReaBodyForce(bodyforce_old+bodyforce,k);

  return;
} // ScaTraEleCalcPoroReac<distype>::MatScaTra

/*----------------------------------------------------------------------*
 | extract element based or nodal values                     vuong 04/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
const std::vector<double> DRT::ELEMENTS::ScaTraEleCalcPoroReac<distype>::ExtractElementAndNodeValues(
    DRT::Element*                 ele,
    Teuchos::ParameterList&       params,
    DRT::Discretization&          discretization,
    DRT::Element::LocationArray&  la
)
{
  return poro::ExtractElementAndNodeValues(ele,params,discretization,la);
}

// template classes

// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::line2>;
template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::line3>;

// 2D elements
//template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::tri3>;
//template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::quad4>;
//template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::quad9>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::hex8>;
//template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::hex27>;
template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::tet4>;
template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::tet10>;
//template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::wedge6>;
template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::pyramid5>;
template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::nurbs9>;
//template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::nurbs27>;


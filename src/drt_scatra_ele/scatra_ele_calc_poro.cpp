/*----------------------------------------------------------------------*/
/*!
 \file scatra_ele_calc_poro.cpp

 \brief

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
 </pre>
 *----------------------------------------------------------------------*/


#include "scatra_ele_calc_poro.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_element.H"

#include "../drt_mat/structporo.H"
#include "../drt_mat/scatra_mat.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcPoro<distype> * DRT::ELEMENTS::ScaTraEleCalcPoro<distype>::Instance(
  const int numdofpernode,
  const int numscal,
  bool create )
{
  static ScaTraEleCalcPoro<distype> * instance;
  if ( create )
  {
    if ( instance==NULL )
    {
      instance = new ScaTraEleCalcPoro<distype>(numdofpernode,numscal);
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
void DRT::ELEMENTS::ScaTraEleCalcPoro<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance( 0, 0, false );
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcPoro<distype>::ScaTraEleCalcPoro(const int numdofpernode,const int numscal)
  : DRT::ELEMENTS::ScaTraEleCalc<distype>::ScaTraEleCalc(numdofpernode,numscal)
{

}

/*----------------------------------------------------------------------*
 |  Material ScaTra                                          ehrl 11/13 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcPoro<distype>::MatScaTra(
  const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
  const int                               k,        //!< id of current scalar
  double&                                 densn,    //!< density at t_(n)
  double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
  double&                                 densam,   //!< density at t_(n+alpha_M)
  Teuchos::RCP<ScaTraEleDiffManager>      diffmanager,  //!< diffusion manager handling diffusivity / diffusivities (in case of systems) or (thermal conductivity/specific heat) in case of loma
  Teuchos::RCP<ScaTraEleReaManager>       reamanager,   //!< reaction manager
  double&                                 visc,     //!< fluid viscosity
  const int                               iquad     //!< id of current gauss point
  )
{
  if(iquad==-1)
    dserror("no gauss point given for evaluation of scatra material. Check your input file.");

  //access structure discretization
  Teuchos::RCP<DRT::Discretization> structdis = Teuchos::null;
  structdis = DRT::Problem::Instance()->GetDis("structure");
  //get corresponding fluid element (it has the same global ID as the scatra element)
  DRT::Element* structele = structdis->gElement(my::eid_);
  if (structele == NULL)
    dserror("Structure element %i not on local processor", my::eid_);

  const Teuchos::RCP<const MAT::StructPoro>& structmat
            = Teuchos::rcp_dynamic_cast<const MAT::StructPoro>(structele->Material());
  if(structmat == Teuchos::null)
    dserror("invalid structure material for poroelasticity");

  const double           porosity   = structmat->GetPorosityAtGP(iquad);

  const Teuchos::RCP<const MAT::ScatraMat>& actmat
    = Teuchos::rcp_dynamic_cast<const MAT::ScatraMat>(material);

  dsassert(my::numdofpernode_==1,"more than 1 dof per node for SCATRA material");

  // set diffusivity (scaled with porosity)
  diffmanager->SetIsotropicDiff(actmat->Diffusivity()*porosity,k);

  // set reaction coefficient (scaled with porosity)
  reamanager->SetReaCoeff(actmat->ReaCoeff()*porosity,k);

  // set densities (scaled with porosity)
  densn = porosity;
  densnp = porosity;
  densam = porosity;


  return;
} // ScaTraEleCalcPoro<distype>::MatScaTra



// template classes

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalcPoro<DRT::Element::line2>;
//template class DRT::ELEMENTS::ScaTraEleCalcPoro<DRT::Element::tri3>;
//template class DRT::ELEMENTS::ScaTraEleCalcPoro<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleCalcPoro<DRT::Element::quad4>;
//template class DRT::ELEMENTS::ScaTraEleCalcPoro<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcPoro<DRT::Element::quad9>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcPoro<DRT::Element::hex8>;
//template class DRT::ELEMENTS::ScaTraEleCalcPoro<DRT::Element::hex20>;
//template class DRT::ELEMENTS::ScaTraEleCalcPoro<DRT::Element::hex27>;
template class DRT::ELEMENTS::ScaTraEleCalcPoro<DRT::Element::tet4>;
template class DRT::ELEMENTS::ScaTraEleCalcPoro<DRT::Element::tet10>;
//template class DRT::ELEMENTS::ScaTraEleCalcPoro<DRT::Element::wedge6>;
//template class DRT::ELEMENTS::ScaTraEleCalcPoro<DRT::Element::pyramid5>;
template class DRT::ELEMENTS::ScaTraEleCalcPoro<DRT::Element::nurbs9>;
//template class DRT::ELEMENTS::ScaTraEleCalcPoro<DRT::Element::nurbs27>;


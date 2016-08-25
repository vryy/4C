/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_calc_poro_reac.cpp

 \brief evaluation class containing routines for calculation of scalar transport
        within porous medium including advanced reactions

\level 2

<pre>
  \maintainer Anh-Tu Vuong
              vuong@lnm.mw.tum.de
              http://www.lnm.mw.tum.de
              089 - 289-15264
</pre>
 *----------------------------------------------------------------------*/
#include "scatra_ele_calc_poro_reac.H"
#include "scatra_ele_parameter_std.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_element.H"

#include "../drt_mat/structporo.H"
#include "../drt_mat/structporo_reaction_ecm.H"
#include "../drt_mat/scatra_mat.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcPoroReac<distype>::ScaTraEleCalcPoroReac(const int numdofpernode,const int numscal,const std::string& disname)
  : DRT::ELEMENTS::ScaTraEleCalc<distype>::ScaTraEleCalc(numdofpernode,numscal,disname),
    DRT::ELEMENTS::ScaTraEleCalcPoro<distype>::ScaTraEleCalcPoro(numdofpernode,numscal,disname),
    DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype>::ScaTraEleCalcAdvReac(numdofpernode,numscal,disname)
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
  const std::string& disname,
  const ScaTraEleCalcPoroReac* delete_me )
{
  static std::map<std::string,ScaTraEleCalcPoroReac<distype>* >  instances;

  if(delete_me == NULL)
  {
    if(instances.find(disname) == instances.end())
      instances[disname] = new ScaTraEleCalcPoroReac<distype>(numdofpernode,numscal,disname);
  }

  else
  {
    for( typename std::map<std::string,ScaTraEleCalcPoroReac<distype>* >::iterator i=instances.begin(); i!=instances.end(); ++i )
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
void DRT::ELEMENTS::ScaTraEleCalcPoroReac<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance( 0, 0, "", this );
}

/*----------------------------------------------------------------------*
 |  evaluate single material  (protected)                    vuong 10/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcPoroReac<distype>::GetMaterialParams(
  const DRT::Element*   ele,       //!< the element we are dealing with
  std::vector<double>&  densn,     //!< density at t_(n)
  std::vector<double>&  densnp,    //!< density at t_(n+1) or t_(n+alpha_F)
  std::vector<double>&  densam,    //!< density at t_(n+alpha_M)
  double&               visc,      //!< fluid viscosity
  const int             iquad      //!< id of current gauss point
  )
{
  //call poro base class to compute porosity
  poro::ComputePorosity(ele);

  //call advreac base class
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
 | extract element based or nodal values                     vuong 04/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcPoroReac<distype>::ExtractElementAndNodeValues(
    DRT::Element*                 ele,
    Teuchos::ParameterList&       params,
    DRT::Discretization&          discretization,
    DRT::Element::LocationArray&  la
)
{
  // call base class routine
  poro::ExtractElementAndNodeValues(ele,params,discretization,la);

  return;
}

// template classes

// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::line2>;
template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::line3>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::tri3>;
template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::tri6>;
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

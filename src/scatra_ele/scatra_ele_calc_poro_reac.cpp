/*----------------------------------------------------------------------*/
/*! \file

 \brief evaluation class containing routines for calculation of scalar transport
        within porous medium including advanced reactions

\level 2

 *----------------------------------------------------------------------*/
#include "scatra_ele_calc_poro_reac.H"
#include "scatra_ele_parameter_std.H"

#include "globalproblem.H"
#include "discret.H"
#include "element.H"

#include "structporo.H"
#include "structporo_reaction_ecm.H"
#include "scatra_mat.H"
#include "singleton_owner.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcPoroReac<distype>::ScaTraEleCalcPoroReac(
    const int numdofpernode, const int numscal, const std::string& disname)
    : DRT::ELEMENTS::ScaTraEleCalc<distype>::ScaTraEleCalc(numdofpernode, numscal, disname),
      DRT::ELEMENTS::ScaTraEleCalcPoro<distype>::ScaTraEleCalcPoro(numdofpernode, numscal, disname),
      DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype>::ScaTraEleCalcAdvReac(
          numdofpernode, numscal, disname)
{
  // safety check
  if (not my::scatrapara_->TauGP())
    dserror("For poro reactions, tau needs to be evaluated by integration-point evaluations!");

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcPoroReac<distype>*
DRT::ELEMENTS::ScaTraEleCalcPoroReac<distype>::Instance(
    const int numdofpernode, const int numscal, const std::string& disname)
{
  static auto singleton_map = ::UTILS::MakeSingletonMap<std::string>(
      [](const int numdofpernode, const int numscal, const std::string& disname)
      {
        return std::unique_ptr<ScaTraEleCalcPoroReac<distype>>(
            new ScaTraEleCalcPoroReac<distype>(numdofpernode, numscal, disname));
      });

  return singleton_map[disname].Instance(
      ::UTILS::SingletonAction::create, numdofpernode, numscal, disname);
}


/*----------------------------------------------------------------------*
 |  evaluate single material  (protected)                    vuong 10/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcPoroReac<distype>::GetMaterialParams(
    const DRT::Element* ele,      //!< the element we are dealing with
    std::vector<double>& densn,   //!< density at t_(n)
    std::vector<double>& densnp,  //!< density at t_(n+1) or t_(n+alpha_F)
    std::vector<double>& densam,  //!< density at t_(n+alpha_M)
    double& visc,                 //!< fluid viscosity
    const int iquad               //!< id of current gauss point
)
{
  // call poro base class to compute porosity
  poro::ComputePorosity(ele);

  // call advreac base class
  advreac::GetMaterialParams(ele, densn, densnp, densam, visc, iquad);

  return;
}

/*----------------------------------------------------------------------*
 |  evaluate single material  (protected)                    thon 02/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcPoroReac<distype>::Materials(
    const Teuchos::RCP<const MAT::Material> material,  //!< pointer to current material
    const int k,                                       //!< id of current scalar
    double& densn,                                     //!< density at t_(n)
    double& densnp,                                    //!< density at t_(n+1) or t_(n+alpha_F)
    double& densam,                                    //!< density at t_(n+alpha_M)
    double& visc,                                      //!< fluid viscosity
    const int iquad                                    //!< id of current gauss point
)
{
  switch (material->MaterialType())
  {
    case INPAR::MAT::m_scatra:
      MatScaTra(material, k, densn, densnp, densam, visc, iquad);
      break;
    default:
      dserror("Material type %i is not supported", material->MaterialType());
      break;
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Material ScaTra                                          thon 02/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcPoroReac<distype>::MatScaTra(
    const Teuchos::RCP<const MAT::Material> material,  //!< pointer to current material
    const int k,                                       //!< id of current scalar
    double& densn,                                     //!< density at t_(n)
    double& densnp,                                    //!< density at t_(n+1) or t_(n+alpha_F)
    double& densam,                                    //!< density at t_(n+alpha_M)
    double& visc,                                      //!< fluid viscosity
    const int iquad                                    //!< id of current gauss point
)
{
  poro::MatScaTra(material, k, densn, densnp, densam, visc, iquad);

  return;
}  // ScaTraEleCalcPoroReac<distype>::MatScaTra

/*----------------------------------------------------------------------*
 | extract element based or nodal values                     vuong 04/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcPoroReac<distype>::ExtractElementAndNodeValues(DRT::Element* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization,
    DRT::Element::LocationArray& la)
{
  // call base class routine
  poro::ExtractElementAndNodeValues(ele, params, discretization, la);

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
// template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::quad9>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::hex8>;
// template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::hex27>;
template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::tet4>;
template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::tet10>;
// template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::wedge6>;
template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::pyramid5>;
template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::nurbs9>;
// template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::nurbs27>;

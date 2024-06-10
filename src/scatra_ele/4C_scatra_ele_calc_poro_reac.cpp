/*----------------------------------------------------------------------*/
/*! \file

 \brief evaluation class containing routines for calculation of scalar transport
        within porous medium including advanced reactions

\level 2

 *----------------------------------------------------------------------*/
#include "4C_scatra_ele_calc_poro_reac.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_fem_general_element.hpp"
#include "4C_global_data.hpp"
#include "4C_mat_scatra.hpp"
#include "4C_mat_structporo.hpp"
#include "4C_mat_structporo_reaction_ecm.hpp"
#include "4C_scatra_ele_parameter_std.hpp"
#include "4C_utils_singleton_owner.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Discret::ELEMENTS::ScaTraEleCalcPoroReac<distype>::ScaTraEleCalcPoroReac(
    const int numdofpernode, const int numscal, const std::string& disname)
    : Discret::ELEMENTS::ScaTraEleCalc<distype>::ScaTraEleCalc(numdofpernode, numscal, disname),
      Discret::ELEMENTS::ScaTraEleCalcPoro<distype>::ScaTraEleCalcPoro(
          numdofpernode, numscal, disname),
      Discret::ELEMENTS::ScaTraEleCalcAdvReac<distype>::ScaTraEleCalcAdvReac(
          numdofpernode, numscal, disname)
{
  // safety check
  if (not my::scatrapara_->TauGP())
    FOUR_C_THROW("For poro reactions, tau needs to be evaluated by integration-point evaluations!");

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Discret::ELEMENTS::ScaTraEleCalcPoroReac<distype>*
Discret::ELEMENTS::ScaTraEleCalcPoroReac<distype>::Instance(
    const int numdofpernode, const int numscal, const std::string& disname)
{
  static auto singleton_map = Core::UTILS::MakeSingletonMap<std::string>(
      [](const int numdofpernode, const int numscal, const std::string& disname)
      {
        return std::unique_ptr<ScaTraEleCalcPoroReac<distype>>(
            new ScaTraEleCalcPoroReac<distype>(numdofpernode, numscal, disname));
      });

  return singleton_map[disname].Instance(
      Core::UTILS::SingletonAction::create, numdofpernode, numscal, disname);
}


/*----------------------------------------------------------------------*
 |  evaluate single material  (protected)                    vuong 10/14 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::ScaTraEleCalcPoroReac<distype>::get_material_params(
    const Core::Elements::Element* ele,  //!< the element we are dealing with
    std::vector<double>& densn,          //!< density at t_(n)
    std::vector<double>& densnp,         //!< density at t_(n+1) or t_(n+alpha_F)
    std::vector<double>& densam,         //!< density at t_(n+alpha_M)
    double& visc,                        //!< fluid viscosity
    const int iquad                      //!< id of current gauss point
)
{
  // call poro base class to compute porosity
  poro::compute_porosity(ele);

  // call advreac base class
  advreac::get_material_params(ele, densn, densnp, densam, visc, iquad);

  return;
}

/*----------------------------------------------------------------------*
 |  evaluate single material  (protected)                    thon 02/14 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::ScaTraEleCalcPoroReac<distype>::materials(
    const Teuchos::RCP<const Core::Mat::Material> material,  //!< pointer to current material
    const int k,                                             //!< id of current scalar
    double& densn,                                           //!< density at t_(n)
    double& densnp,  //!< density at t_(n+1) or t_(n+alpha_F)
    double& densam,  //!< density at t_(n+alpha_M)
    double& visc,    //!< fluid viscosity
    const int iquad  //!< id of current gauss point
)
{
  switch (material->MaterialType())
  {
    case Core::Materials::m_scatra:
      mat_scatra(material, k, densn, densnp, densam, visc, iquad);
      break;
    default:
      FOUR_C_THROW("Material type %i is not supported", material->MaterialType());
      break;
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Material ScaTra                                          thon 02/14 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::ScaTraEleCalcPoroReac<distype>::mat_scatra(
    const Teuchos::RCP<const Core::Mat::Material> material,  //!< pointer to current material
    const int k,                                             //!< id of current scalar
    double& densn,                                           //!< density at t_(n)
    double& densnp,  //!< density at t_(n+1) or t_(n+alpha_F)
    double& densam,  //!< density at t_(n+alpha_M)
    double& visc,    //!< fluid viscosity
    const int iquad  //!< id of current gauss point
)
{
  poro::mat_scatra(material, k, densn, densnp, densam, visc, iquad);

  return;
}  // ScaTraEleCalcPoroReac<distype>::MatScaTra

/*----------------------------------------------------------------------*
 | extract element based or nodal values                     vuong 04/15 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::ScaTraEleCalcPoroReac<distype>::extract_element_and_node_values(
    Core::Elements::Element* ele, Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, Core::Elements::Element::LocationArray& la)
{
  // call base class routine
  poro::extract_element_and_node_values(ele, params, discretization, la);

  return;
}

// template classes

// 1D elements
template class Discret::ELEMENTS::ScaTraEleCalcPoroReac<Core::FE::CellType::line2>;
template class Discret::ELEMENTS::ScaTraEleCalcPoroReac<Core::FE::CellType::line3>;

// 2D elements
template class Discret::ELEMENTS::ScaTraEleCalcPoroReac<Core::FE::CellType::tri3>;
template class Discret::ELEMENTS::ScaTraEleCalcPoroReac<Core::FE::CellType::tri6>;
template class Discret::ELEMENTS::ScaTraEleCalcPoroReac<Core::FE::CellType::quad4>;
// template class Discret::ELEMENTS::ScaTraEleCalcPoroReac<Core::FE::CellType::quad8>;
template class Discret::ELEMENTS::ScaTraEleCalcPoroReac<Core::FE::CellType::quad9>;

// 3D elements
template class Discret::ELEMENTS::ScaTraEleCalcPoroReac<Core::FE::CellType::hex8>;
// template class Discret::ELEMENTS::ScaTraEleCalcPoroReac<Core::FE::CellType::hex20>;
template class Discret::ELEMENTS::ScaTraEleCalcPoroReac<Core::FE::CellType::hex27>;
template class Discret::ELEMENTS::ScaTraEleCalcPoroReac<Core::FE::CellType::tet4>;
template class Discret::ELEMENTS::ScaTraEleCalcPoroReac<Core::FE::CellType::tet10>;
// template class Discret::ELEMENTS::ScaTraEleCalcPoroReac<Core::FE::CellType::wedge6>;
template class Discret::ELEMENTS::ScaTraEleCalcPoroReac<Core::FE::CellType::pyramid5>;
template class Discret::ELEMENTS::ScaTraEleCalcPoroReac<Core::FE::CellType::nurbs9>;
// template class Discret::ELEMENTS::ScaTraEleCalcPoroReac<Core::FE::CellType::nurbs27>;

FOUR_C_NAMESPACE_CLOSE

/*--------------------------------------------------------------------------*/
/*! \file

\brief Factory of artery elements

\level 3

*/
/*--------------------------------------------------------------------------*/
#include "baci_art_net_artery_ele_factory.hpp"

#include "baci_art_net_artery_ele_calc_lin_exp.hpp"
#include "baci_art_net_artery_ele_calc_pres_based.hpp"
#include "baci_art_net_artery_ele_interface.hpp"

BACI_NAMESPACE_OPEN

/*--------------------------------------------------------------------------*
 | (public) kremheller                                                03/18 |
 *--------------------------------------------------------------------------*/
DRT::ELEMENTS::ArteryEleInterface* DRT::ELEMENTS::ArtNetFactory::ProvideImpl(
    CORE::FE::CellType distype, INPAR::ARTDYN::ImplType problem, const std::string& disname)
{
  switch (distype)
  {
    case CORE::FE::CellType::line2:
    {
      return DefineProblemType<CORE::FE::CellType::line2>(problem, disname);

      break;
    }
      // note by J Kremheller:
      // The current implementation relies on the fact that we only use linear elements on several
      // occasions, for instance, when calculating element volumetric flow and element length
      // I currently do not see any application of higher order elements, since the formulation
      // essentially depends on a linear pressure drop prescribed in each element (Hagen-Poiseuille
      // equation)
      // but if this is ever desired the implementation should be checked carefully
    default:
      dserror("Only line2 elements available so far");
      break;
  }
  return nullptr;
}


/*--------------------------------------------------------------------------*
 | (public) kremheller                                                03/18 |
 *--------------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
DRT::ELEMENTS::ArteryEleInterface* DRT::ELEMENTS::ArtNetFactory::DefineProblemType(
    INPAR::ARTDYN::ImplType problem, const std::string& disname)
{
  switch (problem)
  {
    case INPAR::ARTDYN::ImplType::impltype_lin_exp:
    {
      // 2 dofs per node
      return DRT::ELEMENTS::ArteryEleCalcLinExp<distype>::Instance(2, disname);
      break;
    }
    case INPAR::ARTDYN::ImplType::impltype_pressure_based:
    {
      // 1 dof per node (only pressure)
      return DRT::ELEMENTS::ArteryEleCalcPresBased<distype>::Instance(1, disname);
      break;
    }
    default:
    {
      dserror("Defined problem type %d does not exist!!", problem);
      break;
    }
  }

  return nullptr;
}

BACI_NAMESPACE_CLOSE

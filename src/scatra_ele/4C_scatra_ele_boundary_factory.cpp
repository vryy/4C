/*----------------------------------------------------------------------*/
/*! \file

\brief factory for scatra boundary evaluation

\level 2

 */
/*----------------------------------------------------------------------*/
#include "4C_scatra_ele_boundary_factory.hpp"

#include "4C_global_data.hpp"
#include "4C_mat_elchmat.hpp"
#include "4C_scatra_ele.hpp"
#include "4C_scatra_ele_boundary_calc_elch_diffcond.hpp"
#include "4C_scatra_ele_boundary_calc_elch_electrode.hpp"
#include "4C_scatra_ele_boundary_calc_elch_electrode_growth.hpp"
#include "4C_scatra_ele_boundary_calc_elch_electrode_sti_thermo.hpp"
#include "4C_scatra_ele_boundary_calc_elch_NP.hpp"
#include "4C_scatra_ele_boundary_calc_loma.hpp"
#include "4C_scatra_ele_boundary_calc_poro.hpp"
#include "4C_scatra_ele_boundary_calc_refconc_reac.hpp"
#include "4C_scatra_ele_boundary_calc_std.hpp"
#include "4C_scatra_ele_boundary_calc_sti_electrode.hpp"
#include "4C_scatra_ele_boundary_interface.hpp"
#include "4C_scatra_ele_calc.hpp"

FOUR_C_NAMESPACE_OPEN


/*--------------------------------------------------------------------------*
 |                                                 (public) rasthofer 11/13 |
 *--------------------------------------------------------------------------*/
Discret::ELEMENTS::ScaTraBoundaryInterface* Discret::ELEMENTS::ScaTraBoundaryFactory::ProvideImpl(
    const Core::Elements::Element* ele, const enum Inpar::ScaTra::ImplType impltype,
    const int numdofpernode, const int numscal, const std::string& disname)
{
  // number of space dimensions
  const int ndim = disname != "scatra_micro" ? Global::Problem::Instance()->NDim() : 1;

  switch (ele->Shape())
  {
    case Core::FE::CellType::quad4:
    {
      if (ndim == 3)
        return define_problem_type<Core::FE::CellType::quad4, 3>(
            impltype, numdofpernode, numscal, disname);
      else
        FOUR_C_THROW("invalid problem dimension!");
    }
    case Core::FE::CellType::quad8:
    {
      if (ndim == 3)
        return define_problem_type<Core::FE::CellType::quad8, 3>(
            impltype, numdofpernode, numscal, disname);
      else
        FOUR_C_THROW("invalid problem dimension!");
    }
    case Core::FE::CellType::quad9:
    {
      if (ndim == 3)
        return define_problem_type<Core::FE::CellType::quad9, 3>(
            impltype, numdofpernode, numscal, disname);
      else
        FOUR_C_THROW("invalid problem dimension!");
    }
    case Core::FE::CellType::tri3:
    {
      if (ndim == 3)
        return define_problem_type<Core::FE::CellType::tri3, 3>(
            impltype, numdofpernode, numscal, disname);
      else
        FOUR_C_THROW("invalid problem dimension!");
    }
    case Core::FE::CellType::tri6:
    {
      if (ndim == 3)
        return define_problem_type<Core::FE::CellType::tri6, 3>(
            impltype, numdofpernode, numscal, disname);
      else
        FOUR_C_THROW("invalid problem dimension!");
    }
    case Core::FE::CellType::line2:
    {
      if (ndim == 2)
        return define_problem_type<Core::FE::CellType::line2, 2>(
            impltype, numdofpernode, numscal, disname);
      else if (ndim == 3)
        return define_problem_type<Core::FE::CellType::line2, 3>(
            impltype, numdofpernode, numscal, disname);
      else
        FOUR_C_THROW("invalid problem dimension!");
    }
    case Core::FE::CellType::line3:
    {
      if (ndim == 2)
        return define_problem_type<Core::FE::CellType::line3, 2>(
            impltype, numdofpernode, numscal, disname);
      else
        FOUR_C_THROW("invalid problem dimension!");
    }
    case Core::FE::CellType::nurbs3:  // 1D nurbs boundary element
    {
      return define_problem_type<Core::FE::CellType::nurbs3, 2>(
          impltype, numdofpernode, numscal, disname);
    }
    case Core::FE::CellType::nurbs9:  // 2D nurbs boundary element
    {
      return define_problem_type<Core::FE::CellType::nurbs9, 3>(
          impltype, numdofpernode, numscal, disname);
    }
    default:
    {
      FOUR_C_THROW(
          "Element shape %d (%d nodes) not activated. Just do it.", ele->Shape(), ele->num_node());
      break;
    }
  }

  return nullptr;
}


/*-------------------------------------------------------------------------------------------*
 | return instance of element evaluation class depending on implementation type   fang 02/15 |
 *-------------------------------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
Discret::ELEMENTS::ScaTraBoundaryInterface*
Discret::ELEMENTS::ScaTraBoundaryFactory::define_problem_type(
    const enum Inpar::ScaTra::ImplType impltype, const int numdofpernode, const int numscal,
    const std::string& disname)
{
  switch (impltype)
  {
    case Inpar::ScaTra::impltype_advreac:
    case Inpar::ScaTra::impltype_aniso:
    case Inpar::ScaTra::impltype_cardiac_monodomain:
    case Inpar::ScaTra::impltype_chemo:
    case Inpar::ScaTra::impltype_chemoreac:
    case Inpar::ScaTra::impltype_levelset:
    case Inpar::ScaTra::impltype_std:
    case Inpar::ScaTra::impltype_thermo_elch_diffcond:
    case Inpar::ScaTra::impltype_multipororeac:
    {
      return Discret::ELEMENTS::ScaTraEleBoundaryCalcStd<distype, probdim>::Instance(
          numdofpernode, numscal, disname);
      break;
    }
    case Inpar::ScaTra::impltype_loma:
    {
      return Discret::ELEMENTS::ScaTraEleBoundaryCalcLoma<distype, probdim>::Instance(
          numdofpernode, numscal, disname);
      break;
    }
    case Inpar::ScaTra::impltype_elch_electrode:
    {
      return Discret::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<distype, probdim>::Instance(
          numdofpernode, numscal, disname);
      break;
    }
    case Inpar::ScaTra::impltype_elch_electrode_thermo:
    {
      return Discret::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<distype,
          probdim>::Instance(numdofpernode, numscal, disname);
      break;
    }
    case Inpar::ScaTra::impltype_elch_diffcond:
    case Inpar::ScaTra::impltype_elch_diffcond_thermo:
    {
      return Discret::ELEMENTS::ScaTraEleBoundaryCalcElchDiffCond<distype, probdim>::Instance(
          numdofpernode, numscal, disname);
      break;
    }
    case Inpar::ScaTra::impltype_elch_NP:
    {
      return Discret::ELEMENTS::ScaTraEleBoundaryCalcElchNP<distype, probdim>::Instance(
          numdofpernode, numscal, disname);
      break;
    }
    case Inpar::ScaTra::impltype_poro:
    case Inpar::ScaTra::impltype_pororeac:
    {
      return Discret::ELEMENTS::ScaTraEleBoundaryCalcPoro<distype, probdim>::Instance(
          numdofpernode, numscal, disname);
      break;
    }
    case Inpar::ScaTra::impltype_refconcreac:
    {
      return Discret::ELEMENTS::ScaTraEleBoundaryCalcRefConcReac<distype, probdim>::Instance(
          numdofpernode, numscal, disname);
      break;
    }
    case Inpar::ScaTra::impltype_thermo_elch_electrode:
    {
      return Discret::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<distype, probdim>::Instance(
          numdofpernode, numscal, disname);
      break;
    }
    case Inpar::ScaTra::impltype_elch_electrode_growth:
    {
      return Discret::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<distype,
          probdim>::Instance(numdofpernode, numscal, disname);
      break;
    }
    default:
    {
      FOUR_C_THROW("Defined implementation type does not exist!");
      break;
    }
  }  // switch(impltype)

  return nullptr;
}

FOUR_C_NAMESPACE_CLOSE

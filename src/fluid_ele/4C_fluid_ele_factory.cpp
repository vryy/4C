/*----------------------------------------------------------------------*/
/*! \file

\brief Factory class going from the generic evaluation routines to the ones
  templated by the element shape and specialization

\level 1


*/
/*----------------------------------------------------------------------*/


#include "4C_fluid_ele_factory.hpp"

#include "4C_fluid_ele_calc_hdg.hpp"
#include "4C_fluid_ele_calc_hdg_weak_comp.hpp"
#include "4C_fluid_ele_calc_immersed.hpp"
#include "4C_fluid_ele_calc_loma.hpp"
#include "4C_fluid_ele_calc_poro.hpp"
#include "4C_fluid_ele_calc_poro_p1.hpp"
#include "4C_fluid_ele_calc_std.hpp"
#include "4C_fluid_ele_calc_xfem.hpp"
#include "4C_fluid_ele_calc_xwall.hpp"
#include "4C_fluid_ele_interface.hpp"

FOUR_C_NAMESPACE_OPEN

/*--------------------------------------------------------------------------*
 |                                                 (public) rasthofer Jan13 |
 *--------------------------------------------------------------------------*/
Discret::ELEMENTS::FluidEleInterface* Discret::ELEMENTS::FluidFactory::ProvideImpl(
    Core::FE::CellType distype, std::string problem)
{
  switch (distype)
  {
    case Core::FE::CellType::hex8:
    {
      return define_problem_type<Core::FE::CellType::hex8>(problem);
    }
    case Core::FE::CellType::hex20:
    {
      return define_problem_type<Core::FE::CellType::hex20>(problem);
    }
    case Core::FE::CellType::hex27:
    {
      return define_problem_type<Core::FE::CellType::hex27>(problem);
    }
    case Core::FE::CellType::tet4:
    {
      return define_problem_type<Core::FE::CellType::tet4>(problem);
    }
    case Core::FE::CellType::tet10:
    {
      return define_problem_type<Core::FE::CellType::tet10>(problem);
    }
    case Core::FE::CellType::wedge6:
    {
      return define_problem_type<Core::FE::CellType::wedge6>(problem);
    }
    case Core::FE::CellType::wedge15:
    {
      return define_problem_type<Core::FE::CellType::wedge15>(problem);
    }
    case Core::FE::CellType::pyramid5:
    {
      return define_problem_type<Core::FE::CellType::pyramid5>(problem);
    }
    case Core::FE::CellType::quad4:
    {
      return define_problem_type<Core::FE::CellType::quad4>(problem);
    }
    case Core::FE::CellType::quad8:
    {
      return define_problem_type<Core::FE::CellType::quad8>(problem);
    }
    case Core::FE::CellType::quad9:
    {
      return define_problem_type<Core::FE::CellType::quad9>(problem);
    }
    case Core::FE::CellType::tri3:
    {
      return define_problem_type<Core::FE::CellType::tri3>(problem);
    }
    case Core::FE::CellType::tri6:
    {
      return define_problem_type<Core::FE::CellType::tri6>(problem);
    }
    // Nurbs support
    case Core::FE::CellType::nurbs9:
    {
      return define_problem_type<Core::FE::CellType::nurbs9>(problem);
    }
    case Core::FE::CellType::nurbs27:
    {
      return define_problem_type<Core::FE::CellType::nurbs27>(problem);
    }
    // no 1D elements
    default:
      FOUR_C_THROW("Element shape %s not activated. Just do it.",
          Core::FE::CellTypeToString(distype).c_str());
      break;
  }
  return nullptr;
}

/*--------------------------------------------------------------------------*
 |                                                 (public) rasthofer Jan13 |
 *--------------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Discret::ELEMENTS::FluidEleInterface* Discret::ELEMENTS::FluidFactory::define_problem_type(
    std::string problem)
{
  if (problem == "std")
    return Discret::ELEMENTS::FluidEleCalcStd<distype>::Instance();
  else if (problem == "loma")
    return Discret::ELEMENTS::FluidEleCalcLoma<distype>::Instance();
  else if (problem == "std_immersed")
    return Discret::ELEMENTS::FluidEleCalcImmersed<distype>::Instance();
  else if (problem == "poro")
    return Discret::ELEMENTS::FluidEleCalcPoro<distype>::Instance();
  else if (problem == "poro_p1")
    return Discret::ELEMENTS::FluidEleCalcPoroP1<distype>::Instance();
  else if (problem == "hdg")
    return Discret::ELEMENTS::FluidEleCalcHDG<distype>::Instance();
  else if (problem == "hdgweakcomp")
    return Discret::ELEMENTS::FluidEleCalcHDGWeakComp<distype>::Instance();
  else if (problem == "xw")
  {
    // for now we only build the hex8 and tet4 elements for xwall
    // later we might consider other kinds of elements
    if (distype == Core::FE::CellType::hex8)
      return Discret::ELEMENTS::FluidEleCalcXWall<Core::FE::CellType::hex8,
          Discret::ELEMENTS::Fluid::xwall>::Instance();
    else if (distype == Core::FE::CellType::tet4)
      return Discret::ELEMENTS::FluidEleCalcXWall<Core::FE::CellType::tet4,
          Discret::ELEMENTS::Fluid::xwall>::Instance();
    else
      FOUR_C_THROW("only hex8 and tet4 elements compiled for xwall");
  }
  else
    FOUR_C_THROW("Defined problem type does not exist!!");

  return nullptr;
}

/*--------------------------------------------------------------------------*
 |  special implementation of ProvideImpl for XFEM problems                 |
 |  to reduce created template combination         (public) rasthofer Jan13 |
 *--------------------------------------------------------------------------*/
Discret::ELEMENTS::FluidEleInterface* Discret::ELEMENTS::FluidFactory::ProvideImplXFEM(
    Core::FE::CellType distype, std::string problem)
{
  if (problem != "xfem") FOUR_C_THROW("Call ProvideImplXFEM just for xfem problems!");

  switch (distype)
  {
    // only 3D elements
    case Core::FE::CellType::hex8:
    {
      return define_problem_type_xfem<Core::FE::CellType::hex8>(problem);
    }
    case Core::FE::CellType::hex20:
    {
      return define_problem_type_xfem<Core::FE::CellType::hex20>(problem);
    }
    case Core::FE::CellType::hex27:
    {
      return define_problem_type_xfem<Core::FE::CellType::hex27>(problem);
    }
    case Core::FE::CellType::tet4:
    {
      return define_problem_type_xfem<Core::FE::CellType::tet4>(problem);
    }
    case Core::FE::CellType::tet10:
    {
      return define_problem_type_xfem<Core::FE::CellType::tet10>(problem);
    }
    case Core::FE::CellType::wedge6:
    {
      return define_problem_type_xfem<Core::FE::CellType::wedge6>(problem);
    }
    case Core::FE::CellType::wedge15:
    {
      return define_problem_type_xfem<Core::FE::CellType::wedge15>(problem);
    }
      //    case Core::FE::CellType::pyramid5:
      //    {
      //      return define_problem_type_xfem<Core::FE::CellType::pyramid5>(problem);
      //    }
    default:
      FOUR_C_THROW("Element shape %s not activated for XFEM problems. Just do it.",
          Core::FE::CellTypeToString(distype).c_str());
      break;
  }
  return nullptr;
}

/*--------------------------------------------------------------------------*
 |  special implementation of DefineProblemTypeX for XFEM problems          |
 |  to reduce created template combination         (public) rasthofer Jan13 |
 *--------------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Discret::ELEMENTS::FluidEleInterface* Discret::ELEMENTS::FluidFactory::define_problem_type_xfem(
    std::string problem)
{
  if (problem == "xfem")
    return Discret::ELEMENTS::FluidEleCalcXFEM<distype>::Instance();
  else
    FOUR_C_THROW("Defined problem type does not exist!!");

  return nullptr;
}

FOUR_C_NAMESPACE_CLOSE

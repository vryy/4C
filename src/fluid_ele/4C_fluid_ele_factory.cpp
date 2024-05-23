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
DRT::ELEMENTS::FluidEleInterface* DRT::ELEMENTS::FluidFactory::ProvideImpl(
    CORE::FE::CellType distype, std::string problem)
{
  switch (distype)
  {
    case CORE::FE::CellType::hex8:
    {
      return DefineProblemType<CORE::FE::CellType::hex8>(problem);
    }
    case CORE::FE::CellType::hex20:
    {
      return DefineProblemType<CORE::FE::CellType::hex20>(problem);
    }
    case CORE::FE::CellType::hex27:
    {
      return DefineProblemType<CORE::FE::CellType::hex27>(problem);
    }
    case CORE::FE::CellType::tet4:
    {
      return DefineProblemType<CORE::FE::CellType::tet4>(problem);
    }
    case CORE::FE::CellType::tet10:
    {
      return DefineProblemType<CORE::FE::CellType::tet10>(problem);
    }
    case CORE::FE::CellType::wedge6:
    {
      return DefineProblemType<CORE::FE::CellType::wedge6>(problem);
    }
    case CORE::FE::CellType::wedge15:
    {
      return DefineProblemType<CORE::FE::CellType::wedge15>(problem);
    }
    case CORE::FE::CellType::pyramid5:
    {
      return DefineProblemType<CORE::FE::CellType::pyramid5>(problem);
    }
    case CORE::FE::CellType::quad4:
    {
      return DefineProblemType<CORE::FE::CellType::quad4>(problem);
    }
    case CORE::FE::CellType::quad8:
    {
      return DefineProblemType<CORE::FE::CellType::quad8>(problem);
    }
    case CORE::FE::CellType::quad9:
    {
      return DefineProblemType<CORE::FE::CellType::quad9>(problem);
    }
    case CORE::FE::CellType::tri3:
    {
      return DefineProblemType<CORE::FE::CellType::tri3>(problem);
    }
    case CORE::FE::CellType::tri6:
    {
      return DefineProblemType<CORE::FE::CellType::tri6>(problem);
    }
    // Nurbs support
    case CORE::FE::CellType::nurbs9:
    {
      return DefineProblemType<CORE::FE::CellType::nurbs9>(problem);
    }
    case CORE::FE::CellType::nurbs27:
    {
      return DefineProblemType<CORE::FE::CellType::nurbs27>(problem);
    }
    // no 1D elements
    default:
      FOUR_C_THROW("Element shape %s not activated. Just do it.",
          CORE::FE::CellTypeToString(distype).c_str());
      break;
  }
  return nullptr;
}

/*--------------------------------------------------------------------------*
 |                                                 (public) rasthofer Jan13 |
 *--------------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
DRT::ELEMENTS::FluidEleInterface* DRT::ELEMENTS::FluidFactory::DefineProblemType(
    std::string problem)
{
  if (problem == "std")
    return DRT::ELEMENTS::FluidEleCalcStd<distype>::Instance();
  else if (problem == "loma")
    return DRT::ELEMENTS::FluidEleCalcLoma<distype>::Instance();
  else if (problem == "std_immersed")
    return DRT::ELEMENTS::FluidEleCalcImmersed<distype>::Instance();
  else if (problem == "poro")
    return DRT::ELEMENTS::FluidEleCalcPoro<distype>::Instance();
  else if (problem == "poro_p1")
    return DRT::ELEMENTS::FluidEleCalcPoroP1<distype>::Instance();
  else if (problem == "hdg")
    return DRT::ELEMENTS::FluidEleCalcHDG<distype>::Instance();
  else if (problem == "hdgweakcomp")
    return DRT::ELEMENTS::FluidEleCalcHDGWeakComp<distype>::Instance();
  else if (problem == "xw")
  {
    // for now we only build the hex8 and tet4 elements for xwall
    // later we might consider other kinds of elements
    if (distype == CORE::FE::CellType::hex8)
      return DRT::ELEMENTS::FluidEleCalcXWall<CORE::FE::CellType::hex8,
          DRT::ELEMENTS::Fluid::xwall>::Instance();
    else if (distype == CORE::FE::CellType::tet4)
      return DRT::ELEMENTS::FluidEleCalcXWall<CORE::FE::CellType::tet4,
          DRT::ELEMENTS::Fluid::xwall>::Instance();
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
DRT::ELEMENTS::FluidEleInterface* DRT::ELEMENTS::FluidFactory::ProvideImplXFEM(
    CORE::FE::CellType distype, std::string problem)
{
  if (problem != "xfem") FOUR_C_THROW("Call ProvideImplXFEM just for xfem problems!");

  switch (distype)
  {
    // only 3D elements
    case CORE::FE::CellType::hex8:
    {
      return define_problem_type_xfem<CORE::FE::CellType::hex8>(problem);
    }
    case CORE::FE::CellType::hex20:
    {
      return define_problem_type_xfem<CORE::FE::CellType::hex20>(problem);
    }
    case CORE::FE::CellType::hex27:
    {
      return define_problem_type_xfem<CORE::FE::CellType::hex27>(problem);
    }
    case CORE::FE::CellType::tet4:
    {
      return define_problem_type_xfem<CORE::FE::CellType::tet4>(problem);
    }
    case CORE::FE::CellType::tet10:
    {
      return define_problem_type_xfem<CORE::FE::CellType::tet10>(problem);
    }
    case CORE::FE::CellType::wedge6:
    {
      return define_problem_type_xfem<CORE::FE::CellType::wedge6>(problem);
    }
    case CORE::FE::CellType::wedge15:
    {
      return define_problem_type_xfem<CORE::FE::CellType::wedge15>(problem);
    }
      //    case CORE::FE::CellType::pyramid5:
      //    {
      //      return define_problem_type_xfem<CORE::FE::CellType::pyramid5>(problem);
      //    }
    default:
      FOUR_C_THROW("Element shape %s not activated for XFEM problems. Just do it.",
          CORE::FE::CellTypeToString(distype).c_str());
      break;
  }
  return nullptr;
}

/*--------------------------------------------------------------------------*
 |  special implementation of DefineProblemTypeX for XFEM problems          |
 |  to reduce created template combination         (public) rasthofer Jan13 |
 *--------------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
DRT::ELEMENTS::FluidEleInterface* DRT::ELEMENTS::FluidFactory::define_problem_type_xfem(
    std::string problem)
{
  if (problem == "xfem")
    return DRT::ELEMENTS::FluidEleCalcXFEM<distype>::Instance();
  else
    FOUR_C_THROW("Defined problem type does not exist!!");

  return nullptr;
}

FOUR_C_NAMESPACE_CLOSE

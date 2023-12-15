/*----------------------------------------------------------------------*/
/*! \file

\brief Factory class going from the generic evaluation routines to the ones
  templated by the element shape and specialization

\level 1


*/
/*----------------------------------------------------------------------*/


#include "baci_fluid_ele_factory.H"

#include "baci_fluid_ele_calc_hdg.H"
#include "baci_fluid_ele_calc_hdg_weak_comp.H"
#include "baci_fluid_ele_calc_immersed.H"
#include "baci_fluid_ele_calc_loma.H"
#include "baci_fluid_ele_calc_poro.H"
#include "baci_fluid_ele_calc_poro_p1.H"
#include "baci_fluid_ele_calc_std.H"
#include "baci_fluid_ele_calc_xfem.H"
#include "baci_fluid_ele_calc_xwall.H"
#include "baci_fluid_ele_interface.H"

BACI_NAMESPACE_OPEN

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
      dserror("Element shape %s not activated. Just do it.",
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
      dserror("only hex8 and tet4 elements compiled for xwall");
  }
  else
    dserror("Defined problem type does not exist!!");

  return nullptr;
}

/*--------------------------------------------------------------------------*
 |  special implementation of ProvideImpl for XFEM problems                 |
 |  to reduce created template combination         (public) rasthofer Jan13 |
 *--------------------------------------------------------------------------*/
DRT::ELEMENTS::FluidEleInterface* DRT::ELEMENTS::FluidFactory::ProvideImplXFEM(
    CORE::FE::CellType distype, std::string problem)
{
  if (problem != "xfem") dserror("Call ProvideImplXFEM just for xfem problems!");

  switch (distype)
  {
    // only 3D elements
    case CORE::FE::CellType::hex8:
    {
      return DefineProblemTypeXFEM<CORE::FE::CellType::hex8>(problem);
    }
    case CORE::FE::CellType::hex20:
    {
      return DefineProblemTypeXFEM<CORE::FE::CellType::hex20>(problem);
    }
    case CORE::FE::CellType::hex27:
    {
      return DefineProblemTypeXFEM<CORE::FE::CellType::hex27>(problem);
    }
    case CORE::FE::CellType::tet4:
    {
      return DefineProblemTypeXFEM<CORE::FE::CellType::tet4>(problem);
    }
    case CORE::FE::CellType::tet10:
    {
      return DefineProblemTypeXFEM<CORE::FE::CellType::tet10>(problem);
    }
    case CORE::FE::CellType::wedge6:
    {
      return DefineProblemTypeXFEM<CORE::FE::CellType::wedge6>(problem);
    }
    case CORE::FE::CellType::wedge15:
    {
      return DefineProblemTypeXFEM<CORE::FE::CellType::wedge15>(problem);
    }
      //    case CORE::FE::CellType::pyramid5:
      //    {
      //      return DefineProblemTypeXFEM<CORE::FE::CellType::pyramid5>(problem);
      //    }
    default:
      dserror("Element shape %s not activated for XFEM problems. Just do it.",
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
DRT::ELEMENTS::FluidEleInterface* DRT::ELEMENTS::FluidFactory::DefineProblemTypeXFEM(
    std::string problem)
{
  if (problem == "xfem")
    return DRT::ELEMENTS::FluidEleCalcXFEM<distype>::Instance();
  else
    dserror("Defined problem type does not exist!!");

  return nullptr;
}

BACI_NAMESPACE_CLOSE

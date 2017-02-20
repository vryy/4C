/*--------------------------------------------------------------------------*/
/*!
\file fluid_ele_factory.cpp

\brief Factory of fluid elements

<pre>
\maintainer Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15236/245
</pre>

\level 1
*/
/*--------------------------------------------------------------------------*/

#include "fluid_ele_factory.H"
#include "fluid_ele_interface.H"

#include "fluid_ele_calc_std.H"
#include "fluid_ele_calc_loma.H"
#include "fluid_ele_calc_immersed.H"
#include "fluid_ele_calc_poro.H"
#include "fluid_ele_calc_poro_p1.H"
#include "fluid_ele_calc_poro_p1_immersed.H"
#include "fluid_ele_calc_xfem.H"
#include "fluid_ele_calc_xwall.H"
#include "fluid_ele_calc_hdg.H"

#include "../drt_meshfree_discret/meshfree_fluid_cell_interface.H"

#include "../drt_meshfree_discret/meshfree_fluid_cell_calc_std.H"

/*--------------------------------------------------------------------------*
 |                                                 (public) rasthofer Jan13 |
 *--------------------------------------------------------------------------*/
DRT::ELEMENTS::FluidEleInterface* DRT::ELEMENTS::FluidFactory::ProvideImpl(DRT::Element::DiscretizationType distype, std::string problem)
{
  switch(distype)
  {
    case DRT::Element::hex8:
    {
      return DefineProblemType<DRT::Element::hex8>(problem);
    }
    case DRT::Element::hex20:
    {
      return DefineProblemType<DRT::Element::hex20>(problem);
    }
    case DRT::Element::hex27:
    {
      return DefineProblemType<DRT::Element::hex27>(problem);
    }
    case DRT::Element::tet4:
    {
      return DefineProblemType<DRT::Element::tet4>(problem);
    }
    case DRT::Element::tet10:
    {
      return DefineProblemType<DRT::Element::tet10>(problem);
    }
    case DRT::Element::wedge6:
    {
      return DefineProblemType<DRT::Element::wedge6>(problem);
    }
    case DRT::Element::wedge15:
    {
      return DefineProblemType<DRT::Element::wedge15>(problem);
    }
    case DRT::Element::pyramid5:
    {
      return DefineProblemType<DRT::Element::pyramid5>(problem);
    }
    case DRT::Element::quad4:
    {
      return DefineProblemType<DRT::Element::quad4>(problem);
    }
    case DRT::Element::quad8:
    {
      return DefineProblemType<DRT::Element::quad8>(problem);
    }
    case DRT::Element::quad9:
    {
      return DefineProblemType<DRT::Element::quad9>(problem);
    }
    case DRT::Element::tri3:
    {
      return DefineProblemType<DRT::Element::tri3>(problem);
    }
    case DRT::Element::tri6:
    {
      return DefineProblemType<DRT::Element::tri6>(problem);
    }
    // Nurbs support
    case DRT::Element::nurbs9:
    {
      return DefineProblemType<DRT::Element::nurbs9>(problem);
    }
    case DRT::Element::nurbs27:
    {
      return DefineProblemType<DRT::Element::nurbs27>(problem);
    }
    // no 1D elements
    default:
      dserror("Element shape %s not activated. Just do it.",DRT::DistypeToString(distype).c_str());
      break;
    }
  return NULL;
}

/*--------------------------------------------------------------------------*
 |                                                 (public) rasthofer Jan13 |
 *--------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::FluidEleInterface* DRT::ELEMENTS::FluidFactory::DefineProblemType(std::string problem)
{
  if (problem == "std")
    return DRT::ELEMENTS::FluidEleCalcStd<distype>::Instance();
  else if (problem == "loma")
    return DRT::ELEMENTS::FluidEleCalcLoma<distype>::Instance();
  else if (problem == "std_immersed")
    return DRT::ELEMENTS::FluidEleCalcImmersed<distype>::Instance();
  else if (problem == "poro_p1_immersed")
    return DRT::ELEMENTS::FluidEleCalcPoroP1Immersed<distype>::Instance();
  else if (problem == "poro")
    return DRT::ELEMENTS::FluidEleCalcPoro<distype>::Instance();
  else if (problem == "poro_p1")
    return DRT::ELEMENTS::FluidEleCalcPoroP1<distype>::Instance();
  else if (problem == "hdg")
    return DRT::ELEMENTS::FluidEleCalcHDG<distype>::Instance();
  else if (problem == "xw")
  {
    // for now we only build the hex8 and tet4 elements for xwall
    // later we might consider other kinds of elements
    if (distype == DRT::Element::hex8)
      return DRT::ELEMENTS::FluidEleCalcXWall<DRT::Element::hex8,DRT::ELEMENTS::Fluid::xwall>::Instance();
    else if(distype == DRT::Element::tet4)
      return DRT::ELEMENTS::FluidEleCalcXWall<DRT::Element::tet4,DRT::ELEMENTS::Fluid::xwall>::Instance();
    else dserror("only hex8 and tet4 elements compiled for xwall");
  }
  else
    dserror("Defined problem type does not exist!!");

  return NULL;
}

/*--------------------------------------------------------------------------*
 |  special implementation of ProvideImpl for XFEM problems                 |
 |  to reduce created template combination         (public) rasthofer Jan13 |
 *--------------------------------------------------------------------------*/
DRT::ELEMENTS::FluidEleInterface* DRT::ELEMENTS::FluidFactory::ProvideImplXFEM(DRT::Element::DiscretizationType distype, std::string problem)
{
  if(problem != "xfem") dserror("Call ProvideImplXFEM just for xfem problems!");

  switch(distype)
  {
    // only 3D elements
    case DRT::Element::hex8:
    {
      return DefineProblemTypeXFEM<DRT::Element::hex8>(problem);
    }
    case DRT::Element::hex20:
    {
      return DefineProblemTypeXFEM<DRT::Element::hex20>(problem);
    }
    case DRT::Element::hex27:
    {
      return DefineProblemTypeXFEM<DRT::Element::hex27>(problem);
    }
    case DRT::Element::tet4:
    {
      return DefineProblemTypeXFEM<DRT::Element::tet4>(problem);
    }
    case DRT::Element::tet10:
    {
      return DefineProblemTypeXFEM<DRT::Element::tet10>(problem);
    }
    case DRT::Element::wedge6:
    {
      return DefineProblemTypeXFEM<DRT::Element::wedge6>(problem);
    }
//    case DRT::Element::pyramid5:
//    {
//      return DefineProblemTypeXFEM<DRT::Element::pyramid5>(problem);
//    }
    default:
      dserror("Element shape %s not activated for XFEM problems. Just do it.",DRT::DistypeToString(distype).c_str());
      break;
    }
  return NULL;
}

/*--------------------------------------------------------------------------*
 |  special implementation of DefineProblemTypeX for XFEM problems          |
 |  to reduce created template combination         (public) rasthofer Jan13 |
 *--------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::FluidEleInterface* DRT::ELEMENTS::FluidFactory::DefineProblemTypeXFEM(std::string problem)
{
  if (problem == "xfem")
    return DRT::ELEMENTS::FluidEleCalcXFEM<distype>::Instance();
  else
    dserror("Defined problem type does not exist!!");

  return NULL;
}

/*--------------------------------------------------------------------------*
 |  special implementation of ProvideImpl for meshfree problems             |
 |  to reduce created template combination               (public) nis Nov13 |
 *--------------------------------------------------------------------------*/
DRT::ELEMENTS::MeshfreeFluidCellInterface* DRT::ELEMENTS::FluidFactory::ProvideImplMeshfree(DRT::Element::DiscretizationType distype, std::string problem)
{
  switch(distype)
  {
    case DRT::Element::hex8:
    {
      return DefineProblemTypeMeshfree<DRT::Element::hex8>(problem);
    }
    case DRT::Element::tet4:
    {
      return DefineProblemTypeMeshfree<DRT::Element::tet4>(problem);
    }
    case DRT::Element::quad4:
    {
      return DefineProblemTypeMeshfree<DRT::Element::quad4>(problem);
    }
    case DRT::Element::tri3:
    {
      return DefineProblemTypeMeshfree<DRT::Element::tri3>(problem);
    }
    default:
      dserror("Element shape %s not activated for meshfree problems. ",DRT::DistypeToString(distype).c_str());
      break;
    }
  return NULL;
}

/*--------------------------------------------------------------------------*
 |                                                       (public) nis Nov13 |
 *--------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::MeshfreeFluidCellInterface* DRT::ELEMENTS::FluidFactory::DefineProblemTypeMeshfree(std::string problem)
{
  if (problem == "std_meshfree")
    return DRT::ELEMENTS::MeshfreeFluidCellCalcStd<distype>::Instance();
  else
    dserror("Defined problem type does not exist for meshfree problems!!");

  return NULL;
}



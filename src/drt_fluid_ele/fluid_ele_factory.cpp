/*--------------------------------------------------------------------------*/
/*!
\file fluid_ele_factory.cpp

\brief Factory of fluid elements

<pre>
Maintainer: Ursula Rasthofer & Volker Gravemeier
            {rasthofer,vgravem}@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15236/245
</pre>
*/
/*--------------------------------------------------------------------------*/

#include "fluid_ele_factory.H"

#include "fluid_ele_calc_std.H"
#include "fluid_ele_calc_loma.H"
#include "fluid_ele_calc_poro.H"
#include "fluid_ele_calc_xfem.H"


/*--------------------------------------------------------------------------*
 |                                                 (public) rasthofer Jan13 |
 *--------------------------------------------------------------------------*/
DRT::ELEMENTS::FluidEleInterface* DRT::ELEMENTS::FluidFactory::ProvideImpl(DRT::Element::DiscretizationType distype, string problem)
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
    /* wedge15 cannot be used since no mesh generator exists
    case DRT::Element::wedge15:
    {
      return DefineProblemType<DRT::Element::wedge15>(problem);
    }
    */
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
    }
  return NULL;
}

/*--------------------------------------------------------------------------*
 |                                                 (public) rasthofer Jan13 |
 *--------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::FluidEleInterface* DRT::ELEMENTS::FluidFactory::DefineProblemType(string problem)
{
  if (problem == "std")
    return DRT::ELEMENTS::FluidEleCalcStd<distype>::Instance();
  else if (problem == "loma")
    return DRT::ELEMENTS::FluidEleCalcLoma<distype>::Instance();
  else if (problem == "poro")
    return DRT::ELEMENTS::FluidEleCalcPoro<distype>::Instance();
  else
    dserror("Defined problem type does not exist!!");

  return NULL;
}

/*--------------------------------------------------------------------------*
 |  special implementation of ProvideImpl for XFEM problems                 |
 |  to reduce created template combination         (public) rasthofer Jan13 |
 *--------------------------------------------------------------------------*/
DRT::ELEMENTS::FluidEleInterface* DRT::ELEMENTS::FluidFactory::ProvideImplXFEM(DRT::Element::DiscretizationType distype, string problem)
{
  if(problem != "xfem") dserror("Call ProvideImplXFEM just for xfem problems!");

  switch(distype)
  {
    // only 3D elements
    case DRT::Element::hex8:
    {
      return DefineProblemType<DRT::Element::hex8>(problem);
    }
    case DRT::Element::hex20:
    {
      return DefineProblemType<DRT::Element::hex20>(problem);
    }
    default:
      dserror("Element shape %s not activated for XFEM problems. Just do it.",DRT::DistypeToString(distype).c_str());
    }
  return NULL;
}

/*--------------------------------------------------------------------------*
 |                                                 (public) rasthofer Jan13 |
 *--------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::FluidEleInterface* DRT::ELEMENTS::FluidFactory::DefineProblemTypeXFEM(string problem)
{
  if (problem == "xfem")
    return DRT::ELEMENTS::FluidEleCalcXFEM<distype>::Instance();
  else
    dserror("Defined problem type does not exist!!");

  return NULL;
}

/*--------------------------------------------------------------------------*
 |  special implementation of ProvideImpl for meshfree problems             |
 |  to reduce created template combination         (public) rasthofer Jan13 |
 *--------------------------------------------------------------------------*/
DRT::ELEMENTS::FluidEleInterface* DRT::ELEMENTS::FluidFactory::ProvideImplMeshfree(DRT::Element::DiscretizationType distype, string problem)
{
  if(problem != "std_meshfree") dserror("Call ProvideImplMeshfree just for meshfree problems!");

  switch(distype)
  {
    case DRT::Element::hex8:
    {
      return DefineProblemType<DRT::Element::hex8>(problem);
    }
    case DRT::Element::tet4:
    {
      return DefineProblemType<DRT::Element::tet4>(problem);
    }
    case DRT::Element::quad4:
    {
      return DefineProblemType<DRT::Element::quad4>(problem);
    }
    case DRT::Element::tri3:
    {
      return DefineProblemType<DRT::Element::tri3>(problem);
    }
    default:
      dserror("Element shape %s not activated for meshfree problems. ",DRT::DistypeToString(distype).c_str());
    }
  return NULL;
}

/*--------------------------------------------------------------------------*
 |                                                 (public) rasthofer Jan13 |
 *--------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::FluidEleInterface* DRT::ELEMENTS::FluidFactory::DefineProblemTypeMeshfree(string problem)
{
  if (problem == "std_meshfree")
    return NULL;//DRT::ELEMENTS::FluidEleCalcXFEM<distype>::Instance();
  else
    dserror("Defined problem type does not exist!!");

  return NULL;
}


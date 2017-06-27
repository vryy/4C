/*----------------------------------------------------------------------*/
/*!
 \file porofluidmultiphase_ele_factory.cpp

 \brief factory class providing the implementations of the porofluidmultiphase
        element evaluation routines

   \level 3

   \maintainer  Lena Yoshihara
                yoshihara@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
 *----------------------------------------------------------------------*/



#include "porofluidmultiphase_ele_factory.H"

#include "../drt_lib/drt_globalproblem.H"
#include "porofluidmultiphase_ele_calc.H"

/*--------------------------------------------------------------------------*
 | provide the implementation of evaluation class      (public) vuong 08/16 |
 *--------------------------------------------------------------------------*/
DRT::ELEMENTS::PoroFluidMultiPhaseEleInterface* DRT::ELEMENTS::PoroFluidMultiPhaseFactory::ProvideImpl(
    DRT::Element::DiscretizationType distype,
    const int numdofpernode,
    const std::string& disname)
{
  // -------------------------------------- number of degrees of freedom
  // number of degrees of freedom
  static const int ndim = DRT::Problem::Instance()->NDim();

  switch(distype)
  {
  case DRT::Element::quad4:
  {
    if(ndim==2)
      return DefineProblemType<DRT::Element::quad4>(numdofpernode,disname);
    else
      dserror("invalid problem dimension for quad4 porofluidmultiphase element!");
    break;
  }
//  case DRT::Element::quad8:
//  {
//    if(ndim==2)
//      return DefineProblemType<DRT::Element::quad8>(numdofpernode,disname);
//    else
//      dserror("invalid problem dimension for quad8 porofluidmultiphase element!");
//    break;
//  }
  case DRT::Element::quad9:
  {
    if(ndim==2)
      return DefineProblemType<DRT::Element::quad9>(numdofpernode,disname);
    else
      dserror("invalid problem dimension for quad9 porofluidmultiphase element!");
    break;
  }
  case DRT::Element::tri3:
  {
    if(ndim==2)
      return DefineProblemType<DRT::Element::tri3>(numdofpernode,disname);
    else
      dserror("invalid problem dimension for tri3 porofluidmultiphase element!");
    break;
  }
//  case DRT::Element::tri6:
//  {
//    if(ndim==2)
//      return DefineProblemType<DRT::Element::tri6>(numdofpernode,disname);
//    else
//      dserror("invalid problem dimension for tri6 porofluidmultiphase element!");
//    break;
//  }
  case DRT::Element::line2:
  {
    if(ndim==1)
      return DefineProblemType<DRT::Element::line2>(numdofpernode,disname);
    else
      dserror("invalid problem dimension for line2 porofluidmultiphase element!");
    break;
  }
  case DRT::Element::line3:
  {
    if(ndim==1)
      return DefineProblemType<DRT::Element::line3>(numdofpernode,disname);
    else
      dserror("invalid problem dimension for line3 porofluidmultiphase element!");
    break;
  }
  case DRT::Element::hex8:
  {
    if(ndim==3)
      return DefineProblemType<DRT::Element::hex8>(numdofpernode,disname);
    else
      dserror("invalid problem dimension for hex8 porofluidmultiphase element!");
    break;
  }
  case DRT::Element::hex27:
  {
    if(ndim==3)
      return DefineProblemType<DRT::Element::hex27>(numdofpernode,disname);
    else
      dserror("invalid problem dimension for hex27 porofluidmultiphase element!");
    break;
  }
  case DRT::Element::tet4:
  {
    if(ndim==3)
      return DefineProblemType<DRT::Element::tet4>(numdofpernode,disname);
    else
      dserror("invalid problem dimension for tet4 porofluidmultiphase element!");
    break;
  }
  default:
    dserror("Element shape %s not activated. Just do it.",DRT::DistypeToString(distype).c_str());
    break;
  }
  return NULL;
}

/*--------------------------------------------------------------------------*
 | provide the implementation of evaluation class      (public) vuong 08/16 |
 *--------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::PoroFluidMultiPhaseEleInterface* DRT::ELEMENTS::PoroFluidMultiPhaseFactory::DefineProblemType(
    const int numdofpernode,
    const std::string& disname)
{

  return DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<distype>::Instance(numdofpernode,disname);

}


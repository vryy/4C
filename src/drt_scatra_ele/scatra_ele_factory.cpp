/*--------------------------------------------------------------------------*/
/*!
\file scatra_ele_factory.cpp

\brief Factory of scatra elements

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15252
</pre>
*/
/*--------------------------------------------------------------------------*/

#include "scatra_ele_factory.H"

#include "scatra_ele_calc_std.H"
#include "scatra_ele_calc_ls.H"
#include "scatra_ele_calc_lsreinit.H"
#include "scatra_ele_calc_elch_NP.H"
#include "scatra_ele_calc_elch_diffcond.H"
#include "scatra_ele_calc_loma.H"
#include "scatra_ele_calc_poro.H"
#include "scatra_ele_calc_advanced_reaction.H"
#include "scatra_ele_calc_poro_reac.H"

#include "../drt_meshfree_discret/meshfree_scatra_cell_calc_std.H"

#include "../drt_lib/drt_element.H"


/*--------------------------------------------------------------------------*
 |                                               (public) ehrl        Dec13 |
 *--------------------------------------------------------------------------*/
DRT::ELEMENTS::ScaTraEleInterface* DRT::ELEMENTS::ScaTraFactory::ProvideImpl(
  DRT::Element::DiscretizationType distype,
  INPAR::SCATRA::ImplType problem,
  const int numdofpernode,
  const int numscal)
{
  switch(distype)
  {
  case DRT::Element::hex8:
  {
    return DefineProblemType<DRT::Element::hex8>(problem,numdofpernode,numscal);
  }
//  case DRT::Element::hex20:
//  {
//    return ScaTraImpl<DRT::Element::hex20>::Instance(problem,numdofpernode,numscal);
//  } */
  case DRT::Element::hex27:
  {
    return DefineProblemType<DRT::Element::hex27>(problem,numdofpernode,numscal);
  }
//  case DRT::Element::nurbs8:
//  {
//    return ScaTraImpl<DRT::Element::nurbs8>::Instance(problem,numdofpernode,numscal);
//  } */
//  case DRT::Element::nurbs27:
//  {
//    return ScaTraImpl<DRT::Element::nurbs27>::Instance(problem,numdofpernode,numscal);
//  }
  case DRT::Element::tet4:
  {
    return DefineProblemType<DRT::Element::tet4>(problem,numdofpernode,numscal);
  }
  case DRT::Element::tet10:
  {
    return DefineProblemType<DRT::Element::tet10>(problem,numdofpernode,numscal);
  }
//  case DRT::Element::wedge6:
//  {
//    return ScaTraImpl<DRT::Element::wedge6>::Instance(problem,numdofpernode,numscal);
//  } /*
//  case DRT::Element::wedge15:
//  {
//    return ScaTraImpl<DRT::Element::wedge15>::Instance(problem,numdofpernode,numscal);
//  } */
//  case DRT::Element::pyramid5:
//  {
//    return ScaTraImpl<DRT::Element::pyramid5>::Instance(problem,numdofpernode,numscal);
//  }
  case DRT::Element::quad4:
  {
    return DefineProblemType<DRT::Element::quad4>(problem,numdofpernode,numscal);
  } /*
//  case DRT::Element::quad8:
//  {
//    return ScaTraImpl<DRT::Element::quad8>::Instance(problem,numdofpernode,numscal);
//  } */
  case DRT::Element::quad9:
  {
    return DefineProblemType<DRT::Element::quad9>(problem,numdofpernode,numscal);
  } /*
//  case DRT::Element::nurbs4:
//  {
//    return ScaTraImpl<DRT::Element::nurbs4>::Instance(problem,numdofpernode,numscal);
//  } */
  case DRT::Element::nurbs9:
  {
    return DefineProblemType<DRT::Element::nurbs9>(problem,numdofpernode,numscal);
  }
//  case DRT::Element::tri3:
//  {
//    return ScaTraImpl<DRT::Element::tri3>::Instance(problem,numdofpernode,numscal);
//  } /*
//  case DRT::Element::tri6:
//  {
//    return ScaTraImpl<DRT::Element::tri6>::Instance(problem,numdofpernode,numscal);
//  }*/
  case DRT::Element::line2:
  {
    return DefineProblemType<DRT::Element::line2>(problem,numdofpernode,numscal);
  }
//  case DRT::Element::line3:
//  {
//    return ScaTraImpl<DRT::Element::line3>::Instance(problem,numdofpernode,numscal);
//  }
  default:
    dserror("Element shape %s not activated. Just do it.",DRT::DistypeToString(distype).c_str());
    break;
  }
  return NULL;
}

/*--------------------------------------------------------------------------*
 |                                                       (public) nis Jan14 |
 *--------------------------------------------------------------------------*/
DRT::ELEMENTS::MeshfreeScaTraCellInterface* DRT::ELEMENTS::ScaTraFactory::ProvideMeshfreeImpl(
  DRT::Element::DiscretizationType distype,
  INPAR::SCATRA::ImplType problem,
  const int numdofpernode,
  const int numscal)
{
  switch(distype)
  {
  case DRT::Element::hex8:
  {
    return DefineMeshfreeProblemType<DRT::Element::hex8>(problem,numdofpernode,numscal);
  }
  case DRT::Element::tet4:
  {
    return DefineMeshfreeProblemType<DRT::Element::tet4>(problem,numdofpernode,numscal);
  }
  case DRT::Element::quad4:
  {
    return DefineMeshfreeProblemType<DRT::Element::quad4>(problem,numdofpernode,numscal);
  }
  case DRT::Element::tri3:
  {
    return DefineMeshfreeProblemType<DRT::Element::tri3>(problem,numdofpernode,numscal);
  }
  case DRT::Element::line2:
  {
    return DefineMeshfreeProblemType<DRT::Element::line2>(problem,numdofpernode,numscal);
  }
  default:
    dserror("Element shape %s not activated. Just do it.",DRT::DistypeToString(distype).c_str());
    break;
  }
  return NULL;
}

/*--------------------------------------------------------------------------*
 |                                               (public) ehrl        Dec13 |
 *--------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleInterface* DRT::ELEMENTS::ScaTraFactory::DefineProblemType(
  INPAR::SCATRA::ImplType problem,
  const int numdofpernode,
  const int numscal)
{
  if (problem == INPAR::SCATRA::impltype_std)
    return DRT::ELEMENTS::ScaTraEleCalcStd<distype>::Instance(numdofpernode,numscal);
  else if (problem == INPAR::SCATRA::impltype_levelset)
    return DRT::ELEMENTS::ScaTraEleCalcLS<distype>::Instance(numdofpernode,numscal);
  else if (problem == INPAR::SCATRA::impltype_lsreinit)
    return DRT::ELEMENTS::ScaTraEleCalcLsReinit<distype>::Instance(numdofpernode,numscal);
  else if (problem == INPAR::SCATRA::impltype_loma)
    return DRT::ELEMENTS::ScaTraEleCalcLoma<distype>::Instance(numdofpernode,numscal);
  else if (problem == INPAR::SCATRA::impltype_elch_NP)
    return DRT::ELEMENTS::ScaTraEleCalcElchNP<distype>::Instance(numdofpernode,numscal);
  else if (problem == INPAR::SCATRA::impltype_elch_diffcond)
    return DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::Instance(numdofpernode,numscal);
  else if (problem == INPAR::SCATRA::impltype_poro)
    return DRT::ELEMENTS::ScaTraEleCalcPoro<distype>::Instance(numdofpernode,numscal);
  else if (problem == INPAR::SCATRA::impltype_advreac)
    return DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype>::Instance(numdofpernode,numscal);
  else if (problem == INPAR::SCATRA::impltype_pororeac)
      return DRT::ELEMENTS::ScaTraEleCalcPoroReac<distype>::Instance(numdofpernode,numscal);
  else
    dserror("Defined problem type does not exist!!");

  return NULL;
}

/*--------------------------------------------------------------------------*
 |                                                       (public) nis Jan14 |
 *--------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::MeshfreeScaTraCellInterface* DRT::ELEMENTS::ScaTraFactory::DefineMeshfreeProblemType(
    INPAR::SCATRA::ImplType problem,
  const int numdofpernode,
  const int numscal)
{
  if (problem == INPAR::SCATRA::impltype_std_meshfree)
    return DRT::ELEMENTS::MeshfreeScaTraCellCalcStd<distype>::Instance(numdofpernode,numscal);
  else
    dserror("Defined problem type does not exist!!");

  return NULL;
}



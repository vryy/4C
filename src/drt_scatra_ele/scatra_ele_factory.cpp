/*--------------------------------------------------------------------------*/
/*!
\file scatra_ele_factory.cpp

\brief Factory of scatra elements

<pre>
\level 1

\maintainer Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15251
</pre>
*/
/*--------------------------------------------------------------------------*/

#include "scatra_ele_calc_utils.H"
#include "scatra_ele_calc_std.H"
#include "scatra_ele_calc_ls.H"
#include "scatra_ele_calc_lsreinit.H"
#include "scatra_ele_calc_elch_NP.H"
#include "scatra_ele_calc_elch_diffcond.H"
#include "scatra_ele_calc_elch_diffcond_multiscale.H"
#include "scatra_ele_calc_loma.H"
#include "scatra_ele_calc_poro.H"
#include "scatra_ele_calc_advanced_reaction.H"
#include "scatra_ele_calc_refconc_reac.H"
#include "scatra_ele_calc_chemo.H"
#include "scatra_ele_calc_chemo_reac.H"
#include "scatra_ele_calc_bondreac.H"
#include "scatra_ele_calc_poro_reac_ECM.H"
#include "scatra_ele_calc_poro_reac.H"
#include "scatra_ele_calc_multiporo_reac.H"
#include "scatra_ele_calc_aniso.H"
#include "scatra_ele_calc_cardiac_monodomain.H"
#include "scatra_ele_calc_elch_diffcond_sti_thermo.H"
#include "scatra_ele_calc_elch_electrode_sti_thermo.H"
#include "scatra_ele_calc_sti_diffcond.H"
#include "scatra_ele_calc_sti_electrode.H"
#include "scatra_ele_calc_hdg.H"
#include "scatra_ele_calc_variational.H"
#include "scatra_ele_parameter_std.H"

#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_element.H"

#include "scatra_ele_factory.H"

#include "scatra_ele_calc_cardiac_monodomain_hdg.H"


/*--------------------------------------------------------------------------*
 |                                               (public) ehrl        Dec13 |
 *--------------------------------------------------------------------------*/
DRT::ELEMENTS::ScaTraEleInterface* DRT::ELEMENTS::ScaTraFactory::ProvideImpl(
  DRT::Element::DiscretizationType distype,
  INPAR::SCATRA::ImplType problem,
  const int numdofpernode,
  const int numscal,
  const std::string& disname)
{
  // number of space dimensions
  const int ndim = DRT::Problem::Instance(DRT::ELEMENTS::ScaTraEleParameterStd::Instance(disname)->ProbNum())->NDim();

  switch(distype)
  {
  case DRT::Element::hex8:
  {
    if(ndim==3)
      return DefineProblemType<DRT::Element::hex8,3>(problem,numdofpernode,numscal,disname);
    else
      dserror("invalid problem dimension for HEX8 transport element!");
    break;
  }
//  case DRT::Element::hex20:
//  {
//    return ScaTraImpl<DRT::Element::hex20>::Instance(problem,numdofpernode,numscal,disname);
//  } */
  case DRT::Element::hex27:
  {
    if(ndim==3)
      return DefineProblemType<DRT::Element::hex27,3>(problem,numdofpernode,numscal,disname);
    else
      dserror("invalid problem dimension for HEX27 transport element!");
    break;
  }
//  case DRT::Element::nurbs8:
//  {
//    return ScaTraImpl<DRT::Element::nurbs8,3>::Instance(problem,numdofpernode,numscal,disname);
//  } */
//  case DRT::Element::nurbs27:
//  {
//    return ScaTraImpl<DRT::Element::nurbs27,3>::Instance(problem,numdofpernode,numscal,disname);
//  }
  case DRT::Element::tet4:
  {
    if(ndim==3)
      return DefineProblemType<DRT::Element::tet4,3>(problem,numdofpernode,numscal,disname);
    else
      dserror("invalid problem dimension for TET4 transport element!");
    break;
  }
  case DRT::Element::tet10:
  {
    if(ndim==3)
      return DefineProblemType<DRT::Element::tet10,3>(problem,numdofpernode,numscal,disname);
    else
      dserror("invalid problem dimension for TET10 transport element!");
    break;
  }
//  case DRT::Element::wedge6:
//  {
//    return ScaTraImpl<DRT::Element::wedge6,3>::Instance(problem,numdofpernode,numscal,disname);
//  } /*
//  case DRT::Element::wedge15:
//  {
//    return ScaTraImpl<DRT::Element::wedge15,3>::Instance(problem,numdofpernode,numscal,disname);
//  } */
  case DRT::Element::pyramid5:
  {
    if(ndim==3)
      return DefineProblemType<DRT::Element::pyramid5,3>(problem,numdofpernode,numscal,disname);
    else
      dserror("invalid problem dimension for PYRAMID5 transport element!");
    break;
  }
  case DRT::Element::quad4:
  {
    if(ndim==2)
      return DefineProblemType<DRT::Element::quad4,2>(problem,numdofpernode,numscal,disname);
    else if(ndim==3)
      return DefineProblemType<DRT::Element::quad4,3>(problem,numdofpernode,numscal,disname);
    else
      dserror("invalid problem dimension for quad4 transport element!");
    break;
  } /*
//  case DRT::Element::quad8:
//  {
//    return ScaTraImpl<DRT::Element::quad8,2>::Instance(problem,numdofpernode,numscal,disname);
//  } */
  case DRT::Element::quad9:
  {
    if(ndim==2)
      return DefineProblemType<DRT::Element::quad9,2>(problem,numdofpernode,numscal,disname);
//    else if(ndim==3)
//      return DefineProblemType<DRT::Element::quad9,3>(problem,numdofpernode,numscal,disname);
    else
      dserror("QUAD9 transport element not implemented as part of %i-dimensional problem. Just do it",ndim);
    break;
  } /*
//  case DRT::Element::nurbs4:
//  {
//    return ScaTraImpl<DRT::Element::nurbs4,2>::Instance(problem,numdofpernode,numscal,disname);
//  } */
  case DRT::Element::nurbs9:
  {
    if(ndim==2)
      return DefineProblemType<DRT::Element::nurbs9,2>(problem,numdofpernode,numscal,disname);
//    else if(ndim==3)
//      return DefineProblemType<DRT::Element::nurbs9,3>(problem,numdofpernode,numscal,disname);
    else
      dserror("NURBS9 transport element not implemented as part of %i-dimensional problem. Just do it",ndim);
    break;
  }
  case DRT::Element::tri3:
  {
    if(ndim==2)
      return DefineProblemType<DRT::Element::tri3,2>(problem,numdofpernode,numscal,disname);
    else if(ndim==3)
      return DefineProblemType<DRT::Element::tri3,3>(problem,numdofpernode,numscal,disname);
    else
      dserror("invalid problem dimension for tri3 transport element!");
    break;
  }
  case DRT::Element::tri6:
  {
    if(ndim==2)
      return DefineProblemType<DRT::Element::tri6,2>(problem,numdofpernode,numscal,disname);
    else
      dserror("TRI6 transport element not implemented as part of %i-dimensional problem. Just do it",ndim);
    break;
  }
  case DRT::Element::line2:
  {
    if(ndim==1)
      return DefineProblemType<DRT::Element::line2,1>(problem,numdofpernode,numscal,disname);
    else if(ndim==2)
      return DefineProblemType<DRT::Element::line2,2>(problem,numdofpernode,numscal,disname);
    else if(ndim==3)
      return DefineProblemType<DRT::Element::line2,3>(problem,numdofpernode,numscal,disname);
    else
      dserror("invalid problem dimension for LINE2 transport element!");
    break;
  }
  case DRT::Element::line3:
  {
    if(ndim==1)
      return DefineProblemType<DRT::Element::line3,1>(problem,numdofpernode,numscal,disname);
//    else if(ndim==2)
//      return DefineProblemType<DRT::Element::line2,2>(problem,numdofpernode,numscal,disname);
//    else if(ndim==3)
//      return DefineProblemType<DRT::Element::line2,3>(problem,numdofpernode,numscal,disname);
    else
      dserror("LINE3 transport element not implemented as part of %i-dimensional problem. Just do it",ndim);
    break;
  }
  default:
    dserror("Element shape %s not activated. Just do it.",DRT::DistypeToString(distype).c_str());
    break;
  }
  return NULL;
}

/*--------------------------------------------------------------------------*
 |                                                (public) hoermann   11/15 |
 *--------------------------------------------------------------------------*/
DRT::ELEMENTS::ScaTraEleInterface* DRT::ELEMENTS::ScaTraFactory::ProvideImplHDG(
  DRT::Element::DiscretizationType distype,
  INPAR::SCATRA::ImplType problem,
  const int numdofpernode,
  const int numscal,
  const std::string& disname)
{
  // -------------------------------------- number of degrees of freedom
  // number of degrees of freedom
  static const int ndim = DRT::Problem::Instance()->NDim();

  switch(distype)
  {
  case DRT::Element::hex8:
  {
    if(ndim==3)
      return DefineProblemTypeHDG<DRT::Element::hex8,3>(problem,numdofpernode,numscal,disname);
    else
      dserror("invalid problem dimension for HEX8 transport element!");
    break;
  }
//  case DRT::Element::hex20:
//  {
//    return ScaTraImpl<DRT::Element::hex20>::Instance(problem,numdofpernode,numscal,disname);
//  } */
//  case DRT::Element::hex27:
//  {
//    if(ndim==3)
//      return DefineProblemTypeHDG<DRT::Element::hex27,3>(problem,numdofpernode,numscal,disname);
//    else
//      dserror("invalid problem dimension for HEX27 transport element!");
//    break;
//  }
//  case DRT::Element::nurbs8:
//  {
//    return ScaTraImpl<DRT::Element::nurbs8,3>::Instance(problem,numdofpernode,numscal,disname);
//  } */
//  case DRT::Element::nurbs27:
//  {
//    return ScaTraImpl<DRT::Element::nurbs27,3>::Instance(problem,numdofpernode,numscal,disname);
//  }
  case DRT::Element::tet4:
  {
    if(ndim==3)
      return DefineProblemTypeHDG<DRT::Element::tet4,3>(problem,numdofpernode,numscal,disname);
    else
      dserror("invalid problem dimension for TET4 transport element!");
    break;
  }
//  case DRT::Element::tet10:
//  {
//    if(ndim==3)
//      return DefineProblemTypeHDG<DRT::Element::tet10,3>(problem,numdofpernode,numscal,disname);
//    else
//      dserror("invalid problem dimension for TET10 transport element!");
//    break;
//  }
//  case DRT::Element::wedge6:
//  {
//    return ScaTraImpl<DRT::Element::wedge6,3>::Instance(problem,numdofpernode,numscal,disname);
//  } /*
//  case DRT::Element::wedge15:
//  {
//    return ScaTraImpl<DRT::Element::wedge15,3>::Instance(problem,numdofpernode,numscal,disname);
//  } */
  case DRT::Element::pyramid5:
  {
    if(ndim==3)
      return DefineProblemTypeHDG<DRT::Element::pyramid5,3>(problem,numdofpernode,numscal,disname);
    else
      dserror("invalid problem dimension for PYRAMID5 transport element!");
    break;
  }
  case DRT::Element::quad4:
  {
    if(ndim==2)
      return DefineProblemTypeHDG<DRT::Element::quad4,2>(problem,numdofpernode,numscal,disname);
    else if(ndim==3)
      return DefineProblemTypeHDG<DRT::Element::quad4,3>(problem,numdofpernode,numscal,disname);
    else
      dserror("invalid problem dimension for quad4 transport element!");
    break;
  } /*
//  case DRT::Element::quad8:
//  {
//    return ScaTraImpl<DRT::Element::quad8,2>::Instance(problem,numdofpernode,numscal,disname);
//  } */
//  case DRT::Element::quad9:
//  {
//    if(ndim==2)
//      return DefineProblemTypeHDG<DRT::Element::quad9,2>(problem,numdofpernode,numscal,disname);
////    else if(ndim==3)
////      return DefineProblemType<DRT::Element::quad9,3>(problem,numdofpernode,numscal,disname);
//    else
//      dserror("QUAD9 transport element not implemented as part of %i-dimensional problem. Just do it",ndim);
//    break;
//  } /*
//  case DRT::Element::nurbs4:
//  {
//    return ScaTraImpl<DRT::Element::nurbs4,2>::Instance(problem,numdofpernode,numscal,disname);
//  } */
//  case DRT::Element::nurbs9:
//  {
//    if(ndim==2)
//      return DefineProblemTypeHDG<DRT::Element::nurbs9,2>(problem,numdofpernode,numscal,disname);
////    else if(ndim==3)
////      return DefineProblemType<DRT::Element::nurbs9,3>(problem,numdofpernode,numscal,disname);
//    else
//      dserror("NURBS9 transport element not implemented as part of %i-dimensional problem. Just do it",ndim);
//    break;
//  }
  case DRT::Element::tri3:
  {
    return DefineProblemTypeHDG<DRT::Element::tri3,2>(problem,numdofpernode,numscal,disname);
  } /*
  case DRT::Element::tri6:
  {
    return ScaTraImpl<DRT::Element::tri6,2>::Instance(problem,numdofpernode,numscal,disname);
  }*/
//  case DRT::Element::line2:
//  {
//    if(ndim==1)
//      return DefineProblemType<DRT::Element::line2,1>(problem,numdofpernode,numscal,disname);
//    else if(ndim==2)
//      return DefineProblemType<DRT::Element::line2,2>(problem,numdofpernode,numscal,disname);
//    else if(ndim==3)
//      return DefineProblemType<DRT::Element::line2,3>(problem,numdofpernode,numscal,disname);
//    else
//      dserror("invalid problem dimension for LINE2 transport element!");
//    break;
//  }
//  case DRT::Element::line3:
//  {
//    if(ndim==1)
//      return DefineProblemType<DRT::Element::line3,1>(problem,numdofpernode,numscal,disname);
////    else if(ndim==2)
////      return DefineProblemType<DRT::Element::line2,2>(problem,numdofpernode,numscal,disname);
////    else if(ndim==3)
////      return DefineProblemType<DRT::Element::line2,3>(problem,numdofpernode,numscal,disname);
//    else
//      dserror("LINE3 transport element not implemented as part of %i-dimensional problem. Just do it",ndim);
//    break;
//  }
  default:
    dserror("Element shape %s not activated. Just do it.",DRT::DistypeToString(distype).c_str());
    break;
  }
  return NULL;
}


/*--------------------------------------------------------------------------*
 |                                               (public) ehrl        Dec13 |
 *--------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, int probdim>
DRT::ELEMENTS::ScaTraEleInterface* DRT::ELEMENTS::ScaTraFactory::DefineProblemType(
  INPAR::SCATRA::ImplType problem,
  const int numdofpernode,
  const int numscal,
  const std::string& disname)
{
  if(DRT::UTILS::DisTypeToDim<distype>::dim != probdim)
    if( problem != INPAR::SCATRA::impltype_std and
        problem != INPAR::SCATRA::impltype_cardiac_monodomain and
        problem != INPAR::SCATRA::impltype_advreac and
        problem != INPAR::SCATRA::impltype_lsreinit and
        problem != INPAR::SCATRA::impltype_bondreac)
      dserror("ImplType '%s' not implemented for transport on manifolds!",SCATRA::ImplTypeToString(problem).c_str());

  switch(problem)
  {
  case INPAR::SCATRA::impltype_std:
  {
    return DRT::ELEMENTS::ScaTraEleCalcStd<distype,probdim>::Instance(numdofpernode,numscal,disname);
    break;
  }
  case INPAR::SCATRA::impltype_thermo_elch_electrode:
  {
    return DRT::ELEMENTS::ScaTraEleCalcSTIElectrode<distype>::Instance(numdofpernode,numscal,disname);
    break;
  }
  case INPAR::SCATRA::impltype_thermo_elch_diffcond:
  {
    return DRT::ELEMENTS::ScaTraEleCalcSTIDiffCond<distype>::Instance(numdofpernode,numscal,disname);
    break;
  }
  case INPAR::SCATRA::impltype_levelset:
  {
    return DRT::ELEMENTS::ScaTraEleCalcLS<distype>::Instance(numdofpernode,numscal,disname);
    break;
  }
  case INPAR::SCATRA::impltype_lsreinit:
  {
    return DRT::ELEMENTS::ScaTraEleCalcLsReinit<distype,probdim>::Instance(numdofpernode,numscal,disname);
    break;
  }
  case INPAR::SCATRA::impltype_loma:
  {
    return DRT::ELEMENTS::ScaTraEleCalcLoma<distype>::Instance(numdofpernode,numscal,disname);
    break;
  }
  case INPAR::SCATRA::impltype_elch_NP:
  {
    return DRT::ELEMENTS::ScaTraEleCalcElchNP<distype>::Instance(numdofpernode,numscal,disname);
    break;
  }
  case INPAR::SCATRA::impltype_elch_electrode:
  case INPAR::SCATRA::impltype_elch_electrode_growth:
  {
    return DRT::ELEMENTS::ScaTraEleCalcElchElectrode<distype>::Instance(numdofpernode,numscal,disname);
    break;
  }
  case INPAR::SCATRA::impltype_elch_electrode_thermo:
  {
    return DRT::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<distype>::Instance(numdofpernode,numscal,disname);
    break;
  }
  case INPAR::SCATRA::impltype_variational_diffusion:
  {
    return DRT::ELEMENTS::ScaTraEleCalVariational<distype>::Instance(numdofpernode,numscal,disname);
    break;
  }
  case INPAR::SCATRA::impltype_elch_diffcond:
  {
    return DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::Instance(numdofpernode,numscal,disname);
    break;
  }
  case INPAR::SCATRA::impltype_elch_diffcond_multiscale:
  {
    return DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<distype>::Instance(numdofpernode,numscal,disname);
    break;
  }
  case INPAR::SCATRA::impltype_elch_diffcond_thermo:
  {
    return DRT::ELEMENTS::ScaTraEleCalcElchDiffCondSTIThermo<distype>::Instance(numdofpernode,numscal,disname);
    break;
  }
  case INPAR::SCATRA::impltype_poro:
  {
    return DRT::ELEMENTS::ScaTraEleCalcPoro<distype>::Instance(numdofpernode,numscal,disname);
    break;
  }
  case INPAR::SCATRA::impltype_advreac:
  {
    return DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype,probdim>::Instance(numdofpernode,numscal,disname);
    break;
  }
  case INPAR::SCATRA::impltype_refconcreac:
  {
    return DRT::ELEMENTS::ScaTraEleCalcRefConcReac<distype>::Instance(numdofpernode,numscal,disname);
    break;
  }
  case INPAR::SCATRA::impltype_chemo:
  {
    return DRT::ELEMENTS::ScaTraEleCalcChemo<distype>::Instance(numdofpernode,numscal,disname);
    break;
  }
  case INPAR::SCATRA::impltype_chemoreac:
  {
    return DRT::ELEMENTS::ScaTraEleCalcChemoReac<distype>::Instance(numdofpernode,numscal,disname);
    break;
  }
  case INPAR::SCATRA::impltype_multipororeac:
  {
    return DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>::Instance(numdofpernode,numscal,disname);
    break;
  }
  case INPAR::SCATRA::impltype_pororeac:
  {
    return DRT::ELEMENTS::ScaTraEleCalcPoroReac<distype>::Instance(numdofpernode,numscal,disname);
    break;
  }
  case INPAR::SCATRA::impltype_pororeacECM:
  {
    return DRT::ELEMENTS::ScaTraEleCalcPoroReacECM<distype>::Instance(numdofpernode,numscal,disname);
    break;
  }
  case INPAR::SCATRA::impltype_aniso:
  {
    return DRT::ELEMENTS::ScaTraEleCalcAniso<distype,probdim>::Instance(numdofpernode,numscal,disname);
    break;
  }
  case INPAR::SCATRA::impltype_cardiac_monodomain:
  {
    return DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<distype,probdim>::Instance(numdofpernode,numscal,disname);
    break;
  }
  case INPAR::SCATRA::impltype_bondreac:
  {
    return DRT::ELEMENTS::ScaTraEleCalcBondReac<distype,probdim>::Instance(numdofpernode,numscal,disname);
        break;
  }
  default:
  {
    dserror("Defined problem type does not exist!!");
    break;
  }
  }

  return NULL;
}

/*--------------------------------------------------------------------------*
 |                                               (public) hoermann    09/15 |
 *--------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, int probdim>
DRT::ELEMENTS::ScaTraEleInterface* DRT::ELEMENTS::ScaTraFactory::DefineProblemTypeHDG(
  INPAR::SCATRA::ImplType problem,
  const int numdofpernode,
  const int numscal,
  const std::string& disname)
{
  if(DRT::UTILS::DisTypeToDim<distype>::dim != probdim)
    if( problem != INPAR::SCATRA::impltype_std and
        problem != INPAR::SCATRA::impltype_cardiac_monodomain)
      dserror("ImplType '%s' not implemented for transport on manifolds!",SCATRA::ImplTypeToString(problem).c_str());

  switch(problem)
  {
  case INPAR::SCATRA::impltype_std_hdg:
  {
    return DRT::ELEMENTS::ScaTraEleCalcHDG<distype,probdim>::Instance(numdofpernode,numscal,disname);
    break;
  }
  case INPAR::SCATRA::impltype_cardiac_monodomain_hdg:
  {
    return DRT::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<distype,probdim>::Instance(numdofpernode,numscal,disname);
    break;
  }
  default:
  {
    dserror("Defined problem type does not exist!!");
    break;
  }
  }

  return NULL;
}


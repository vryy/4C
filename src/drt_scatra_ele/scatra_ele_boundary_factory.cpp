/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_boundary_factory.cpp

\brief factory for scatra boundary evaluation

\level 2

<pre>
\maintainer Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15251
</pre>
 */
/*----------------------------------------------------------------------*/
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/elchmat.H"

#include "scatra_ele.H"
#include "scatra_ele_boundary_calc_elch_diffcond.H"
#include "scatra_ele_boundary_calc_elch_electrode.H"
#include "scatra_ele_boundary_calc_elch_electrode_growth.H"
#include "scatra_ele_boundary_calc_elch_electrode_sti_thermo.H"
#include "scatra_ele_boundary_calc_elch_NP.H"
#include "scatra_ele_boundary_calc_loma.H"
#include "scatra_ele_boundary_calc_poro.H"
#include "scatra_ele_boundary_calc_refconc_reac.H"
#include "scatra_ele_boundary_calc_std.H"
#include "scatra_ele_boundary_calc_sti_electrode.H"
#include "scatra_ele_boundary_interface.H"
#include "scatra_ele_calc.H"
#include "scatra_ele_boundary_factory.H"


/*--------------------------------------------------------------------------*
 |                                                 (public) rasthofer 11/13 |
 *--------------------------------------------------------------------------*/
DRT::ELEMENTS::ScaTraBoundaryInterface* DRT::ELEMENTS::ScaTraBoundaryFactory::ProvideImpl(const DRT::Element* ele,
                                                                                          const enum INPAR::SCATRA::ImplType impltype,
                                                                                          const int numdofpernode,
                                                                                          const int numscal,
                                                                                          const std::string& disname)
{
  switch(ele->Shape())
  {
  case DRT::Element::quad4:
  {
    return DefineProblemType<DRT::Element::quad4>(impltype,numdofpernode,numscal,disname);
  }
  case DRT::Element::quad8:
  {
    return DefineProblemType<DRT::Element::quad8>(impltype,numdofpernode,numscal,disname);
  }
  case DRT::Element::quad9:
  {
    return DefineProblemType<DRT::Element::quad9>(impltype,numdofpernode,numscal,disname);
  }
  case DRT::Element::tri3:
  {
    return DefineProblemType<DRT::Element::tri3>(impltype,numdofpernode,numscal,disname);
  }
  case DRT::Element::tri6:
  {
    return DefineProblemType<DRT::Element::tri6>(impltype,numdofpernode,numscal,disname);
  }
  case DRT::Element::line2:
  {
    return DefineProblemType<DRT::Element::line2>(impltype,numdofpernode,numscal,disname);
  }
  case DRT::Element::line3:
  {
    return DefineProblemType<DRT::Element::line3>(impltype,numdofpernode,numscal,disname);
  }
  case DRT::Element::nurbs3:    // 1D nurbs boundary element
  {
    return DefineProblemType<DRT::Element::nurbs3>(impltype,numdofpernode,numscal,disname);
  }
  case DRT::Element::nurbs9:    // 2D nurbs boundary element
  {
    return DefineProblemType<DRT::Element::nurbs9>(impltype,numdofpernode,numscal,disname);
  }
  default:
  {
    dserror("Element shape %d (%d nodes) not activated. Just do it.", ele->Shape(), ele->NumNode());
    break;
  }
  }

  return NULL;
}


/*-------------------------------------------------------------------------------------------*
 | return instance of element evaluation class depending on implementation type   fang 02/15 |
 *-------------------------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraBoundaryInterface* DRT::ELEMENTS::ScaTraBoundaryFactory::DefineProblemType(
    const enum INPAR::SCATRA::ImplType impltype,
    const int numdofpernode,
    const int numscal,
    const std::string& disname
    )
{
  switch(impltype)
  {
  case INPAR::SCATRA::impltype_advreac:
  case INPAR::SCATRA::impltype_aniso:
  case INPAR::SCATRA::impltype_cardiac_monodomain:
  case INPAR::SCATRA::impltype_chemo:
  case INPAR::SCATRA::impltype_chemoreac:
  case INPAR::SCATRA::impltype_levelset:
  case INPAR::SCATRA::impltype_std:
  case INPAR::SCATRA::impltype_thermo_elch_diffcond:
  case INPAR::SCATRA::impltype_bondreac:
  {
    return DRT::ELEMENTS::ScaTraEleBoundaryCalcStd<distype>::Instance(numdofpernode,numscal,disname);
    break;
  }
  case INPAR::SCATRA::impltype_loma:
  {
    return DRT::ELEMENTS::ScaTraEleBoundaryCalcLoma<distype>::Instance(numdofpernode,numscal,disname);
    break;
  }
  case INPAR::SCATRA::impltype_elch_electrode:
  {
    return DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<distype>::Instance(numdofpernode,numscal,disname);
    break;
  }
  case INPAR::SCATRA::impltype_elch_electrode_thermo:
  {
    return DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<distype>::Instance(numdofpernode,numscal,disname);
    break;
  }
  case INPAR::SCATRA::impltype_elch_diffcond:
  case INPAR::SCATRA::impltype_elch_diffcond_thermo:
  {
    return DRT::ELEMENTS::ScaTraEleBoundaryCalcElchDiffCond<distype>::Instance(numdofpernode,numscal,disname);
    break;
  }
  case INPAR::SCATRA::impltype_elch_NP:
  {
    return DRT::ELEMENTS::ScaTraEleBoundaryCalcElchNP<distype>::Instance(numdofpernode,numscal,disname);
    break;
  }
  case INPAR::SCATRA::impltype_poro:
  case INPAR::SCATRA::impltype_pororeac:
  {
    return DRT::ELEMENTS::ScaTraEleBoundaryCalcPoro<distype>::Instance(numdofpernode,numscal,disname);
    break;
  }
  case INPAR::SCATRA::impltype_refconcreac:
  {
    return DRT::ELEMENTS::ScaTraEleBoundaryCalcRefConcReac<distype>::Instance(numdofpernode,numscal,disname);
    break;
  }
  case INPAR::SCATRA::impltype_thermo_elch_electrode:
  {
    return DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<distype>::Instance(numdofpernode,numscal,disname);
    break;
  }
  case INPAR::SCATRA::impltype_elch_electrode_growth:
  {
    return DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<distype>::Instance(numdofpernode,numscal,disname);
    break;
  }
  default:
  {
    dserror("Defined implementation type does not exist!");
    break;
  }
  } // switch(impltype)

  return NULL;
}

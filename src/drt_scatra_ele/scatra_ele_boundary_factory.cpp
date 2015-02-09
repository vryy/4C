/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_boundary_factory.cpp

\brief factory for fluid scatra boundary evaluation

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
 */
/*----------------------------------------------------------------------*/
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/elchmat.H"

#include "scatra_ele.H"
#include "scatra_ele_boundary_calc_elch_diffcond.H"
#include "scatra_ele_boundary_calc_elch_electrode.H"
#include "scatra_ele_boundary_calc_elch_NP.H"
#include "scatra_ele_boundary_calc_poro.H"
#include "scatra_ele_boundary_calc_std.H"
#include "scatra_ele_boundary_interface.H"
#include "scatra_ele_calc.H"
#include "scatra_ele_boundary_factory.H"


/*--------------------------------------------------------------------------*
 |                                                 (public) rasthofer 11/13 |
 *--------------------------------------------------------------------------*/
DRT::ELEMENTS::ScaTraBoundaryInterface* DRT::ELEMENTS::ScaTraBoundaryFactory::ProvideImpl(const DRT::Element* ele,
                                                                                          const enum INPAR::SCATRA::ImplType impltype,
                                                                                          const int numdofpernode,
                                                                                          const int numscal)
{
  switch(ele->Shape())
  {
  case DRT::Element::quad4:
  {
    return DefineProblemType<DRT::Element::quad4>(impltype,numdofpernode,numscal);
  }
  case DRT::Element::quad8:
  {
    return DefineProblemType<DRT::Element::quad8>(impltype,numdofpernode,numscal);
  }
  case DRT::Element::quad9:
  {
    return DefineProblemType<DRT::Element::quad9>(impltype,numdofpernode,numscal);
  }
  case DRT::Element::tri3:
  {
    return DefineProblemType<DRT::Element::tri3>(impltype,numdofpernode,numscal);
  }
  case DRT::Element::tri6:
  {
    return DefineProblemType<DRT::Element::tri6>(impltype,numdofpernode,numscal);
  }
  case DRT::Element::line2:
  {
    return DefineProblemType<DRT::Element::line2>(impltype,numdofpernode,numscal);
  }
  case DRT::Element::line3:
  {
    return DefineProblemType<DRT::Element::line3>(impltype,numdofpernode,numscal);
  }
  case DRT::Element::nurbs3:    // 1D nurbs boundary element
  {
    return DefineProblemType<DRT::Element::nurbs3>(impltype,numdofpernode,numscal);
  }
  case DRT::Element::nurbs9:    // 2D nurbs boundary element
  {
    return DefineProblemType<DRT::Element::nurbs9>(impltype,numdofpernode,numscal);
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
    const int numscal
    )
{
  switch(impltype)
  {
  case INPAR::SCATRA::impltype_std:
  case INPAR::SCATRA::impltype_advreac:
  case INPAR::SCATRA::impltype_aniso:
  case INPAR::SCATRA::impltype_cardiac_monodomain:
  case INPAR::SCATRA::impltype_levelset:
  case INPAR::SCATRA::impltype_loma:
  {
    return DRT::ELEMENTS::ScaTraEleBoundaryCalcStd<distype>::Instance(numdofpernode,numscal);
    break;
  }
  case INPAR::SCATRA::impltype_elch_electrode:
  {
    return DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<distype>::Instance(numdofpernode,numscal);
    break;
  }
  case INPAR::SCATRA::impltype_elch_diffcond:
  {
    return DRT::ELEMENTS::ScaTraEleBoundaryCalcElchDiffCond<distype>::Instance(numdofpernode,numscal);
    break;
  }
  case INPAR::SCATRA::impltype_elch_NP:
  {
    return DRT::ELEMENTS::ScaTraEleBoundaryCalcElchNP<distype>::Instance(numdofpernode,numscal);
    break;
  }
  case INPAR::SCATRA::impltype_poro:
  case INPAR::SCATRA::impltype_pororeac:
    {
    return DRT::ELEMENTS::ScaTraEleBoundaryCalcPoro<distype>::Instance(numdofpernode,numscal);
    break;
    }
  default:
  {
    dserror("Defined problem type does not exist!!");
    break;
  }
  } // switch(impltype)

  return NULL;
}

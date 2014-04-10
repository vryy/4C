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

#include "scatra_ele_boundary_factory.H"
#include "scatra_ele_boundary_interface.H"
#include "scatra_ele_boundary_calc_std.H"
#include "scatra_ele_boundary_calc_poro.H"
#include "scatra_ele_boundary_calc_elch.H"
#include "scatra_ele_calc.H"
#include "scatra_ele.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/elchmat.H"

/*--------------------------------------------------------------------------*
 |                                                 (public) rasthofer 11/13 |
 *--------------------------------------------------------------------------*/
DRT::ELEMENTS::ScaTraBoundaryInterface* DRT::ELEMENTS::ScaTraBoundaryFactory::ProvideImpl(const DRT::Element* ele,
                                                                                          const enum INPAR::SCATRA::ScaTraType scatratype,
                                                                                          const int numdofpernode,
                                                                                          const int numscal)
{
  switch (ele->Shape())
  {
  case DRT::Element::quad4:
  {
    return DefineProblemType<DRT::Element::quad4>(scatratype,numdofpernode,numscal);
  }
  case DRT::Element::quad8:
  {
    return DefineProblemType<DRT::Element::quad8>(scatratype,numdofpernode,numscal);
  }
  case DRT::Element::quad9:
  {
    return DefineProblemType<DRT::Element::quad9>(scatratype,numdofpernode,numscal);
  }
  case DRT::Element::tri3:
  {
    return DefineProblemType<DRT::Element::tri3>(scatratype,numdofpernode,numscal);
  }
  case DRT::Element::tri6:
  {
    return DefineProblemType<DRT::Element::tri6>(scatratype,numdofpernode,numscal);
  }
  case DRT::Element::line2:
  {
    return DefineProblemType<DRT::Element::line2>(scatratype,numdofpernode,numscal);
  }
  case DRT::Element::line3:
  {
    return DefineProblemType<DRT::Element::line3>(scatratype,numdofpernode,numscal);
  }
  case DRT::Element::nurbs3:    // 1D nurbs boundary element
  {
    return DefineProblemType<DRT::Element::nurbs3>(scatratype,numdofpernode,numscal);
  }
  case DRT::Element::nurbs9:    // 2D nurbs boundary element
  {
    return DefineProblemType<DRT::Element::nurbs9>(scatratype,numdofpernode,numscal);
  }
  default:
  {
    dserror("Element shape %d (%d nodes) not activated. Just do it.", ele->Shape(), ele->NumNode());
    break;
  }
  break;
  }
  return NULL;
}


/*--------------------------------------------------------------------------*
 |                                                 (public) rasthofer 11/13 |
 *--------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraBoundaryInterface* DRT::ELEMENTS::ScaTraBoundaryFactory::DefineProblemType(const enum INPAR::SCATRA::ScaTraType scatratype,
                                                                                                const int numdofpernode,
                                                                                                const int numscal
                                                                                                )
{
  switch (scatratype)
  {
  case INPAR::SCATRA::scatratype_undefined:
  case INPAR::SCATRA::scatratype_condif:
  case INPAR::SCATRA::scatratype_loma:
  case INPAR::SCATRA::scatratype_levelset:
  case INPAR::SCATRA::scatratype_advreac:
  case INPAR::SCATRA::scatratype_cardio_monodomain:
    {
      return DRT::ELEMENTS::ScaTraEleBoundaryCalcStd<distype>::Instance(numdofpernode,numscal);
      break;
    }
  case INPAR::SCATRA::scatratype_elch:
  {
    return DRT::ELEMENTS::ScaTraEleBoundaryCalcElch<distype>::Instance(numdofpernode,numscal);
    break;
  }
  case INPAR::SCATRA::scatratype_pororeac:
  case INPAR::SCATRA::scatratype_poro:
    {
    return DRT::ELEMENTS::ScaTraEleBoundaryCalcPoro<distype>::Instance(numdofpernode,numscal);
    break;
    }
  default:
  {
    dserror("Defined problem type does not exist!!");
    break;
  }
  break;
  } // switch scatratype

  return NULL;
}


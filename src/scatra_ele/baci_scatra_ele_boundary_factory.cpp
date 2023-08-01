/*----------------------------------------------------------------------*/
/*! \file

\brief factory for scatra boundary evaluation

\level 2

 */
/*----------------------------------------------------------------------*/
#include "baci_lib_globalproblem.H"
#include "baci_mat_elchmat.H"

#include "baci_scatra_ele.H"
#include "baci_scatra_ele_boundary_calc_elch_diffcond.H"
#include "baci_scatra_ele_boundary_calc_elch_electrode.H"
#include "baci_scatra_ele_boundary_calc_elch_electrode_growth.H"
#include "baci_scatra_ele_boundary_calc_elch_electrode_sti_thermo.H"
#include "baci_scatra_ele_boundary_calc_elch_NP.H"
#include "baci_scatra_ele_boundary_calc_loma.H"
#include "baci_scatra_ele_boundary_calc_poro.H"
#include "baci_scatra_ele_boundary_calc_refconc_reac.H"
#include "baci_scatra_ele_boundary_calc_std.H"
#include "baci_scatra_ele_boundary_calc_sti_electrode.H"
#include "baci_scatra_ele_boundary_interface.H"
#include "baci_scatra_ele_calc.H"
#include "baci_scatra_ele_boundary_factory.H"


/*--------------------------------------------------------------------------*
 |                                                 (public) rasthofer 11/13 |
 *--------------------------------------------------------------------------*/
DRT::ELEMENTS::ScaTraBoundaryInterface* DRT::ELEMENTS::ScaTraBoundaryFactory::ProvideImpl(
    const DRT::Element* ele, const enum INPAR::SCATRA::ImplType impltype, const int numdofpernode,
    const int numscal, const std::string& disname)
{
  // number of space dimensions
  const int ndim = disname != "scatra_micro" ? DRT::Problem::Instance()->NDim() : 1;

  switch (ele->Shape())
  {
    case DRT::Element::quad4:
    {
      if (ndim == 3)
        return DefineProblemType<DRT::Element::quad4, 3>(impltype, numdofpernode, numscal, disname);
      else
        dserror("invalid problem dimension!");
    }
    case DRT::Element::quad8:
    {
      if (ndim == 3)
        return DefineProblemType<DRT::Element::quad8, 3>(impltype, numdofpernode, numscal, disname);
      else
        dserror("invalid problem dimension!");
    }
    case DRT::Element::quad9:
    {
      if (ndim == 3)
        return DefineProblemType<DRT::Element::quad9, 3>(impltype, numdofpernode, numscal, disname);
      else
        dserror("invalid problem dimension!");
    }
    case DRT::Element::tri3:
    {
      if (ndim == 3)
        return DefineProblemType<DRT::Element::tri3, 3>(impltype, numdofpernode, numscal, disname);
      else
        dserror("invalid problem dimension!");
    }
    case DRT::Element::tri6:
    {
      if (ndim == 3)
        return DefineProblemType<DRT::Element::tri6, 3>(impltype, numdofpernode, numscal, disname);
      else
        dserror("invalid problem dimension!");
    }
    case DRT::Element::line2:
    {
      if (ndim == 2)
        return DefineProblemType<DRT::Element::line2, 2>(impltype, numdofpernode, numscal, disname);
      else if (ndim == 3)
        return DefineProblemType<DRT::Element::line2, 3>(impltype, numdofpernode, numscal, disname);
      else
        dserror("invalid problem dimension!");
    }
    case DRT::Element::line3:
    {
      if (ndim == 2)
        return DefineProblemType<DRT::Element::line3, 2>(impltype, numdofpernode, numscal, disname);
      else
        dserror("invalid problem dimension!");
    }
    case DRT::Element::nurbs3:  // 1D nurbs boundary element
    {
      return DefineProblemType<DRT::Element::nurbs3, 2>(impltype, numdofpernode, numscal, disname);
    }
    case DRT::Element::nurbs9:  // 2D nurbs boundary element
    {
      return DefineProblemType<DRT::Element::nurbs9, 3>(impltype, numdofpernode, numscal, disname);
    }
    default:
    {
      dserror(
          "Element shape %d (%d nodes) not activated. Just do it.", ele->Shape(), ele->NumNode());
      break;
    }
  }

  return NULL;
}


/*-------------------------------------------------------------------------------------------*
 | return instance of element evaluation class depending on implementation type   fang 02/15 |
 *-------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
DRT::ELEMENTS::ScaTraBoundaryInterface* DRT::ELEMENTS::ScaTraBoundaryFactory::DefineProblemType(
    const enum INPAR::SCATRA::ImplType impltype, const int numdofpernode, const int numscal,
    const std::string& disname)
{
  switch (impltype)
  {
    case INPAR::SCATRA::impltype_advreac:
    case INPAR::SCATRA::impltype_aniso:
    case INPAR::SCATRA::impltype_cardiac_monodomain:
    case INPAR::SCATRA::impltype_chemo:
    case INPAR::SCATRA::impltype_chemoreac:
    case INPAR::SCATRA::impltype_levelset:
    case INPAR::SCATRA::impltype_std:
    case INPAR::SCATRA::impltype_thermo_elch_diffcond:
    case INPAR::SCATRA::impltype_multipororeac:
    {
      return DRT::ELEMENTS::ScaTraEleBoundaryCalcStd<distype, probdim>::Instance(
          numdofpernode, numscal, disname);
      break;
    }
    case INPAR::SCATRA::impltype_loma:
    {
      return DRT::ELEMENTS::ScaTraEleBoundaryCalcLoma<distype, probdim>::Instance(
          numdofpernode, numscal, disname);
      break;
    }
    case INPAR::SCATRA::impltype_elch_electrode:
    {
      return DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<distype, probdim>::Instance(
          numdofpernode, numscal, disname);
      break;
    }
    case INPAR::SCATRA::impltype_elch_electrode_thermo:
    {
      return DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<distype, probdim>::Instance(
          numdofpernode, numscal, disname);
      break;
    }
    case INPAR::SCATRA::impltype_elch_diffcond:
    case INPAR::SCATRA::impltype_elch_diffcond_thermo:
    {
      return DRT::ELEMENTS::ScaTraEleBoundaryCalcElchDiffCond<distype, probdim>::Instance(
          numdofpernode, numscal, disname);
      break;
    }
    case INPAR::SCATRA::impltype_elch_NP:
    {
      return DRT::ELEMENTS::ScaTraEleBoundaryCalcElchNP<distype, probdim>::Instance(
          numdofpernode, numscal, disname);
      break;
    }
    case INPAR::SCATRA::impltype_poro:
    case INPAR::SCATRA::impltype_pororeac:
    {
      return DRT::ELEMENTS::ScaTraEleBoundaryCalcPoro<distype, probdim>::Instance(
          numdofpernode, numscal, disname);
      break;
    }
    case INPAR::SCATRA::impltype_refconcreac:
    {
      return DRT::ELEMENTS::ScaTraEleBoundaryCalcRefConcReac<distype, probdim>::Instance(
          numdofpernode, numscal, disname);
      break;
    }
    case INPAR::SCATRA::impltype_thermo_elch_electrode:
    {
      return DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<distype, probdim>::Instance(
          numdofpernode, numscal, disname);
      break;
    }
    case INPAR::SCATRA::impltype_elch_electrode_growth:
    {
      return DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<distype, probdim>::Instance(
          numdofpernode, numscal, disname);
      break;
    }
    default:
    {
      dserror("Defined implementation type does not exist!");
      break;
    }
  }  // switch(impltype)

  return NULL;
}

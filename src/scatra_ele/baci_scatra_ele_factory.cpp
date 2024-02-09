/*--------------------------------------------------------------------------*/
/*! \file

\brief Factory of scatra elements

\level 1

*/
/*--------------------------------------------------------------------------*/

#include "baci_scatra_ele_factory.hpp"

#include "baci_global_data.hpp"
#include "baci_scatra_ele_calc_advanced_reaction.hpp"
#include "baci_scatra_ele_calc_aniso.hpp"
#include "baci_scatra_ele_calc_artery.hpp"
#include "baci_scatra_ele_calc_cardiac_monodomain.hpp"
#include "baci_scatra_ele_calc_cardiac_monodomain_hdg.hpp"
#include "baci_scatra_ele_calc_chemo.hpp"
#include "baci_scatra_ele_calc_chemo_reac.hpp"
#include "baci_scatra_ele_calc_elch_diffcond.hpp"
#include "baci_scatra_ele_calc_elch_diffcond_multiscale.hpp"
#include "baci_scatra_ele_calc_elch_diffcond_sti_thermo.hpp"
#include "baci_scatra_ele_calc_elch_electrode.hpp"
#include "baci_scatra_ele_calc_elch_electrode_sti_thermo.hpp"
#include "baci_scatra_ele_calc_elch_NP.hpp"
#include "baci_scatra_ele_calc_elch_scl.hpp"
#include "baci_scatra_ele_calc_hdg.hpp"
#include "baci_scatra_ele_calc_loma.hpp"
#include "baci_scatra_ele_calc_ls.hpp"
#include "baci_scatra_ele_calc_lsreinit.hpp"
#include "baci_scatra_ele_calc_multiporo_reac.hpp"
#include "baci_scatra_ele_calc_no_physics.hpp"
#include "baci_scatra_ele_calc_poro.hpp"
#include "baci_scatra_ele_calc_poro_reac.hpp"
#include "baci_scatra_ele_calc_poro_reac_ECM.hpp"
#include "baci_scatra_ele_calc_refconc_reac.hpp"
#include "baci_scatra_ele_calc_std.hpp"
#include "baci_scatra_ele_calc_sti_diffcond.hpp"
#include "baci_scatra_ele_calc_sti_electrode.hpp"
#include "baci_scatra_ele_calc_utils.hpp"
#include "baci_scatra_ele_parameter_std.hpp"

BACI_NAMESPACE_OPEN


/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
DRT::ELEMENTS::ScaTraEleInterface* DRT::ELEMENTS::ScaTraFactory::ProvideImpl(
    CORE::FE::CellType distype, INPAR::SCATRA::ImplType problem, const int numdofpernode,
    const int numscal, const std::string& disname)
{
  // number of space dimensions
  const int ndim = disname != "scatra_micro"
                       ? GLOBAL::Problem::Instance(
                             DRT::ELEMENTS::ScaTraEleParameterStd::Instance(disname)->ProbNum())
                             ->NDim()
                       : 1;

  switch (distype)
  {
    case CORE::FE::CellType::hex8:
    {
      if (ndim == 3)
        return DefineProblemType<CORE::FE::CellType::hex8, 3>(
            problem, numdofpernode, numscal, disname);
      else
        dserror("invalid problem dimension for HEX8 transport element!");
      break;
    }
    case CORE::FE::CellType::hex27:
    {
      if (ndim == 3)
        return DefineProblemType<CORE::FE::CellType::hex27, 3>(
            problem, numdofpernode, numscal, disname);
      else
        dserror("invalid problem dimension for HEX27 transport element!");
      break;
    }
    case CORE::FE::CellType::tet4:
    {
      if (ndim == 3)
        return DefineProblemType<CORE::FE::CellType::tet4, 3>(
            problem, numdofpernode, numscal, disname);
      else
        dserror("invalid problem dimension for TET4 transport element!");
      break;
    }
    case CORE::FE::CellType::tet10:
    {
      if (ndim == 3)
        return DefineProblemType<CORE::FE::CellType::tet10, 3>(
            problem, numdofpernode, numscal, disname);
      else
        dserror("invalid problem dimension for TET10 transport element!");
      break;
    }
    case CORE::FE::CellType::pyramid5:
    {
      if (ndim == 3)
      {
        return DefineProblemType<CORE::FE::CellType::pyramid5, 3>(
            problem, numdofpernode, numscal, disname);
      }
      else
        dserror("invalid problem dimension for PYRAMID5 transport element!");
      break;
    }
    case CORE::FE::CellType::quad4:
    {
      if (ndim == 2)
        return DefineProblemType<CORE::FE::CellType::quad4, 2>(
            problem, numdofpernode, numscal, disname);
      else if (ndim == 3)
        return DefineProblemType<CORE::FE::CellType::quad4, 3>(
            problem, numdofpernode, numscal, disname);
      else
        dserror("invalid problem dimension for quad4 transport element!");
      break;
    }
    case CORE::FE::CellType::quad9:
    {
      if (ndim == 2)
        return DefineProblemType<CORE::FE::CellType::quad9, 2>(
            problem, numdofpernode, numscal, disname);
      else
      {
        dserror(
            "QUAD9 transport element not implemented as part of %i-dimensional problem. Just do it",
            ndim);
      }
      break;
    }
    case CORE::FE::CellType::nurbs9:
    {
      if (ndim == 2)
        return DefineProblemType<CORE::FE::CellType::nurbs9, 2>(
            problem, numdofpernode, numscal, disname);
      else
      {
        dserror(
            "NURBS9 transport element not implemented as part of %i-dimensional problem. Just do "
            "it",
            ndim);
      }
      break;
    }
    case CORE::FE::CellType::tri3:
    {
      if (ndim == 2)
        return DefineProblemType<CORE::FE::CellType::tri3, 2>(
            problem, numdofpernode, numscal, disname);
      else if (ndim == 3)
        return DefineProblemType<CORE::FE::CellType::tri3, 3>(
            problem, numdofpernode, numscal, disname);
      else
        dserror("invalid problem dimension for tri3 transport element!");
      break;
    }
    case CORE::FE::CellType::tri6:
    {
      if (ndim == 2)
        return DefineProblemType<CORE::FE::CellType::tri6, 2>(
            problem, numdofpernode, numscal, disname);
      else
      {
        dserror(
            "TRI6 transport element not implemented as part of %i-dimensional problem. Just do it",
            ndim);
      }
      break;
    }
    case CORE::FE::CellType::line2:
    {
      if (ndim == 1)
        return DefineProblemType<CORE::FE::CellType::line2, 1>(
            problem, numdofpernode, numscal, disname);
      else if (ndim == 2)
        return DefineProblemType<CORE::FE::CellType::line2, 2>(
            problem, numdofpernode, numscal, disname);
      else if (ndim == 3)
        return DefineProblemType<CORE::FE::CellType::line2, 3>(
            problem, numdofpernode, numscal, disname);
      else
        dserror("invalid problem dimension for LINE2 transport element!");
      break;
    }
    case CORE::FE::CellType::line3:
    {
      if (ndim != 1)
      {
        dserror(
            "LINE3 transport element not implemented as part of %i-dimensional problem. Just do it",
            ndim);
      }
      else
        return DefineProblemType<CORE::FE::CellType::line3, 1>(
            problem, numdofpernode, numscal, disname);
      break;
    }
    default:
      dserror("Element shape %s not activated. Just do it.",
          CORE::FE::CellTypeToString(distype).c_str());
      break;
  }
  return nullptr;
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
DRT::ELEMENTS::ScaTraEleInterface* DRT::ELEMENTS::ScaTraFactory::ProvideImplHDG(
    CORE::FE::CellType distype, INPAR::SCATRA::ImplType problem, const int numdofpernode,
    const int numscal, const std::string& disname)
{
  // -------------------------------------- number of degrees of freedom
  // number of degrees of freedom
  static const int ndim = GLOBAL::Problem::Instance()->NDim();

  switch (distype)
  {
    case CORE::FE::CellType::hex8:
    {
      if (ndim == 3)
      {
        return DefineProblemTypeHDG<CORE::FE::CellType::hex8, 3>(
            problem, numdofpernode, numscal, disname);
      }
      else
        dserror("invalid problem dimension for HEX8 transport element!");
      break;
    }
    case CORE::FE::CellType::tet4:
    {
      if (ndim == 3)
      {
        return DefineProblemTypeHDG<CORE::FE::CellType::tet4, 3>(
            problem, numdofpernode, numscal, disname);
      }
      else
        dserror("invalid problem dimension for TET4 transport element!");
      break;
    }
    case CORE::FE::CellType::tet10:
    {
      if (ndim == 3)
      {
        return DefineProblemTypeHDG<CORE::FE::CellType::tet10, 3>(
            problem, numdofpernode, numscal, disname);
      }
      else
        dserror("invalid problem dimension for TET10 transport element!");
      break;
    }
    case CORE::FE::CellType::pyramid5:
    {
      if (ndim == 3)
      {
        return DefineProblemTypeHDG<CORE::FE::CellType::pyramid5, 3>(
            problem, numdofpernode, numscal, disname);
      }
      else
        dserror("invalid problem dimension for PYRAMID5 transport element!");
      break;
    }
    case CORE::FE::CellType::quad4:
    {
      if (ndim == 2)
      {
        return DefineProblemTypeHDG<CORE::FE::CellType::quad4, 2>(
            problem, numdofpernode, numscal, disname);
      }
      else if (ndim == 3)
      {
        return DefineProblemTypeHDG<CORE::FE::CellType::quad4, 3>(
            problem, numdofpernode, numscal, disname);
      }
      else
        dserror("invalid problem dimension for quad4 transport element!");
      break;
    }
    case CORE::FE::CellType::tri3:
    {
      return DefineProblemTypeHDG<CORE::FE::CellType::tri3, 2>(
          problem, numdofpernode, numscal, disname);
    }
    default:
      dserror("Element shape %s not activated. Just do it.",
          CORE::FE::CellTypeToString(distype).c_str());
      break;
  }
  return nullptr;
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
DRT::ELEMENTS::ScaTraEleInterface* DRT::ELEMENTS::ScaTraFactory::DefineProblemType(
    INPAR::SCATRA::ImplType problem, const int numdofpernode, const int numscal,
    const std::string& disname)
{
  if ((probdim - CORE::FE::dim<distype>) == 1)
  {
    if (problem != INPAR::SCATRA::impltype_std and
        problem != INPAR::SCATRA::impltype_cardiac_monodomain and
        problem != INPAR::SCATRA::impltype_advreac and
        problem != INPAR::SCATRA::impltype_lsreinit and
        problem != INPAR::SCATRA::impltype_one_d_artery and
        problem != INPAR::SCATRA::impltype_no_physics and
        problem != INPAR::SCATRA::impltype_elch_electrode and
        problem != INPAR::SCATRA::impltype_elch_diffcond)
      dserror("ImplType '%s' not implemented for transport on manifolds!",
          SCATRA::ImplTypeToString(problem).c_str());
  }

  switch (problem)
  {
    case INPAR::SCATRA::impltype_std:
    {
      return DRT::ELEMENTS::ScaTraEleCalcStd<distype, probdim>::Instance(
          numdofpernode, numscal, disname);
    }
    case INPAR::SCATRA::impltype_thermo_elch_electrode:
    {
      return DRT::ELEMENTS::ScaTraEleCalcSTIElectrode<distype>::Instance(
          numdofpernode, numscal, disname);
    }
    case INPAR::SCATRA::impltype_thermo_elch_diffcond:
    {
      return DRT::ELEMENTS::ScaTraEleCalcSTIDiffCond<distype>::Instance(
          numdofpernode, numscal, disname);
    }
    case INPAR::SCATRA::impltype_levelset:
    {
      return DRT::ELEMENTS::ScaTraEleCalcLS<distype>::Instance(numdofpernode, numscal, disname);
    }
    case INPAR::SCATRA::impltype_lsreinit:
    {
      return DRT::ELEMENTS::ScaTraEleCalcLsReinit<distype, probdim>::Instance(
          numdofpernode, numscal, disname);
    }
    case INPAR::SCATRA::impltype_loma:
    {
      return DRT::ELEMENTS::ScaTraEleCalcLoma<distype>::Instance(numdofpernode, numscal, disname);
    }
    case INPAR::SCATRA::impltype_elch_NP:
    {
      return DRT::ELEMENTS::ScaTraEleCalcElchNP<distype>::Instance(numdofpernode, numscal, disname);
    }
    case INPAR::SCATRA::impltype_elch_electrode:
    case INPAR::SCATRA::impltype_elch_electrode_growth:
    {
      return DRT::ELEMENTS::ScaTraEleCalcElchElectrode<distype, probdim>::Instance(
          numdofpernode, numscal, disname);
    }
    case INPAR::SCATRA::impltype_elch_electrode_thermo:
    {
      return DRT::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<distype>::Instance(
          numdofpernode, numscal, disname);
    }
    case INPAR::SCATRA::impltype_elch_diffcond:
    {
      return DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype, probdim>::Instance(
          numdofpernode, numscal, disname);
    }
    case INPAR::SCATRA::impltype_elch_diffcond_multiscale:
    {
      return DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<distype, probdim>::Instance(
          numdofpernode, numscal, disname);
    }
    case INPAR::SCATRA::impltype_elch_diffcond_thermo:
    {
      return DRT::ELEMENTS::ScaTraEleCalcElchDiffCondSTIThermo<distype>::Instance(
          numdofpernode, numscal, disname);
    }
    case INPAR::SCATRA::impltype_elch_scl:
    {
      return DRT::ELEMENTS::ScaTraEleCalcElchScl<distype, probdim>::Instance(
          numdofpernode, numscal, disname);
    }
    case INPAR::SCATRA::impltype_poro:
    {
      return DRT::ELEMENTS::ScaTraEleCalcPoro<distype>::Instance(numdofpernode, numscal, disname);
    }
    case INPAR::SCATRA::impltype_advreac:
    {
      return DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype, probdim>::Instance(
          numdofpernode, numscal, disname);
    }
    case INPAR::SCATRA::impltype_refconcreac:
    {
      return DRT::ELEMENTS::ScaTraEleCalcRefConcReac<distype>::Instance(
          numdofpernode, numscal, disname);
    }
    case INPAR::SCATRA::impltype_chemo:
    {
      return DRT::ELEMENTS::ScaTraEleCalcChemo<distype>::Instance(numdofpernode, numscal, disname);
    }
    case INPAR::SCATRA::impltype_chemoreac:
    {
      return DRT::ELEMENTS::ScaTraEleCalcChemoReac<distype>::Instance(
          numdofpernode, numscal, disname);
    }
    case INPAR::SCATRA::impltype_multipororeac:
    {
      return DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>::Instance(
          numdofpernode, numscal, disname);
    }
    case INPAR::SCATRA::impltype_pororeac:
    {
      return DRT::ELEMENTS::ScaTraEleCalcPoroReac<distype>::Instance(
          numdofpernode, numscal, disname);
    }
    case INPAR::SCATRA::impltype_pororeacECM:
    {
      return DRT::ELEMENTS::ScaTraEleCalcPoroReacECM<distype>::Instance(
          numdofpernode, numscal, disname);
    }
    case INPAR::SCATRA::impltype_aniso:
    {
      return DRT::ELEMENTS::ScaTraEleCalcAniso<distype, probdim>::Instance(
          numdofpernode, numscal, disname);
    }
    case INPAR::SCATRA::impltype_cardiac_monodomain:
    {
      return DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<distype, probdim>::Instance(
          numdofpernode, numscal, disname);
    }
    case INPAR::SCATRA::impltype_one_d_artery:
    {
      return DRT::ELEMENTS::ScaTraEleCalcArtery<distype, probdim>::Instance(
          numdofpernode, numscal, disname);
    }
    case INPAR::SCATRA::impltype_no_physics:
      return DRT::ELEMENTS::ScaTraEleCalcNoPhysics<distype, probdim>::Instance(
          numdofpernode, numscal, disname);

    default:
    {
      dserror("Defined problem type does not exist!!");
      break;
    }
  }

  return nullptr;
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
DRT::ELEMENTS::ScaTraEleInterface* DRT::ELEMENTS::ScaTraFactory::DefineProblemTypeHDG(
    INPAR::SCATRA::ImplType problem, const int numdofpernode, const int numscal,
    const std::string& disname)
{
  if (CORE::FE::dim<distype> != probdim)
  {
    if (problem != INPAR::SCATRA::impltype_std and
        problem != INPAR::SCATRA::impltype_cardiac_monodomain)
      dserror("ImplType '%s' not implemented for transport on manifolds!",
          SCATRA::ImplTypeToString(problem).c_str());
  }

  switch (problem)
  {
    case INPAR::SCATRA::impltype_std_hdg:
    {
      return DRT::ELEMENTS::ScaTraEleCalcHDG<distype, probdim>::Instance(
          numdofpernode, numscal, disname);
    }
    case INPAR::SCATRA::impltype_cardiac_monodomain_hdg:
    {
      return DRT::ELEMENTS::ScaTraEleCalcHDGCardiacMonodomain<distype, probdim>::Instance(
          numdofpernode, numscal, disname);
    }
    default:
    {
      dserror("Defined problem type does not exist!!");
      break;
    }
  }

  return nullptr;
}

BACI_NAMESPACE_CLOSE

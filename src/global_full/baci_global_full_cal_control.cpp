/*----------------------------------------------------------------------*/
/*! \file

\brief routine to control execution phase

\level 1


*----------------------------------------------------------------------*/

#include "baci_ale_dyn.H"
#include "baci_art_net_dyn_drt.H"
#include "baci_ehl_dyn.H"
#include "baci_elch_dyn.H"
#include "baci_elemag_dyn.H"
#include "baci_fluid_dyn_nln_drt.H"
#include "baci_fpsi_dyn.H"
#include "baci_fs3i_dyn.H"
#include "baci_fsi_dyn.H"
#include "baci_immersed_problem_dyn.H"
#include "baci_levelset_dyn.H"
#include "baci_lib_globalproblem.H"
#include "baci_loma_dyn.H"
#include "baci_lubrication_dyn.H"
#include "baci_particle_algorithm_sim.H"
#include "baci_pasi_dyn.H"
#include "baci_poroelast_dyn.H"
#include "baci_poroelast_scatra_dyn.H"
#include "baci_porofluidmultiphase_dyn.H"
#include "baci_poromultiphase_dyn.H"
#include "baci_poromultiphase_scatra_dyn.H"
#include "baci_red_airways_dyn_drt.H"
#include "baci_scatra_cardiac_monodomain_dyn.H"
#include "baci_scatra_dyn.H"
#include "baci_ssi_dyn.H"
#include "baci_ssti_dyn.H"
#include "baci_sti_dyn.H"
#include "baci_stru_multi_microstatic_npsupport.H"
#include "baci_structure_dyn_nln_drt.H"
#include "baci_thermo_dyn.H"
#include "baci_tsi_dyn.H"
#include "baci_wear_dyn.H"

/*----------------------------------------------------------------------*
 |  routine to control execution phase                   m.gee 6/01     |
 *----------------------------------------------------------------------*/
void ntacal()
{
  using namespace BACI;

  int restart = DRT::Problem::Instance()->Restart();

  // choose the entry-routine depending on the problem type
  switch (DRT::Problem::Instance()->GetProblemType())
  {
    case ProblemType::structure:
    case ProblemType::polymernetwork:
      caldyn_drt();
      break;
    case ProblemType::fluid:
    case ProblemType::fluid_redmodels:
      dyn_fluid_drt(restart);
      break;
    case ProblemType::lubrication:
      lubrication_dyn(restart);
      break;
    case ProblemType::ehl:
      ehl_dyn();
      break;
    case ProblemType::scatra:
      scatra_dyn(restart);
      break;
    case ProblemType::cardiac_monodomain:
      scatra_cardiac_monodomain_dyn(restart);
      break;
    case ProblemType::sti:
      sti_dyn(restart);
      break;
    case ProblemType::fluid_xfem:
      fluid_xfem_drt();
      break;
      break;
    case ProblemType::fluid_ale:
      fluid_ale_drt();
      break;
    case ProblemType::freesurf:
      fluid_freesurf_drt();
      break;

    case ProblemType::fsi:
    case ProblemType::fsi_redmodels:
    case ProblemType::fsi_lung:
      fsi_ale_drt();
      break;
    case ProblemType::fsi_xfem:
      xfsi_drt();
      break;
    case ProblemType::fpsi_xfem:
      xfpsi_drt();
      break;
    case ProblemType::gas_fsi:
    case ProblemType::ac_fsi:
    case ProblemType::biofilm_fsi:
    case ProblemType::thermo_fsi:
    case ProblemType::fps3i:
      fs3i_dyn();
      break;
    case ProblemType::fbi:
      fsi_immersed_drt();
      break;

    case ProblemType::ale:
      dyn_ale_drt();
      break;

    case ProblemType::thermo:
      thr_dyn_drt();
      break;

    case ProblemType::tsi:
      tsi_dyn_drt();
      break;

    case ProblemType::loma:
      loma_dyn(restart);
      break;

    case ProblemType::elch:
      elch_dyn(restart);
      break;

    case ProblemType::art_net:
      dyn_art_net_drt();
      break;

    case ProblemType::red_airways:
      dyn_red_airways_drt();
      break;

    case ProblemType::struct_ale:
      wear_dyn_drt(restart);
      break;

    case ProblemType::immersed_fsi:
      immersed_problem_drt();
      break;

    case ProblemType::poroelast:
      poroelast_drt();
      break;
    case ProblemType::poroscatra:
      poro_scatra_drt();
      break;
    case ProblemType::porofluidmultiphase:
      porofluidmultiphase_dyn(restart);
      break;
    case ProblemType::poromultiphase:
      poromultiphase_dyn(restart);
      break;
    case ProblemType::poromultiphasescatra:
      poromultiphasescatra_dyn(restart);
      break;
    case ProblemType::fpsi:
      fpsi_drt();
      break;
    case ProblemType::ssi:
      ssi_drt();
      break;
    case ProblemType::ssti:
      ssti_drt();
      break;
    case ProblemType::redairways_tissue:
      redairway_tissue_dyn();
      break;

    case ProblemType::particle:
      particle_drt();
      break;

    case ProblemType::pasi:
      pasi_dyn();
      break;

    case ProblemType::level_set:
      levelset_dyn(restart);
      break;

    case ProblemType::np_support:
      STRUMULTI::np_support_drt();
      break;

    case ProblemType::elemag:
      electromagnetics_drt();
      break;

    default:
      dserror("solution of unknown problemtyp %d requested",
          DRT::Problem::Instance()->GetProblemType());
      break;
  }
}

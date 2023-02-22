/*----------------------------------------------------------------------*/
/*! \file

\brief routine to control execution phase

\level 1


*----------------------------------------------------------------------*/

#include "lib_globalproblem.H"

#include "ale_dyn.H"
#include "art_net_dyn_drt.H"
#include "ehl_dyn.H"
#include "elch_dyn.H"
#include "elemag_dyn.H"
#include "fluid_dyn_nln_drt.H"
#include "fpsi_dyn.H"
#include "fs3i_dyn.H"
#include "fsi_dyn.H"
#include "immersed_problem_dyn.H"
#include "levelset_dyn.H"
#include "loma_dyn.H"
#include "lubrication_dyn.H"
#include "particle_algorithm_sim.H"
#include "pasi_dyn.H"
#include "poroelast_dyn.H"
#include "poroelast_scatra_dyn.H"
#include "poromultiphase_dyn.H"
#include "poromultiphase_scatra_dyn.H"
#include "porofluidmultiphase_dyn.H"
#include "red_airways_dyn_drt.H"
#include "scatra_cardiac_monodomain_dyn.H"
#include "scatra_dyn.H"
#include "ssi_dyn.H"
#include "ssti_dyn.H"
#include "sti_dyn.H"
#include "stru_multi_microstatic_npsupport.H"
#include "structure_dyn_nln_drt.H"
#include "thermo_dyn.H"
#include "tsi_dyn.H"
#include "two_phase_flow_dyn.H"
#include "wear_dyn.H"
#include "inv_analysis_cal_drt.H"
#include "tutorial_dyn.H"

/*----------------------------------------------------------------------*
 |  routine to control execution phase                   m.gee 6/01     |
 *----------------------------------------------------------------------*/
void ntacal()
{
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

    case ProblemType::two_phase_flow:
      two_phase_dyn(restart);
      break;
    case ProblemType::fluid_xfem_ls:
      fluid_xfem_ls_drt(restart);  // Exists in two_phase_flow subfolder
      break;

    case ProblemType::invana:
      invana_cal();
      break;

    case ProblemType::tutorial:
      tutorial_drt();
      break;

    default:
      dserror("solution of unknown problemtyp %d requested",
          DRT::Problem::Instance()->GetProblemType());
      break;
  }
}

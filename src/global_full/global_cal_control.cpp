/*----------------------------------------------------------------------*/
/*! \file

\brief routine to control execution phase

\level 1


*----------------------------------------------------------------------*/

#include "../drt_lib/drt_globalproblem.H"

#include "../drt_ale/ale_dyn.H"
#include "../drt_art_net/art_net_dyn_drt.H"
#include "../drt_ehl/ehl_dyn.H"
#include "../drt_elch/elch_dyn.H"
#include "../drt_elemag/elemag_dyn.H"
#include "../drt_fluid/fluid_dyn_nln_drt.H"
#include "../drt_fpsi/fpsi_dyn.H"
#include "../drt_fs3i/fs3i_dyn.H"
#include "../drt_fsi/fsi_dyn.H"
#include "../drt_immersed_problem/immersed_problem_dyn.H"
#include "../drt_levelset/levelset_dyn.H"
#include "../drt_loma/loma_dyn.H"
#include "../drt_lubrication/lubrication_dyn.H"
#include "../drt_opti/topopt_dyn.H"
#include "../drt_particle_algorithm/particle_sim.H"
#include "../drt_pasi/pasi_dyn.H"
#include "../drt_poroelast/poro_dyn.H"
#include "../drt_poromultiphase/poromultiphase_dyn.H"
#include "../drt_poromultiphase_scatra/poromultiphase_scatra_dyn.H"
#include "../drt_porofluidmultiphase/porofluidmultiphase_dyn.H"
#include "../drt_red_airways/red_airways_dyn_drt.H"
#include "../drt_scatra/scatra_cardiac_monodomain_dyn.H"
#include "../drt_scatra/scatra_dyn.H"
#include "../drt_ssi/ssi_dyn.H"
#include "../drt_ssti/ssti_dyn.H"
#include "../drt_sti/sti_dyn.H"
#include "../drt_stru_multi/microstatic_npsupport.H"
#include "../drt_structure/stru_dyn_nln_drt.H"
#include "../drt_thermo/thr_dyn.H"
#include "../drt_tsi/tsi_dyn.H"
#include "../drt_two_phase_flow/two_phase_dyn.H"
#include "../drt_wear/wear_dyn.H"
#include "../drt_inv_analysis/invana_cal_drt.H"
#include "../drt_tutorial/tutorial_dyn.H"

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

    case ProblemType::fluid_topopt:
      fluid_topopt_dyn();
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
      fluid_xfem_ls_drt(restart);  // Exists in drt_two_phase_flow subfolder
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

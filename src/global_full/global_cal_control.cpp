/*----------------------------------------------------------------------*/
/*! \file

\brief routine to control execution phase

\level 1

\maintainer Martin Kronbichler

*----------------------------------------------------------------------*/

#include "../drt_lib/drt_globalproblem.H"

#include "../drt_acou/acou_dyn.H"
#include "../drt_ale/ale_dyn.H"
#include "../drt_art_net/art_net_dyn_drt.H"
#include "../drt_contact_xcontact/xcontact_dyn.H"
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
#include "../drt_particle/particle_dyn.H"
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
#include "../drt_sti/sti_dyn.H"
#include "../drt_stru_multi/microstatic_npsupport.H"
#include "../drt_structure/stru_dyn_nln_drt.H"
#include "../drt_thermo/thr_dyn.H"
#include "../drt_tsi/tsi_dyn.H"
#include "../drt_two_phase_flow/two_phase_dyn.H"
#include "../drt_variational/var_chemdiff_dyn.H"
#include "../drt_wear/wear_dyn.H"
#ifdef HAVE_FFTW
#include "../drt_mlmc/drt_uq_dyn.H"
#endif
#include "../drt_inv_analysis/invana_cal_drt.H"
#include "../drt_tutorial/tutorial_dyn.H"

/*----------------------------------------------------------------------*
 |  routine to control execution phase                   m.gee 6/01     |
 *----------------------------------------------------------------------*/
void ntacal()
{
  int restart = DRT::Problem::Instance()->Restart();

  // choose the entry-routine depending on the problem type
  switch (DRT::Problem::Instance()->ProblemType())
  {
    case prb_structure:
    case prb_polymernetwork:
      caldyn_drt();
      break;
    case prb_xcontact:
      xcontact_dyn(restart);
      break;
    case prb_fluid:
    case prb_fluid_redmodels:
      dyn_fluid_drt(restart);
      break;
    case prb_lubrication:
      lubrication_dyn(restart);
      break;
    case prb_ehl:
      ehl_dyn();
      break;
    case prb_scatra_endoexocytosis:
    case prb_scatra:
      scatra_dyn(restart);
      break;
    case prb_cardiac_monodomain:
      scatra_cardiac_monodomain_dyn(restart);
      break;
    case prb_sti:
      sti_dyn(restart);
      break;
    case prb_fluid_xfem:
      fluid_xfem_drt();
      break;
      break;
    case prb_fluid_ale:
      fluid_ale_drt();
      break;
    case prb_freesurf:
      fluid_freesurf_drt();
      break;

    case prb_fsi:
    case prb_fsi_redmodels:
    case prb_fsi_lung:
      fsi_ale_drt();
      break;
    case prb_fsi_xfem:
      xfsi_drt();
      break;
    case prb_fpsi_xfem:
      xfpsi_drt();
      break;
    case prb_gas_fsi:
    case prb_ac_fsi:
    case prb_biofilm_fsi:
    case prb_thermo_fsi:
    case prb_fps3i:
      fs3i_dyn();
      break;

    case prb_ale:
      dyn_ale_drt();
      break;

    case prb_thermo:
      thr_dyn_drt();
      break;

    case prb_tsi:
      tsi_dyn_drt();
      break;

    case prb_loma:
      loma_dyn(restart);
      break;

    case prb_elch:
      elch_dyn(restart);
      break;

    case prb_var_chemdiff:
      var_chemdiff_dyn(restart);
      break;

    case prb_fluid_topopt:
      fluid_topopt_dyn();
      break;

    case prb_art_net:
      dyn_art_net_drt();
      break;

    case prb_red_airways:
      dyn_red_airways_drt();
      break;

    case prb_struct_ale:
      wear_dyn_drt(restart);
      break;

    case prb_immersed_fsi:
    case prb_immersed_cell:
    case prb_immersed_ale_fsi:
    case prb_immersed_membrane_fsi:
      immersed_problem_drt();
      break;

    case prb_poroelast:
      poroelast_drt();
      break;
    case prb_poroscatra:
      poro_scatra_drt();
      break;
    case prb_porofluidmultiphase:
      porofluidmultiphase_dyn(restart);
      break;
    case prb_poromultiphase:
      poromultiphase_dyn(restart);
      break;
    case prb_poromultiphasescatra:
      poromultiphasescatra_dyn(restart);
      break;
    case prb_fpsi:
      fpsi_drt();
      break;
    case prb_ssi:
      ssi_drt();
      break;
    case prb_redairways_tissue:
      redairway_tissue_dyn();
      break;

    case prb_particle:
      particle_drt();
      break;

    case prb_particle_old:
      particle_old_drt();
      break;

    case prb_pasi:
      pasi_dyn();
      break;

    case prb_level_set:
      levelset_dyn(restart);
      break;

    case prb_np_support:
      STRUMULTI::np_support_drt();
      break;

    case prb_acou:
      acoustics_drt();
      break;

    case prb_elemag:
      electromagnetics_drt();
      break;

    case prb_two_phase_flow:
      two_phase_dyn(restart);
      break;
    case prb_fluid_xfem_ls:
      fluid_xfem_ls_drt(restart);  // Exists in drt_two_phase_flow subfolder
      break;

    case prb_uq:
#ifdef HAVE_FFTW
      dyn_uq();
#else
      dserror("Uncertainty Quantification only works with FFTW ");
#endif
      break;

    case prb_invana:
      invana_cal();
      break;

    case prb_tutorial:
      tutorial_drt();
      break;

    default:
      dserror(
          "solution of unknown problemtyp %d requested", DRT::Problem::Instance()->ProblemType());
      break;
  }
}

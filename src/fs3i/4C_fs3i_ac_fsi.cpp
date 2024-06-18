/*----------------------------------------------------------------------*/
/*! \file


\brief cpp-file associated with algorithmic routines for two-way coupled partitioned
       solution approaches to fluid-structure-scalar-scalar interaction
       (FS3I). Specifically related version for multiscale approches. This file thereby holds
       all functions related with the small time scales simulation and all the basic control
structures.

\level 3


*----------------------------------------------------------------------*/


#include "4C_fs3i_ac_fsi.hpp"

#include "4C_adapter_ale_fsi.hpp"
#include "4C_adapter_fld_fluid_ac_fsi.hpp"
#include "4C_adapter_str_fsiwrapper.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fluid_utils_mapextractor.hpp"
#include "4C_fsi_monolithic.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_fs3i.hpp"
#include "4C_inpar_fsi.hpp"
#include "4C_inpar_scatra.hpp"
#include "4C_io.hpp"
#include "4C_io_control.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_scatra_algorithm.hpp"
#include "4C_scatra_timint_implicit.hpp"
#include "4C_structure_aux.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | constructor                                               Thon 12/14 |
 *----------------------------------------------------------------------*/
FS3I::ACFSI::ACFSI(const Epetra_Comm& comm)
    : PartFS3I(comm),
      structureincrement_(Teuchos::null),
      fluidincrement_(Teuchos::null),
      aleincrement_(Teuchos::null),
      fluidphinp_lp_(Teuchos::null),
      structurephinp_blts_(Teuchos::null),
      growth_updates_counter_(0),
      wall_shear_stress_lp_(Teuchos::null),
      fsiperiod_(Global::Problem::Instance()->FS3IDynamicParams().sublist("AC").get<double>(
          "PERIODICITY")),
      dt_large_(Global::Problem::Instance()->FS3IDynamicParams().sublist("AC").get<double>(
          "LARGE_TIMESCALE_TIMESTEP")),
      fsiisperiodic_(false),
      scatraisperiodic_(false),
      fsineedsupdate_(false),
      meanmanager_(Teuchos::null)
{
}

/*----------------------------------------------------------------------*
 | Init                                                     rauch 09/16 |
 *----------------------------------------------------------------------*/
void FS3I::ACFSI::init()
{
  FS3I::PartFS3I::init();

  // Some AC FSI specific testings:

  const Teuchos::ParameterList& fs3idyn = Global::Problem::Instance()->FS3IDynamicParams();
  const Teuchos::ParameterList& fs3idynpart = fs3idyn.sublist("PARTITIONED");
  const Teuchos::ParameterList& fs3idynac = fs3idyn.sublist("AC");

  // if ( fs3idyn.get<int>("RESULTSEVRY") != 1 )
  //   FOUR_C_THROW("If you want the fsi problem to be periodically repeated from some point, you
  //   have to have RESULTSEVRY set to 1!");
  // if ( fs3idyn.get<int>("RESTARTEVRY") != 1 )
  //   FOUR_C_THROW("If you want the fsi problem to be periodically repeated from some point, you
  //   have to have RESTARTEVRY set to 1!");

  if (fsiperiod_ <= 1e-14)
    FOUR_C_THROW(
        "You need to specify the PERIODICITY. If you don't want that change your PROBLEMTYPE to "
        "gas_fluid_structure_interaction!");

  if ((dt_large_ + 1e-14) <= dt_)
    FOUR_C_THROW(
        "You need to specify a LARGE_TIMESCALE_TIMESTEP and it must be larger than the 'normal' "
        "fs3i TIMESTEP!");

  std::vector<Core::Conditions::Condition*> ImpCond;
  Global::Problem::Instance()->GetDis("fluid")->GetCondition("ImpedanceCond", ImpCond);
  for (auto& i : ImpCond)
  {
    const double thisperiod = i->parameters().get<double>("TIMEPERIOD");

    if (thisperiod != fsiperiod_)
    {
      FOUR_C_THROW("your impedance TIMEPERIOD and your fs3i PERIODICITY do not match!");
    }
  }

  if (not modulo_is_realtive_zero(fsiperiod_, dt_, fsiperiod_))
    FOUR_C_THROW("Choose a time step such that TIMESTEP = PERIODIC / n with n being an integer!");

  if (not modulo_is_realtive_zero(dt_large_, fsiperiod_, dt_large_))
    FOUR_C_THROW("Choose LARGE_TIMESCALE_TIMESTEP as a multiple of PERIODICITY!");

  if (infperm_)
    FOUR_C_THROW("AC-FS3I does have a finite interface permeability. So set INF_PERM to NO!");

  const int fsiperssisteps = fs3idynac.get<int>("FSI_STEPS_PER_SCATRA_STEP");
  const Inpar::FS3I::SolutionSchemeOverFields couplingalgo =
      Core::UTILS::IntegralValue<Inpar::FS3I::SolutionSchemeOverFields>(fs3idynpart, "COUPALGO");
  const double wk_rel_tol = fs3idynac.get<double>("WINDKESSEL_REL_TOL");

  if (couplingalgo == Inpar::FS3I::fs3i_IterStagg and fsiperssisteps != 1)
    FOUR_C_THROW("In an iteratively staggered FS3I scheme subcycling is not supported!");

  if (couplingalgo == Inpar::FS3I::fs3i_IterStagg and not IsRealtiveEqualTo(wk_rel_tol, -1.0, 1.0))
    FOUR_C_THROW(
        "In an iteratively staggered FS3I scheme periodical repetition of the fsi problem is not "
        "supported!");

  if (not IsRealtiveEqualTo(dt_, fsiperssisteps * fsi_->fluid_field()->Dt(), 1.0))
    FOUR_C_THROW("Your fluid time step does not match!");
  if (not IsRealtiveEqualTo(dt_, fsiperssisteps * fsi_->structure_field()->Dt(), 1.0))
    FOUR_C_THROW("Your structure time step does not match!");
  if (not IsRealtiveEqualTo(dt_, fsiperssisteps * fsi_->ale_field()->Dt(), 1.0))
    FOUR_C_THROW("Your ale time step does not match!");
  if (not IsRealtiveEqualTo(dt_, scatravec_[0]->ScaTraField()->Dt(), 1.0))
    FOUR_C_THROW("Your fluid scatra time step does not match!");
  if (not IsRealtiveEqualTo(dt_, scatravec_[1]->ScaTraField()->Dt(), 1.0))
    FOUR_C_THROW("Your structure scatra time step does not match!");

  if (not(Inpar::FLUID::wss_standard ==
          Core::UTILS::IntegralValue<Inpar::FLUID::WSSType>(
              Global::Problem::Instance()->FluidDynamicParams(), "WSS_TYPE")))
    FOUR_C_THROW(
        "WSS_TYPE must be 'Standard', we will mean the WSS by using the fs3i mean manager!");

  return;
}

/*----------------------------------------------------------------------*
 | Setup                                                    rauch 09/16 |
 *----------------------------------------------------------------------*/
void FS3I::ACFSI::setup()
{
  FS3I::PartFS3I::setup();

  meanmanager_ = Teuchos::rcp(new FS3I::MeanManager(*fsi_->fluid_field()->dof_row_map(0),
      *scatravec_[0]->ScaTraField()->dof_row_map(), *fsi_->fluid_field()->PressureRowMap()));

  structureincrement_ = Core::LinAlg::CreateVector(*fsi_->structure_field()->dof_row_map(0), true);
  fluidincrement_ = Core::LinAlg::CreateVector(*fsi_->fluid_field()->dof_row_map(0), true);
  aleincrement_ = Core::LinAlg::CreateVector(*fsi_->ale_field()->dof_row_map(), true);
  fluidphinp_lp_ = Core::LinAlg::CreateVector(*scatravec_[0]->ScaTraField()->dof_row_map(), true);
  structurephinp_blts_ =
      Core::LinAlg::CreateVector(*scatravec_[1]->ScaTraField()->dof_row_map(), true);
  wall_shear_stress_lp_ = Core::LinAlg::CreateVector(*fsi_->fluid_field()->dof_row_map(0), true);

  extractjthstructscalar_ = BuildMapExtractor();

  return;
}

/*----------------------------------------------------------------------*
 | Read restart                                              Thon 12/14 |
 *----------------------------------------------------------------------*/
void FS3I::ACFSI::read_restart()
{
  // standard restart business
  PartFS3I::read_restart();

  // AC specific restart stuff
  const int restart = Global::Problem::Instance()->restart();

  if (restart)
  {
    const Teuchos::ParameterList& fs3idynac = Global::Problem::Instance()->FS3IDynamicParams();
    const bool restartfrompartfsi =
        Core::UTILS::IntegralValue<int>(fs3idynac, "RESTART_FROM_PART_FSI");

    auto input_control_file = Global::Problem::Instance()->InputControlFile();

    if (not restartfrompartfsi)  // standard restart
    {
      Core::IO::DiscretizationReader fluidreader = Core::IO::DiscretizationReader(
          fsi_->fluid_field()->discretization(), input_control_file, restart);
      meanmanager_->read_restart(fluidreader);

      fsiisperiodic_ = (bool)fluidreader.read_int("fsi_periodic_flag");
      scatraisperiodic_ = (bool)fluidreader.read_int("scatra_periodic_flag");

      fluidreader.read_vector(wall_shear_stress_lp_, "wss_mean");

      if (not(fsiisperiodic_ and scatraisperiodic_))  // restart while in a small time scale loop
      {
        // reconstruct WallShearStress_lp_
        const int beginnperiodstep = get_step_of_beginn_of_this_period_and_prepare_reading(
            fsi_->fluid_field()->Step(), fsi_->fluid_field()->Time(), fsi_->fluid_field()->Dt());
        Core::IO::DiscretizationReader fluidreaderbeginnperiod = Core::IO::DiscretizationReader(
            fsi_->fluid_field()->discretization(), input_control_file, beginnperiodstep);

        // some safety check:
        Teuchos::RCP<Epetra_Vector> WallShearStress_lp_new =
            Core::LinAlg::CreateVector(*fsi_->fluid_field()->dof_row_map(0), true);

        fluidreaderbeginnperiod.read_vector(WallShearStress_lp_new, "SumWss");
        const double SumDtWss = fluidreaderbeginnperiod.read_double("SumDtWss");
        if (abs(SumDtWss - fsiperiod_) > 1e-14)
          FOUR_C_THROW(
              "SumWss and SumDtWss must be read from a step, which was written at the end of a fsi "
              "period!");

        WallShearStress_lp_new->Scale(1 / SumDtWss);
        WallShearStress_lp_new->Update(1.0, *wall_shear_stress_lp_, -1.0);
        double diff_norm(0.0);
        WallShearStress_lp_new->Norm2(&diff_norm);
        if (diff_norm > 1e-10) FOUR_C_THROW("WallShearStress_lp_ is not written/read correctly!");

        // reconstruct fluidphinp_lp_
        fluidreaderbeginnperiod.read_vector(fluidphinp_lp_, "SumPhi");
        fluidphinp_lp_->Scale(1 / SumDtWss);
      }
      else  // restart while in a large time scale loop
      {
      }
    }
    else  // we do not want to read the scatras values and the lagrange multiplyer, since we start
          // from a partitioned FSI
    {
      // AC-FSI specific input
      Core::IO::DiscretizationReader reader = Core::IO::DiscretizationReader(
          fsi_->fluid_field()->discretization(), input_control_file, restart);
      reader.read_vector(wall_shear_stress_lp_, "wss");
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 | Timeloop                                                  Thon 12/14 |
 *----------------------------------------------------------------------*/
void FS3I::ACFSI::Timeloop()
{
  check_is_init();
  check_is_setup();

  // prepare time loop
  fsi_->PrepareTimeloop();
  SetFSISolution();

  // calculate inital time derivative, when restart was done from a part. FSI simulation
  if (Global::Problem::Instance()->restart() and
      Core::UTILS::IntegralValue<int>(
          Global::Problem::Instance()->FS3IDynamicParams(), "RESTART_FROM_PART_FSI"))
  {
    scatravec_[0]->ScaTraField()->prepare_first_time_step();
    scatravec_[1]->ScaTraField()->prepare_first_time_step();
  }

  // output of initial state
  //  if (step_ == 0)
  {
    constexpr bool force_prepare = true;
    fsi_->prepare_output(force_prepare);
    FsiOutput();
    ScatraOutput();
  }

  // do time loop
  while (NotFinished())
  {
    SmallTimeScaleLoop();

    if (NotFinished()) LargeTimeScaleLoop();
  }
}

/*----------------------------------------------------------------------*
 | timeloop for small time scales                            Thon 07/15 |
 *----------------------------------------------------------------------*/
void FS3I::ACFSI::SmallTimeScaleLoop()
{
  // print info
  if (Comm().MyPID() == 0)
  {
    std::cout << "\n\n************************************************************************"
                 "\n                         SMALL TIME SCALE LOOP"
                 "\n************************************************************************"
              << std::endl;
  }

  while (small_time_scale_loop_not_finished())
  {
    // Do a time step
    small_time_scale_prepare_time_step();

    small_time_scale_outer_loop();

    small_time_scale_update_and_output();
  }
}

/*--------------------------------------------------------------------------*
 | flag whether small time scale time loop is finished           Thon 07/15 |
 *--------------------------------------------------------------------------*/
bool FS3I::ACFSI::small_time_scale_loop_not_finished()
{
  return (NotFinished() and not(fsiisperiodic_ and scatraisperiodic_));
}

/*----------------------------------------------------------------------*
 | Prepare time step                                         Thon 12/14 |
 *----------------------------------------------------------------------*/
void FS3I::ACFSI::small_time_scale_prepare_time_step()
{
  increment_time_and_step();

  // Print to screen
  if (Comm().MyPID() == 0)
  {
    std::cout << "\n"
              << "TIME:  " << std::scientific << std::setprecision(12) << time_ << "/"
              << std::setprecision(4) << timemax_ << "     DT = " << std::scientific << dt_
              << "     STEP = " << std::setw(4) << step_ << "/" << std::setw(4) << numstep_ << "\n"
              << std::endl;
  }

  // iff this is the beginning of a new fsi cycle
  if (step_ > 1 and modulo_is_realtive_zero(time_ - dt_, fsiperiod_, time_))
  {
    meanmanager_->Reset();  // Reset mean Manager
    if (Comm().MyPID() == 0)
    {
      std::cout << "Reseting mean manager\n" << std::endl;
    }
  }

  // NOTE: We have to set the concentrations here, since the structure does call an
  // evaluate inside prepare_time_step() in order to calculate the predictor error. Hence first:
  set_struct_scatra_solution();
  // prepare time step for fsi subproblem
  fsi_->prepare_time_step();

  // prepare time step for both fluid- and structure-based scatra field
  // SetFSISolution();
  for (unsigned i = 0; i < scatravec_.size(); ++i)
  {
    scatravec_[i]->ScaTraField()->prepare_time_step();
  }
}

/*----------------------------------------------------------------------*
 | outer_loop                                                 Thon 12/14 |
 *----------------------------------------------------------------------*/
void FS3I::ACFSI::small_time_scale_outer_loop()
{
  const Teuchos::ParameterList& fs3idynpart =
      Global::Problem::Instance()->FS3IDynamicParams().sublist("PARTITIONED");
  // get coupling algorithm from input file
  const Inpar::FS3I::SolutionSchemeOverFields couplingalgo =
      Core::UTILS::IntegralValue<Inpar::FS3I::SolutionSchemeOverFields>(fs3idynpart, "COUPALGO");

  switch (couplingalgo)
  {
    case Inpar::FS3I::fs3i_SequStagg:
      small_time_scale_outer_loop_sequ_stagg();
      break;
    case Inpar::FS3I::fs3i_IterStagg:
      small_time_scale_outer_loop_iter_stagg();
      break;
    default:
      FOUR_C_THROW("partitioned FS3I coupling scheme not implemented!");
      break;
  }
}

/*----------------------------------------------------------------------*
 | outer_loop for sequentially staggered FS3I scheme          Thon 12/14 |
 *----------------------------------------------------------------------*/
void FS3I::ACFSI::small_time_scale_outer_loop_sequ_stagg()
{
  set_struct_scatra_solution();

  DoFSIStep();

  SetFSISolution();

  small_time_scale_do_scatra_step();
}

/*----------------------------------------------------------------------*
 | outer_loop for iterative staggered FS3I scheme             Thon 12/14 |
 *----------------------------------------------------------------------*/
void FS3I::ACFSI::small_time_scale_outer_loop_iter_stagg()
{
  int itnum = 0;

  bool stopnonliniter = false;

  if (Comm().MyPID() == 0)
  {
    std::cout << "\n************************************************************************\n"
                 "                         OUTER ITERATION START"
              << "\n************************************************************************"
              << std::endl;
  }

  while (stopnonliniter == false)
  {
    itnum++;

    structureincrement_->Update(1.0, *fsi_->structure_field()->Dispnp(), 0.0);
    fluidincrement_->Update(1.0, *fsi_->fluid_field()->Velnp(), 0.0);
    aleincrement_->Update(1.0, *fsi_->ale_field()->Dispnp(), 0.0);

    set_struct_scatra_solution();

    DoFSIStep();

    SetFSISolution();

    small_time_scale_do_scatra_step();

    stopnonliniter = part_fs3i_convergence_ckeck(itnum);
  }
}

/*----------------------------------------------------------------------*
 | Do a single fsi step                                                 |
 | (including subcycling and periodic repetition )           Thon 12/14 |
 *----------------------------------------------------------------------*/
void FS3I::ACFSI::DoFSIStep()
{
  if (not fsiisperiodic_)
  {
    const int fsiperssisteps =
        Global::Problem::Instance()->FS3IDynamicParams().sublist("AC").get<int>(
            "FSI_STEPS_PER_SCATRA_STEP");
    if (fsiperssisteps == 1)  // no subcycling
    {
      // this is the normal case, here we also check before
      check_if_times_and_steps_and_dts_match();

      DoFSIStepStandard();
    }
    else  // subcycling
    {
      DoFSIStepSubcycled(fsiperssisteps);
    }
  }
  else  // fsi problem is periodic
  {
    DoFSIStepPeriodic();
  }

  // just for safety reasons we check if all marching quantities match
  // NOTE: we do this after the actual solving since in case of subcycling the may differ in between
  check_if_times_and_steps_and_dts_match();
}

void FS3I::ACFSI::is_small_time_scale_periodic()
{
  if (modulo_is_realtive_zero(time_, fsiperiod_, time_))  // iff a new period is about to begin
  {
    // Check fsi periodicity
    IsFsiPeriodic();

    // Check if scatra is periodic
    IsScatraPeriodic();
  }
}

/*----------------------------------------------------------------------*
 | Decide if fsi problem is already periodic                 Thon 12/14 |
 *----------------------------------------------------------------------*/
void FS3I::ACFSI::IsFsiPeriodic()
{
  fsiisperiodic_ = true;

  //------------- first check convergence of Windkessels ------------------------------
  {
    Teuchos::RCP<Adapter::FluidACFSI> fluid =
        Teuchos::rcp_dynamic_cast<Adapter::FluidACFSI>(fsi_->fluid_field());
    std::vector<double> wk_rel_errors = fluid->GetWindkesselErrors();

    const double wk_rel_tol =
        Global::Problem::Instance()->FS3IDynamicParams().sublist("AC").get<double>(
            "WINDKESSEL_REL_TOL");

    if (Comm().MyPID() == 0) std::cout << std::endl;

    for (unsigned int i = 0; i < wk_rel_errors.size(); i++)
    {
      if (Comm().MyPID() == 0)
        std::cout << std::scientific << std::setprecision(2) << "The " << i + 1
                  << "-th Windkessel has a relative error of      " << abs(wk_rel_errors[i])
                  << "  (tol " << wk_rel_tol << ")" << std::endl;

      if (abs(wk_rel_errors[i]) > wk_rel_tol) fsiisperiodic_ = false;
    }
  }

  //------------- second check convergence of WSS ------------------------------
  {
    const double wss_rel_tol =
        Global::Problem::Instance()->FS3IDynamicParams().sublist("AC").get<double>("WSS_REL_TOL");

    // Note: we compare the difference of mean wssnp and mean wssn only on the FS3I interface, since
    // this is the interesting part for the multiscale problem Therefore we have to do nothing
    // special since the rest is zero anyways
    const Teuchos::RCP<const Epetra_Vector> wss_bar_boundary =
        meanmanager_->GetMeanValue("mean_wss");  // mean wall shear stresses at fs3i interface
    const Teuchos::RCP<Epetra_Vector> wssdiff_bar_boundary =
        Core::LinAlg::CreateVector(*(fsi_->fluid_field()->dof_row_map()));
    wssdiff_bar_boundary->Update(1.0, *wss_bar_boundary, -1.0, *wall_shear_stress_lp_, 0.0);

    double wss_bar_boundary_norm(0.0);
    wall_shear_stress_lp_->Norm2(&wss_bar_boundary_norm);
    if (wss_bar_boundary_norm < 1e-10)  // e.g. at the first time step this will happen...
    {
      wss_bar_boundary_norm = 1.0;  //... so we simply check the absolut error.
    }

    double wssdiff_bar_boundary_norm(0.0);
    wssdiff_bar_boundary->Norm2(&wssdiff_bar_boundary_norm);

    // compute the actual relative error
    double wss_rel_error = wssdiff_bar_boundary_norm / wss_bar_boundary_norm;

    if (Comm().MyPID() == 0)
    {
      std::cout << std::scientific << std::setprecision(2)
                << "The wall shear stresses have a relative error of " << wss_rel_error << "  (tol "
                << wss_rel_tol << ")" << std::endl;
    }

    wall_shear_stress_lp_->Update(1.0, *(meanmanager_->GetMeanValue("mean_wss")), 0.0);  // Update

    if (wss_rel_error > wss_rel_tol) fsiisperiodic_ = false;
  }

  return;
}

///*----------------------------------------------------------------------*
// |  Extract wall shear stresses                              Thon 11/14 |
// *----------------------------------------------------------------------*/
void FS3I::ACFSI::set_wall_shear_stresses() const
{
  // set current WSS
  PartFS3I::set_wall_shear_stresses();

  // Set MeanWSS
  // set_mean_wall_shear_stresses();
}

/*----------------------------------------------------------------------*
 | Decide if fluid scatra problem is periodic                Thon 07/15 |
 *----------------------------------------------------------------------*/
void FS3I::ACFSI::IsScatraPeriodic()
{
  scatraisperiodic_ = true;

  const double fluid_scatra_rel_tol =
      Global::Problem::Instance()->FS3IDynamicParams().sublist("AC").get<double>(
          "FLUID_SCATRA_REL_TOL");
  const int numscal = scatravec_[0]->ScaTraField()->NumScal();  // numscal of fluid scatra field

  // we test all scalars individually
  for (int i = 0; i < numscal; i++)
  {
    // Note: we compare the difference of mean phinp and mean phin only on the FS3I interface, since
    // this is the interesting part for the multiscale problem

    // extract interface concentrations
    const Teuchos::RCP<const Epetra_Vector> phinp_bar_boundary =
        scatrafieldexvec_[0]->ExtractVector(
            meanmanager_->GetMeanValue("mean_phi"), 1);  // mean fluidscatra at fs3i interface
    // move to struct interface (with full struct scatra dof map)
    const Teuchos::RCP<const Epetra_Vector> phinp_bar_boundary_on_struct =
        scatrafieldexvec_[1]->InsertVector(Scatra1ToScatra2(phinp_bar_boundary), 1);
    // extract i-th scalar
    const Teuchos::RCP<const Epetra_Vector> ith_phinp_bar_boundary_on_struct =
        extractjthstructscalar_[i]->ExtractCondVector(phinp_bar_boundary_on_struct);

    // extract interface concentrations
    const Teuchos::RCP<const Epetra_Vector> phin_bar_boundary = scatrafieldexvec_[0]->ExtractVector(
        fluidphinp_lp_, 1);  // mean fluidscatra at fs3i interface
    // move to struct interface (with full struct scatra dof map)
    const Teuchos::RCP<const Epetra_Vector> phin_bar_boundary_on_struct =
        scatrafieldexvec_[1]->InsertVector(Scatra1ToScatra2(phin_bar_boundary), 1);
    // extract i-th scalar
    const Teuchos::RCP<const Epetra_Vector> ith_phin_bar_boundary_on_struct =
        extractjthstructscalar_[i]->ExtractCondVector(phin_bar_boundary_on_struct);

    // calculate the difference vector
    const Teuchos::RCP<Epetra_Vector> ith_phi_diff_bar_boundary =
        extractjthstructscalar_[i]->ExtractCondVector(
            Core::LinAlg::CreateVector(*scatravec_[1]->ScaTraField()->dof_row_map(0), true));
    ith_phi_diff_bar_boundary->Update(
        1.0, *ith_phinp_bar_boundary_on_struct, -1.0, *ith_phin_bar_boundary_on_struct, 0.0);

    double ith_phin_bar_boundary_on_struct_norm(0.0);
    ith_phin_bar_boundary_on_struct->Norm2(&ith_phin_bar_boundary_on_struct_norm);
    if (ith_phin_bar_boundary_on_struct_norm <
        1e-10)  // e.g. at the first time step this will happen...
    {
      ith_phin_bar_boundary_on_struct_norm = 1.0;  //... so we simply check the absolut error.
    }

    double ith_phi_diff_bar_boundary_norm(0.0);
    ith_phi_diff_bar_boundary->Norm2(&ith_phi_diff_bar_boundary_norm);

    // compute the actual relative error
    double ith_fluid_scatra_rel_error =
        ith_phi_diff_bar_boundary_norm / ith_phin_bar_boundary_on_struct_norm;

    if (Comm().MyPID() == 0)
    {
      std::cout << std::scientific << std::setprecision(2) << "The " << i + 1
                << "-th fluid-scatra has a relative error of    " << ith_fluid_scatra_rel_error
                << "  (tol " << fluid_scatra_rel_tol << ")" << std::endl;
    }

    if (ith_fluid_scatra_rel_error > fluid_scatra_rel_tol) scatraisperiodic_ = false;
  }

  fluidphinp_lp_->Update(1.0, *meanmanager_->GetMeanValue("mean_phi"), 0.0);

  return;
}

/*----------------------------------------------------------------------*
 | Do a standard fsi step                                    Thon 12/14 |
 *----------------------------------------------------------------------*/
void FS3I::ACFSI::DoFSIStepStandard()
{
  if (Comm().MyPID() == 0)
  {
    std::cout << "\n************************************************************************"
                 "\n                               FSI SOLVER "
                 "\n************************************************************************"
              << std::endl;
  }
  fsi_->TimeStep(fsi_);
}

/*----------------------------------------------------------------------*
 | Do a fsi step with subcycling                             Thon 12/14 |
 *----------------------------------------------------------------------*/
void FS3I::ACFSI::DoFSIStepSubcycled(const int subcyclingsteps)
{
  int subcyclingiter = 0;

  while (subcyclingiter < subcyclingsteps)
  {
    subcyclingiter++;

    if (subcyclingiter != 1)  // for the first subcycling step we...
    {
      constexpr bool force_prepare = false;
      fsi_->prepare_output(force_prepare);  //... will do this in update_and_output()
      fsi_->update();                       //... will do this in update_and_output()
      fsi_->prepare_time_step();            //... have already done this in prepare_time_step()
      // now fix the step_ counter. When subcycling the fsi subproblem we do not want to proceed the
      // step_ AND the time_, but just the time_.
      SetTimeAndStepInFSI(fsi_->fluid_field()->Time(), step_);
    }

    if (Comm().MyPID() == 0)
    {
      std::cout << "\n************************************************************************"
                   "\n                     FSI SUBCYCLING SOLVER "
                << subcyclingiter << "/" << subcyclingsteps
                << "\n************************************************************************"
                << std::endl;
    }

    fsi_->TimeStep(fsi_);  // all necessary changes for the fsi problem (i.e. adapting dt) has
                           // already been done in PartFS3I::ManipulateDt()
  }

  // do time and step in fsi (including subfields) and fs3i match after we will have called
  // small_time_scale_update_and_output()
  if (not IsRealtiveEqualTo(fsi_->fluid_field()->Time(), time_, time_))
    FOUR_C_THROW("After the subcycling the fsi time and fs3i time do not match anymore!");

  if (not IsRealtiveEqualTo(fsi_->fluid_field()->Step(), step_, step_))
    FOUR_C_THROW("After the subcycling the fsi step and fs3i step do not match anymore!");
}

/*----------------------------------------------------------------------*
 | Get fsi solution from one period before                 Thon 12/14 |
 *----------------------------------------------------------------------*/
void FS3I::ACFSI::DoFSIStepPeriodic()
{
  // instead of solving the FSI problem for the present time step, we take the
  // solution from the FSI period before. We do this by replacing all values
  // in FSI via read_restart(). Afterwards we just have to repair everything we
  // destroyed by calling read_restart().

  // this is the related step of the last period
  int previousperiodstep = get_step_of_one_period_ago_and_prepare_reading(step_, time_);

  if (Comm().MyPID() == 0)
  {
    std::cout << "\n************************************************************************"
                 "\n                         PERIODICAL FSI STEP"
                 "\n************************************************************************\n"
                 "\nUsing results from timestep "
              << previousperiodstep << " as solution of the current FSI step" << std::endl;
  }

  // do the reading
  fsi_->read_restart(previousperiodstep);  // ReadRestartfromTime()

  // we first fix the grid velocity of the fluid. This calculation is normally done
  // in Adapter::FluidFSI::apply_mesh_displacement(), but since we never call this function
  // we have to call it ourself
  fsi_->fluid_field()->UpdateGridv();  // calculate grid velocity via FD approximation

  // update time and step in FSI and all subproblems
  SetTimeAndStepInFSI(time_, step_);
}

/*--------------------------------------------------------------------------------------*
 | Get step number of on cycle ago                                           Thon 07/15 |
 *--------------------------------------------------------------------------------------*/
double FS3I::ACFSI::get_step_of_one_period_ago_and_prepare_reading(
    const int actstep, const double acttime)
{
  if (not modulo_is_realtive_zero(fsiperiod_, dt_, fsiperiod_))
    FOUR_C_THROW("PERIODICITY should be an multiple of TIMESTEP!");

  const int previousperiodstep =
      actstep - round(fsiperiod_ / dt_);  // here we assume a constant timestep over the last period

  // Is this the right step? Let's check:
  {
    // we have to clean the mapstack, otherwise it would fill with each iteration
    Global::Problem::Instance()->GetDis("structure")->Writer()->clear_map_cache();

    // get filename in which the equivalent step of the last period is written
    std::string filename = GetFileName(previousperiodstep);
    // we always have to recreate the InputControl() since our Inputfile (=Outputfile) has changed
    // in since the last reading (new timestep written)
    Teuchos::RCP<Core::IO::InputControl> inputreader =
        Teuchos::rcp(new Core::IO::InputControl(filename, Comm()));
    // overwrite existing InputControl()
    Global::Problem::Instance()->SetInputControlFile(inputreader);

    // AC-FSI specific input
    Core::IO::DiscretizationReader reader =
        Core::IO::DiscretizationReader(fsi_->fluid_field()->discretization(),
            Global::Problem::Instance()->InputControlFile(), previousperiodstep);

    double previousperiodtime = reader.read_double("time");

    // Now check if the candidate is right
    if (not IsRealtiveEqualTo(previousperiodtime + fsiperiod_, acttime, acttime))
      FOUR_C_THROW("You can't change your TIMESTEP when you are within the FSI PERIODIC cycle!");
  }

  return previousperiodstep;
}

/*--------------------------------------------------------------------------------------*
 | Get step number of the beginning of this cycle                            Thon /15 |
 *--------------------------------------------------------------------------------------*/
double FS3I::ACFSI::get_step_of_beginn_of_this_period_and_prepare_reading(
    const int actstep, const double acttime, const double dt)
{
  if (not modulo_is_realtive_zero(fsiperiod_, dt, fsiperiod_))
    FOUR_C_THROW("PERIODICITY should be an multiple of TIMESTEP!");

  const double beginnperiodtime = acttime - (fmod(acttime + 1e-14, fsiperiod_) - 1e-14);
  const int teststep =
      actstep - round((fmod(acttime + 1e-14, fsiperiod_) - 1e-14) /
                      dt);  // here we assume a constant timestep over the last period

  // Is this the right step? Let's check:
  {
    // we have to clean the mapstack, otherwise it would fill with each iteration
    Global::Problem::Instance()->GetDis("structure")->Writer()->clear_map_cache();

    // get filename in which the equivalent step of the last period is written
    std::string filename = GetFileName(teststep);
    // we always have to recreate the InputControl() since our Inputfile (=Outputfile) has changed
    // in since the last reading (new timestep written)
    Teuchos::RCP<Core::IO::InputControl> inputreader =
        Teuchos::rcp(new Core::IO::InputControl(filename, Comm()));
    // overwrite existing InputControl()
    Global::Problem::Instance()->SetInputControlFile(inputreader);

    // AC-FSI specific input
    Core::IO::DiscretizationReader reader =
        Core::IO::DiscretizationReader(fsi_->fluid_field()->discretization(),
            Global::Problem::Instance()->InputControlFile(), teststep);

    double testtime = reader.read_double("time");

    if (testtime - beginnperiodtime > 1e-10)
    {
      return get_step_of_beginn_of_this_period_and_prepare_reading(teststep, testtime, 0.5 * dt);
    }
    else if (beginnperiodtime - testtime > 1e-10)
    {
      return get_step_of_beginn_of_this_period_and_prepare_reading(actstep, acttime, 2 * dt);
    }

    // Now check if the candidate is right
    if (not modulo_is_realtive_zero(testtime, fsiperiod_, testtime))
      FOUR_C_THROW("Why is the time not a multiple of FSI PERIODIC cycle??");
  }

  return teststep;
}

/*--------------------------------------------------------------------------------------*
 | Get filename in which the equivalent step of the last period is written   Thon 12/14 |
 *--------------------------------------------------------------------------------------*/
std::string FS3I::ACFSI::GetFileName(const int step)
{
  // we have to be careful since in the fist steps after a restart the last fsi period is written
  // in the input and not in the output file

  std::string filename;

  const int restart = Global::Problem::Instance()->restart();
  if (restart)
  {
    const int crit_step = (int)(restart + fsiperiod_ / dt_ + 1e-14);

    if (step <= crit_step)  // the last period is written in the file we have restarted from
    {
      filename = Global::Problem::Instance()->InputControlFile()->file_name();
    }
    else  // the last period is written in the newly written output file
    {
      filename = Global::Problem::Instance()->OutputControlFile()->file_name();
    }
  }
  else
  {
    filename = Global::Problem::Instance()->OutputControlFile()->file_name();
  }
  return filename;
}

/*----------------------------------------------------------------------*
 | Set time and step in FSI and all subfields                Thon 12/14 |
 *----------------------------------------------------------------------*/
void FS3I::ACFSI::SetTimeAndStepInFSI(const double time, const int step)
{
  // Set time and step in FSI. This looks a bit strange but in case of subcyling
  // we want a 'proper' screen output. And since the fsi time and step are only used
  // for this purpose, it's not really important what we do here
  fsi_->SetTimeStep(time, step);
  // Note: The last function did not touch the subfields, so we have to do it yourself

  // Set time and step in structure field
  fsi_->structure_field()->set_time(time - fsi_->Dt());
  fsi_->structure_field()->SetTimen(time);
  fsi_->structure_field()->SetStep(step - 1);
  fsi_->structure_field()->SetStepn(step);

  // Set time and step in fluid field
  fsi_->fluid_field()->SetTimeStep(time, step);

  // Set time and step in ale field
  fsi_->ale_field()->SetTimeStep(time, step);
}

/*----------------------------------------------------------------------*
 | Do a single scatra step                                   Thon 12/14 |
 *----------------------------------------------------------------------*/
void FS3I::ACFSI::small_time_scale_do_scatra_step()
{
  if (Comm().MyPID() == 0)
  {
    std::cout << "\n************************************************************************"
                 "\n                        AC COUPLED SCATRA SOLVER"
                 "\n************************************************************************\n"
              << std::endl;

    std::cout << "+- step/max -+-- scal-res/ abs-tol [norm] -+-- scal-inc/ rel-tol [norm] -+"
              << std::endl;
  }

  bool stopnonliniter = false;
  int itnum = 0;

  while (stopnonliniter == false)
  {
    scatra_evaluate_solve_iter_update();
    itnum++;
    if (scatra_convergence_check(itnum)) break;
  }
}

/*----------------------------------------------------------------------*
 | Update and output the small time scale                    Thon 12/14 |
 *----------------------------------------------------------------------*/
void FS3I::ACFSI::small_time_scale_update_and_output()
{
  // Update mean manager
  meanmanager_->AddValue("phi", scatravec_[0]->ScaTraField()->Phinp(), dt_);
  meanmanager_->AddValue(
      "pressure", fsi_->fluid_field()->ExtractPressurePart(fsi_->fluid_field()->Velnp()), dt_);
  if (not fsiisperiodic_)
    meanmanager_->AddValue("wss", fsi_->fluid_field()->calculate_wall_shear_stresses(),
        dt_);  // add instantaneous wss vector
  else
  {
    // fsi is periodic, hence wss do not change any more. But to satisfy the manager we have to add
    // something..
    meanmanager_->AddValue("wss", wall_shear_stress_lp_, dt_);
  }

  // NOTE: we can not reset the mean manager here, since we first need to write its data, to be able
  // to restart. Hence the correct order is: 1. Calculate mean values; 2. Write output; 3. Reset
  // mean Manager (in next time step)

  // NOTE: it is important to update the periodic flags AFTER doing the windkessel
  // time update and BEFORE writing the output. So be careful with the order here!!

  // Update field variables
  constexpr bool force_prepare = false;
  fsi_->prepare_output(force_prepare);
  fsi_->update();
  UpdateScatraFields();

  // Update the isperiodic_ flags
  is_small_time_scale_periodic();

  // write outputs
  FsiOutput();
  ScatraOutput();
}

/*----------------------------------------------------------------------*
 | Write FSI output                                          Thon 03/15 |
 *----------------------------------------------------------------------*/
void FS3I::ACFSI::FsiOutput()
{
  int coupling =
      Teuchos::getIntegralValue<int>(Global::Problem::Instance()->FSIDynamicParams(), "COUPALGO");

  const Teuchos::ParameterList& fs3idyn = Global::Problem::Instance()->FS3IDynamicParams();
  const int uprestart = fs3idyn.get<int>("RESTARTEVRY");
  const int upresults = fs3idyn.get<int>("RESULTSEVRY");

  //    /* Note: The order is important here! In here control file entries are
  //     * written. And these entries define the order in which the filters handle
  //     * the Discretizations, which in turn defines the dof number ordering of the
  //     * Discretizations.

  // structure output
  fsi_->structure_field()->output();
  if (coupling == fsi_iter_monolithicstructuresplit) fsi_->OutputLambda();

  // fluid output
  fsi_->fluid_field()->output();
  if (coupling == fsi_iter_monolithicfluidsplit) fsi_->OutputLambda();

  if ((step_ % upresults == 0) or (uprestart != 0 && step_ % uprestart == 0) or
      modulo_is_realtive_zero(time_, fsiperiod_, time_))
  {
    Teuchos::RCP<Core::IO::DiscretizationWriter> fluiddiskwriter =
        fsi_->fluid_field()->DiscWriter();
    fluiddiskwriter->write_vector("wss_mean", wall_shear_stress_lp_);
  }
  // AC specific fluid output iff it is a restart step or we are at the end of a fsi circle
  if ((uprestart != 0 && step_ % uprestart == 0) or
      modulo_is_realtive_zero(time_, fsiperiod_, time_))
  {
    Teuchos::RCP<Core::IO::DiscretizationWriter> fluiddiskwriter =
        fsi_->fluid_field()->DiscWriter();

    // fluiddiskwriter->write_vector("wss", fsi_->fluid_field()->calculate_wall_shear_stresses());
    meanmanager_->write_restart(fluiddiskwriter);

    fluiddiskwriter->write_int("fsi_periodic_flag", (int)fsiisperiodic_);
    fluiddiskwriter->write_int("scatra_periodic_flag", (int)scatraisperiodic_);
  }

  // ale output
  fsi_->ale_field()->output();
}

/*----------------------------------------------------------------------*
 | check convergence of scatra fields                        Thon 12/14 |
 *----------------------------------------------------------------------*/
bool FS3I::ACFSI::scatra_convergence_check(const int itnum)
{
  // some input parameters for the scatra fields
  const Teuchos::ParameterList& scatradyn =
      Global::Problem::Instance()->scalar_transport_dynamic_params();
  const int scatraitemax = scatradyn.sublist("NONLINEAR").get<int>("ITEMAX");
  const double scatraittol = scatradyn.sublist("NONLINEAR").get<double>("CONVTOL");
  const double scatraabstolres = scatradyn.sublist("NONLINEAR").get<double>("ABSTOLRES");


  double conresnorm(0.0);
  scatrarhs_->Norm2(&conresnorm);
  double incconnorm(0.0);
  scatraincrement_->Norm2(&incconnorm);

  // set up vector of absolute concentrations
  Teuchos::RCP<Epetra_Vector> con = Teuchos::rcp(new Epetra_Vector(scatraincrement_->Map()));
  Teuchos::RCP<const Epetra_Vector> scatra1 = scatravec_[0]->ScaTraField()->Phinp();  // fluidscatra
  Teuchos::RCP<const Epetra_Vector> scatra2 =
      scatravec_[1]->ScaTraField()->Phinp();  // structurescatra
  setup_coupled_scatra_vector(con, scatra1, scatra2);

  double connorm(0.0);
  con->Norm2(&connorm);

  // care for the case that nothing really happens in the concentration field
  if (connorm < 1e-5) connorm = 1.0;

  // print the screen info
  if (Comm().MyPID() == 0)
  {
    printf("|   %3d/%3d  |  %1.3E/ %1.1E [L_2 ]  |  %1.3E/ %1.1E [L_2 ]  |\n", itnum, scatraitemax,
        conresnorm, scatraabstolres, incconnorm / connorm, scatraittol);
  }

  // this is the convergence check
  // We always require at least one solve. We test the L_2-norm of the
  // current residual. Norm of residual is just printed for information
  if (conresnorm <= scatraabstolres and incconnorm / connorm <= scatraittol)
  {
    if (Comm().MyPID() == 0)
    {
      // print 'finish line'
      printf("+------------+-----------------------------+-----------------------------+\n\n");
    }
    return true;
  }
  // if itemax is reached without convergence stop the simulation
  else if (itnum == scatraitemax)
  {
    if (Comm().MyPID() == 0)
    {
      printf("+---------------------------------------------------------------+\n");
      printf("|    scalar-scalar field did not converge in itemax steps!     |\n");
      printf("+---------------------------------------------------------------+\n");
    }
    // yes, we stop!
    //    FOUR_C_THROW("Scatra not converged in itemax steps!");
    return true;
  }
  else
    return false;
}

/*----------------------------------------------------------------------*
 | Convergence check for iterative staggered FS3I scheme     Thon 12/14 |
 *----------------------------------------------------------------------*/
bool FS3I::ACFSI::part_fs3i_convergence_ckeck(const int itnum)
{
  const Teuchos::ParameterList& fs3idynpart =
      Global::Problem::Instance()->FS3IDynamicParams().sublist("PARTITIONED");
  // get control parameters from input file
  double ittol = fs3idynpart.get<double>("CONVTOL");
  int itmax = fs3idynpart.get<int>("ITEMAX");

  // convergence check based on the scalar increment
  bool stopnonliniter = false;

  // calculate fsi increments. scatra increment is already done in the scatra fields convergence
  // check
  structureincrement_->Update(1.0, *fsi_->structure_field()->Dispnp(), -1.0);
  fluidincrement_->Update(1.0, *fsi_->fluid_field()->Velnp(), -1.0);
  aleincrement_->Update(1.0, *fsi_->ale_field()->Dispnp(), -1.0);

  // L2-norms of increment vectors
  double scatraincconnorm_L2(0.0);
  scatraincrement_->Norm2(&scatraincconnorm_L2);
  double structureincconnorm_L2(0.0);
  structureincrement_->Norm2(&structureincconnorm_L2);
  double fluidincconnorm_L2(0.0);
  fluidincrement_->Norm2(&fluidincconnorm_L2);
  double aleincconnorm_L2(0.0);
  aleincrement_->Norm2(&aleincconnorm_L2);

  // set up vector of absolute concentrations
  Teuchos::RCP<Epetra_Vector> scatra = Teuchos::rcp(new Epetra_Vector(scatraincrement_->Map()));
  Teuchos::RCP<const Epetra_Vector> scatra1 = scatravec_[0]->ScaTraField()->Phinp();  // fluidscatra
  Teuchos::RCP<const Epetra_Vector> scatra2 =
      scatravec_[1]->ScaTraField()->Phinp();  // structurescatra
  setup_coupled_scatra_vector(scatra, scatra1, scatra2);

  // norms of solution vectors
  double scatranorm_L2(0.0);
  scatra->Norm2(&scatranorm_L2);
  if (scatranorm_L2 < 1e-05) scatranorm_L2 = 1.0;
  double structurenorm_L2(0.0);
  fsi_->structure_field()->Dispnp()->Norm2(&structurenorm_L2);
  if (structurenorm_L2 < 1e-05) structurenorm_L2 = 1.0;
  double fluidnorm_L2(0.0);
  fsi_->fluid_field()->Velnp()->Norm2(&fluidnorm_L2);
  if (fluidnorm_L2 < 1e-05) fluidnorm_L2 = 1.0;
  double alenorm_L2(0.0);
  fsi_->ale_field()->Dispnp()->Norm2(&alenorm_L2);
  if (alenorm_L2 < 1e-05) alenorm_L2 = 1.0;

  // print the incremental based convergence check to the screen
  if (Comm().MyPID() == 0)
  {
    std::cout << "\n************************************************************************\n"
                 "                       OUTER ITERATION STEP "
              << itnum << "/" << itmax << std::endl;

    printf(
        "+--------------+---------------------+----------------+---------------+---------------+---"
        "-----------+\n");
    printf(
        "|   step/max   |   tol      [norm]   |   scalar-inc   |   disp-inc    |   vel-inc     |   "
        "ale-inc    |\n");
    printf("|   %3d/%3d    |  %10.3E[L_2 ]   |  %10.3E    |  %10.3E   |  %10.3E   | %10.3E   |\n",
        itnum, itmax, ittol, scatraincconnorm_L2 / scatranorm_L2,
        structureincconnorm_L2 / structurenorm_L2, fluidincconnorm_L2 / fluidnorm_L2,
        aleincconnorm_L2 / alenorm_L2);
    printf(
        "+--------------+---------------------+----------------+---------------+---------------+---"
        "-----------+\n");
    std::cout << "************************************************************************\n";
  }

  if (scatraincconnorm_L2 / scatranorm_L2 <= ittol and
      structureincconnorm_L2 / structurenorm_L2 <= ittol and
      fluidincconnorm_L2 / fluidnorm_L2 <= ittol and
      aleincconnorm_L2 / alenorm_L2 <= ittol)  // converged!
  {
    stopnonliniter = true;
    if ((Comm().MyPID() == 0))
    {
      std::cout << "************************************************************************\n"
                   "                     OUTER ITERATION STEP CONVERGED"
                   "\n************************************************************************\n"
                << std::endl;
    }
  }
  else if (itnum == itmax)
  {
    stopnonliniter = true;
    if ((Comm().MyPID() == 0))
    {
      std::cout << "************************************************************************\n"
                   "           OUTER ITERATION STEP NOT CONVERGED IN ITEMAX STEPS"
                   "\n************************************************************************\n"
                << std::endl;
    }
    //    FOUR_C_THROW("The partitioned FS3I solver did not converge in ITEMAX steps!");
  }

  return stopnonliniter;
}

/*----------------------------------------------------------------------*
 | Compare if two doubles are relatively equal               Thon 08/15 |
 *----------------------------------------------------------------------*/
void FS3I::ACFSI::check_if_times_and_steps_and_dts_match()
{
  // NOTE: this check should pass each time after all prepare_time_step() have been called,
  // i.e. directly before a actual computation begins. I wish you luck :-)

  // check times
  const double fluidtime = fsi_->fluid_field()->Time();
  const double structuretime = fsi_->structure_field()->Time();
  const double aletime = fsi_->ale_field()->Time();
  const double fsitime = fsi_->Time();
  const double fluidscatratime = scatravec_[0]->ScaTraField()->Time();
  const double structurescatratime = scatravec_[1]->ScaTraField()->Time();

  if (not IsRealtiveEqualTo(fluidtime, time_, time_))
    FOUR_C_THROW("Your fluid time %f does not match the fs3i time %f!", fluidtime, time_);
  if (not IsRealtiveEqualTo(structuretime, time_, time_))
    FOUR_C_THROW("Your structure time %f does not match the fs3i time %f!", structuretime, time_);
  if (not IsRealtiveEqualTo(aletime, time_, time_))
    FOUR_C_THROW("Your ale time %f does not match the fs3i time %f!", aletime, time_);
  if (not IsRealtiveEqualTo(fsitime, time_, time_))
    FOUR_C_THROW("Your fsi time %f does not match the fs3i time %f!", fsitime, time_);
  if (not IsRealtiveEqualTo(fluidscatratime, time_, time_))
    FOUR_C_THROW(
        "Your fluid-scalar time %f does not match the fs3i time %f!", fluidscatratime, time_);
  if (not IsRealtiveEqualTo(structurescatratime, time_, time_))
    FOUR_C_THROW("Your structure-scalar time %f does not match the fs3i time %f!",
        structurescatratime, time_);

  // check steps
  const int fluidstep = fsi_->fluid_field()->Step();
  const int structurestep = fsi_->structure_field()->Step();
  const int alestep = fsi_->ale_field()->Step();
  const int fsistep = fsi_->Step();
  const int fluidscatrastep = scatravec_[0]->ScaTraField()->Step();
  const int structurescatrastep = scatravec_[1]->ScaTraField()->Step();

  if (not IsRealtiveEqualTo(fluidstep, step_, step_))
    FOUR_C_THROW("Your fluid step %i does not match the fs3i step %i!", fluidstep, step_);
  if (not IsRealtiveEqualTo(structurestep, step_, step_))
    FOUR_C_THROW("Your structure step %i does not match the fs3i step %i!", structurestep, step_);
  if (not IsRealtiveEqualTo(alestep, step_, step_))
    FOUR_C_THROW("Your ale step %i does not match the fs3i step %i!", alestep, step_);
  if (not IsRealtiveEqualTo(fsistep, step_, step_))
    FOUR_C_THROW("Your fsi step %i does not match the fs3i step %i!", fsistep, step_);
  if (not IsRealtiveEqualTo(fluidscatrastep, step_, step_))
    FOUR_C_THROW(
        "Your fluid-scalar step %i does not match the fs3i step %i!", fluidscatrastep, step_);
  if (not IsRealtiveEqualTo(structurescatrastep, step_, step_))
    FOUR_C_THROW("Your structure-scalar step %i does not match the fs3i step %i!",
        structurescatrastep, step_);

  // check dts
  const double fsiperssisteps =
      (double)(Global::Problem::Instance()->FS3IDynamicParams().sublist("AC").get<int>(
          "FSI_STEPS_PER_SCATRA_STEP"));
  const double fluiddt = fsi_->fluid_field()->Dt() * fsiperssisteps;
  const double structuredt = fsi_->structure_field()->Dt() * fsiperssisteps;
  const double aledt = fsi_->ale_field()->Dt() * fsiperssisteps;
  const double fsidt = fsi_->Dt() * fsiperssisteps;
  const double fluidscatradt = scatravec_[0]->ScaTraField()->Dt();
  const double structurescatradt = scatravec_[1]->ScaTraField()->Dt();

  if (not IsRealtiveEqualTo(fluiddt, dt_, 1.0))
    FOUR_C_THROW("Your fluid dt %f does not match the fs3i time %f!", fluiddt, dt_);
  if (not IsRealtiveEqualTo(structuredt, dt_, 1.0))
    FOUR_C_THROW("Your structure dt %f does not match the fs3i time %f!", structuredt, dt_);
  if (not IsRealtiveEqualTo(aledt, dt_, 1.0))
    FOUR_C_THROW("Your ale dt %f does not match the fs3i time %f!", aledt, dt_);
  if (not IsRealtiveEqualTo(fsidt, dt_, 1.0))
    FOUR_C_THROW("Your fsi dt %f does not match the fs3i time %f!", fsidt, dt_);
  if (not IsRealtiveEqualTo(fluidscatradt, dt_, 1.0))
    FOUR_C_THROW("Your fluid-scalar dt %f does not match the fs3i time %f!", fluidscatradt, dt_);
  if (not IsRealtiveEqualTo(structurescatradt, dt_, 1.0))
    FOUR_C_THROW(
        "Your structure-scalar dt %f does not match the fs3i time %f!", structurescatradt, dt_);
}

FOUR_C_NAMESPACE_CLOSE

/*----------------------------------------------------------------------*/
/*! \file


\brief cpp-file associated with algorithmic routines for two-way coupled partitioned
       solution approaches to fluid-structure-scalar-scalar interaction
       (FS3I). Specifically related version for multiscale approches. This file thereby holds
       all functions related with the large time scale simulation and
       the small to large to small time scale 'communication'.

\level 3


----------------------------------------------------------------------*/


#include "4C_adapter_ale_fsi.hpp"
#include "4C_adapter_fld_fluid_ac_fsi.hpp"
#include "4C_adapter_str_fsiwrapper.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fem_general_element.hpp"
#include "4C_fluid_utils_mapextractor.hpp"
#include "4C_fs3i_ac_fsi.hpp"
#include "4C_fsi_monolithic.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_material.hpp"
#include "4C_io.hpp"
#include "4C_io_control.hpp"
#include "4C_linalg_mapextractor.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_linear_solver_method_linalg.hpp"
#include "4C_mat_growth.hpp"
#include "4C_mat_growth_law.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_material_base.hpp"
#include "4C_scatra_algorithm.hpp"
#include "4C_scatra_timint_implicit.hpp"
#include "4C_structure_aux.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | timeloop for small time scales                            Thon 07/15 |
 *----------------------------------------------------------------------*/
void FS3I::ACFSI::LargeTimeScaleLoop()
{
  prepare_large_time_scale_loop();

  while (large_time_scale_loop_not_finished())
  {
    large_time_scale_prepare_time_step();

    large_time_scale_outer_loop();

    large_time_scale_update_and_output();
  }

  finish_large_time_scale_loop();
}

/*----------------------------------------------------------------------*
 | Prepare the large time scale loop                         Thon 08/15 |
 *----------------------------------------------------------------------*/
void FS3I::ACFSI::prepare_large_time_scale_loop()
{
  // print info
  if (Comm().MyPID() == 0)
  {
    std::cout << "\n************************************************************************"
                 "\n                         LARGE TIME SCALE LOOP"
                 "\n************************************************************************"
              << std::endl;
  }
  // Set large time scale time step in both scatra fields
  scatravec_[0]->ScaTraField()->set_dt(dt_large_);
  scatravec_[1]->ScaTraField()->set_dt(dt_large_);

  // set mean values in scatra fields
  large_time_scale_set_fsi_solution();

  // set back large time scale flags
  fsineedsupdate_ = false;
  growth_updates_counter_ = 0;

  // Save the phinp vector at the beginning of the large time scale loop in
  // in order to estimate the so far induced growth
  *structurephinp_blts_ = *scatravec_[1]->ScaTraField()->Phinp();
}

/*----------------------------------------------------------------------*
 |  Set mean wall shear stresses in scatra fields            Thon 11/15 |
 *----------------------------------------------------------------------*/
void FS3I::ACFSI::set_mean_wall_shear_stresses() const
{
  std::vector<Teuchos::RCP<const Epetra_Vector>> wss;

  // ############ Fluid Field ###############
  scatravec_[0]->ScaTraField()->set_wall_shear_stresses(FluidToFluidScalar(wall_shear_stress_lp_));

  // ############ Structure Field ###############

  // extract FSI-Interface from fluid field
  Teuchos::RCP<Epetra_Vector> WallShearStress =
      fsi_->fluid_field()->Interface()->ExtractFSICondVector(wall_shear_stress_lp_);

  // replace global fluid interface dofs through structure interface dofs
  WallShearStress = fsi_->fluid_to_struct(WallShearStress);

  // insert structure interface entries into vector with full structure length
  Teuchos::RCP<Epetra_Vector> structurewss =
      Core::LinAlg::CreateVector(*(fsi_->structure_field()->Interface()->FullMap()), true);

  // Parameter int block of function InsertVector: (0: inner dofs of structure, 1: interface dofs of
  // structure, 2: inner dofs of porofluid, 3: interface dofs of porofluid )
  fsi_->structure_field()->Interface()->InsertVector(WallShearStress, 1, structurewss);
  scatravec_[1]->ScaTraField()->set_wall_shear_stresses(
      structure_to_structure_scalar(structurewss));
}

/*----------------------------------------------------------------------*
 |  Set mean concentration of the fluid scatra field         Thon 11/15 |
 *----------------------------------------------------------------------*/
void FS3I::ACFSI::set_mean_fluid_scatra_concentration()
{
  Teuchos::RCP<const Epetra_Vector> MeanFluidConc = meanmanager_->GetMeanValue("mean_phi");

  scatravec_[0]->ScaTraField()->set_mean_concentration(MeanFluidConc);
}

/*----------------------------------------------------------------------*
 |  Set zero velocity field in scatra fields                 Thon 11/14 |
 *----------------------------------------------------------------------*/
void FS3I::ACFSI::set_zero_velocity_field()
{
  Teuchos::RCP<Epetra_Vector> zeros =
      Teuchos::rcp(new Epetra_Vector(fsi_->fluid_field()->Velnp()->Map(), true));
  scatravec_[0]->ScaTraField()->set_velocity_field(
      FluidToFluidScalar(zeros), Teuchos::null, FluidToFluidScalar(zeros), Teuchos::null);
  Teuchos::RCP<Epetra_Vector> zeros2 =
      Teuchos::rcp(new Epetra_Vector(fsi_->structure_field()->Velnp()->Map(), true));
  scatravec_[1]->ScaTraField()->set_velocity_field(structure_to_structure_scalar(zeros2),
      Teuchos::null, structure_to_structure_scalar(zeros2), Teuchos::null);
}

/*-------------------------------------------------------------------------------*
 | Evaluate surface permeability condition for struct scatra field    Thon 08/15 |
 *-------------------------------------------------------------------------------*/
void FS3I::ACFSI::evaluateith_scatra_surface_permeability(const int i  // id of scalar to evaluate
)
{
  // Note: 0 corresponds to fluid-scatra
  //      1 corresponds to structure-scatra

  //----------------------------------------------------------------------
  // set membrane concentrations
  //----------------------------------------------------------------------
  set_membrane_concentration();

  //----------------------------------------------------------------------
  // evaluate simplified kedem-katchalsy condtion
  //----------------------------------------------------------------------
  Teuchos::RCP<Epetra_Vector> rhs_scal = scatracoupforce_[i];
  Teuchos::RCP<Core::LinAlg::SparseMatrix> mat_scal = scatracoupmat_[i];

  rhs_scal->PutScalar(0.0);
  mat_scal->Zero();

  scatravec_[i]->ScaTraField()->SurfacePermeability(mat_scal, rhs_scal);

  // apply Dirichlet boundary conditions to coupling matrix and vector
  const Teuchos::RCP<const Epetra_Map> dbcmap =
      scatravec_[i]->ScaTraField()->DirichMaps()->CondMap();
  mat_scal->ApplyDirichlet(*dbcmap, false);
  Core::LinAlg::apply_dirichlet_to_system(*rhs_scal, *scatrazeros_[i], *dbcmap);
}

/*----------------------------------------------------------------------*
 | Finish the large time scale loop                          Thon 08/15 |
 *----------------------------------------------------------------------*/
void FS3I::ACFSI::finish_large_time_scale_loop()
{
  // Set small time scale time step size
  scatravec_[0]->ScaTraField()->set_dt(dt_);
  scatravec_[1]->ScaTraField()->set_dt(dt_);

  // Fix time and step in fsi and fluid scatra field

  // We start the small time scale with a new cycle. But since dt_large is a
  // multiple of fsiperiod_ we are already at the this time_
  // We do not modify the step_ counter; we just keep counting..

  double tmp = fmod(time_, dt_large_);
  tmp = tmp - fmod(tmp, 10.0 * fsiperiod_) + 10.0 * fsiperiod_;
  time_ = tmp;

  SetTimeAndStepInFSI(time_, step_);
  scatravec_[0]->ScaTraField()->SetTimeStep(time_, step_);
  scatravec_[1]->ScaTraField()->SetTimeStep(time_, step_);

  // we now have to fix the time_ and step_ of the structure field, since this is not shifted
  // in prepare_time_step(), but in update(), which we here will not call. So..
  fsi_->structure_field()->set_time(time_);
  fsi_->structure_field()->SetTimen(time_ + fsi_->fluid_field()->Dt());
  fsi_->structure_field()->SetStep(step_);
  fsi_->structure_field()->SetStepn(step_ + 1);

  // we start with a clean small time scale loop
  fsiisperiodic_ = false;
  scatraisperiodic_ = false;


  // NOTE: we start a new output file since paraview does only read floating point numbers.
  // Hence the upcoming small time scale calculation may not be displayable in paraview due
  // to the large variety in the time scales. Bad thing :(
  //*-------------------------------------------------------------------------------*
  // | create new output file
  //*-------------------------------------------------------------------------------*/
  Teuchos::RCP<Core::IO::DiscretizationWriter> output_writer =
      Global::Problem::Instance()->GetDis("structure")->Writer();
  output_writer->new_result_file(step_);
  // and write all meshes
  output_writer->create_new_result_and_mesh_file();
  output_writer->write_mesh(0, 0.0);
  output_writer = Global::Problem::Instance()->GetDis("fluid")->Writer();
  output_writer->create_new_result_and_mesh_file();
  output_writer->write_mesh(0, 0.0);
  output_writer = Global::Problem::Instance()->GetDis("ale")->Writer();
  output_writer->create_new_result_and_mesh_file();
  output_writer->write_mesh(0, 0.0);
  output_writer = Global::Problem::Instance()->GetDis("scatra1")->Writer();
  output_writer->create_new_result_and_mesh_file();
  output_writer->write_mesh(0, 0.0);
  output_writer = Global::Problem::Instance()->GetDis("scatra2")->Writer();
  output_writer->create_new_result_and_mesh_file();
  output_writer->write_mesh(0, 0.0);

  // write outputs in new file
  constexpr bool force_prepare = false;
  fsi_->prepare_output(force_prepare);

  FsiOutput();
  ScatraOutput();
}

/*----------------------------------------------------------------------*
 | timeloop for large time scales                            Thon 07/15 |
 *----------------------------------------------------------------------*/
bool FS3I::ACFSI::large_time_scale_loop_not_finished()
{
  return NotFinished() and not fsineedsupdate_;
}

/*----------------------------------------------------------------------*
 | Prepare small time scale time step                        Thon 07/15 |
 *----------------------------------------------------------------------*/
void FS3I::ACFSI::large_time_scale_prepare_time_step()
{
  // Set large time scale time step in both scatra fields
  scatravec_[0]->ScaTraField()->set_dt(dt_large_);
  scatravec_[1]->ScaTraField()->set_dt(dt_large_);

  // Increment time and step
  step_ += 1;
  time_ += dt_large_;

  // Print to screen
  if (Comm().MyPID() == 0)
  {
    std::cout << "\n\n"
              << "TIME:  " << std::scientific << std::setprecision(12) << time_ << "/"
              << std::setprecision(4) << timemax_ << "     DT = " << std::scientific << dt_large_
              << "     STEP = " << std::setw(4) << step_ << "/" << std::setw(4) << numstep_ << "\n";
  }

  // prepare structure scatra field
  scatravec_[1]->ScaTraField()->prepare_time_step();
}

/*----------------------------------------------------------------------*
 | outer_loop for sequentially staggered FS3I scheme          Thon 08/15 |
 *----------------------------------------------------------------------*/
void FS3I::ACFSI::large_time_scale_outer_loop()
{
  DoStructScatraStep();

  if (does_growth_needs_update())  // includes the check for fsineedsupdate_
  {
    large_time_scale_do_growth_update();
  }
}

/*----------------------------------------------------------------------*
 | Do a large time scale structe scatra step                 Thon 08/15 |
 *----------------------------------------------------------------------*/
void FS3I::ACFSI::DoStructScatraStep()
{
  if (Comm().MyPID() == 0)
  {
    std::cout << "\n************************************************************************"
                 "\n                       AC STRUCTURE SCATRA SOLVER"
                 "\n************************************************************************\n"
              << std::endl;

    std::cout << "+- step/max -+-- scal-res/ abs-tol [norm] -+-- scal-inc/ rel-tol [norm] -+"
              << std::endl;
  }

  bool stopnonliniter = false;
  int itnum = 0;

  while (stopnonliniter == false)
  {
    struct_scatra_evaluate_solve_iter_update();
    itnum++;
    if (struct_scatra_convergence_check(itnum)) break;
  }
}

/*--------------------------------------------------------------------------------*
 | evaluate, solver and iteratively update structure scalar problem    Thon 08/15 |
 *--------------------------------------------------------------------------------*/
void FS3I::ACFSI::struct_scatra_evaluate_solve_iter_update()
{
  if (infperm_) FOUR_C_THROW("This not a valid option!");  // just for safety

  const Teuchos::RCP<ScaTra::ScaTraTimIntImpl> scatra =
      scatravec_[1]->ScaTraField();  // structure scatra

  //----------------------------------------------------------------------
  // evaluate the structure scatra field
  //----------------------------------------------------------------------
  scatra->PrepareLinearSolve();

  //----------------------------------------------------------------------
  // calculate contributions due to finite interface permeability
  //----------------------------------------------------------------------
  evaluateith_scatra_surface_permeability(1);

  //----------------------------------------------------------------------
  // recalculate fluid scatra contributions due to possible changed time step size
  // and the using of mean wss and mean phi for the fluid scatra field
  //----------------------------------------------------------------------
  evaluateith_scatra_surface_permeability(0);

  //----------------------------------------------------------------------
  // add coupling to the resiudal
  //----------------------------------------------------------------------
  const Teuchos::RCP<Epetra_Vector> rhs_struct_scal = scatracoupforce_[1];
  const Teuchos::RCP<Core::LinAlg::SparseMatrix> mat_struct_scal = scatracoupmat_[1];
  const Teuchos::RCP<Epetra_Vector> residual = scatra->Residual();

  residual->Update(1.0, *rhs_struct_scal, 1.0);

  // add contribution of the fluid field
  Teuchos::RCP<Epetra_Vector> rhs_fluid_scal_boundary =
      scatrafieldexvec_[0]->ExtractVector(scatracoupforce_[0], 1);
  Teuchos::RCP<Epetra_Vector> rhs_fluid_scal =
      scatrafieldexvec_[1]->InsertVector(Scatra1ToScatra2(rhs_fluid_scal_boundary), 1);

  residual->Update(-1.0, *rhs_fluid_scal, 1.0);

  //----------------------------------------------------------------------
  // add coupling to the sysmat
  //----------------------------------------------------------------------
  const Teuchos::RCP<Core::LinAlg::SparseMatrix> sysmat = scatra->SystemMatrix();
  sysmat->Add(*mat_struct_scal, false, 1.0, 1.0);

  //----------------------------------------------------------------------
  // solve the scatra problem
  //----------------------------------------------------------------------
  const Teuchos::RCP<Epetra_Vector> structurescatraincrement =
      Core::LinAlg::CreateVector(*scatra->dof_row_map(), true);

  Core::LinAlg::SolverParams solver_params;
  solver_params.refactor = true;
  solver_params.reset = true;
  scatra->Solver()->Solve(
      sysmat->EpetraOperator(), structurescatraincrement, residual, solver_params);

  //----------------------------------------------------------------------
  // update the strucutre scatra increment
  //----------------------------------------------------------------------
  scatra->UpdateIter(structurescatraincrement);
}

/*----------------------------------------------------------------------*
 | check convergence of structure scatra field               Thon 08/15 |
 *----------------------------------------------------------------------*/
bool FS3I::ACFSI::struct_scatra_convergence_check(const int itnum)
{
  const Teuchos::RCP<ScaTra::ScaTraTimIntImpl> scatra =
      scatravec_[1]->ScaTraField();  // structure scatra

  // some input parameters for the scatra fields
  const Teuchos::ParameterList& scatradyn =
      Global::Problem::Instance()->scalar_transport_dynamic_params();
  const int scatraitemax = scatradyn.sublist("NONLINEAR").get<int>("ITEMAX");
  const double scatraittol = scatradyn.sublist("NONLINEAR").get<double>("CONVTOL");
  const double scatraabstolres = scatradyn.sublist("NONLINEAR").get<double>("ABSTOLRES");


  double conresnorm(0.0);
  scatra->Residual()->Norm2(&conresnorm);
  double incconnorm(0.0);
  scatra->Increment()->Norm2(&incconnorm);
  double phinpnorm(0.0);
  scatra->Phinp()->Norm2(&phinpnorm);

  // care for the case that nothing really happens in the concentration field
  if (phinpnorm < 1e-5) phinpnorm = 1.0;

  // print the screen info
  if (Comm().MyPID() == 0)
  {
    printf("|   %3d/%3d  |  %1.3E/ %1.1E [L_2 ]  |  %1.3E/ %1.1E [L_2 ]  |\n", itnum, scatraitemax,
        conresnorm, scatraabstolres, incconnorm / phinpnorm, scatraittol);
  }

  // this is the convergence check
  // We always require at least one solve. We test the L_2-norm of the
  // current residual. Norm of residual is just printed for information
  if (conresnorm <= scatraabstolres and incconnorm / phinpnorm <= scatraittol)
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
    //    FOUR_C_THROW("Structure scatra not converged in itemax steps!");
    return true;
  }
  else
    return false;
}

/*----------------------------------------------------------------------*
 | Do we need to update the structure scatra displacments               |
 | due to growth                                             Thon 08/15 |
 *----------------------------------------------------------------------*/
bool FS3I::ACFSI::does_growth_needs_update()
{
  bool growthneedsupdate = false;

  // check if the structure material is a growth material. We assume here
  // that the structure has the same material for the whole discretiazation.
  // Hence we check only the first element:
  Teuchos::RCP<Core::FE::Discretization> structuredis = fsi_->structure_field()->discretization();
  const int GID = structuredis->ElementColMap()->GID(0);  // global element ID

  Teuchos::RCP<Core::Mat::Material> structurematerial = structuredis->gElement(GID)->Material();

  if (structurematerial->MaterialType() != Core::Materials::m_growth_volumetric)
  {
    FOUR_C_THROW("In AC-FS3I we want growth, so use a growth material like MAT_GrowthVolumetric!");
  }
  else
  {
    //----------------------------------------------------------------------------------------------------
    // get alpha and growth inducing scalar
    //----------------------------------------------------------------------------------------------------
    double alpha = 0.0;
    int sc1 = 1;

    Teuchos::RCP<Mat::GrowthVolumetric> growthmaterial =
        Teuchos::rcp_dynamic_cast<Mat::GrowthVolumetric>(structurematerial);

    if (growthmaterial == Teuchos::null)
      FOUR_C_THROW("Dynamic cast to Mat::GrowthVolumetric failed!");

    Teuchos::RCP<Mat::GrowthLaw> growthlaw = growthmaterial->Parameter()->growthlaw_;

    switch (growthlaw->MaterialType())
    {
      case Core::Materials::m_growth_ac:
      {
        Teuchos::RCP<Mat::GrowthLawAC> growthlawac =
            Teuchos::rcp_dynamic_cast<Mat::GrowthLawAC>(growthlaw);
        if (growthmaterial == Teuchos::null)
          FOUR_C_THROW("Dynamic cast to Mat::GrowthLawAC failed!");
        alpha = growthlawac->Parameter()->alpha_;
        sc1 = growthlawac->Parameter()->Sc1_;
        break;
      }
      case Core::Materials::m_growth_ac_radial:
      {
        Teuchos::RCP<Mat::GrowthLawACRadial> growthlawacradial =
            Teuchos::rcp_dynamic_cast<Mat::GrowthLawACRadial>(growthlaw);
        if (growthlawacradial == Teuchos::null)
          FOUR_C_THROW("Dynamic cast to Mat::GrowthLawACRadial failed!");
        alpha = growthlawacradial->Parameter()->alpha_;
        sc1 = growthlawacradial->Parameter()->Sc1_;
        break;
      }
      case Core::Materials::m_growth_ac_radial_refconc:
      {
        Teuchos::RCP<Mat::GrowthLawACRadialRefConc> growthlawacradialrefconc =
            Teuchos::rcp_dynamic_cast<Mat::GrowthLawACRadialRefConc>(growthlaw);
        if (growthlawacradialrefconc == Teuchos::null)
          FOUR_C_THROW("Dynamic cast to Mat::GrowthLawACRadialRefConc failed!");
        alpha = growthlawacradialrefconc->Parameter()->alpha_;
        sc1 = growthlawacradialrefconc->Parameter()->Sc1_;
        break;
      }
      default:
      {
        FOUR_C_THROW("Growth law not supported in AC-FS3I!");
        break;
      }
    }
    // Puh! That was exhausting. But we have to keep going.

    //----------------------------------------------------------------------------------------------------
    // get the approx. increase of volume due to growth since the beginning of the large time scale
    // loop
    //----------------------------------------------------------------------------------------------------
    const Teuchos::RCP<ScaTra::ScaTraTimIntImpl> scatra =
        scatravec_[1]->ScaTraField();                                 // structure scatra
    const Teuchos::RCP<const Epetra_Vector> phinp = scatra->Phinp();  // fluidscatra

    // build difference vector with the reference
    const Teuchos::RCP<Epetra_Vector> phidiff_bltsl_ =
        Core::LinAlg::CreateVector(*scatra->dof_row_map(), true);
    phidiff_bltsl_->Update(1.0, *phinp, -1.0, *structurephinp_blts_, 0.0);

    // Extract the dof of interest
    Teuchos::RCP<Epetra_Vector> phidiff_bltsl_j =
        extractjthstructscalar_[sc1 - 1]->ExtractCondVector(phidiff_bltsl_);

    // get the maximum
    double max_phidiff_bltsl = 0.0;
    phidiff_bltsl_j->MaxValue(&max_phidiff_bltsl);

    //----------------------------------------------------------------------------------------------------
    // screen output
    //----------------------------------------------------------------------------------------------------
    const int growth_updates =
        Global::Problem::Instance()->FS3IDynamicParams().sublist("AC").get<int>("GROWTH_UPDATES");
    const double fsi_update_tol =
        Global::Problem::Instance()->FS3IDynamicParams().sublist("AC").get<double>(
            "FSI_UPDATE_TOL");

    if (Comm().MyPID() == 0)
      std::cout << std::scientific << std::setprecision(3)
                << "The maximal relative local growth since the small time scale is "
                << alpha * max_phidiff_bltsl << " (tol "
                << ((double)growth_updates_counter_ + 1.0) / (double)growth_updates * fsi_update_tol
                << ", iter " << growth_updates_counter_ << "/" << growth_updates << ")"
                << std::endl;

    // some safety check
    if (growth_updates_counter_ > growth_updates)
      FOUR_C_THROW("It should not be possible to have done so much growth updates. Sorry!");

    //----------------------------------------------------------------------------------------------------
    // now the actual comparison
    //----------------------------------------------------------------------------------------------------
    // do we need a growth update?
    if (max_phidiff_bltsl * alpha >=
        ((double)growth_updates_counter_ + 1.0) / (double)growth_updates * fsi_update_tol)
    {
      growthneedsupdate = true;
    }

    // are we done with the current large time scale loop?
    if (max_phidiff_bltsl * alpha >= fsi_update_tol)
    {
      fsineedsupdate_ = true;
    }
  }

  return growthneedsupdate;
}

/*-------------------------------------------------------------------------*
 | update the structure scatra displacments due to growth       Thon 08/15 |
 *-------------------------------------------------------------------------*/
void FS3I::ACFSI::large_time_scale_do_growth_update()
{
  const int growth_updates =
      Global::Problem::Instance()->FS3IDynamicParams().sublist("AC").get<int>("GROWTH_UPDATES");

  const Teuchos::RCP<ScaTra::ScaTraTimIntImpl> fluidscatra = scatravec_[0]->ScaTraField();
  const Teuchos::RCP<ScaTra::ScaTraTimIntImpl> structurescatra = scatravec_[1]->ScaTraField();

  // Note: we never do never proceed with time_ and step_, so this really just about updating the
  // growth, i.e. the displacements of the structure scatra fields

  //----------------------------------------------------------------------
  // print to screen
  //----------------------------------------------------------------------
  if (Comm().MyPID() == 0)
  {
    std::cout << "\n************************************************************************"
                 "\n                         AC GROWTH UPDATE "
              << growth_updates_counter_ + 1 << "/" << growth_updates
              << "\n************************************************************************"
              << std::endl;
  }


  //----------------------------------------------------------------------
  // finish present structure scatra time step (no output)
  //----------------------------------------------------------------------
  structurescatra->update();

  //----------------------------------------------------------------------
  // Switch time step of scatra fields
  //----------------------------------------------------------------------
  // Switch back the time step to do the update with the same (small) timestep as the fsi
  // (subcycling time step possible!)
  fluidscatra->set_dt(dt_);
  structurescatra->set_dt(dt_);

  //----------------------------------------------------------------------
  // Fix time_ and step_ counters
  //----------------------------------------------------------------------
  // time_+=dt_;

  SetTimeAndStepInFSI(time_ - dt_, step_ - 1);
  fluidscatra->SetTimeStep(time_ - dt_, step_ - 1);
  structurescatra->SetTimeStep(time_ - dt_, step_ - 1);

  // we now have to fix the time_ and step_ of the structure field, since this is not shifted
  // in prepare_time_step(), but in update(), which we here will not call. So..
  fsi_->structure_field()->set_time(time_ - dt_);
  fsi_->structure_field()->SetTimen(time_);
  fsi_->structure_field()->SetStep(step_ - 1);
  fsi_->structure_field()->SetStepn(step_);

  //----------------------------------------------------------------------
  // Prepare time steps
  //----------------------------------------------------------------------
  // fsi problem
  set_struct_scatra_solution();
  fsi_->prepare_time_step();
  // scatra fields
  fluidscatra->prepare_time_step();
  structurescatra->prepare_time_step();

  //----------------------------------------------------------------------
  // do the growth update
  //----------------------------------------------------------------------
  // Safety check:
  check_if_times_and_steps_and_dts_match();

  // the actual calculations
  large_time_scale_outer_loop_iter_stagg();

  //----------------------------------------------------------------------
  // write the output
  //----------------------------------------------------------------------
  // write fsi output. Scatra outputs are done later
  // fsi output
  constexpr bool force_prepare = false;
  fsi_->prepare_output(force_prepare);
  // NOTE: we have to call this functions, otherwise the structure displacements are not applied
  fsi_->update();
  FsiOutput();
  // fluid scatra update. Structure scatra is done later
  fluidscatra->update();
  fluidscatra->check_and_write_output_and_restart();

  //----------------------------------------------------------------------
  // Switch back time steps and set mean values in scatra fields
  //----------------------------------------------------------------------
  // Now set the time step back:
  fluidscatra->set_dt(dt_large_);
  structurescatra->set_dt(dt_large_);

  // set mean values in scatra fields
  large_time_scale_set_fsi_solution();

  //----------------------------------------------------------------------
  // higher growth counter
  //----------------------------------------------------------------------
  growth_updates_counter_++;
}

/*-------------------------------------------------------------------------------*
 | outer_loop for large time scale iterative staggered FS3I scheme     Thon 11/15 |
 *-------------------------------------------------------------------------------*/
void FS3I::ACFSI::large_time_scale_outer_loop_iter_stagg()
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

    DoFSIStepStandard();
    // subcycling is not allowed, since we use this function for the growth update. Nevertheless it
    // should work.. periodical repetition is not allowed, since we want to converge the problems

    large_time_scale_set_fsi_solution();

    small_time_scale_do_scatra_step();

    stopnonliniter = part_fs3i_convergence_ckeck(itnum);
  }
}

/*-----------------------------------------------------------------------*
 | set mean FSI values in scatra fields                       Thon 11/15 |
 *---------------------------------------------------- ------------------*/
void FS3I::ACFSI::large_time_scale_set_fsi_solution()
{
  // we clear every state, including the states of the secondary dof sets
  for (unsigned i = 0; i < scatravec_.size(); ++i)
  {
    scatravec_[i]->ScaTraField()->discretization()->ClearState(true);
    // we have to manually clear this since this can not be saved directly in the
    // primary dof set (because it is cleared in between)
    scatravec_[i]->ScaTraField()->clear_external_concentrations();
  }

  set_mesh_disp();
  set_mean_wall_shear_stresses();
  set_mean_fluid_scatra_concentration();
  set_membrane_concentration();
  // Set zeros velocities since we assume that the large time scale can not see the deformation of
  // the small time scale
  set_zero_velocity_field();
}

/*----------------------------------------------------------------------*
 | Update and output the large time scale                    Thon 08/15 |
 *----------------------------------------------------------------------*/
void FS3I::ACFSI::large_time_scale_update_and_output()
{
  // keep fsi time and fluid scatra field up to date
  SetTimeAndStepInFSI(time_, step_);
  scatravec_[0]->ScaTraField()->SetTimeStep(time_, step_);

  // NOTE: fsi output is already updated and written in large_time_scale_do_growth_update()
  // NOTE: fluid scatra is already updated and written in large_time_scale_do_growth_update()

  // now update and output the structure scatra field
  scatravec_[1]->ScaTraField()->update();
  scatravec_[1]->ScaTraField()->check_and_write_output_and_restart();
}

/*----------------------------------------------------------------------*
 | Build map extractor which extracts the j-th dof           Thon 08/15 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<Core::LinAlg::MapExtractor>> FS3I::ACFSI::BuildMapExtractor()
{
  std::vector<Teuchos::RCP<Core::LinAlg::MapExtractor>> extractjthscalar;

  const Teuchos::RCP<ScaTra::ScaTraTimIntImpl> scatra =
      scatravec_[1]->ScaTraField();  // structure scatra
  const int numscal = scatra->NumScal();
  const Teuchos::RCP<const Core::FE::Discretization> dis = scatra->discretization();

  for (int k = 0; k < numscal; k++)
  {
    std::set<int> conddofset;
    std::set<int> otherdofset;

    int numrownodes = dis->NumMyRowNodes();
    for (int i = 0; i < numrownodes; ++i)
    {
      Core::Nodes::Node* node = dis->lRowNode(i);

      std::vector<int> dof = dis->Dof(0, node);
      if (dof.size() != (unsigned)scatravec_[1]->ScaTraField()->NumScal())
        FOUR_C_THROW("There was some error building the Map Extractor!");
      for (unsigned j = 0; j < dof.size(); ++j)
      {
        // test for dof position
        if (j != static_cast<unsigned>(k))
        {
          otherdofset.insert(dof[j]);
        }
        else
        {
          conddofset.insert(dof[j]);
        }
      }
    }
    std::vector<int> conddofmapvec;
    conddofmapvec.reserve(conddofset.size());
    conddofmapvec.assign(conddofset.begin(), conddofset.end());
    conddofset.clear();
    Teuchos::RCP<Epetra_Map> conddofmap = Teuchos::rcp(
        new Epetra_Map(-1, conddofmapvec.size(), conddofmapvec.data(), 0, dis->Comm()));
    conddofmapvec.clear();

    std::vector<int> otherdofmapvec;
    otherdofmapvec.reserve(otherdofset.size());
    otherdofmapvec.assign(otherdofset.begin(), otherdofset.end());
    otherdofset.clear();
    Teuchos::RCP<Epetra_Map> otherdofmap = Teuchos::rcp(
        new Epetra_Map(-1, otherdofmapvec.size(), otherdofmapvec.data(), 0, dis->Comm()));
    otherdofmapvec.clear();

    Teuchos::RCP<Core::LinAlg::MapExtractor> getjdof = Teuchos::rcp(new Core::LinAlg::MapExtractor);
    getjdof->setup(*dis->dof_row_map(), conddofmap, otherdofmap);
    extractjthscalar.push_back(getjdof);
  }

  return extractjthscalar;
}

/*----------------------------------------------------------------------*
 | Compare if two doubles are relatively equal               Thon 08/15 |
 *----------------------------------------------------------------------*/
bool FS3I::ACFSI::IsRealtiveEqualTo(const double A, const double B, const double Ref)
{
  return ((fabs(A - B) / Ref) < 1e-12);
}

/*----------------------------------------------------------------------*
 | Compare if A mod B is relatively equal to zero            Thon 08/15 |
 *----------------------------------------------------------------------*/
bool FS3I::ACFSI::modulo_is_realtive_zero(const double value, const double modulo, const double Ref)
{
  return IsRealtiveEqualTo(fmod(value + modulo / 2, modulo) - modulo / 2, 0.0, Ref);
}

/*----------------------------------------------------------------------*
 | Compare if A mod B is relatively equal to zero            Thon 10/15 |
 *----------------------------------------------------------------------*/
FS3I::MeanManager::MeanManager(
    const Epetra_Map& wssmap, const Epetra_Map& phimap, const Epetra_Map& pressuremap)
    : sum_wss_(Core::LinAlg::CreateVector(wssmap, true)),
      sum_phi_(Core::LinAlg::CreateVector(phimap, true)),
      sum_pres_(Core::LinAlg::CreateVector(pressuremap, true)),
      sum_dt_wss_(0.0),
      sum_dt_phi_(0.0),
      sum_dt_pres_(0.0)
{
}


/*----------------------------------------------------------------------*
 | add value into the mean manager                           Thon 10/15 |
 *----------------------------------------------------------------------*/
void FS3I::MeanManager::AddValue(
    const std::string type, const Teuchos::RCP<const Epetra_Vector> value, const double dt)
{
  if (type == "wss")
  {
#ifdef FOUR_C_DEBUG
    // check, whether maps are the same
    if (not value->Map().PointSameAs(sum_wss_->Map()))
    {
      FOUR_C_THROW("Maps do not match, but they have to.");
    }
#endif

    sum_wss_->Update(dt, *value, 1.0);  // weighted sum of all prior stresses
    sum_dt_wss_ += dt;
  }
  else if (type == "phi")
  {
#ifdef FOUR_C_DEBUG
    // check, whether maps are the same
    if (not value->Map().PointSameAs(sum_phi_->Map()))
    {
      FOUR_C_THROW("Maps do not match, but they have to.");
    }
#endif

    sum_phi_->Update(dt, *value, 1.0);  // weighted sum of all prior stresses
    sum_dt_phi_ += dt;
  }
  else if (type == "pressure")
  {
#ifdef FOUR_C_DEBUG
    // check, whether maps are the same
    if (not value->Map().PointSameAs(sum_pres_->Map()))
    {
      FOUR_C_THROW("Maps do not match, but they have to.");
    }
#endif

    sum_pres_->Update(dt, *value, 1.0);  // weighted sum of all prior stresses
    sum_dt_pres_ += dt;
  }
  else
    FOUR_C_THROW("Mean Manager does not support the given value '%s'.", type.c_str());

  return;
}

/*----------------------------------------------------------------------*
 | reset mean manager                                        Thon 10/15 |
 *----------------------------------------------------------------------*/
void FS3I::MeanManager::reset()
{
  // first some checking
  if (abs(sum_dt_wss_ - sum_dt_phi_) > 1e-14 or abs(sum_dt_wss_ - sum_dt_pres_) > 1e-14)
    FOUR_C_THROW("The time ranges you did mean over do not match!");

  sum_wss_->PutScalar(0.0);
  sum_dt_wss_ = 0.0;
  sum_phi_->PutScalar(0.0);
  sum_dt_phi_ = 0.0;
  sum_pres_->PutScalar(0.0);
  sum_dt_pres_ = 0.0;
}

/*----------------------------------------------------------------------*
 | get some mean value                                       Thon 10/15 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> FS3I::MeanManager::GetMeanValue(const std::string type) const
{
  Teuchos::RCP<Epetra_Vector> meanvector;

  if (type == "mean_wss")
  {
    meanvector = Teuchos::rcp(new Epetra_Vector(sum_wss_->Map(), true));

    if (sum_dt_wss_ > 1e-12)  // iff we have actually calculated some mean wss
      meanvector->Update(1.0 / sum_dt_wss_, *sum_wss_, 0.0);  // weighted sum of all prior stresses
    else
    {
      double norm = 0.0;
      meanvector->NormInf(&norm);
      if (norm > 1e-12)
        FOUR_C_THROW("SumDtWss_ is zero, but SumWss_ not.. Something is terribly wrong!");
    }
  }
  else if (type == "osi")
  {
    FOUR_C_THROW("Oscillatory shear index is yet not supported!");
  }
  else if (type == "mean_phi")
  {
    meanvector = Teuchos::rcp(new Epetra_Vector(sum_phi_->Map(), true));

    if (sum_dt_phi_ > 1e-12)  // iff we have actually calculated some mean wss
      meanvector->Update(1.0 / sum_dt_phi_, *sum_phi_, 0.0);  // weighted sum of all prior stresses
    else
    {
      double norm = 0.0;
      meanvector->NormInf(&norm);
      if (norm > 1e-12)
        FOUR_C_THROW("SumDtPhi_ is zero, but SumPhi_ not.. Something is terribly wrong!");
    }
  }
  else if (type == "mean_pressure")
  {
    meanvector = Teuchos::rcp(new Epetra_Vector(sum_pres_->Map(), true));

    if (sum_dt_pres_ > 1e-12)  // iff we have actually calculated some mean wss
      meanvector->Update(
          1.0 / sum_dt_pres_, *sum_pres_, 0.0);  // weighted sum of all prior stresses
    else
    {
      double norm = 0.0;
      meanvector->NormInf(&norm);
      if (norm > 1e-12)
        FOUR_C_THROW("SumDtPres_ is zero, but SumPres_ not.. Something is terribly wrong!");
    }
  }
  else
    FOUR_C_THROW("Mean Manager does not support the given value '%s'.", type.c_str());

  return meanvector;
}

/*----------------------------------------------------------------------*
 | Write restart of mean manager                             Thon 10/15 |
 *----------------------------------------------------------------------*/
void FS3I::MeanManager::write_restart(
    Teuchos::RCP<Core::IO::DiscretizationWriter> fluidwriter) const
{
  // first some checking
  if (abs(sum_dt_wss_ - sum_dt_phi_) > 1e-14 or abs(sum_dt_wss_ - sum_dt_pres_) > 1e-14)
    FOUR_C_THROW("The time ranges you did mean over do not match!");

  // write all values
  fluidwriter->write_vector("SumWss", sum_wss_);
  fluidwriter->write_vector("SumPhi", sum_phi_);
  //  fluidwriter->write_vector("SumPres", SumPres_);
  // we need only one SumDt since they are all the same
  fluidwriter->write_double("SumDtWss", sum_dt_wss_);
}

/*----------------------------------------------------------------------*
 | Read restart of mean manager                             Thon 10/15 |
 *----------------------------------------------------------------------*/
void FS3I::MeanManager::read_restart(Core::IO::DiscretizationReader& fluidreader)
{
  // read all values...
  fluidreader.read_vector(sum_wss_, "SumWss");
  fluidreader.read_vector(sum_phi_, "SumPhi");
  //  fluidreader.read_vector(SumPres_, "SumPres");
  sum_dt_wss_ = fluidreader.read_double("SumDtWss");
  //...and recover the rest
  sum_dt_phi_ = sum_dt_wss_;
  sum_dt_pres_ = sum_dt_wss_;
}

FOUR_C_NAMESPACE_CLOSE

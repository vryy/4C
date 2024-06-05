/*----------------------------------------------------------------------*/
/*! \file
 \brief one way coupled partitioned scalar structure interaction

 \level 2


 *------------------------------------------------------------------------------------------------*/

#include "4C_ssi_partitioned_1wc.hpp"

#include "4C_adapter_scatra_base_algorithm.hpp"
#include "4C_adapter_str_ssiwrapper.hpp"
#include "4C_adapter_str_wrapper.hpp"
#include "4C_global_data.hpp"
#include "4C_io.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_scatra_timint_cardiac_monodomain.hpp"
#include "4C_scatra_timint_implicit.hpp"

FOUR_C_NAMESPACE_OPEN

SSI::SSIPart1WC::SSIPart1WC(const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams)
    : SSIPart(comm, globaltimeparams), isscatrafromfile_(false)
{
  // Keep this constructor empty!
  // First do everything on the more basic objects like the discretizations, like e.g.
  // redistribution of elements. Only then call the setup to this class. This will call he setup to
  // all classes in the inheritance hierarchy. This way, this class may also override a method that
  // is called during Setup() in a base class.
}

/*----------------------------------------------------------------------*
 | Setup this class                                         rauch 08/16 |
 *----------------------------------------------------------------------*/
void SSI::SSIPart1WC::Init(const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams,
    const Teuchos::ParameterList& scatraparams, const Teuchos::ParameterList& structparams,
    const std::string& struct_disname, const std::string& scatra_disname, bool isAle)
{
  // call setup of base class
  SSI::SSIPart::Init(
      comm, globaltimeparams, scatraparams, structparams, struct_disname, scatra_disname, isAle);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSIPart1WC::do_struct_step()
{
  if (Comm().MyPID() == 0)
  {
    std::cout << "\n***********************\n STRUCTURE SOLVER \n***********************\n";
  }

  // Newton-Raphson iteration
  structure_field()->Solve();
  // calculate stresses, strains, energies
  constexpr bool force_prepare = false;
  structure_field()->prepare_output(force_prepare);
  // update all single field solvers
  structure_field()->Update();
  // write output to files
  structure_field()->Output();
  // write output to screen
  structure_field()->PrintStep();
  // clean up
  structure_field()->discretization()->ClearState(true);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSIPart1WC::do_scatra_step()
{
  if (Comm().MyPID() == 0)
  {
    std::cout << "\n***********************\n TRANSPORT SOLVER \n***********************\n";
  }

  // -------------------------------------------------------------------
  //      load solution from previously performed scatra simulation
  // -------------------------------------------------------------------
  if (isscatrafromfile_)
  {
    int diffsteps = structure_field()->Dt() / ScaTraField()->Dt();
    if (ScaTraField()->Step() % diffsteps == 0)
    {
      Teuchos::RCP<CORE::IO::DiscretizationReader> reader =
          Teuchos::rcp(new CORE::IO::DiscretizationReader(ScaTraField()->discretization(),
              GLOBAL::Problem::Instance()->InputControlFile(), ScaTraField()->Step()));

      // check if this is a cardiac monodomain problem
      Teuchos::RCP<SCATRA::TimIntCardiacMonodomain> cardmono =
          Teuchos::rcp_dynamic_cast<SCATRA::TimIntCardiacMonodomain>(ScaTraField());

      if (cardmono == Teuchos::null)
      {
        // read phinp from restart file
        Teuchos::RCP<Epetra_MultiVector> phinptemp = reader->ReadVector("phinp");

        // replace old scatra map with new map since ssi map has more dofs
        int err = phinptemp->ReplaceMap(*ScaTraField()->dof_row_map());
        if (err) FOUR_C_THROW("Replacing old scatra map with new scatra map in ssi failed!");

        // update phinp
        ScaTraField()->Phinp()->Update(1.0, *phinptemp, 0.0);
      }
      else
      {
        // create vector with noderowmap from previously performed scatra calculation
        Teuchos::RCP<Epetra_Vector> phinptemp =
            CORE::LINALG::CreateVector(*cardmono->discretization()->NodeRowMap());

        // read phinp from restart file
        reader->ReadVector(phinptemp, "phinp");

        // replace old scatra map with new map since ssi map has more dofs
        int err = phinptemp->ReplaceMap(*ScaTraField()->dof_row_map());
        if (err) FOUR_C_THROW("Replacing old scatra map with new scatra map in ssi failed!");

        // update phinp
        ScaTraField()->Phinp()->Update(1.0, *phinptemp, 0.0);
      }
    }
  }
  // -------------------------------------------------------------------
  //                  solve nonlinear / linear equation
  // -------------------------------------------------------------------
  else
    ScaTraField()->Solve();


  // -------------------------------------------------------------------
  //                         update solution
  //        current solution becomes old solution of next timestep
  // -------------------------------------------------------------------
  ScaTraField()->Update();

  // -------------------------------------------------------------------
  // evaluate error for problems with analytical solution
  // -------------------------------------------------------------------
  ScaTraField()->evaluate_error_compared_to_analytical_sol();

  // -------------------------------------------------------------------
  //                         output of solution
  // -------------------------------------------------------------------
  ScaTraField()->check_and_write_output_and_restart();

  // cleanup
  ScaTraField()->discretization()->ClearState();
}

/*----------------------------------------------------------------------*/
// prepare time step
/*----------------------------------------------------------------------*/
void SSI::SSIPart1WCSolidToScatra::prepare_time_step(bool printheader)
{
  increment_time_and_step();

  if (printheader) print_header();

  // if adaptive time stepping: calculate time step in scatra (prepare_time_step() of Scatra) and
  // pass to structure
  if (ScaTraField()->TimeStepAdapted()) set_dt_from_sca_tra_to_structure();

  structure_field()->prepare_time_step();

  const int diffsteps = ScaTraField()->Dt() / structure_field()->Dt();

  if (structure_field()->Step() % diffsteps == 0)
  {
    if (is_s2_i_kinetics_with_pseudo_contact()) structure_field()->determine_stress_strain();
    set_struct_solution(structure_field()->Dispn(), structure_field()->Veln(),
        is_s2_i_kinetics_with_pseudo_contact());
    ScaTraField()->prepare_time_step();
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
SSI::SSIPart1WCSolidToScatra::SSIPart1WCSolidToScatra(
    const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams)
    : SSIPart1WC(comm, globaltimeparams)
{
  // Keep this constructor empty!
  // First do everything on the more basic objects like the discretizations, like e.g.
  // redistribution of elements. Only then call the setup to this class. This will call he setup to
  // all classes in the inheritance hierarchy. This way, this class may also override a method that
  // is called during Setup() in a base class.
}

/*----------------------------------------------------------------------*
 | Setup this class                                         rauch 08/16 |
 *----------------------------------------------------------------------*/
void SSI::SSIPart1WCSolidToScatra::Init(const Epetra_Comm& comm,
    const Teuchos::ParameterList& globaltimeparams, const Teuchos::ParameterList& scatraparams,
    const Teuchos::ParameterList& structparams, const std::string& struct_disname,
    const std::string& scatra_disname, bool isAle)
{
  // call setup of base class
  SSI::SSIPart1WC::Init(
      comm, globaltimeparams, scatraparams, structparams, struct_disname, scatra_disname, isAle);

  // do some checks
  {
    auto convform = CORE::UTILS::IntegralValue<INPAR::SCATRA::ConvForm>(scatraparams, "CONVFORM");
    if (convform != INPAR::SCATRA::convform_conservative)
    {
      FOUR_C_THROW(
          "If the scalar tranport problem is solved on the deforming domain, the conservative form "
          "must be "
          "used to include volume changes! Set 'CONVFORM' to 'conservative' in the SCALAR "
          "TRANSPORT DYNAMIC section!");
    }
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSIPart1WCSolidToScatra::Timeloop()
{
  // safety checks
  check_is_init();
  check_is_setup();

  if (structure_field()->Dt() > ScaTraField()->Dt())
  {
    FOUR_C_THROW(
        "Timestepsize of scatra should be equal or bigger than solid timestep in solid to scatra "
        "interaction");
  }

  const int diffsteps = ScaTraField()->Dt() / structure_field()->Dt();

  while (NotFinished())
  {
    prepare_time_step(false);
    do_struct_step();  // It has its own time and timestep variables, and it increments them by
                       // itself.
    if (structure_field()->Step() % diffsteps == 0)
    {
      if (is_s2_i_kinetics_with_pseudo_contact()) structure_field()->determine_stress_strain();
      set_struct_solution(structure_field()->Dispnp(), structure_field()->Velnp(),
          is_s2_i_kinetics_with_pseudo_contact());
      do_scatra_step();  // It has its own time and timestep variables, and it increments them by
                         // itself.
    }
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
SSI::SSIPart1WCScatraToSolid::SSIPart1WCScatraToSolid(
    const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams)
    : SSIPart1WC(comm, globaltimeparams)
{
  // Keep this constructor empty!
  // First do everything on the more basic objects like the discretizations, like e.g.
  // redistribution of elements. Only then call the setup to this class. This will call he setup to
  // all classes in the inheritance hierarchy. This way, this class may also override a method that
  // is called during Setup() in a base class.
}

/*----------------------------------------------------------------------*
 | Setup this class                                         rauch 08/16 |
 *----------------------------------------------------------------------*/
void SSI::SSIPart1WCScatraToSolid::Init(const Epetra_Comm& comm,
    const Teuchos::ParameterList& globaltimeparams, const Teuchos::ParameterList& scatraparams,
    const Teuchos::ParameterList& structparams, const std::string& struct_disname,
    const std::string& scatra_disname, bool isAle)
{
  // call setup of base class
  SSI::SSIPart1WC::Init(
      comm, globaltimeparams, scatraparams, structparams, struct_disname, scatra_disname, isAle);

  // Flag for reading scatra result from restart file instead of computing it
  isscatrafromfile_ = CORE::UTILS::IntegralValue<bool>(
      GLOBAL::Problem::Instance()->SSIControlParams(), "SCATRA_FROM_RESTART_FILE");
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSIPart1WCScatraToSolid::Timeloop()
{
  if (structure_field()->Dt() < ScaTraField()->Dt())
  {
    FOUR_C_THROW(
        "Timestepsize of solid should be equal or bigger than scatra timestep in scatra to solid "
        "interaction");
  }

  // set zero velocity and displacement field for scatra
  auto zeros_structure = CORE::LINALG::CreateVector(*structure_field()->dof_row_map(), true);
  set_struct_solution(zeros_structure, zeros_structure, false);

  ScaTraField()->prepare_time_loop();

  const int diffsteps = structure_field()->Dt() / ScaTraField()->Dt();
  while (!Finished())
  {
    prepare_time_step();
    do_scatra_step();  // It has its own time and timestep variables, and it increments them by
                       // itself.
    if (ScaTraField()->Step() % diffsteps == 0)
    {
      SetScatraSolution(ScaTraField()->Phinp());

      // set micro scale value (projected to macro scale) to structure field
      if (macro_scale()) set_micro_scatra_solution(ScaTraField()->PhinpMicro());

      // evaluate temperature from function and set to structural discretization
      evaluate_and_set_temperature_field();

      // prepare_time_step() is called after solving the scalar transport, because then the
      // predictor will include the new scalar solution
      structure_field()->prepare_time_step();
      do_struct_step();  // It has its own time and timestep variables, and it increments them by
                         // itself.
    }
  }
}

/*----------------------------------------------------------------------*/
// prepare time step
/*----------------------------------------------------------------------*/
void SSI::SSIPart1WCScatraToSolid::prepare_time_step(bool printheader)
{
  increment_time_and_step();
  print_header();

  ScaTraField()->prepare_time_step();
  // prepare_time_step of structure field is called later

  // copy time step to SSI problem, in case it was modified in ScaTra
  if (ScaTraField()->TimeStepAdapted()) set_dt_from_sca_tra_to_ssi();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool SSI::SSIPart1WCScatraToSolid::Finished() const
{
  if (diff_time_step_size())
    return !NotFinished();
  else
    return !(NotFinished() and ScaTraField()->NotFinished());
}

FOUR_C_NAMESPACE_CLOSE

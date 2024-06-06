/*----------------------------------------------------------------------*/
/*! \file
 \brief base class for all porous multiphase flow through elastic medium problems

   \level 3

 *----------------------------------------------------------------------*/

#include "4C_poromultiphase_base.hpp"

#include "4C_adapter_porofluidmultiphase_wrapper.hpp"
#include "4C_adapter_str_factory.hpp"
#include "4C_adapter_str_structure_new.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_porofluidmultiphase.hpp"
#include "4C_io.hpp"
#include "4C_io_control.hpp"
#include "4C_lib_discret.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_porofluidmultiphase_utils.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | constructor                                              vuong 08/16  |
 *----------------------------------------------------------------------*/
POROMULTIPHASE::PoroMultiPhaseBase::PoroMultiPhaseBase(
    const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams)
    : AlgorithmBase(comm, globaltimeparams),
      structure_(Teuchos::null),
      fluid_(Teuchos::null),
      struct_zeros_(Teuchos::null),
      solve_structure_(true),
      artery_coupl_(Core::UTILS::IntegralValue<int>(globaltimeparams, "ARTERY_COUPLING"))
{
}

/*----------------------------------------------------------------------*
 | initialize algorithm                                    vuong 08/16  |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseBase::Init(const Teuchos::ParameterList& globaltimeparams,
    const Teuchos::ParameterList& algoparams, const Teuchos::ParameterList& structparams,
    const Teuchos::ParameterList& fluidparams, const std::string& struct_disname,
    const std::string& fluid_disname, bool isale, int nds_disp, int nds_vel, int nds_solidpressure,
    int ndsporofluid_scatra, const std::map<int, std::set<int>>* nearbyelepairs)
{
  // access the global problem
  Global::Problem* problem = Global::Problem::Instance();

  // Create the two uncoupled subproblems.
  // access the structural discretization
  Teuchos::RCP<Discret::Discretization> structdis = problem->GetDis(struct_disname);

  // build structural time integrator
  Teuchos::RCP<Adapter::StructureBaseAlgorithmNew> adapterbase =
      Adapter::build_structure_algorithm(structparams);
  adapterbase->Init(globaltimeparams, const_cast<Teuchos::ParameterList&>(structparams), structdis);
  adapterbase->Setup();
  structure_ = adapterbase->structure_field();

  // initialize zero vector for convenience
  struct_zeros_ = Core::LinAlg::CreateVector(*structure_->dof_row_map(), true);
  // do we also solve the structure, this is helpful in case of fluid-scatra coupling without mesh
  // deformation
  solve_structure_ = Core::UTILS::IntegralValue<int>(algoparams, "SOLVE_STRUCTURE");

  // get the solver number used for ScalarTransport solver
  const int linsolvernumber = fluidparams.get<int>("LINEAR_SOLVER");


  // -------------------------------------------------------------------
  // access the fluid discretization
  // -------------------------------------------------------------------
  Teuchos::RCP<Discret::Discretization> fluiddis =
      Global::Problem::Instance()->GetDis(fluid_disname);

  // -------------------------------------------------------------------
  // set degrees of freedom in the discretization
  // -------------------------------------------------------------------
  if (!fluiddis->Filled()) fluiddis->fill_complete();

  // -------------------------------------------------------------------
  // context for output and restart
  // -------------------------------------------------------------------
  Teuchos::RCP<Core::IO::DiscretizationWriter> output = fluiddis->Writer();
  output->WriteMesh(0, 0.0);

  // -------------------------------------------------------------------
  // algorithm construction depending on
  // time-integration (or stationary) scheme
  // -------------------------------------------------------------------
  Inpar::POROFLUIDMULTIPHASE::TimeIntegrationScheme timintscheme =
      Core::UTILS::IntegralValue<Inpar::POROFLUIDMULTIPHASE::TimeIntegrationScheme>(
          fluidparams, "TIMEINTEGR");

  // build poro fluid time integrator
  Teuchos::RCP<Adapter::PoroFluidMultiphase> porofluid =
      POROFLUIDMULTIPHASE::UTILS::CreateAlgorithm(
          timintscheme, fluiddis, linsolvernumber, globaltimeparams, fluidparams, output);

  // wrap it
  fluid_ = Teuchos::rcp(new Adapter::PoroFluidMultiphaseWrapper(porofluid));
  // initialize it
  fluid_->Init(isale, nds_disp, nds_vel, nds_solidpressure, ndsporofluid_scatra, nearbyelepairs);

  // done.
  return;
}

/*----------------------------------------------------------------------*
 | read restart information for given time step (public)   vuong 08/16  |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseBase::read_restart(int restart)
{
  if (restart)
  {
    // read restart data for structure field (will set time and step internally)
    structure_->read_restart(restart);

    // read restart data for fluid field (will set time and step internally)
    fluid_->read_restart(restart);

    // reset time and step for the global algorithm
    SetTimeStep(structure_->TimeOld(), restart);
  }


  return;
}

/*----------------------------------------------------------------------*
 | time loop                                            kremheller 03/17 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseBase::Timeloop()
{
  // prepare the loop
  prepare_time_loop();

  // time loop
  while (NotFinished())
  {
    prepare_time_step();

    TimeStep();

    update_and_output();
  }

  return;
}

/*----------------------------------------------------------------------*
 | prepare the time loop                                     vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseBase::prepare_time_loop()
{
  // initial output
  if (solve_structure_)
  {
    constexpr bool force_prepare = true;
    structure_field()->prepare_output(force_prepare);
    structure_field()->Output();
    set_struct_solution(structure_field()->Dispnp(), structure_field()->Velnp());
  }
  else
  {
    // Inform user that structure field has been disabled
    print_structure_disabled_info();
    // just set displacements and velocities to zero
    set_struct_solution(struct_zeros_, struct_zeros_);
  }
  fluid_field()->prepare_time_loop();

  return;
}

/*----------------------------------------------------------------------*
 | prepare one time step                                     vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseBase::prepare_time_step()
{
  increment_time_and_step();

  structure_field()->discretization()->set_state(1, "porofluid", fluid_field()->Phinp());

  if (solve_structure_)
  {
    // NOTE: the predictor of the structure is called in here
    structure_field()->prepare_time_step();
    set_struct_solution(structure_field()->Dispnp(), structure_field()->Velnp());
  }
  else
    set_struct_solution(struct_zeros_, struct_zeros_);

  fluid_field()->prepare_time_step();
}

/*----------------------------------------------------------------------*
 | Test the results of all subproblems                       vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseBase::CreateFieldTest()
{
  Global::Problem* problem = Global::Problem::Instance();

  problem->AddFieldTest(structure_->CreateFieldTest());
  problem->AddFieldTest(fluid_->CreateFieldTest());
}

/*------------------------------------------------------------------------*
 | communicate the solution of the structure to the fluid    vuong 08/16  |
 *------------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseBase::set_struct_solution(
    Teuchos::RCP<const Epetra_Vector> disp, Teuchos::RCP<const Epetra_Vector> vel)
{
  set_mesh_disp(disp);
  set_velocity_fields(vel);
}

/*------------------------------------------------------------------------*
 | communicate the structure velocity  to the fluid           vuong 08/16  |
 *------------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseBase::set_velocity_fields(Teuchos::RCP<const Epetra_Vector> vel)
{
  fluid_->set_velocity_field(vel);
}

/*------------------------------------------------------------------------*
 | communicate the scatra solution to the fluid             vuong 08/16  |
 *------------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseBase::SetScatraSolution(
    unsigned nds, Teuchos::RCP<const Epetra_Vector> scalars)
{
  fluid_->SetScatraSolution(nds, scalars);
}


/*------------------------------------------------------------------------*
 | communicate the structure displacement to the fluid        vuong 08/16  |
 *------------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseBase::set_mesh_disp(Teuchos::RCP<const Epetra_Vector> disp)
{
  fluid_->ApplyMeshMovement(disp);
}

/*----------------------------------------------------------------------*
 | update fields and output results                         vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseBase::update_and_output()
{
  // prepare the output
  constexpr bool force_prepare = false;
  structure_field()->prepare_output(force_prepare);

  // update single fields
  structure_field()->Update();
  fluid_field()->Update();

  // evaluate error if desired
  fluid_field()->evaluate_error_compared_to_analytical_sol();

  // set structure on fluid (necessary for possible domain integrals)
  set_struct_solution(structure_field()->Dispnp(), structure_field()->Velnp());

  // output single fields
  structure_field()->Output();
  fluid_field()->Output();
}

/*------------------------------------------------------------------------*
 | dof map of vector of unknowns of structure field           vuong 08/16  |
 *------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> POROMULTIPHASE::PoroMultiPhaseBase::StructDofRowMap() const
{
  return structure_->dof_row_map();
}

/*------------------------------------------------------------------------*
 | dof map of vector of unknowns of fluid field           vuong 08/16  |
 *------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> POROMULTIPHASE::PoroMultiPhaseBase::FluidDofRowMap() const
{
  return fluid_->dof_row_map();
}

/*------------------------------------------------------------------------*
 | coupled artery-porofluid system matrix                kremheller 05/18 |
 *------------------------------------------------------------------------*/
Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase>
POROMULTIPHASE::PoroMultiPhaseBase::artery_porofluid_sysmat() const
{
  return fluid_->artery_porofluid_sysmat();
}

/*------------------------------------------------------------------------*
 | dof map of vector of unknowns of artery field         kremheller 05/18 |
 *------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> POROMULTIPHASE::PoroMultiPhaseBase::ArteryDofRowMap() const
{
  return fluid_->ArteryDofRowMap();
}

/*------------------------------------------------------------------------*
 | return structure displacements                             vuong 08/16  |
 *------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> POROMULTIPHASE::PoroMultiPhaseBase::StructDispnp() const
{
  return structure_->Dispnp();
}

/*------------------------------------------------------------------------*
 | return structure velocities                               vuong 08/16  |
 *------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> POROMULTIPHASE::PoroMultiPhaseBase::StructVelnp() const
{
  return structure_->Velnp();
}

/*------------------------------------------------------------------------*
 | return fluid Flux                                         vuong 08/16  |
 *------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_MultiVector> POROMULTIPHASE::PoroMultiPhaseBase::FluidFlux() const
{
  return fluid_->Flux();
}

/*------------------------------------------------------------------------*
 | return fluid Flux                                         vuong 08/16  |
 *------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> POROMULTIPHASE::PoroMultiPhaseBase::FluidPhinp() const
{
  return fluid_->Phinp();
}

/*------------------------------------------------------------------------*
 | return fluid Flux                                         vuong 08/16  |
 *------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> POROMULTIPHASE::PoroMultiPhaseBase::FluidSaturation() const
{
  return fluid_->Saturation();
}

/*------------------------------------------------------------------------*
 | return fluid Flux                                         vuong 08/16  |
 *------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> POROMULTIPHASE::PoroMultiPhaseBase::FluidPressure() const
{
  return fluid_->Pressure();
}

/*------------------------------------------------------------------------*
 | return fluid Flux                                         vuong 08/16  |
 *------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> POROMULTIPHASE::PoroMultiPhaseBase::SolidPressure() const
{
  return fluid_->SolidPressure();
}

/*----------------------------------------------------------------------*
 | inform user that structure is not solved            kremheller 04/17 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseBase::print_structure_disabled_info()
{
  // print out Info
  if (Comm().MyPID() == 0)
  {
    std::cout << "\n";
    std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                 "++++++++++++++++++++++++++++++++\n";
    std::cout << " INFO:    STRUCTURE FIELD IS NOT SOLVED   \n";
    std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                 "++++++++++++++++++++++++++++++++\n";
  }
}

FOUR_C_NAMESPACE_CLOSE

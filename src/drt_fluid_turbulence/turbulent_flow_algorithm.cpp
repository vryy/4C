/*----------------------------------------------------------------------*/
/*!

\brief algorithm for turbulent flows with separate inflow section

\maintainer Martin Kronbichler

\level 2

*/
/*----------------------------------------------------------------------*/

#include "turbulent_flow_algorithm.H"
#include "../drt_fluid/fluid_discret_extractor.H"

#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_dserror.H"
#include <Teuchos_TimeMonitor.hpp>


/*----------------------------------------------------------------------*
 | Destructor (public)                                   rasthofer 06/11|
 *----------------------------------------------------------------------*/
FLD::TurbulentFlowAlgorithm::~TurbulentFlowAlgorithm() { return; }


/*----------------------------------------------------------------------*
 | Constructor (public)                                  rasthofer 06/11|
 *----------------------------------------------------------------------*/
FLD::TurbulentFlowAlgorithm::TurbulentFlowAlgorithm(
    const Epetra_Comm& comm, const Teuchos::ParameterList& fdyn)
    : step_(0)
{
  if (comm.MyPID() == 0)
  {
    std::cout << "#-----------------------------------------------#" << std::endl;
    std::cout << "#       INITIALIZE BASIC FLUID ALGORITHM        #" << std::endl;
    std::cout << "#-----------------------------------------------#" << std::endl;
  }
  // initialize fluid algorithm
  // this is the first and main fluid algorithm
  fluidalgo_ = Teuchos::rcp(new ADAPTER::FluidBaseAlgorithm(fdyn, fdyn, "fluid", false));

  // get the compete fluid discretization
  fluiddis_ = fluidalgo_->FluidField()->Discretization();
  if (comm.MyPID() == 0)
  {
    std::cout << "#-----------------------------------------------#" << std::endl;
    std::cout << "#         EXTRACT INFLOW DISCRETIZATION         #" << std::endl;
    std::cout << "#-----------------------------------------------#" << std::endl;
  }
  // build extra discretization for turbulent inflow generation
  inflowgenerator_ =
      Teuchos::rcp(new FluidDiscretExtractor(fluiddis_, "TurbulentInflowSection", true));
  // and get this discretization
  inflowdis_ = inflowgenerator_->GetChildDiscretization();

  // set number of time steps
  numtimesteps_ = fdyn.sublist("TURBULENT INFLOW").get<int>("NUMINFLOWSTEP");

  if (comm.MyPID() == 0)
  {
    std::cout << "#-----------------------------------------------#" << std::endl;
    std::cout << "#       INITIALIZE INFLOW FLUID ALGORITHM       #" << std::endl;
    std::cout << "#-----------------------------------------------#" << std::endl;
  }

  // initialize fluid inflow algorithm
  // this is a second fluid algorithm
  inflowfluidalgo_ = Teuchos::rcp(new ADAPTER::FluidBaseAlgorithm(fdyn, inflowdis_));

  return;
}


/*--------------------------------------------------------------------------------*
 | Algorithm for development of turbulent flow in inflow section   rasthofer 06/11|
 *--------------------------------------------------------------------------------*/
void FLD::TurbulentFlowAlgorithm::TimeLoop()
{
  if (fluiddis_->Comm().MyPID() == 0)
  {
    std::cout << "#-----------------------------------------------#" << std::endl;
    std::cout << "#       START TURBULENT INFLOW COMPUTATION      #" << std::endl;
    std::cout << "#-----------------------------------------------#\n" << std::endl;
  }

  while (step_ < numtimesteps_)
  {
    step_++;

    // prepare time integration
    inflowfluidalgo_->FluidField()->PrepareTimeStep();
    if (fluiddis_->Comm().MyPID() == 0)
      printf("#   STEP = %4d/%4d     TIME: %11.4E  DT = %11.4E \n", step_, numtimesteps_,
          inflowfluidalgo_->FluidField()->Time(), inflowfluidalgo_->FluidField()->Dt());
    // slove nonlinear problem
    inflowfluidalgo_->FluidField()->Solve();
    // update time integration
    inflowfluidalgo_->FluidField()->Update();
    // write output of statistics only
    // remark: does also gmsh-output if required
    inflowfluidalgo_->FluidField()->StatisticsOutput();

    // transfer solution of inflow section to fluid discretization
    TransferInflowVelocity();

    // increase time and step only
    fluidalgo_->FluidField()->IncrementTimeAndStep();
    // velnp is set manually instead of being computed in Solve()
    // replaces Solve
    fluidalgo_->FluidField()->SetVelocityField(velnp_);
    // update time integration with given velocity field
    fluidalgo_->FluidField()->Update();
    // write output
    fluidalgo_->FluidField()->Output();
  }

  if (fluiddis_->Comm().MyPID() == 0)
  {
    std::cout << "#-----------------------------------------------#" << std::endl;
    std::cout << "#     FINISHED TURBULENT INFLOW COMPUTATION     #" << std::endl;
    std::cout << "#     -> problem ready for restart              #" << std::endl;
    std::cout << "#-----------------------------------------------#\n" << std::endl;
  }

  // summarize the performance measurements
  Teuchos::TimeMonitor::summarize();

  return;
}


/*-------------------------------------------------------------------------------------------*
 | transfer solution of inflow section to the complete fluid discretization   rasthofer 06/11|
 *-------------------------------------------------------------------------------------------*/
void FLD::TurbulentFlowAlgorithm::TransferInflowVelocity()
{
  if (fluiddis_->Comm().MyPID() == 0)
    std::cout << "#   transfer solution of inflow section ..." << std::flush;

  // velocity/pressure at time n+1 of inflow section
  Teuchos::RCP<const Epetra_Vector> inflowvelnp = inflowfluidalgo_->FluidField()->Velnp();

  // velocity/pressure at time n+1 to be transferred to the complete fluid field
  // get a vector layout from the complete discretization
  velnp_ = LINALG::CreateVector(*fluiddis_->DofRowMap(), true);

  // get exporter for transfer of dofs from inflow discretization to complete fluid discretization
  Epetra_Export exporter(inflowvelnp->Map(), velnp_->Map());
  // export inflow velocity
  int err = velnp_->Export(*inflowvelnp, exporter, Insert);
  if (err != 0) dserror("Export using exporter returned err=%d", err);

  if (fluiddis_->Comm().MyPID() == 0) std::cout << "done\n" << std::endl;

  return;
}


/*---------------------------------------------------------------------------*
 | read restart                                               rasthofer 06/11|
 *---------------------------------------------------------------------------*/
void FLD::TurbulentFlowAlgorithm::ReadRestart(const int restart)
{
  if (fluiddis_->Comm().MyPID() == 0)
  {
    std::cout << "#-----------------------------------------------#" << std::endl;
    std::cout << "#                 READ RESTART                  #" << std::endl;
    std::cout << "#-----------------------------------------------#\n" << std::endl;
  }
  // As we don't write a separate output for the inflow section, we first read
  // the values of the complete discretization, then extract the values belonging
  // to the inflow section and, finally, set them manually as restart values
  // in the fluid time integration.

  // set step
  step_ = restart;

  // read restart for complete discretization
  fluidalgo_->FluidField()->ReadRestart(restart);

  // vectors to be transferred to the inflow field
  // get a vector layout from the inflow discretization
  Teuchos::RCP<Epetra_Vector> velnp;
  velnp = LINALG::CreateVector(*inflowdis_->DofRowMap(), true);
  Teuchos::RCP<Epetra_Vector> veln;
  veln = LINALG::CreateVector(*inflowdis_->DofRowMap(), true);
  Teuchos::RCP<Epetra_Vector> velnm;
  velnm = LINALG::CreateVector(*inflowdis_->DofRowMap(), true);
  Teuchos::RCP<Epetra_Vector> accnp;
  accnp = LINALG::CreateVector(*inflowdis_->DofRowMap(), true);
  Teuchos::RCP<Epetra_Vector> accn;
  accn = LINALG::CreateVector(*inflowdis_->DofRowMap(), true);

  // get all vectors of restart
  Teuchos::RCP<const Epetra_Vector> fluidvelnp = fluidalgo_->FluidField()->Velnp();
  Teuchos::RCP<const Epetra_Vector> fluidveln = fluidalgo_->FluidField()->Veln();
  Teuchos::RCP<const Epetra_Vector> fluidvelnm = fluidalgo_->FluidField()->Velnm();
  Teuchos::RCP<const Epetra_Vector> fluidaccnp = fluidalgo_->FluidField()->Accnp();
  Teuchos::RCP<const Epetra_Vector> fluidaccn = fluidalgo_->FluidField()->Accn();

  // export vectors to inflow discretization
  int err = 0;
  Epetra_Export exportvelnp(fluidvelnp->Map(), velnp->Map());
  err = velnp->Export(*fluidvelnp, exportvelnp, Insert);
  if (err != 0) dserror("Export using exporter returned err=%d", err);
  Epetra_Export exportveln(fluidveln->Map(), veln->Map());
  err = veln->Export(*fluidveln, exportveln, Insert);
  if (err != 0) dserror("Export using exporter returned err=%d", err);
  Epetra_Export exportvelnm(fluidvelnm->Map(), velnm->Map());
  err = velnm->Export(*fluidvelnm, exportvelnm, Insert);
  if (err != 0) dserror("Export using exporter returned err=%d", err);
  Epetra_Export exportaccnp(fluidaccnp->Map(), accnp->Map());
  err = accnp->Export(*fluidaccnp, exportaccnp, Insert);
  if (err != 0) dserror("Export using exporter returned err=%d", err);
  Epetra_Export exportaccn(fluidaccn->Map(), accn->Map());
  err = accn->Export(*fluidaccn, exportaccn, Insert);
  if (err != 0) dserror("Export using exporter returned err=%d", err);

  // set values in the inflow field
  inflowfluidalgo_->FluidField()->SetRestart(
      restart, fluidalgo_->FluidField()->Time(), velnp, veln, velnm, accnp, accn);

  if (fluiddis_->Comm().MyPID() == 0) std::cout << "#   ... done \n" << std::endl;

  return;
}

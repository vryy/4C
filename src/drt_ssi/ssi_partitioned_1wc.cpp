/*----------------------------------------------------------------------*/
/*! \file
 \brief one way coupled partitioned scalar structure interaction

 \level 2

 \maintainer Christoph Schmidt

 *------------------------------------------------------------------------------------------------*/

#include "ssi_partitioned_1wc.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../linalg/linalg_utils_sparse_algebra_create.H"

#include "../drt_adapter/ad_str_wrapper.H"
#include "../drt_adapter/ad_str_ssiwrapper.H"
#include "../drt_adapter/adapter_scatra_base_algorithm.H"

#include "../drt_scatra/scatra_timint_implicit.H"
#include "../drt_scatra/scatra_timint_cardiac_monodomain.H"

#include "../drt_io/io.H"

SSI::SSI_Part1WC::SSI_Part1WC(
    const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams)
    : SSI_Part(comm, globaltimeparams), isscatrafromfile_(false)
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
int SSI::SSI_Part1WC::Init(const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams,
    const Teuchos::ParameterList& scatraparams, const Teuchos::ParameterList& structparams,
    const std::string struct_disname, const std::string scatra_disname, bool isAle)
{
  int returnvar = 0;

  // call setup of base class
  returnvar = SSI::SSI_Part::Init(
      comm, globaltimeparams, scatraparams, structparams, struct_disname, scatra_disname, isAle);

  return returnvar;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSI_Part1WC::DoStructStep()
{
  if (Comm().MyPID() == 0)
  {
    std::cout << "\n***********************\n STRUCTURE SOLVER \n***********************\n";
  }

  // Newton-Raphson iteration
  structure_->Solve();
  // calculate stresses, strains, energies
  structure_->PrepareOutput();
  // update all single field solvers
  structure_->Update();
  // write output to files
  structure_->Output();
  // write output to screen
  structure_->PrintStep();
  // clean up
  structure_->Discretization()->ClearState(true);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSI_Part1WC::DoScatraStep()
{
  if (Comm().MyPID() == 0)
  {
    std::cout << "\n***********************\n TRANSPORT SOLVER \n***********************\n";
  }

  // -------------------------------------------------------------------
  //                  solve nonlinear / linear equation
  // -------------------------------------------------------------------
  if (isscatrafromfile_)
  {
    int diffsteps = structure_->Dt() / scatra_->ScaTraField()->Dt();
    if (scatra_->ScaTraField()->Step() % diffsteps == 0)
    {
      // check if this is a cardiac monodomain problem
      Teuchos::RCP<SCATRA::TimIntCardiacMonodomain> cardmono =
          Teuchos::rcp_dynamic_cast<SCATRA::TimIntCardiacMonodomain>(scatra_->ScaTraField());
      if (cardmono == Teuchos::null)
        dserror("SCATRA_FROM_RESTART_FILE works only for Cardiac Monodoamin problems");
      // create vector with dofrowmap from previously performed scarta calculation
      Teuchos::RCP<Epetra_Vector> phinptemp = LINALG::CreateVector(*cardmono->DofRowMapScatra());
      Teuchos::RCP<IO::DiscretizationReader> reader = Teuchos::rcp(new IO::DiscretizationReader(
          scatra_->ScaTraField()->Discretization(), scatra_->ScaTraField()->Step()));
      // read phinp from restart file
      reader->ReadVector(phinptemp, "phinp");
      // replace old scatra map with new map since ssi map has more dofs
      phinptemp->ReplaceMap(*scatra_->ScaTraField()->DofRowMap());
      // update phinp
      scatra_->ScaTraField()->Phinp()->Update(1.0, *phinptemp, 0.0);
    }
  }
  else
    scatra_->ScaTraField()->Solve();  // really solve scatra problem


  // -------------------------------------------------------------------
  //                         update solution
  //        current solution becomes old solution of next timestep
  // -------------------------------------------------------------------
  scatra_->ScaTraField()->Update();

  // -------------------------------------------------------------------
  // evaluate error for problems with analytical solution
  // -------------------------------------------------------------------
  scatra_->ScaTraField()->EvaluateErrorComparedToAnalyticalSol();

  // -------------------------------------------------------------------
  //                         output of solution
  // -------------------------------------------------------------------
  scatra_->ScaTraField()->Output();

  // cleanup
  scatra_->ScaTraField()->Discretization()->ClearState();
}

/*----------------------------------------------------------------------*/
// prepare time step
/*----------------------------------------------------------------------*/
void SSI::SSI_Part1WC_SolidToScatra::PrepareTimeStep(bool printheader)
{
  IncrementTimeAndStep();

  if (printheader) PrintHeader();

  // if adaptive time stepping: calculate time step in scatra (PrepareTimeStep() of Scatra) and pass
  // to structure
  if (AdaptiveTimeStepping())
  {
    StructureField()->SetDt(scatra_->ScaTraField()->Dt());
    StructureField()->SetTimen(scatra_->ScaTraField()->Time());
  }
  structure_->PrepareTimeStep();

  const int diffsteps = scatra_->ScaTraField()->Dt() / structure_->Dt();

  if (structure_->Step() % diffsteps == 0)
  {
    SetStructSolution(structure_->Dispn(), structure_->Veln());
    scatra_->ScaTraField()->PrepareTimeStep();
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
SSI::SSI_Part1WC_SolidToScatra::SSI_Part1WC_SolidToScatra(
    const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams)
    : SSI_Part1WC(comm, globaltimeparams)
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
int SSI::SSI_Part1WC_SolidToScatra::Init(const Epetra_Comm& comm,
    const Teuchos::ParameterList& globaltimeparams, const Teuchos::ParameterList& scatraparams,
    const Teuchos::ParameterList& structparams, const std::string struct_disname,
    const std::string scatra_disname, bool isAle)
{
  int returnvar = 0;

  // call setup of base class
  returnvar = SSI::SSI_Part1WC::Init(
      comm, globaltimeparams, scatraparams, structparams, struct_disname, scatra_disname, isAle);

  // do some checks
  {
    INPAR::SCATRA::ConvForm convform =
        DRT::INPUT::IntegralValue<INPAR::SCATRA::ConvForm>(scatraparams, "CONVFORM");
    if (convform != INPAR::SCATRA::convform_conservative)
      dserror(
          "If the scalar tranport problem is solved on the deforming domain, the conservative form "
          "must be "
          "used to include volume changes! Set 'CONVFORM' to 'conservative' in the SCALAR "
          "TRANSPORT DYNAMIC section!");
  }

  return returnvar;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSI_Part1WC_SolidToScatra::Timeloop()
{
  // safety checks
  CheckIsInit();
  CheckIsSetup();

  if (structure_->Dt() > scatra_->ScaTraField()->Dt())
    dserror(
        "Timestepsize of scatra should be equal or bigger than solid timestep in solid to scatra "
        "interaction");

  const int diffsteps = scatra_->ScaTraField()->Dt() / structure_->Dt();

  while (NotFinished())
  {
    PrepareTimeStep(false);
    DoStructStep();  // It has its own time and timestep variables, and it increments them by
                     // itself.
    if (structure_->Step() % diffsteps == 0)
    {
      SetStructSolution(structure_->Dispnp(), structure_->Velnp());
      DoScatraStep();  // It has its own time and timestep variables, and it increments them by
                       // itself.
    }
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
SSI::SSI_Part1WC_ScatraToSolid::SSI_Part1WC_ScatraToSolid(
    const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams)
    : SSI_Part1WC(comm, globaltimeparams)
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
int SSI::SSI_Part1WC_ScatraToSolid::Init(const Epetra_Comm& comm,
    const Teuchos::ParameterList& globaltimeparams, const Teuchos::ParameterList& scatraparams,
    const Teuchos::ParameterList& structparams, const std::string struct_disname,
    const std::string scatra_disname, bool isAle)
{
  int returnvar = 0;

  // call setup of base class
  returnvar = SSI::SSI_Part1WC::Init(
      comm, globaltimeparams, scatraparams, structparams, struct_disname, scatra_disname, isAle);

  // Flag for reading scatra result from restart file instead of computing it
  isscatrafromfile_ = DRT::INPUT::IntegralValue<bool>(
      DRT::Problem::Instance()->SSIControlParams(), "SCATRA_FROM_RESTART_FILE");

  return returnvar;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSI_Part1WC_ScatraToSolid::Timeloop()
{
  if (structure_->Dt() < scatra_->ScaTraField()->Dt())
  {
    dserror(
        "Timestepsize of solid should be equal or bigger than scatra timestep in scatra to solid "
        "interaction");
  }

  // set zero velocity field for scatra
  scatra_->ScaTraField()->SetVelocityField(1);

  const int diffsteps = structure_->Dt() / scatra_->ScaTraField()->Dt();
  while (NotFinished())
  {
    PrepareTimeStep();
    DoScatraStep();  // It has its own time and timestep variables, and it increments them by
                     // itself.
    if (scatra_->ScaTraField()->Step() % diffsteps == 0)
    {
      SetScatraSolution(scatra_->ScaTraField()->Phinp());
      // PrepareTimeStep() is called after solving the scalar transport, because then the predictor
      // will include the new scalar solution
      structure_->PrepareTimeStep();
      DoStructStep();  // It has its own time and timestep variables, and it increments them by
                       // itself.
    }
  }
}

/*----------------------------------------------------------------------*/
// prepare time step
/*----------------------------------------------------------------------*/
void SSI::SSI_Part1WC_ScatraToSolid::PrepareTimeStep(bool printheader)
{
  IncrementTimeAndStep();
  // PrintHeader();

  scatra_->ScaTraField()->PrepareTimeStep();
  // PrepareTimeStep of structure field is called later
}

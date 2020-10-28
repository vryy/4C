/*----------------------------------------------------------------------*/
/*! \file
 \brief one way coupled partitioned scalar structure interaction

 \level 2


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
int SSI::SSIPart1WC::Init(const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams,
    const Teuchos::ParameterList& scatraparams, const Teuchos::ParameterList& structparams,
    const std::string struct_disname, const std::string scatra_disname, bool isAle)
{
  // call setup of base class
  int returnvar = SSI::SSIPart::Init(
      comm, globaltimeparams, scatraparams, structparams, struct_disname, scatra_disname, isAle);

  return returnvar;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSIPart1WC::DoStructStep()
{
  if (Comm().MyPID() == 0)
  {
    std::cout << "\n***********************\n STRUCTURE SOLVER \n***********************\n";
  }

  // Newton-Raphson iteration
  StructureField()->Solve();
  // calculate stresses, strains, energies
  StructureField()->PrepareOutput();
  // update all single field solvers
  StructureField()->Update();
  // write output to files
  StructureField()->Output();
  // write output to screen
  StructureField()->PrintStep();
  // clean up
  StructureField()->Discretization()->ClearState(true);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSIPart1WC::DoScatraStep()
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
    int diffsteps = StructureField()->Dt() / ScaTraField()->Dt();
    if (ScaTraField()->Step() % diffsteps == 0)
    {
      Teuchos::RCP<IO::DiscretizationReader> reader = Teuchos::rcp(
          new IO::DiscretizationReader(ScaTraField()->Discretization(), ScaTraField()->Step()));

      // check if this is a cardiac monodomain problem
      Teuchos::RCP<SCATRA::TimIntCardiacMonodomain> cardmono =
          Teuchos::rcp_dynamic_cast<SCATRA::TimIntCardiacMonodomain>(ScaTraField());

      if (cardmono == Teuchos::null)
      {
        // read phinp from restart file
        Teuchos::RCP<Epetra_MultiVector> phinptemp = reader->ReadVector("phinp");

        // replace old scatra map with new map since ssi map has more dofs
        int err = phinptemp->ReplaceMap(*ScaTraField()->DofRowMap());
        if (err) dserror("Replacing old scatra map with new scatra map in ssi failed!");

        // update phinp
        ScaTraField()->Phinp()->Update(1.0, *phinptemp, 0.0);
      }
      else
      {
        // create vector with dofrowmap from previously performed scatra calculation
        Teuchos::RCP<Epetra_Vector> phinptemp = LINALG::CreateVector(*cardmono->DofRowMapScatra());

        // read phinp from restart file
        reader->ReadVector(phinptemp, "phinp");

        // replace old scatra map with new map since ssi map has more dofs
        int err = phinptemp->ReplaceMap(*ScaTraField()->DofRowMap());
        if (err) dserror("Replacing old scatra map with new scatra map in ssi failed!");

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
  ScaTraField()->EvaluateErrorComparedToAnalyticalSol();

  // -------------------------------------------------------------------
  //                         output of solution
  // -------------------------------------------------------------------
  ScaTraField()->Output();

  // cleanup
  ScaTraField()->Discretization()->ClearState();
}

/*----------------------------------------------------------------------*/
// prepare time step
/*----------------------------------------------------------------------*/
void SSI::SSIPart1WCSolidToScatra::PrepareTimeStep(bool printheader)
{
  IncrementTimeAndStep();

  if (printheader) PrintHeader();

  // if adaptive time stepping: calculate time step in scatra (PrepareTimeStep() of Scatra) and pass
  // to structure
  if (ScaTraField()->TimeStepAdapted()) SetDtFromScaTraToStructure();

  StructureField()->PrepareTimeStep();

  const int diffsteps = ScaTraField()->Dt() / StructureField()->Dt();

  if (StructureField()->Step() % diffsteps == 0)
  {
    SetStructSolution(StructureField()->Dispn(), StructureField()->Veln());
    ScaTraField()->PrepareTimeStep();
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
int SSI::SSIPart1WCSolidToScatra::Init(const Epetra_Comm& comm,
    const Teuchos::ParameterList& globaltimeparams, const Teuchos::ParameterList& scatraparams,
    const Teuchos::ParameterList& structparams, const std::string struct_disname,
    const std::string scatra_disname, bool isAle)
{
  // call setup of base class
  int returnvar = SSI::SSIPart1WC::Init(
      comm, globaltimeparams, scatraparams, structparams, struct_disname, scatra_disname, isAle);

  // do some checks
  {
    auto convform = DRT::INPUT::IntegralValue<INPAR::SCATRA::ConvForm>(scatraparams, "CONVFORM");
    if (convform != INPAR::SCATRA::convform_conservative)
    {
      dserror(
          "If the scalar tranport problem is solved on the deforming domain, the conservative form "
          "must be "
          "used to include volume changes! Set 'CONVFORM' to 'conservative' in the SCALAR "
          "TRANSPORT DYNAMIC section!");
    }
  }

  return returnvar;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSIPart1WCSolidToScatra::Timeloop()
{
  // safety checks
  CheckIsInit();
  CheckIsSetup();

  if (StructureField()->Dt() > ScaTraField()->Dt())
  {
    dserror(
        "Timestepsize of scatra should be equal or bigger than solid timestep in solid to scatra "
        "interaction");
  }

  const int diffsteps = ScaTraField()->Dt() / StructureField()->Dt();

  while (NotFinished())
  {
    PrepareTimeStep(false);
    DoStructStep();  // It has its own time and timestep variables, and it increments them by
                     // itself.
    if (StructureField()->Step() % diffsteps == 0)
    {
      SetStructSolution(StructureField()->Dispnp(), StructureField()->Velnp());
      DoScatraStep();  // It has its own time and timestep variables, and it increments them by
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
int SSI::SSIPart1WCScatraToSolid::Init(const Epetra_Comm& comm,
    const Teuchos::ParameterList& globaltimeparams, const Teuchos::ParameterList& scatraparams,
    const Teuchos::ParameterList& structparams, const std::string struct_disname,
    const std::string scatra_disname, bool isAle)
{
  // call setup of base class
  int returnvar = SSI::SSIPart1WC::Init(
      comm, globaltimeparams, scatraparams, structparams, struct_disname, scatra_disname, isAle);

  // Flag for reading scatra result from restart file instead of computing it
  isscatrafromfile_ = DRT::INPUT::IntegralValue<bool>(
      DRT::Problem::Instance()->SSIControlParams(), "SCATRA_FROM_RESTART_FILE");

  return returnvar;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSIPart1WCScatraToSolid::Timeloop()
{
  if (StructureField()->Dt() < ScaTraField()->Dt())
  {
    dserror(
        "Timestepsize of solid should be equal or bigger than scatra timestep in scatra to solid "
        "interaction");
  }

  // set zero velocity field for scatra
  ScaTraField()->SetVelocityField(1);

  const int diffsteps = StructureField()->Dt() / ScaTraField()->Dt();
  while (NotFinished())
  {
    PrepareTimeStep();
    DoScatraStep();  // It has its own time and timestep variables, and it increments them by
                     // itself.
    if (ScaTraField()->Step() % diffsteps == 0)
    {
      SetScatraSolution(ScaTraField()->Phinp());

      // evaluate temperature from function and set to structural discretization
      EvaluateAndSetTemperatureField();

      // PrepareTimeStep() is called after solving the scalar transport, because then the predictor
      // will include the new scalar solution
      StructureField()->PrepareTimeStep();
      DoStructStep();  // It has its own time and timestep variables, and it increments them by
                       // itself.
    }
  }
}

/*----------------------------------------------------------------------*/
// prepare time step
/*----------------------------------------------------------------------*/
void SSI::SSIPart1WCScatraToSolid::PrepareTimeStep(bool printheader)
{
  IncrementTimeAndStep();
  // PrintHeader();

  ScaTraField()->PrepareTimeStep();
  // PrepareTimeStep of structure field is called later
}

/*!------------------------------------------------------------------------------------------------*
 \brief two way coupled partitioned scalar structure interaction

 \level 2

 \maintainer Christoph Schmidt

 *------------------------------------------------------------------------------------------------*/

#include "ssi_partitioned_2wc.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../linalg/linalg_utils.H"

#include "../drt_adapter/ad_str_ssiwrapper.H"
#include "../drt_adapter/adapter_scatra_base_algorithm.H"

#include "../drt_scatra_ele/scatra_ele.H"
#include "../drt_scatra/scatra_timint_implicit.H"

#include <Teuchos_Time.hpp>

/*----------------------------------------------------------------------*
 | constructor                                               Thon 12/14 |
 *----------------------------------------------------------------------*/
SSI::SSI_Part2WC::SSI_Part2WC(
    const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams)
    : SSI_Part(comm, globaltimeparams),
      scaincnp_(Teuchos::null),
      dispincnp_(Teuchos::null),
      ittol_(-1.0),
      itmax_(-1)
{
  // Keep this constructor empty!
  // First do everything on the more basic objects like the discretizations, like e.g.
  // redistribution of elements. Only then call the setup to this class. This will call the setup to
  // all classes in the inheritance hierarchy. This way, this class may also override a method that
  // is called during Setup() in a base class.
}

/*----------------------------------------------------------------------*
 | Init this class                                          rauch 08/16 |
 *----------------------------------------------------------------------*/
int SSI::SSI_Part2WC::Init(const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams,
    const Teuchos::ParameterList& scatraparams, const Teuchos::ParameterList& structparams,
    const std::string struct_disname, const std::string scatra_disname, bool isAle)
{
  int returnvar = 0;

  // call setup of base class
  returnvar = SSI::SSI_Part::Init(
      comm, globaltimeparams, scatraparams, structparams, struct_disname, scatra_disname, isAle);

  // call the SSI parameter lists
  const Teuchos::ParameterList& ssicontrol = DRT::Problem::Instance()->SSIControlParams();
  const Teuchos::ParameterList& ssicontrolpart =
      DRT::Problem::Instance()->SSIControlParams().sublist("PARTITIONED");

  // do some checks
  {
    INPAR::STR::DynamicType structtimealgo =
        DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(structparams, "DYNAMICTYP");
    if (structtimealgo == INPAR::STR::dyna_statics)
      dserror(
          "If you use statics as the structural time integrator no velocities will be calculated "
          "and hence"
          "the deformations will not be applied to the scalar transport problem!");

    INPAR::SCATRA::ConvForm convform =
        DRT::INPUT::IntegralValue<INPAR::SCATRA::ConvForm>(scatraparams, "CONVFORM");
    if (convform == INPAR::SCATRA::convform_convective)
    {
      // get scatra discretization
      Teuchos::RCP<DRT::Discretization> scatradis =
          DRT::Problem::Instance()->GetDis(scatra_disname);

      // loop over all elements of scatra discretization to check if impltype is correct or not
      for (int i = 0; i < scatradis->NumMyColElements(); ++i)
      {
        if ((dynamic_cast<DRT::ELEMENTS::Transport*>(scatradis->lColElement(i)))->ImplType() !=
            INPAR::SCATRA::impltype_refconcreac)
          dserror(
              "If the scalar transport problem is solved on the deforming domain, the conservative "
              "form must be "
              "used to include volume changes! Set 'CONVFORM' to 'conservative' in the SCALAR "
              "TRANSPORT DYNAMIC section or use RefConcReac as ScatraType!");
      }
    }
  }

  if (DRT::INPUT::IntegralValue<int>(ssicontrol, "DIFFTIMESTEPSIZE"))
  {
    dserror("Different time stepping for two way coupling not implemented yet.");
  }

  // Get the parameters for the ConvergenceCheck
  itmax_ = ssicontrol.get<int>("ITEMAX");          // default: =10
  ittol_ = ssicontrolpart.get<double>("CONVTOL");  // default: =1e-6

  return returnvar;
}

/*----------------------------------------------------------------------*
 | Setup this class                                          rauch 08/16 |
 *----------------------------------------------------------------------*/
void SSI::SSI_Part2WC::Setup()
{
  // call setup in base class
  SSI::SSI_Part::Setup();

  // construct increment vectors
  scaincnp_ = LINALG::CreateVector(*scatra_->ScaTraField()->Discretization()->DofRowMap(0), true);
  dispincnp_ = LINALG::CreateVector(*structure_->DofRowMap(0), true);
  return;
}

/*----------------------------------------------------------------------*
 | Setup the discretizations                                rauch 08/16 |
 *----------------------------------------------------------------------*/
void SSI::SSI_Part2WC::InitDiscretizations(
    const Epetra_Comm& comm, const std::string& struct_disname, const std::string& scatra_disname)
{
  // call init in base class
  SSI::SSI_Base::InitDiscretizations(comm, struct_disname, scatra_disname);
  return;
}

/*----------------------------------------------------------------------*
 | Timeloop for 2WC SSI problems                             Thon 12/14 |
 *----------------------------------------------------------------------*/
void SSI::SSI_Part2WC::Timeloop()
{
  // safety checks
  CheckIsInit();
  CheckIsSetup();

  PrepareTimeLoop();

  // time loop
  while (NotFinished() and scatra_->ScaTraField()->NotFinished())
  {
    PrepareTimeStep();

    // store time before calling nonlinear solver
    double time = Teuchos::Time::wallTime();

    OuterLoop();

    // determine time spent by nonlinear solver and take maximum over all processors via
    // communication
    double mydtnonlinsolve(Teuchos::Time::wallTime() - time), dtnonlinsolve(0.);
    Comm().MaxAll(&mydtnonlinsolve, &dtnonlinsolve, 1);

    // output performance statistics associated with nonlinear solver into *.csv file if applicable
    if (DRT::INPUT::IntegralValue<int>(
            *scatra_->ScaTraField()->ScatraParameterList(), "OUTPUTNONLINSOLVERSTATS"))
      scatra_->ScaTraField()->OutputNonlinSolverStats(iter_, dtnonlinsolve, Step(), Comm());

    UpdateAndOutput();
  }
}

/*----------------------------------------------------------------------*
 | Solve structure filed                                     Thon 12/14 |
 *----------------------------------------------------------------------*/
void SSI::SSI_Part2WC::DoStructStep()
{
  if (Comm().MyPID() == 0)
  {
    std::cout << "\n***********************\n STRUCTURE SOLVER \n***********************\n";
  }

  // Newton-Raphson iteration
  structure_->Solve();

  // set mesh displacement and velocity fields
  return SetStructSolution(structure_->Dispnp(), structure_->Velnp());
}

/*----------------------------------------------------------------------*
 | Solve Scatra field                                        Thon 12/14 |
 *----------------------------------------------------------------------*/
void SSI::SSI_Part2WC::DoScatraStep()
{
  if (Comm().MyPID() == 0)
  {
    std::cout << "\n***********************\n  TRANSPORT SOLVER \n***********************\n";
  }

  // -------------------------------------------------------------------
  //                  solve nonlinear / linear equation
  // -------------------------------------------------------------------
  scatra_->ScaTraField()->Solve();

  // set structure-based scalar transport values
  return SetScatraSolution(scatra_->ScaTraField()->Phinp());
}

/*----------------------------------------------------------------------*
 | Solve Scatra field                                        Thon 12/14 |
 *----------------------------------------------------------------------*/
void SSI::SSI_Part2WC::PreOperator1()
{
  if (iter_ != 1 and use_old_structure_)
  {
    // NOTE: the predictor is NOT called in here. Just the screen output is not correct.
    // we only get norm of the evaluation of the structure problem
    structure_->PreparePartitionStep();
  }
}

/*----------------------------------------------------------------------*/
// prepare time step
/*----------------------------------------------------------------------*/
void SSI::SSI_Part2WC::PrepareTimeLoop()
{
  // initial output
  structure_->PrepareOutput();
  structure_->Output();
  SetStructSolution(structure_->Dispnp(), structure_->Velnp());
  scatra_->ScaTraField()->PrepareTimeLoop();
}

/*----------------------------------------------------------------------*/
// prepare time step
/*----------------------------------------------------------------------*/
void SSI::SSI_Part2WC::PrepareTimeStep(bool printheader)
{
  IncrementTimeAndStep();
  if (printheader) PrintHeader();

  SetStructSolution(structure_->Dispnp(), structure_->Velnp());
  scatra_->ScaTraField()->PrepareTimeStep();

  SetScatraSolution(scatra_->ScaTraField()->Phinp());
  // NOTE: the predictor of the structure is called in here
  structure_->PrepareTimeStep();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSI_Part2WC::UpdateAndOutput()
{
  structure_->PrepareOutput();

  structure_->Update();
  scatra_->ScaTraField()->Update();

  scatra_->ScaTraField()->EvaluateErrorComparedToAnalyticalSol();

  structure_->Output();
  scatra_->ScaTraField()->Output();
}


/*----------------------------------------------------------------------*
 | update the current states in every iteration             rauch 05/16 |
 *----------------------------------------------------------------------*/
void SSI::SSI_Part2WC::IterUpdateStates()
{
  // store last solutions (current states).
  // will be compared in ConvergenceCheck to the solutions,
  // obtained from the next Struct and Scatra steps.
  scaincnp_->Update(1.0, *scatra_->ScaTraField()->Phinp(), 0.0);
  dispincnp_->Update(1.0, *structure_->Dispnp(), 0.0);

  return;
}  // IterUpdateStates()


/*----------------------------------------------------------------------*
 | Outer Timeloop for 2WC SSi without relaxation            rauch 06/17 |
 *----------------------------------------------------------------------*/
void SSI::SSI_Part2WC::OuterLoop()
{
  // reset iteration number
  iter_ = 0;
  bool stopnonliniter = false;

  if (Comm().MyPID() == 0)
  {
    std::cout << "\n****************************************\n          OUTER ITERATION "
                 "LOOP\n****************************************\n";
  }

  while (stopnonliniter == false)
  {
    // increment iteration number
    iter_++;

    // update the states to the last solutions obtained
    IterUpdateStates();

    // standard ssi_2wc :
    // 1.) solve structural system
    // 2.) set disp and vel states in scatra field
    PreOperator1();
    Operator1();
    PostOperator1();

    // standard ssi_2wc :
    // 1.) solve scalar transport equation
    // 2.) set phi state in structure field
    PreOperator2();
    Operator2();
    PostOperator2();

    // check convergence for all fields
    // stop iteration loop if converged
    stopnonliniter = ConvergenceCheck(iter_);
  }

  return;
}

/*----------------------------------------------------------------------*
 | convergence check for both fields (scatra & structure) (copied form tsi)
 *----------------------------------------------------------------------*/
bool SSI::SSI_Part2WC::ConvergenceCheck(int itnum)
{
  // convergence check based on the scalar increment
  bool stopnonliniter = false;

  //    | scalar increment |_2
  //  -------------------------------- < Tolerance
  //    | scalar state n+1 |_2
  //
  // AND
  //
  //    | scalar increment |_2
  //  -------------------------------- < Tolerance , with n := global length of vector
  //        dt * sqrt(n)
  //
  // The same is checked for the structural displacements.
  //

  // variables to save different L2 - Norms
  // define L2-norm of incremental scalar and scalar
  double scaincnorm_L2(0.0);
  double scanorm_L2(0.0);
  double dispincnorm_L2(0.0);
  double dispnorm_L2(0.0);

  // build the current scalar increment Inc T^{i+1}
  // \f Delta T^{k+1} = Inc T^{k+1} = T^{k+1} - T^{k}  \f
  scaincnp_->Update(1.0, *(scatra_->ScaTraField()->Phinp()), -1.0);
  dispincnp_->Update(1.0, *(structure_->Dispnp()), -1.0);

  // build the L2-norm of the scalar increment and the scalar
  scaincnp_->Norm2(&scaincnorm_L2);
  scatra_->ScaTraField()->Phinp()->Norm2(&scanorm_L2);
  dispincnp_->Norm2(&dispincnorm_L2);
  structure_->Dispnp()->Norm2(&dispnorm_L2);

  // care for the case that there is (almost) zero scalar
  if (scanorm_L2 < 1e-6) scanorm_L2 = 1.0;
  if (dispnorm_L2 < 1e-6) dispnorm_L2 = 1.0;

  // print the incremental based convergence check to the screen
  if (Comm().MyPID() == 0)
  {
    std::cout << "\n";
    std::cout
        << "***********************************************************************************\n";
    std::cout << "    OUTER ITERATION STEP    \n";
    std::cout
        << "***********************************************************************************\n";
    printf(
        "+--------------+---------------------+----------------+------------------+----------------"
        "----+------------------+\n");
    printf(
        "|-  step/max  -|-  tol      [norm]  -|-  scalar-inc  -|-  disp-inc      -|-  "
        "scalar-rel-inc  -|-  disp-rel-inc  -|\n");
    printf(
        "|   %3d/%3d    |  %10.3E[L_2 ]   |  %10.3E    |  %10.3E      |  %10.3E        |  %10.3E   "
        "   |",
        itnum, itmax_, ittol_, scaincnorm_L2 / Dt() / sqrt(scaincnp_->GlobalLength()),
        dispincnorm_L2 / Dt() / sqrt(dispincnp_->GlobalLength()), scaincnorm_L2 / scanorm_L2,
        dispincnorm_L2 / dispnorm_L2);
    printf("\n");
    printf(
        "+--------------+---------------------+----------------+------------------+----------------"
        "----+------------------+\n");
  }

  // converged
  if (((scaincnorm_L2 / scanorm_L2) <= ittol_) and ((dispincnorm_L2 / dispnorm_L2) <= ittol_) and
      ((dispincnorm_L2 / Dt() / sqrt(dispincnp_->GlobalLength())) <= ittol_) and
      ((scaincnorm_L2 / Dt() / sqrt(scaincnp_->GlobalLength())) <= ittol_))
  {
    stopnonliniter = true;
    if (Comm().MyPID() == 0)
    {
      printf(
          "|  Outer Iteration loop converged after iteration %3d/%3d !                             "
          "                         |\n",
          itnum, itmax_);
      printf(
          "+--------------+---------------------+----------------+------------------+--------------"
          "------+------------------+\n");
    }
  }

  // stop if itemax is reached without convergence
  // timestep
  if ((itnum == itmax_) and
      (((scaincnorm_L2 / scanorm_L2) > ittol_) or ((dispincnorm_L2 / dispnorm_L2) > ittol_) or
          ((dispincnorm_L2 / Dt() / sqrt(dispincnp_->GlobalLength())) > ittol_) or
          (scaincnorm_L2 / Dt() / sqrt(scaincnp_->GlobalLength())) > ittol_))
  {
    stopnonliniter = true;
    if ((Comm().MyPID() == 0))
    {
      printf(
          "|     >>>>>> not converged in itemax steps!                                             "
          "                         |\n");
      printf(
          "+--------------+---------------------+----------------+------------------+--------------"
          "------+------------------+\n");
      printf("\n");
      printf("\n");
    }
    dserror("The partitioned SSI solver did not converge in ITEMAX steps!");
  }

  return stopnonliniter;
}

/*----------------------------------------------------------------------*
 | calculate velocities by a FD approximation                Thon 14/11 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> SSI::SSI_Part2WC::CalcVelocity(Teuchos::RCP<const Epetra_Vector> dispnp)
{
  Teuchos::RCP<Epetra_Vector> vel = Teuchos::null;
  // copy D_n onto V_n+1
  vel = Teuchos::rcp(new Epetra_Vector(*(structure_->Dispn())));
  // calculate velocity with timestep Dt()
  //  V_n+1^k = (D_n+1^k - D_n) / Dt
  vel->Update(1. / Dt(), *dispnp, -1. / Dt());

  return vel;
}  // CalcVelocity()

/*----------------------------------------------------------------------*
 | Constructor                                               Thon 12/14 |
 *----------------------------------------------------------------------*/
SSI::SSI_Part2WC_SolidToScatra_Relax::SSI_Part2WC_SolidToScatra_Relax(
    const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams)
    : SSI_Part2WC(comm, globaltimeparams), omega_(-1.0)
{
  // Keep this constructor empty!
  // First do everything on the more basic objects like the discretizations, like e.g.
  // redistribution of elements. Only then call the setup to this class. This will call the setup to
  // all classes in the inheritance hierarchy. This way, this class may also override a method that
  // is called during Setup() in a base class.
}

/*----------------------------------------------------------------------*
 | Init this class                                          rauch 08/16 |
 *----------------------------------------------------------------------*/
int SSI::SSI_Part2WC_SolidToScatra_Relax::Init(const Epetra_Comm& comm,
    const Teuchos::ParameterList& globaltimeparams, const Teuchos::ParameterList& scatraparams,
    const Teuchos::ParameterList& structparams, const std::string struct_disname,
    const std::string scatra_disname, bool isAle)
{
  int returnvar = 0;

  // call init of base class
  returnvar = SSI::SSI_Part2WC::Init(
      comm, globaltimeparams, scatraparams, structparams, struct_disname, scatra_disname, isAle);

  const Teuchos::ParameterList& ssicontrolpart =
      DRT::Problem::Instance()->SSIControlParams().sublist("PARTITIONED");

  // Get minimal relaxation parameter from input file
  omega_ = ssicontrolpart.get<double>("STARTOMEGA");

  return returnvar;
}

/*----------------------------------------------------------------------*
 | Outer Timeloop for 2WC SSi with relaxed displacements     Thon 12/14 |
 *----------------------------------------------------------------------*/
void SSI::SSI_Part2WC_SolidToScatra_Relax::OuterLoop()
{
  iter_ = 0;
  bool stopnonliniter = false;

  if (Comm().MyPID() == 0)
  {
    std::cout << "\n****************************************\n          OUTER ITERATION "
                 "LOOP\n****************************************\n";
  }

  // these are the relaxed inputs
  Teuchos::RCP<Epetra_Vector> dispnp = LINALG::CreateVector(*(structure_->DofRowMap(0)), true);
  Teuchos::RCP<Epetra_Vector> velnp = LINALG::CreateVector(*(structure_->DofRowMap(0)), true);

  while (stopnonliniter == false)
  {
    ++iter_;

    if (iter_ == 1)
    {
      dispnp->Update(1.0, *(structure_->Dispnp()), 0.0);  // TSI does Dispn()
      velnp->Update(1.0, *(structure_->Velnp()), 0.0);
    }

    // store scalars and displacements for the convergence check later
    scaincnp_->Update(1.0, *scatra_->ScaTraField()->Phinp(), 0.0);
    dispincnp_->Update(1.0, *dispnp, 0.0);

    // begin nonlinear solver / outer iteration ***************************

    // set relaxed mesh displacements and velocity field
    SetStructSolution(dispnp, velnp);

    // solve scalar transport equation
    DoScatraStep();

    // prepare a partitioned structure step
    if (iter_ != 1 and use_old_structure_) structure_->PreparePartitionStep();

    // solve structural system
    DoStructStep();

    // end nonlinear solver / outer iteration *****************************

    // check convergence for all fields and stop iteration loop if
    // convergence is achieved overall
    stopnonliniter = ConvergenceCheck(iter_);

    // get relaxation parameter
    CalcOmega(omega_, iter_);

    // do the relaxation
    // d^{i+1} = omega^{i+1} . d^{i+1} + (1- omega^{i+1}) d^i
    //         = d^i + omega^{i+1} * ( d^{i+1} - d^i )
    dispnp->Update(omega_, *dispincnp_, 1.0);

    // since the velocity field has to fit to the relaxated displacements we also have to relaxate
    // them.
    if (use_old_structure_)
      // since the velocity depends nonlinear on the displacements we can just approximate them via
      // finite differences here.
      velnp = CalcVelocity(dispnp);
    else
    {
      // consistent derivation of velocity field from displacement field
      // SetState(dispnp) will automatically be undone during the next evaluation of the structural
      // field
      structure_->SetState(dispnp);
      velnp->Update(1., *structure_->Velnp(), 0.);
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | Calculate relaxation parameter                            Thon 12/14 |
 *----------------------------------------------------------------------*/
void SSI::SSI_Part2WC_SolidToScatra_Relax::CalcOmega(double& omega, const int itnum)
{
  // nothing to do in here since we have a constant relaxation parameter: omega != startomega_;
  if (Comm().MyPID() == 0)
    std::cout << "Fixed relaxation parameter omega is: " << omega << std::endl;
}

/*----------------------------------------------------------------------*
 | Constructor                                               Thon 12/14 |
 *----------------------------------------------------------------------*/
SSI::SSI_Part2WC_SolidToScatra_Relax_Aitken::SSI_Part2WC_SolidToScatra_Relax_Aitken(
    const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams)
    : SSI_Part2WC_SolidToScatra_Relax(comm, globaltimeparams), dispincnpold_(Teuchos::null)
{
  // Keep this constructor empty!
  // First do everything on the more basic objects like the discretizations, like e.g.
  // redistribution of elements. Only then call the setup to this class. This will call he setup to
  // all classes in the inheritance hierarchy. This way, this class may also override a method that
  // is called during Setup() in a base class.
}

/*----------------------------------------------------------------------*
 | Setup this class                                          fang 01/18 |
 *----------------------------------------------------------------------*/
void SSI::SSI_Part2WC_SolidToScatra_Relax_Aitken::Setup()
{
  // call setup of base class
  SSI::SSI_Part2WC::Setup();

  // setup old scatra increment vector
  dispincnpold_ = LINALG::CreateVector(*structure_->DofRowMap(0), true);

  return;
}

/*----------------------------------------------------------------------*
 | Calculate relaxation parameter via Aitken                 Thon 12/14 |
 *----------------------------------------------------------------------*/
void SSI::SSI_Part2WC_SolidToScatra_Relax_Aitken::CalcOmega(double& omega, const int itnum)
{
  const Teuchos::ParameterList& ssicontrolpart =
      DRT::Problem::Instance()->SSIControlParams().sublist("PARTITIONED");
  ;
  // Get maximal relaxation parameter from input file
  const double maxomega = ssicontrolpart.get<double>("MAXOMEGA");
  // Get minimal relaxation parameter from input file
  const double minomega = ssicontrolpart.get<double>("MINOMEGA");

  // calculate difference of current (i+1) and old (i) residual vector
  // dispincnpdiff = ( r^{i+1}_{n+1} - r^i_{n+1} )
  Teuchos::RCP<Epetra_Vector> dispincnpdiff =
      LINALG::CreateVector(*(structure_->DofRowMap(0)), true);
  dispincnpdiff->Update(
      1.0, *dispincnp_, (-1.0), *dispincnpold_, 0.0);  // update r^{i+1}_{n+1} - r^i_{n+1}

  double dispincnpdiffnorm = 0.0;
  dispincnpdiff->Norm2(&dispincnpdiffnorm);
  if (dispincnpdiffnorm <= 1e-06 and Comm().MyPID() == 0)
    std::cout << "Warning: The structure increment is to small in order to use it for Aitken "
                 "relaxation. Using the previous Omega instead!"
              << std::endl;

  // calculate dot product
  double dispincsdot = 0.0;  // delsdot = ( r^{i+1}_{n+1} - r^i_{n+1} )^T . r^{i+1}_{n+1}
  dispincnpdiff->Dot(*dispincnp_, &dispincsdot);

  if (itnum != 1 and dispincnpdiffnorm > 1e-06)
  {  // relaxation parameter
    // omega^{i+1} = 1- mu^{i+1} and nu^{i+1} = nu^i + (nu^i -1) . (r^{i+1} - r^i)^T . (-r^{i+1}) /
    // |r^{i+1} - r^{i}|^2 results in
    omega = omega *
            (1 - (dispincsdot) / (dispincnpdiffnorm *
                                     dispincnpdiffnorm));  // compare e.g. PhD thesis U. Kuettler

    // we force omega to be in the range defined in the input file
    if (omega < minomega)
    {
      if (Comm().MyPID() == 0)
        std::cout << "Warning: The calculation of the relaxation parameter omega via Aitken did "
                     "lead to a value smaller than MINOMEGA!"
                  << std::endl;
      omega = minomega;
    }
    if (omega > maxomega)
    {
      if (Comm().MyPID() == 0)
        std::cout << "Warning: The calculation of the relaxation parameter omega via Aitken did "
                     "lead to a value bigger than MAXOMEGA!"
                  << std::endl;
      omega = maxomega;
    }
  }

  // else //if itnum==1 nothing is to do here since we want to take the last omega from the previous
  // step
  if (Comm().MyPID() == 0)
    std::cout << "Using Aitken the relaxation parameter omega was estimated to: " << omega
              << std::endl;

  // update history vector old increment r^i_{n+1}
  dispincnpold_->Update(1.0, *dispincnp_, 0.0);
}


/*----------------------------------------------------------------------*
 | Constructor                                               Thon 12/14 |
 *----------------------------------------------------------------------*/
SSI::SSI_Part2WC_ScatraToSolid_Relax::SSI_Part2WC_ScatraToSolid_Relax(
    const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams)
    : SSI_Part2WC(comm, globaltimeparams), omega_(-1.0)
{
  // Keep this constructor empty!
  // First do everything on the more basic objects like the discretizations, like e.g.
  // redistribution of elements. Only then call the setup to this class. This will call the setup to
  // all classes in the inheritance hierarchy. This way, this class may also override a method that
  // is called during Setup() in a base class.
}

/*----------------------------------------------------------------------*
 | Setup this class                                         rauch 08/16 |
 *----------------------------------------------------------------------*/
int SSI::SSI_Part2WC_ScatraToSolid_Relax::Init(const Epetra_Comm& comm,
    const Teuchos::ParameterList& globaltimeparams, const Teuchos::ParameterList& scatraparams,
    const Teuchos::ParameterList& structparams, const std::string struct_disname,
    const std::string scatra_disname, bool isAle)
{
  int returnvar = 0;

  // call setup of base class
  returnvar = SSI::SSI_Part2WC::Init(
      comm, globaltimeparams, scatraparams, structparams, struct_disname, scatra_disname, isAle);

  const Teuchos::ParameterList& ssicontrolpart =
      DRT::Problem::Instance()->SSIControlParams().sublist("PARTITIONED");

  // Get start relaxation parameter from input file
  omega_ = ssicontrolpart.get<double>("STARTOMEGA");

  if (Comm().MyPID() == 0)
  {
    std::cout << "\n#########################################################################\n  "
              << std::endl;
    std::cout << "The ScatraToSolid relaxations are not well tested . Keep your eyes open!  "
              << std::endl;
    std::cout << "\n#########################################################################\n  "
              << std::endl;
  }

  return returnvar;
}

/*----------------------------------------------------------------------*
 | Outer Timeloop for 2WC SSi with relaxed scalar             Thon 12/14 |
 *----------------------------------------------------------------------*/
void SSI::SSI_Part2WC_ScatraToSolid_Relax::OuterLoop()
{
  iter_ = 0;
  bool stopnonliniter = false;

  if (Comm().MyPID() == 0)
  {
    std::cout << "\n****************************************\n          OUTER ITERATION "
                 "LOOP\n****************************************\n";
  }

  // this is the relaxed input
  Teuchos::RCP<Epetra_Vector> phinp =
      LINALG::CreateVector(*scatra_->ScaTraField()->Discretization()->DofRowMap(0), true);

  while (stopnonliniter == false)
  {
    ++iter_;

    if (iter_ == 1)
    {
      phinp->Update(1.0, *(scatra_->ScaTraField()->Phinp()), 0.0);  // TSI does Dispn()
    }

    // store scalars and displacements for the convergence check later
    scaincnp_->Update(1.0, *phinp, 0.0);
    dispincnp_->Update(1.0, *structure_->Dispnp(), 0.0);


    // begin nonlinear solver / outer iteration ***************************

    // set relaxed scalars
    SetScatraSolution(phinp);

    // prepare partitioned structure step
    if (iter_ != 1 and use_old_structure_) structure_->PreparePartitionStep();

    // solve structural system
    DoStructStep();

    // solve scalar transport equation
    DoScatraStep();

    // check convergence for all fields and stop iteration loop if
    // convergence is achieved overall
    stopnonliniter = ConvergenceCheck(iter_);

    // get relaxation parameter
    CalcOmega(omega_, iter_);

    // do the relaxation
    // d^{i+1} = omega^{i+1} . d^{i+1} + (1- omega^{i+1}) d^i
    //         = d^i + omega^{i+1} * ( d^{i+1} - d^i )
    phinp->Update(omega_, *scaincnp_, 1.0);
  }
}

/*----------------------------------------------------------------------*
 | Calculate relaxation parameter                            Thon 12/14 |
 *----------------------------------------------------------------------*/
void SSI::SSI_Part2WC_ScatraToSolid_Relax::CalcOmega(double& omega, const int itnum)
{
  // nothing to do in here since we have a constant relaxation parameter: omega != startomega_;
  if (Comm().MyPID() == 0)
    std::cout << "Fixed relaxation parameter omega is: " << omega << std::endl;
}

/*----------------------------------------------------------------------*
 | Constructor                                               Thon 12/14 |
 *----------------------------------------------------------------------*/
SSI::SSI_Part2WC_ScatraToSolid_Relax_Aitken::SSI_Part2WC_ScatraToSolid_Relax_Aitken(
    const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams)
    : SSI_Part2WC_ScatraToSolid_Relax(comm, globaltimeparams), scaincnpold_(Teuchos::null)
{
  // Keep this constructor empty!
  // First do everything on the more basic objects like the discretizations, like e.g.
  // redistribution of elements. Only then call the setup to this class. This will call the setup to
  // all classes in the inheritance hierarchy. This way, this class may also override a method that
  // is called during Setup() in a base class.
}

/*----------------------------------------------------------------------*
 | Setup this class                                          fang 01/18 |
 *----------------------------------------------------------------------*/
void SSI::SSI_Part2WC_ScatraToSolid_Relax_Aitken::Setup()
{
  // call setup of base class
  SSI::SSI_Part2WC::Setup();

  // setup old scatra increment vector
  scaincnpold_ = LINALG::CreateVector(*scatra_->ScaTraField()->Discretization()->DofRowMap(), true);

  return;
}

/*----------------------------------------------------------------------*
 | Calculate relaxation parameter via Aitken                 Thon 12/14 |
 *----------------------------------------------------------------------*/
void SSI::SSI_Part2WC_ScatraToSolid_Relax_Aitken::CalcOmega(double& omega, const int itnum)
{
  const Teuchos::ParameterList& ssicontrolpart =
      DRT::Problem::Instance()->SSIControlParams().sublist("PARTITIONED");
  ;
  // Get maximal relaxation parameter from input file
  const double maxomega = ssicontrolpart.get<double>("MAXOMEGA");
  // Get minimal relaxation parameter from input file
  const double minomega = ssicontrolpart.get<double>("MINOMEGA");

  // scaincnpdiff =  r^{i+1}_{n+1} - r^i_{n+1}
  Teuchos::RCP<Epetra_Vector> scaincnpdiff =
      LINALG::CreateVector(*scatra_->ScaTraField()->Discretization()->DofRowMap(0), true);
  scaincnpdiff->Update(1.0, *scaincnp_, (-1.0), *scaincnpold_, 0.0);

  double scaincnpdiffnorm = 0.0;
  scaincnpdiff->Norm2(&scaincnpdiffnorm);

  if (scaincnpdiffnorm <= 1e-06 and Comm().MyPID() == 0)
    std::cout << "Warning: The scalar increment is to small in order to use it for Aitken "
                 "relaxation. Using the previous omega instead!"
              << std::endl;

  // calculate dot product
  double scaincsdot = 0.0;  // delsdot = ( r^{i+1}_{n+1} - r^i_{n+1} )^T . r^{i+1}_{n+1}
  scaincnpdiff->Dot(*scaincnp_, &scaincsdot);

  if (itnum != 1 and scaincnpdiffnorm > 1e-06)
  {  // relaxation parameter
    // omega^{i+1} = 1- mu^{i+1} and nu^{i+1} = nu^i + (nu^i -1) . (r^{i+1} - r^i)^T . (-r^{i+1}) /
    // |r^{i+1} - r^{i}|^2 results in
    omega = omega *
            (1 - (scaincsdot) /
                     (scaincnpdiffnorm * scaincnpdiffnorm));  // compare e.g. PhD thesis U. Kuettler

    // we force omega to be in the range defined in the input file
    if (omega < minomega)
    {
      if (Comm().MyPID() == 0)
        std::cout << "Warning: The calculation of the relaxation parameter omega via Aitken did "
                     "lead to a value smaller than MINOMEGA!"
                  << std::endl;
      omega = minomega;
    }
    if (omega > maxomega)
    {
      if (Comm().MyPID() == 0)
        std::cout << "Warning: The calculation of the relaxation parameter omega via Aitken did "
                     "lead to a value bigger than MAXOMEGA!"
                  << std::endl;
      omega = maxomega;
    }
  }

  // else //if itnum==1 nothing is to do here since we want to take the last omega from the previous
  // step
  if (Comm().MyPID() == 0)
    std::cout << "Using Aitken the relaxation parameter omega was estimated to: " << omega
              << std::endl;

  // update history vector old increment r^i_{n+1}
  scaincnpold_->Update(1.0, *scaincnp_, 0.0);
}

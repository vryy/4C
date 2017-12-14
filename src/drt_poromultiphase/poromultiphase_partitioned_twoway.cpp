/*----------------------------------------------------------------------*/
/*!
 \file poromultiphase_partitioned_twoway.cpp

 \brief two-way coupled partitioned solution algorithm
        for porous multiphase flow through elastic medium problems

   \level 3

   \maintainer  Lena Yoshihara
                yoshihara@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
 *----------------------------------------------------------------------*/

#include "poromultiphase_partitioned_twoway.H"

#include "../drt_adapter/ad_porofluidmultiphase_wrapper.H"
#include "../drt_adapter/ad_str_wrapper.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../linalg/linalg_utils.H"
#include "../drt_io/io.H"

/*----------------------------------------------------------------------*
 | constructor                                              vuong 08/16 |
 *----------------------------------------------------------------------*/
POROMULTIPHASE::PoroMultiPhasePartitionedTwoWay::PoroMultiPhasePartitionedTwoWay(
    const Epetra_Comm& comm,
    const Teuchos::ParameterList& globaltimeparams):
    PoroMultiPhasePartitioned(comm, globaltimeparams),
    phiincnp_(Teuchos::null),
    dispincnp_(Teuchos::null),
    fluidphinp_(Teuchos::null),
    fluidphioldnp_(Teuchos::null),
    fluidphiincnp_(Teuchos::null),
    ittol_(0.0),
    omega_(1.0),
    startomega_(1.0),
    omegamin_(1.0),
    omegamax_(1.0),
    itmax_(0),
    itnum_(0),
    writerestartevery_(-1),
    relaxationmethod_(INPAR::POROMULTIPHASE::RelaxationMethods::relaxation_none)
{

}

/*----------------------------------------------------------------------*
 | initialization                                            vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhasePartitionedTwoWay::Init(
    const Teuchos::ParameterList& globaltimeparams,
    const Teuchos::ParameterList& algoparams,
    const Teuchos::ParameterList& structparams,
    const Teuchos::ParameterList& fluidparams,
    const std::string& struct_disname,
    const std::string& fluid_disname,
    bool isale,
    int nds_disp,
    int nds_vel,
    int nds_solidpressure,
    int ndsporofluid_scatra)
{
  //call base class
  POROMULTIPHASE::PoroMultiPhasePartitioned::Init(
      globaltimeparams,
      algoparams,
      structparams,
      fluidparams,
      struct_disname,
      fluid_disname,
      isale,
      nds_disp,
      nds_vel,
      nds_solidpressure,
      ndsporofluid_scatra);

  // initialize increment vectors
  phiincnp_ = LINALG::CreateVector(*FluidField()->DofRowMap(0),true);
  dispincnp_ = LINALG::CreateVector(*StructureField()->DofRowMap(0),true);

  // initialize fluid vectors
  fluidphinp_ = LINALG::CreateVector(*FluidField()->Discretization()->DofRowMap(), true);
  fluidphioldnp_ = LINALG::CreateVector(*FluidField()->Discretization()->DofRowMap(), true);
  fluidphiincnp_ = LINALG::CreateVector(*FluidField()->Discretization()->DofRowMap(), true);
  fluidphiincnpold_ = LINALG::CreateVector(*FluidField()->Discretization()->DofRowMap(), true);

  // Get the parameters for the ConvergenceCheck
  itmax_ = algoparams.get<int>("ITEMAX");
  ittol_ = algoparams.sublist("PARTITIONED").get<double>("CONVTOL");

  // restart
  writerestartevery_ = globaltimeparams.get<int>("RESTARTEVRY");

  // relaxation parameters
  startomega_ = algoparams.sublist("PARTITIONED").get<double>("STARTOMEGA");
  omegamin_ = algoparams.sublist("PARTITIONED").get<double>("MINOMEGA");
  omegamax_ = algoparams.sublist("PARTITIONED").get<double>("MAXOMEGA");

  relaxationmethod_ = DRT::INPUT::IntegralValue<INPAR::POROMULTIPHASE::RelaxationMethods>(algoparams.sublist("PARTITIONED"),"RELAXATION");
}

/*----------------------------------------------------------------------*
 | setup the system if necessary                             vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhasePartitionedTwoWay::SetupSystem()
{
  // Do nothing, just monolithic coupling needs this method
  return;
}

/*----------------------------------------------------------------------*
 | Outer Timeloop  without relaxation                        vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhasePartitionedTwoWay::OuterLoop()
{
  // reset counter
  itnum_ = 0;
  bool stopnonliniter = false;

  if (Comm().MyPID()==0)
  {
    std::cout<<"********************************************************************************" <<
        "***********************************************\n";
    std::cout<<"* PARTITIONED OUTER ITERATION LOOP ----- FLUID  <-------> STRUCTURE         " <<
        "                                                  *\n";
    std::cout<<"* STEP: " << std::setw(5) << std::setprecision(4) << std::scientific << Step() << "/"
        << std::setw(5) << std::setprecision(4) << std::scientific << NStep() << ", Time: "
        << std::setw(11) << std::setprecision(4) << std::scientific << Time() << "/"
        << std::setw(11) << std::setprecision(4) << std::scientific << MaxTime() << ", Dt: "
        << std::setw(11) << std::setprecision(4) << std::scientific << Dt() <<
        "                                                           *"<< std::endl;
  }

  while (stopnonliniter==false)
  {
    // increment number of iteration
    itnum_++;

    // update the states to the last solutions obtained
    IterUpdateStates();

    if(solve_structure_)
    {
      // 1.) solve structural system, note: in first iteration solid pressure has already been set in PrepareTimeStep()
      DoStructStep();
      // 2.) set disp and vel states in porofluid field
      SetStructSolution(StructureField()->Dispnp(),StructureField()->Velnp());
    }
    else
    {
      // Inform user that structure field has been disabled
      PrintStructureDisabledInfo();
      // just set displacements and velocities to zero
      SetStructSolution(struct_zeros_,struct_zeros_);
    }

    // 1.) solve scalar transport equation
    DoFluidStep();

    // perform relaxation
    PerformRelaxation(FluidField()->Phinp(), itnum_);

    // 2.) set fluid solution in structure field
    SetRelaxedFluidSolution();

    // check convergence for all fields
    // stop iteration loop if converged
    stopnonliniter = ConvergenceCheck(itnum_);
  }

  return;
}

/*----------------------------------------------------------------------*
 | convergence check for both fields (copied form tsi)       vuong 08/16 |
 *----------------------------------------------------------------------*/
bool POROMULTIPHASE::PoroMultiPhasePartitionedTwoWay::ConvergenceCheck(int itnum)
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
  // define L2-norm of increments
  double phiincnorm_L2(0.0);
  double phinorm_L2(0.0);
  double dispincnorm_L2(0.0);
  double dispnorm_L2(0.0);

  // build the current scalar increment Inc T^{i+1}
  // \f Delta T^{k+1} = Inc T^{k+1} = T^{k+1} - T^{k}  \f
  phiincnp_->Update(1.0,*(FluidField()->Phinp()),-1.0);
  dispincnp_->Update(1.0,*(StructureField()->Dispnp()),-1.0);

  // build the L2-norm of the scalar increment and the scalar
  phiincnp_->Norm2(&phiincnorm_L2);
  FluidField()->Phinp()->Norm2(&phinorm_L2);
  dispincnp_->Norm2(&dispincnorm_L2);
  StructureField()->Dispnp()->Norm2(&dispnorm_L2);

  // care for the case that there is (almost) zero scalar
  if (phinorm_L2 < 1e-6) phinorm_L2 = 1.0;
  if (dispnorm_L2 < 1e-6) dispnorm_L2 = 1.0;

  // print the incremental based convergence check to the screen
  if (Comm().MyPID()==0 )
  {
    std::cout<<"                                                                                                                              *\n";
    std::cout<<"+----------------------------------------------------------------------------------------------------------------+            *\n";
    std::cout<<"| PARTITIONED OUTER ITERATION STEP ----- FLUID  <-------> STRUCTURE                                              |            *\n";
    printf("+--------------+---------------------+----------------+------------------+--------------------+------------------+            *\n");
    printf("|-  step/max  -|-  tol      [norm]  -|-  fluid-inc   -|-  disp-inc      -|-  fluid-rel-inc   -|-  disp-rel-inc  -|            *\n");
    printf("|   %3d/%3d    |  %10.3E[L_2 ]   |  %10.3E    |  %10.3E      |  %10.3E        |  %10.3E      |",
         itnum,itmax_,ittol_,phiincnorm_L2/Dt()/sqrt(phiincnp_->GlobalLength()),dispincnorm_L2/Dt()/sqrt(dispincnp_->GlobalLength()),phiincnorm_L2/phinorm_L2,dispincnorm_L2/dispnorm_L2);
    printf("            *\n");
    printf("+--------------+---------------------+----------------+------------------+--------------------+------------------+            *\n");
  }

  // converged
  if ( ((phiincnorm_L2/phinorm_L2) <= ittol_) and ((dispincnorm_L2/dispnorm_L2) <= ittol_) and ((dispincnorm_L2/Dt()/sqrt(dispincnp_->GlobalLength()))<=ittol_) and ((phiincnorm_L2/Dt()/sqrt(phiincnp_->GlobalLength()))<=ittol_) )
  {
    stopnonliniter = true;
    if (Comm().MyPID()==0 )
    {
      printf("* FLUID  <-------> STRUCTURE Outer Iteration loop converged after iteration %3d/%3d !                                         *\n", itnum,itmax_);
      printf("*******************************************************************************************************************************\n");
    }
  }

  // stop if itemax is reached without convergence
  // timestep
  if ( (itnum==itmax_) and (((phiincnorm_L2/phinorm_L2) > ittol_) or ((dispincnorm_L2/dispnorm_L2) > ittol_) or ((dispincnorm_L2/Dt()/sqrt(dispincnp_->GlobalLength()))>ittol_) or (phiincnorm_L2/Dt()/sqrt(phiincnp_->GlobalLength()))>ittol_ ) )
  {
    stopnonliniter = true;
    if ((Comm().MyPID()==0) )
    {
      printf("|     >>>>>> not converged in itemax steps!                                                                      |\n");
      printf("+--------------+---------------------+----------------+------------------+--------------------+------------------+\n");
      printf("\n");
      printf("\n");
    }
    dserror("The partitioned solver did not converge in ITEMAX steps!");
  }

  return stopnonliniter;
}

/*----------------------------------------------------------------------*
 | Solve structure filed                                     vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhasePartitionedTwoWay::DoStructStep()
{
  if (Comm().MyPID() == 0)
  {
    std::cout<<"\n";
    std::cout<<"*****************************************************************************************************************\n";
    std::cout<<"STRUCTURE SOLVER   \n";
    std::cout<<"*****************************************************************************************************************\n";

  }

  // Newton-Raphson iteration
  StructureField()-> Solve();

  return;
}

/*----------------------------------------------------------------------*
 | Solve poro fluid field                                    vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhasePartitionedTwoWay::DoFluidStep()
{
  // -------------------------------------------------------------------
  //                  solve nonlinear / linear equation
  // -------------------------------------------------------------------
  FluidField()->Solve();

  return;
}

/*----------------------------------------------------------------------*
 | Set relaxed fluid solution on structure             kremheller 09/16 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhasePartitionedTwoWay::SetRelaxedFluidSolution()
{

  // set fluid solution on structure
  StructureField()->Discretization()->SetState(1,"porofluid",fluidphinp_);

  return;
}

/*----------------------------------------------------------------------*
 | Calculate relaxation parameter omega                kremheller 09/16 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhasePartitionedTwoWay::PerformRelaxation(
    Teuchos::RCP<const Epetra_Vector> phi,
    const int itnum
    )
{
  // get the increment vector
  fluidphiincnp_->Update(1.0,*phi,-1.0,*fluidphioldnp_,0.0);

  // perform relaxation
  switch(relaxationmethod_)
  {
    case INPAR::POROMULTIPHASE::RelaxationMethods::relaxation_none:
    {
      // no relaxation
      omega_ = 1.0;
      break;
    }

    case INPAR::POROMULTIPHASE::RelaxationMethods::relaxation_constant:
    {
      // constant relaxation parameter omega
      omega_ = startomega_;
      if (Comm().MyPID()==0)
        std::cout << "Fixed relaxation parameter omega is: " << omega_ << std::endl;
      break;
    }

    case INPAR::POROMULTIPHASE::RelaxationMethods::relaxation_aitken:
    {
      // Aitken
      AitkenRelaxation(omega_,itnum);
      if (Comm().MyPID()==0)
        std::cout << "Aitken relaxation parameter omega is: " << omega_ << std::endl;
      break;
    }

    default:
    {
      dserror("Relaxation method not yet implemented!");
      break;
    }
  }

  // calculate the relaxed fluid solution as phi,n+1^i+1 = phi,n+1^i + \omega*(phi,n+1^i+1 - phi,n+1^i)
  // note: in first iteration step, omega = 1.0
  fluidphinp_->Update(1.0,*fluidphioldnp_,omega_,*fluidphiincnp_,0.0);

  // save the old fluid solution
  fluidphioldnp_->Update(1.0,*fluidphinp_,0.0);

  return;
}

/*----------------------------------------------------------------------*
 | Perform Aitken relaxation                           kremheller 09/16 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhasePartitionedTwoWay::AitkenRelaxation(
    double&                              omega,
    const int itnum)
{

  // fluidphiincnpdiff =  r^{i+1}_{n+1} - r^i_{n+1}
  Teuchos::RCP<Epetra_Vector> fluidphiincnpdiff = LINALG::CreateVector(*FluidField()->Discretization()->DofRowMap(), true);
  fluidphiincnpdiff->Update(1.0,*fluidphiincnp_,(-1.0),*fluidphiincnpold_,0.0);

  double fluidphiincnpdiffnorm = 0.0;
  fluidphiincnpdiff->Norm2(&fluidphiincnpdiffnorm);

  if (fluidphiincnpdiffnorm <=1e-06 and Comm().MyPID()==0)
    std::cout<<"Warning: The scalar increment is too small in order to use it for Aitken relaxation. Using the previous omega instead!"<<std::endl;

  // calculate dot product
  double fluidphiincsdot = 0.0; //delsdot = ( r^{i+1}_{n+1} - r^i_{n+1} )^T . r^{i+1}_{n+1}
  fluidphiincnpdiff->Dot(*fluidphiincnp_,&fluidphiincsdot);

  if (itnum != 1 and fluidphiincnpdiffnorm > 1e-06)
  {
    // relaxation parameter
    // omega^{i+1} = 1- mu^{i+1} and nu^{i+1} = nu^i + (nu^i -1) . (r^{i+1} - r^i)^T . (-r^{i+1}) / |r^{i+1} - r^{i}|^2 results in
    omega = omega * (1.0 - (fluidphiincsdot)/(fluidphiincnpdiffnorm *fluidphiincnpdiffnorm)); //compare e.g. PhD thesis U. Kuettler

    // we force omega to be in the range defined in the input file
    if (omega < omegamin_)
    {
      if (Comm().MyPID()==0)
        std::cout<<"Warning: The calculation of the relaxation parameter omega via Aitken did lead to a value smaller than MINOMEGA!"<<std::endl;
      omega= omegamin_;
    }
    if (omega > omegamax_)
    {
      if (Comm().MyPID()==0)
        std::cout<<"Warning: The calculation of the relaxation parameter omega via Aitken did lead to a value bigger than MAXOMEGA!"<<std::endl;
      omega= omegamax_;
    }
  }

  // update history vector old increment r^i_{n+1}
  fluidphiincnpold_->Update(1.0,*fluidphiincnp_,0.0);
}

/*----------------------------------------------------------------------*
 | update the current states in every iteration             vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhasePartitionedTwoWay::IterUpdateStates()
{
  // store last solutions (current states).
  // will be compared in ConvergenceCheck to the solutions,
  // obtained from the next Struct and Scatra steps.
  phiincnp_ ->Update(1.0,*FluidField()->Phinp(),0.0);
  dispincnp_->Update(1.0,*StructureField()->Dispnp(),0.0);

  return;
}  // IterUpdateStates()

/*-------------------------------------------------------------------------*
 | read restart information for given time step (public) kremheller 09/17  |
 *-------------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhasePartitionedTwoWay::ReadRestart( int restart )
{
  if (restart)
  {
    // call base class
    POROMULTIPHASE::PoroMultiPhaseBase::ReadRestart(restart);

    IO::DiscretizationReader reader(FluidField()->Discretization(), restart);
    if (restart != reader.ReadInt("step"))
      dserror("Time step on file not equal to given step");

    // get omega_ from restart
    omega_ = reader.ReadDouble("omega_");

    // get omega_ from restart
    reader.ReadVector(fluidphioldnp_,"fluidphioldnp_");
  }


  return;
}

/*----------------------------------------------------------------------*
 | update fields and output results                    kremheller 09/17 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhasePartitionedTwoWay::UpdateAndOutput()
{
  // call base class
  POROMULTIPHASE::PoroMultiPhaseBase::UpdateAndOutput();

  // write interface force and relaxation parameter in restart
  if (writerestartevery_ and Step() % writerestartevery_ == 0)
  {
    FluidField()->Discretization()->Writer()->WriteDouble("omega_",omega_);
    FluidField()->Discretization()->Writer()->WriteVector("fluidphioldnp_",fluidphioldnp_);
  }
}

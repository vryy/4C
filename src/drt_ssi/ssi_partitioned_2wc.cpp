/*!------------------------------------------------------------------------------------------------*
 \file ssi_partitioned_2wc.cpp

 \brief two way coupled partitioned scalar structure interaction

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
 </pre>
 *------------------------------------------------------------------------------------------------*/

#include "ssi_partitioned_2wc.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../linalg/linalg_utils.H"

#include "../drt_adapter/ad_str_wrapper.H"
#include "../drt_adapter/adapter_scatra_base_algorithm.H"

#include "../drt_scatra/scatra_timint_implicit.H"

/*----------------------------------------------------------------------*
 | constructor                                               Thon 12/14 |
 *----------------------------------------------------------------------*/
SSI::SSI_Part2WC::SSI_Part2WC(const Epetra_Comm& comm,
    const Teuchos::ParameterList& globaltimeparams,
    const Teuchos::ParameterList& scatraparams,
    const Teuchos::ParameterList& structparams,
    const std::string struct_disname,
    const std::string scatra_disname)
  : SSI_Part(comm, globaltimeparams, scatraparams, structparams,struct_disname,scatra_disname),
    scaincnp_(LINALG::CreateVector(*scatra_->ScaTraField()->Discretization()->DofRowMap(0),true)),
    dispincnp_(LINALG::CreateVector(*structure_->DofRowMap(0),true)),
    itnum_(0)
{

  // call the SSI parameter lists
  const Teuchos::ParameterList& ssicontrol = DRT::Problem::Instance()->SSIControlParams();
  const Teuchos::ParameterList& ssicontrolpart = DRT::Problem::Instance()->SSIControlParams().sublist("PARTITIONED");

  if (DRT::INPUT::IntegralValue<int>(ssicontrol, "DIFFTIMESTEPSIZE")){
    dserror("Different time stepping for two way coupling not implemented yet.");
  }

  // Get the parameters for the ConvergenceCheck
  itmax_ = ssicontrol.get<int>("ITEMAX"); // default: =10
  ittol_ = ssicontrolpart.get<double>("CONVTOL"); // default: =1e-6

  //do some checks
  {
    INPAR::STR::DynamicType structtimealgo = DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(structparams,"DYNAMICTYP");
    if ( structtimealgo == INPAR::STR::dyna_statics )
      dserror("If you use statics as the structural time integrator no velocities will be calculated and hence"
          "the deformations will not be applied to the scalar transport problem!");

    INPAR::SCATRA::ConvForm convform
    = DRT::INPUT::IntegralValue<INPAR::SCATRA::ConvForm>(scatraparams,"CONVFORM");
    if ( not (convform == INPAR::SCATRA::convform_conservative) )
      dserror("If the scalar transport problem is solved on the deforming domain, the conservative form must be"
          "used to include volume changes! Set 'CONVFORM' to 'conservative' in the SCALAR TRANSPORT DYNAMIC section!");
  }
}

/*----------------------------------------------------------------------*
 | Timeloop for 2WC SSI problems                             Thon 12/14 |
 *----------------------------------------------------------------------*/
void SSI::SSI_Part2WC::Timeloop()
{
  //initial output
  structure_-> PrepareOutput();
  structure_-> Output();
  SetStructSolution(structure_->Dispnp(),structure_->Velnp());
  scatra_->ScaTraField()->Output();

  //time loop
  while (NotFinished())
  {
    PrepareTimeStep();

    OuterLoop();

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
    std::cout
        << "\n***********************\n STRUCTURE SOLVER \n***********************\n";
  }

  // Newton-Raphson iteration
  structure_-> Solve();
}

/*----------------------------------------------------------------------*
 | Solve Scatra field                                        Thon 12/14 |
 *----------------------------------------------------------------------*/
void SSI::SSI_Part2WC::DoScatraStep()
{
  if (Comm().MyPID() == 0)
  {
    std::cout
        << "\n***********************\n  TRANSPORT SOLVER \n***********************\n";
  }

  // -------------------------------------------------------------------
  //                  solve nonlinear / linear equation
  // -------------------------------------------------------------------
  scatra_->ScaTraField()->Solve();

}

/*----------------------------------------------------------------------*/
//prepare time step
/*----------------------------------------------------------------------*/
void SSI::SSI_Part2WC::PrepareTimeStep(bool printheader)
{
  IncrementTimeAndStep();
  if(printheader)
    PrintHeader();

  SetScatraSolution(scatra_->ScaTraField()->Phin());
  structure_-> PrepareTimeStep();

  SetStructSolution(structure_->Dispn(),structure_->Veln());
  scatra_->ScaTraField()->PrepareTimeStep();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSI_Part2WC::UpdateAndOutput()
{
  structure_-> PrepareOutput();

  structure_-> Update();
  scatra_->ScaTraField()->Update();

  scatra_->ScaTraField()->EvaluateErrorComparedToAnalyticalSol();

  structure_-> Output();
  scatra_->ScaTraField()->Output();
}


/*----------------------------------------------------------------------*
 | Outer Timeloop for 2WC SSi without relaxation
 *----------------------------------------------------------------------*/
void SSI::SSI_Part2WC::OuterLoop()
{
  int  itnum_ = 0;
  bool stopnonliniter = false;

  if (Comm().MyPID()==0)
  {
    std::cout<<"\n****************************************\n          OUTER ITERATION LOOP\n****************************************\n";
  }

  while (stopnonliniter==false)
  {
    itnum_++;

    // store scalar from first solution for convergence check (like in
    // elch_algorithm: use current values)
    scaincnp_->Update(1.0,*scatra_->ScaTraField()->Phinp(),0.0);
    dispincnp_->Update(1.0,*structure_->Dispnp(),0.0);

    // set structure-based scalar transport values
    SetScatraSolution(scatra_->ScaTraField()->Phinp());

    if(itnum_!=1)
      structure_->PreparePartitionStep();
    // solve structural system
    DoStructStep();

    // set mesh displacement and velocity fields
    SetStructSolution(structure_->Dispnp(),structure_->Velnp());

    // solve scalar transport equation
    DoScatraStep();
    //ScatraEvaluateSolveIterUpdate();

    // check convergence for all fields and stop iteration loop if
    // convergence is achieved overall
    stopnonliniter = ConvergenceCheck(itnum_);
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
  //     | scalar+1 |_2
  //
  // AND
  //
  //    | scalar increment |_2
  //  -------------------------------- < Tolerance
  //             dt * n

  // variables to save different L2 - Norms
  // define L2-norm of incremental scalar and scalar
  double scaincnorm_L2(0.0);
  double scanorm_L2(0.0);
  double dispincnorm_L2(0.0);
  double dispnorm_L2(0.0);

  // build the current scalar increment Inc T^{i+1}
  // \f Delta T^{k+1} = Inc T^{k+1} = T^{k+1} - T^{k}  \f
  scaincnp_->Update(1.0,*(scatra_->ScaTraField()->Phinp()),-1.0);
  dispincnp_->Update(1.0,*(structure_->Dispnp()),-1.0);

  // build the L2-norm of the scalar increment and the scalar
  scaincnp_->Norm2(&scaincnorm_L2);
  scatra_->ScaTraField()->Phinp()->Norm2(&scanorm_L2);
  dispincnp_->Norm2(&dispincnorm_L2);
  structure_->Dispnp()->Norm2(&dispnorm_L2);

  // care for the case that there is (almost) zero scalar
  if (scanorm_L2 < 1e-6) scanorm_L2 = 1.0;
  if (dispnorm_L2 < 1e-6) dispnorm_L2 = 1.0;

  // print the incremental based convergence check to the screen
  if (Comm().MyPID()==0 )
  {
    std::cout<<"\n";
    std::cout<<"***********************************************************************************\n";
    std::cout<<"    OUTER ITERATION STEP    \n";
    std::cout<<"***********************************************************************************\n";
    printf("+--------------+---------------------+----------------+------------------+--------------------+------------------+\n");
    printf("|-  step/max  -|-  tol      [norm]  -|-  scalar-inc  -|-  disp-inc      -|-  scalar-rel-inc  -|-  disp-rel-inc  -|\n");
    printf("|   %3d/%3d    |  %10.3E[L_2 ]   |  %10.3E    |  %10.3E      |  %10.3E        |  %10.3E      |",
         itnum,itmax_,ittol_,scaincnorm_L2/Dt()/sqrt(scaincnp_->GlobalLength()),dispincnorm_L2/Dt()/sqrt(dispincnp_->GlobalLength()),scaincnorm_L2/scanorm_L2,dispincnorm_L2/dispnorm_L2);
    printf("\n");
    printf("+--------------+---------------------+----------------+------------------+--------------------+------------------+\n");
  }

  // converged
  if ( ((scaincnorm_L2/scanorm_L2) <= ittol_) and ((dispincnorm_L2/dispnorm_L2) <= ittol_) and ((dispincnorm_L2/Dt()/sqrt(dispincnp_->GlobalLength()))<=ittol_) and ((scaincnorm_L2/Dt()/sqrt(scaincnp_->GlobalLength()))<=ittol_) )
  {
    stopnonliniter = true;
    if (Comm().MyPID()==0 )
    {
      printf("|  Outer Iteration loop converged after iteration %3d/%3d !                                                      |\n", itnum,itmax_);
      printf("+--------------+---------------------+----------------+------------------+--------------------+------------------+\n");
    }
  }

  // stop if itemax is reached without convergence
  // timestep
  if ( (itnum==itmax_) and (((scaincnorm_L2/scanorm_L2) > ittol_) or ((dispincnorm_L2/dispnorm_L2) > ittol_) or ((dispincnorm_L2/Dt()/sqrt(dispincnp_->GlobalLength()))>ittol_) or (scaincnorm_L2/Dt()/sqrt(scaincnp_->GlobalLength()))>ittol_ ) )
  {
    stopnonliniter = true;
    if ((Comm().MyPID()==0) )
    {
      printf("|     >>>>>> not converged in itemax steps!                                                                      |\n");
      printf("+--------------+---------------------+----------------+------------------+--------------------+------------------+\n");
      printf("\n");
      printf("\n");
    }
    dserror("The partitioned SSI solver did not converge in ITEMAX steps!");
  }

  return stopnonliniter;
}

/*----------------------------------------------------------------------*
 | calculate velocities by a FD approximation                 Thon 14/11 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> SSI::SSI_Part2WC::CalcVelocity(
  Teuchos::RCP<const Epetra_Vector> dispnp
  )
{
  Teuchos::RCP<Epetra_Vector> vel = Teuchos::null;
  // copy D_n onto V_n+1
  vel = Teuchos::rcp(new Epetra_Vector( *(structure_->Dispn()) ) );
  // calculate velocity with timestep Dt()
  //  V_n+1^k = (D_n+1^k - D_n) / Dt
  vel->Update(1./Dt(), *dispnp, -1./Dt());

  return vel;
}  // CalcVelocity()

/*----------------------------------------------------------------------*
 | Constructor                                               Thon 12/14 |
 *----------------------------------------------------------------------*/
SSI::SSI_Part2WC_SolidToScatra_Relax::SSI_Part2WC_SolidToScatra_Relax(const Epetra_Comm& comm,
    const Teuchos::ParameterList& globaltimeparams,
    const Teuchos::ParameterList& scatraparams,
    const Teuchos::ParameterList& structparams,
    const std::string struct_disname,
    const std::string scatra_disname)
  : SSI_Part2WC(comm, globaltimeparams, scatraparams, structparams, struct_disname, scatra_disname)
{
  const Teuchos::ParameterList& ssicontrolpart = DRT::Problem::Instance()->SSIControlParams().sublist("PARTITIONED");

  //Get minimal relaxation parameter from input file
  omega_ = ssicontrolpart.get<double>("STARTOMEGA");
}

/*----------------------------------------------------------------------*
 | Outer Timeloop for 2WC SSi with relaxed displacements     Thon 12/14 |
 *----------------------------------------------------------------------*/
void SSI::SSI_Part2WC_SolidToScatra_Relax::OuterLoop()
{
  int  itnum = 0;
  bool stopnonliniter = false;

  if (Comm().MyPID()==0)
  {
    std::cout<<"\n****************************************\n          OUTER ITERATION LOOP\n****************************************\n";
  }

  //these are the relaxed inputs
  Teuchos::RCP<Epetra_Vector> dispnp
    = LINALG::CreateVector(*(structure_->DofRowMap(0)), true);
  Teuchos::RCP<Epetra_Vector> velnp
    = LINALG::CreateVector(*(structure_->DofRowMap(0)), true);

  while (stopnonliniter==false)
  {
    itnum++;

    if (itnum == 1)
    {
      dispnp->Update(1.0,*(structure_->Dispnp()),0.0); //TSI does Dispn()
      velnp->Update(1.0,*(structure_->Velnp()),0.0);
    }

    // store scalars and displacements for the convergence check later
    scaincnp_->Update(1.0,*scatra_->ScaTraField()->Phinp(),0.0);
    dispincnp_->Update(1.0,*structure_->Dispnp(),0.0);


    // begin nonlinear solver / outer iteration ***************************

    // set relaxed mesh displacements and velocity field
    SetStructSolution( dispnp , velnp);

    // solve scalar transport equation
    DoScatraStep();

    // set scalar fields
    SetScatraSolution(scatra_->ScaTraField()->Phinp());

    //prepare a partitioned structure step
    if(itnum!=1)
      structure_->PreparePartitionStep();

    //update variable whic is going to be relaxed later
    dispnp->Update(1.0,*structure_()->Dispnp(),0.0);

    // solve structural system
    DoStructStep();

    // end nonlinear solver / outer iteration *****************************

    // check convergence for all fields and stop iteration loop if
    // convergence is achieved overall
    stopnonliniter = ConvergenceCheck(itnum);

    // get relaxation parameter
    CalcOmega(omega_,itnum);

    // do the relaxation
    // d^{i+1} = omega^{i+1} . d^{i+1} + (1- omega^{i+1}) d^i
    //         = d^i + omega^{i+1} * ( d^{i+1} - d^i )
    dispnp->Update(omega_,*dispincnp_,1.0);

    // since the velocity field has to fit to the relaxated displacements we also have to relaxate them.
    // since the velocity depends nonlinear on the displacements we can just approximate them via finite differences here.
    velnp = CalcVelocity(dispnp);
  }

}

/*----------------------------------------------------------------------*
 | Calculate relaxation parameter                            Thon 12/14 |
 *----------------------------------------------------------------------*/
void SSI::SSI_Part2WC_SolidToScatra_Relax::CalcOmega(double& omega, const int itnum)
{
  //nothing to do in here since we have a constant relaxation parameter: omega != startomega_;
  if (Comm().MyPID()==0 )
    std::cout<<"Fixed relaxation parameter omega is: " <<omega<<std::endl;
}

/*----------------------------------------------------------------------*
 | Constructor                                               Thon 12/14 |
 *----------------------------------------------------------------------*/
SSI::SSI_Part2WC_SolidToScatra_Relax_Aitken::SSI_Part2WC_SolidToScatra_Relax_Aitken(const Epetra_Comm& comm,
    const Teuchos::ParameterList& globaltimeparams,
    const Teuchos::ParameterList& scatraparams,
    const Teuchos::ParameterList& structparams,
    const std::string struct_disname,
    const std::string scatra_disname)
  : SSI_Part2WC_SolidToScatra_Relax(comm, globaltimeparams, scatraparams, structparams, struct_disname, scatra_disname),
    del_(LINALG::CreateVector(*structure_->DofRowMap(0),true)),
    delhist_(LINALG::CreateVector(*structure_->DofRowMap(0),true))
{
  const Teuchos::ParameterList& ssicontrolpart = DRT::Problem::Instance()->SSIControlParams().sublist("PARTITIONED");;

  //Get maximal relaxation parameter from input file
  maxomega_ = ssicontrolpart.get<double>("MAXOMEGA");

  //Get minimal relaxation parameter from input file
  minomega_ = ssicontrolpart.get<double>("MINOMEGA");

  del_->PutScalar(1.0e05);
  delhist_->PutScalar(0.0);
}

/*----------------------------------------------------------------------*
 | Calculate relaxation parameter via Aitken                 Thon 12/14 |
 *----------------------------------------------------------------------*/
void SSI::SSI_Part2WC_SolidToScatra_Relax_Aitken::CalcOmega(double& omega, const int itnum)
{
  // calculate difference of current (i+1) and old (i) residual vector
  // delhist = ( r^{i+1}_{n+1} - r^i_{n+1} )
  // update history vector old increment r^i_{n+1}

  delhist_->Update(1.0,*dispincnp_,(-1.0),*del_,0.0);  // update r^{i+1}_{n+1} - r^i_{n+1}

  // del_ = r^{i+1}_{n+1} = Phi^{i+1}_{n+1} - Phi^{i}_{n+1}
  del_->Update(1.0,*dispincnp_,0.0);

  double delhistnorm = 0.0;
  delhist_->Norm2(&delhistnorm);
  if ( delhistnorm <=1e-05 and Comm().MyPID()==0 )
    std::cout<<"Warning: The structure increment is to small in order to use it for Aitken relaxation. Using the previous Omega instead!"<<std::endl;

  // calculate dot product
  double delsdot = 0.0; //delsdot = ( r^{i+1}_{n+1} - r^i_{n+1} )^T . r^{i+1}_{n+1}
  delhist_->Dot(*del_,&delsdot);

  if (itnum != 1 and delhistnorm > 1e-05)
  { // relaxation parameter
    // omega^{i+1} = 1- mu^{i+1} and nu^{i+1} = nu^i + (nu^i -1) . (r^{i+1} - r^i)^T . (-r^{i+1}) / |r^{i+1} - r^{i}|^2 results in
    omega = omega*(1  - (delsdot)/(delhistnorm * delhistnorm)); //compare e.g. PhD thesis U. Kuettler

    //we force omega to be in the range defined in the input file
    if (omega < minomega_)
    {
      if (Comm().MyPID()==0)
        std::cout<<"Warning: The calculation of the relaxation parameter omega via Aitken did lead to a value smaller than MINOMEGA!"<<std::endl;
      omega= minomega_;
    }
    if (omega > maxomega_)
    {
      if (Comm().MyPID()==0)
        std::cout<<"Warning: The calculation of the relaxation parameter omega via Aitken did lead to a value bigger than MAXOMEGA!"<<std::endl;
      omega= maxomega_;
    }
  }

  //else //if itnum==1 nothing is to do here since we want to take the last omega from the previous step
  if (Comm().MyPID()==0 )
    std::cout<<"Using Aitken the relaxation parameter omega was estimated to: " <<omega<<std::endl;
}


/*----------------------------------------------------------------------*
 | Constructor                                               Thon 12/14 |
 *----------------------------------------------------------------------*/
SSI::SSI_Part2WC_ScatraToSolid_Relax::SSI_Part2WC_ScatraToSolid_Relax(const Epetra_Comm& comm,
    const Teuchos::ParameterList& globaltimeparams,
    const Teuchos::ParameterList& scatraparams,
    const Teuchos::ParameterList& structparams,
    const std::string struct_disname,
    const std::string scatra_disname)
  : SSI_Part2WC(comm, globaltimeparams, scatraparams, structparams, struct_disname, scatra_disname)
{
  const Teuchos::ParameterList& ssicontrolpart = DRT::Problem::Instance()->SSIControlParams().sublist("PARTITIONED");

  //Get start relaxation parameter from input file
  omega_ = ssicontrolpart.get<double>("STARTOMEGA");

  if (Comm().MyPID()==0)
  {
    std::cout<<"\n#########################################################################\n  "<<std::endl;
    std::cout<<"The ScatraToSolid relaxations are not well tested . Keep your eyes open!  "<<std::endl;
    std::cout<<"\n#########################################################################\n  "<<std::endl;
  }
}

/*----------------------------------------------------------------------*
 | Outer Timeloop for 2WC SSi with relaxed scalar             Thon 12/14 |
 *----------------------------------------------------------------------*/
void SSI::SSI_Part2WC_ScatraToSolid_Relax::OuterLoop()
{
  int  itnum = 0;
  bool stopnonliniter = false;

  if (Comm().MyPID()==0)
  {
    std::cout<<"\n****************************************\n          OUTER ITERATION LOOP\n****************************************\n";
  }

  //this is the relaxed input
  Teuchos::RCP<Epetra_Vector> phinp
    = LINALG::CreateVector(*scatra_->ScaTraField()->Discretization()->DofRowMap(0), true);

  while (stopnonliniter==false)
  {
    itnum++;

    if (itnum == 1)
    {
      phinp->Update(1.0,*(scatra_->ScaTraField()->Phinp()),0.0); //TSI does Dispn()
    }

    // store scalars and displacements for the convergence check later
    scaincnp_->Update(1.0,*scatra_->ScaTraField()->Phinp(),0.0);
    dispincnp_->Update(1.0,*structure_->Dispnp(),0.0);


    // begin nonlinear solver / outer iteration ***************************

    // set relaxed scalars
    SetScatraSolution(phinp);

    //prepare partioned tructure step
    if(itnum!=1)
      structure_->PreparePartitionStep();

    // solve structural system
    DoStructStep();

    // set mesh displacement and velocity fields
    SetStructSolution( structure_->Dispnp() , structure_->Velnp() );

    //update variable whic is going to be relaxed later
    phinp->Update(1.0,*scatra_->ScaTraField()->Phinp(),0.0);

    // solve scalar transport equation
    DoScatraStep();

    // check convergence for all fields and stop iteration loop if
    // convergence is achieved overall
    stopnonliniter = ConvergenceCheck(itnum);

    // get relaxation parameter
    CalcOmega(omega_,itnum);

    // do the relaxation
    // d^{i+1} = omega^{i+1} . d^{i+1} + (1- omega^{i+1}) d^i
    //         = d^i + omega^{i+1} * ( d^{i+1} - d^i )
    phinp->Update(omega_,*scaincnp_,1.0);
    }
}

/*----------------------------------------------------------------------*
 | Calculate relaxation parameter                            Thon 12/14 |
 *----------------------------------------------------------------------*/
void SSI::SSI_Part2WC_ScatraToSolid_Relax::CalcOmega(double& omega, const int itnum)
{
  //nothing to do in here since we have a constant relaxation parameter: omega != startomega_;
  if (Comm().MyPID()==0 )
    std::cout<<"Fixed relaxation parameter omega is: " <<omega<<std::endl;
}

/*----------------------------------------------------------------------*
 | Constructor                                               Thon 12/14 |
 *----------------------------------------------------------------------*/
SSI::SSI_Part2WC_ScatraToSolid_Relax_Aitken::SSI_Part2WC_ScatraToSolid_Relax_Aitken(const Epetra_Comm& comm,
    const Teuchos::ParameterList& globaltimeparams,
    const Teuchos::ParameterList& scatraparams,
    const Teuchos::ParameterList& structparams,
    const std::string struct_disname,
    const std::string scatra_disname)
  : SSI_Part2WC_ScatraToSolid_Relax(comm, globaltimeparams, scatraparams, structparams, struct_disname, scatra_disname),
    del_(LINALG::CreateVector(*scatra_->ScaTraField()->Discretization()->DofRowMap(0),true)),
    delhist_(LINALG::CreateVector(*scatra_->ScaTraField()->Discretization()->DofRowMap(0),true))
{
  const Teuchos::ParameterList& ssicontrolpart = DRT::Problem::Instance()->SSIControlParams().sublist("PARTITIONED");;

  //Get maximal relaxation parameter from input file
  maxomega_ = ssicontrolpart.get<double>("MAXOMEGA");

  //Get minimal relaxation parameter from input file
  minomega_ = ssicontrolpart.get<double>("MINOMEGA");

  del_->PutScalar(1.0e05);
  delhist_->PutScalar(0.0);
}

/*----------------------------------------------------------------------*
 | Calculate relaxation parameter via Aitken                 Thon 12/14 |
 *----------------------------------------------------------------------*/
void SSI::SSI_Part2WC_ScatraToSolid_Relax_Aitken::CalcOmega(double& omega, const int itnum)
{
  //delhist =  r^{i+1}_{n+1} - r^i_{n+1}
  delhist_->Update(1.0,*scaincnp_,(-1.0),*del_,0.0);

  // del_ = r^{i+1}_{n+1} = Phi^{i+1}_{n+1} - Phi^{i}_{n+1}
  del_->Update(1.0,*scaincnp_,0.0);

  double delhistnorm = 0.0;
  delhist_->Norm2(&delhistnorm);

  if (delhistnorm <=1e-05 and Comm().MyPID()==0 )
    std::cout<<"Warning: The scalar increment is to small in order to use it for Aitken relaxation. Using the previous omega instead!"<<std::endl;

  // calculate dot product
  double delsdot = 0.0; //delsdot = ( r^{i+1}_{n+1} - r^i_{n+1} )^T . r^{i+1}_{n+1}
  delhist_->Dot(*del_,&delsdot);

  if (itnum != 1 and delhistnorm > 1e-05)
  { // relaxation parameter
    // omega^{i+1} = 1- mu^{i+1} and nu^{i+1} = nu^i + (nu^i -1) . (r^{i+1} - r^i)^T . (-r^{i+1}) / |r^{i+1} - r^{i}|^2 results in
    omega = omega*(1  - (delsdot)/(delhistnorm * delhistnorm)); //compare e.g. PhD thesis U. Kuettler

    //we force omega to be in the range defined in the input file
    if (omega < minomega_)
    {
      if (Comm().MyPID()==0)
        std::cout<<"Warning: The calculation of the relaxation parameter omega via Aitken did lead to a value smaller than MINOMEGA!"<<std::endl;
      omega= minomega_;
    }
    if (omega > maxomega_)
    {
      if (Comm().MyPID()==0)
        std::cout<<"Warning: The calculation of the relaxation parameter omega via Aitken did lead to a value bigger than MAXOMEGA!"<<std::endl;
      omega= maxomega_;
    }
  }

  //else //if itnum==1 nothing is to do here since we want to take the last omega from the previous step
  if (Comm().MyPID()==0 )
    std::cout<<"Using Aitken the relaxation parameter omega was estimated to: " <<omega<<std::endl;

}

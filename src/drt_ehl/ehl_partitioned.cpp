/*--------------------------------------------------------------------------*/
/*!
\file ehl_partitioned.cpp

\brief class for partitioned elastohydrodynamic lubrication (lubrication structure interaction)

<pre>
Maintainer: Andy Wirtz
            wirtz@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15270
</pre>
*/
/*--------------------------------------------------------------------------*/


#include "../linalg/linalg_utils.H"

#include "../drt_lib/drt_globalproblem.H"

#include "../drt_adapter/ad_str_wrapper.H"
#include "../drt_adapter/adapter_lubrication.H"

#include "../drt_lubrication/lubrication_timint_implicit.H"

#include "ehl_partitioned.H"

/*----------------------------------------------------------------------*
 | constructor                                              wirtz 12/15 |
 *----------------------------------------------------------------------*/
EHL::Partitioned::Partitioned(const Epetra_Comm& comm,
    const Teuchos::ParameterList& globaltimeparams,
    const Teuchos::ParameterList& lubricationparams,
    const Teuchos::ParameterList& structparams,
    const std::string struct_disname,
    const std::string lubrication_disname)
  : Base(comm, globaltimeparams,lubricationparams,structparams,struct_disname,lubrication_disname),
    preincnp_(LINALG::CreateVector(*lubrication_->LubricationField()->Discretization()->DofRowMap(0),true)),
    dispincnp_(LINALG::CreateVector(*structure_->DofRowMap(0),true))
{

  // call the EHL parameter lists
  const Teuchos::ParameterList& ehlparams = DRT::Problem::Instance()->ElastoHydroDynamicParams();
  const Teuchos::ParameterList& ehlparamspart = DRT::Problem::Instance()->ElastoHydroDynamicParams().sublist("PARTITIONED");

  if (DRT::INPUT::IntegralValue<int>(ehlparams, "DIFFTIMESTEPSIZE")){
    dserror("Different time stepping for two way coupling not implemented yet.");
  }

  // Get the parameters for the ConvergenceCheck
  itmax_ = ehlparams.get<int>("ITEMAX"); // default: =10
  ittol_ = ehlparamspart.get<double>("CONVTOL"); // default: =1e-6

}

/*----------------------------------------------------------------------*
 | Timeloop for EHL problems                                wirtz 12/15 |
 *----------------------------------------------------------------------*/
void EHL::Partitioned::Timeloop()
{
  while (NotFinished())
  {
    PrepareTimeStep();

    OuterLoop();

    UpdateAndOutput();

  }
}


/*----------------------------------------------------------------------*
 | prepare time step                                        wirtz 12/15 |
 *----------------------------------------------------------------------*/
void EHL::Partitioned::PrepareTimeStep()
{
  IncrementTimeAndStep();
  PrintHeader();

  SetStructSolution(structure_->Dispn(),structure_->Veln());
  structure_-> PrepareTimeStep();
//  SetLubricationSolution(lubrication_->LubricationField()->Quantity()); // todo: what quantity
  lubrication_->LubricationField()->PrepareTimeStep();
}


/*----------------------------------------------------------------------*
 | outer Timeloop for EHL without relaxation                wirtz 12/15 |
 *----------------------------------------------------------------------*/
void EHL::Partitioned::OuterLoop()
{
  int  itnum = 0;
  bool stopnonliniter = false;

  if (Comm().MyPID()==0)
  {
    std::cout<<"\n****************************************\n          OUTER ITERATION LOOP\n****************************************\n";
  }

  while (stopnonliniter==false)
  {
    itnum++;

    // store pressure from first solution for convergence check (like in
    // elch_algorithm: use current values)
    preincnp_->Update(1.0,*lubrication_->LubricationField()->Prenp(),0.0);
    dispincnp_->Update(1.0,*structure_->Dispnp(),0.0);

    // set structure-based pressure values
    SetLubricationSolution(lubrication_->LubricationField()->Prenp());

    if(itnum!=1)
      structure_->PreparePartitionStep();
    // solve structural system
    DoStructStep();

    // set mesh displacement and velocity fields
    SetStructSolution(structure_->Dispnp(),structure_->Velnp());

    // solve lubrication equation
    DoLubricationStep();
    //LubricationEvaluateSolveIterUpdate();

    // check convergence for all fields and stop iteration loop if
    // convergence is achieved overall
    stopnonliniter = ConvergenceCheck(itnum);
  }

  return;
}


/*----------------------------------------------------------------------*
 | constructor                                              wirtz 12/15 |
 *----------------------------------------------------------------------*/
void EHL::Partitioned::UpdateAndOutput()
{
  structure_-> PrepareOutput();

  structure_-> Update();
  lubrication_->LubricationField()->Update();

  lubrication_->LubricationField()->EvaluateErrorComparedToAnalyticalSol();

  structure_-> Output();
  lubrication_->LubricationField()->Output();
}


/*----------------------------------------------------------------------*
 | solve structure filed                                    wirtz 12/15 |
 *----------------------------------------------------------------------*/
void EHL::Partitioned::DoStructStep()
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
 | solve Lubrication field                                  wirtz 12/15 |
 *----------------------------------------------------------------------*/
void EHL::Partitioned::DoLubricationStep()
{
  if (Comm().MyPID() == 0)
  {
    std::cout
        << "\n***********************\n  LUBRICATION SOLVER \n***********************\n";
  }

  // -------------------------------------------------------------------
  //                           solve nonlinear
  // -------------------------------------------------------------------
  lubrication_->LubricationField()->Solve();

}


/*----------------------------------------------------------------------*
 | convergence check for both fields (lubrication & structure)          |
 |                                                          wirtz 12/15 |
 *----------------------------------------------------------------------*/
bool EHL::Partitioned::ConvergenceCheck(int itnum)
{

  // convergence check based on the pressure increment
  bool stopnonliniter = false;

  //    | pressure increment |_2
  //  -------------------------------- < Tolerance
  //     | pressure+1 |_2
  //
  // AND
  //
  //    | pressure increment |_2
  //  -------------------------------- < Tolerance
  //             dt * n

  // variables to save different L2 - Norms
  // define L2-norm of incremental pressure and pressure
  double preincnorm_L2(0.0);
  double prenorm_L2(0.0);
  double dispincnorm_L2(0.0);
  double dispnorm_L2(0.0);

  // build the current pressure increment Inc T^{i+1}
  // \f Delta T^{k+1} = Inc T^{k+1} = T^{k+1} - T^{k}  \f
  preincnp_->Update(1.0,*(lubrication_->LubricationField()->Prenp()),-1.0);
  dispincnp_->Update(1.0,*(structure_->Dispnp()),-1.0);

  // build the L2-norm of the pressure increment and the pressure
  preincnp_->Norm2(&preincnorm_L2);
  lubrication_->LubricationField()->Prenp()->Norm2(&prenorm_L2);
  dispincnp_->Norm2(&dispincnorm_L2);
  structure_->Dispnp()->Norm2(&dispnorm_L2);

  // care for the case that there is (almost) zero pressure
  if (prenorm_L2 < 1e-6) prenorm_L2 = 1.0;
  if (dispnorm_L2 < 1e-6) dispnorm_L2 = 1.0;

  // print the incremental based convergence check to the screen
  if (Comm().MyPID()==0 )
  {
    std::cout<<"\n";
    std::cout<<"***********************************************************************************\n";
    std::cout<<"    OUTER ITERATION STEP    \n";
    std::cout<<"***********************************************************************************\n";
    printf("+--------------+---------------------+------------------+-----------------+----------------------+------------------+\n");
    printf("|-  step/max  -|-  tol      [norm]  -|-  pressure-inc  -|  disp-inc      -|-  pressure-rel-inc  -|-  disp-rel-inc  -|\n");
    printf("|   %3d/%3d    |  %10.3E[L_2 ]   |  %10.3E      |  %10.3E     |  %10.3E          |  %10.3E      |",
         itnum,itmax_,ittol_,preincnorm_L2/Dt()/sqrt(preincnp_->GlobalLength()),dispincnorm_L2/Dt()/sqrt(dispincnp_->GlobalLength()),preincnorm_L2/prenorm_L2,dispincnorm_L2/dispnorm_L2);
    printf("\n");
    printf("+--------------+---------------------+------------------+-----------------+----------------------+------------------+\n");
  }

  // converged
  if ( ((preincnorm_L2/prenorm_L2) <= ittol_) and ((dispincnorm_L2/dispnorm_L2) <= ittol_) and ((dispincnorm_L2/Dt()/sqrt(dispincnp_->GlobalLength()))<=ittol_) and ((preincnorm_L2/Dt()/sqrt(preincnp_->GlobalLength()))<=ittol_) )
  {
    stopnonliniter = true;
    if (Comm().MyPID()==0 )
    {
      printf("\n");
      printf("|  Outer Iteration loop converged after iteration %3d/%3d !                                                         |\n", itnum,itmax_);
      printf("+--------------+---------------------+------------------+-----------------+----------------------+------------------+\n");
    }
  }

  // stop if itemax is reached without convergence
  // timestep
  if ( (itnum==itmax_) and (((preincnorm_L2/prenorm_L2) > ittol_) or ((dispincnorm_L2/dispnorm_L2) > ittol_) or ((dispincnorm_L2/Dt()/sqrt(dispincnp_->GlobalLength()))>ittol_) or (preincnorm_L2/Dt()/sqrt(preincnp_->GlobalLength()))>ittol_ ) )
  {
    stopnonliniter = true;
    if ((Comm().MyPID()==0) )
    {
      printf("|     >>>>>> not converged in itemax steps!                                                                         |\n");
      printf("+--------------+---------------------+------------------+-----------------+----------------------+------------------+\n");
      printf("\n");
      printf("\n");
    }
    dserror("The partitioned EHL solver did not converge in ITEMAX steps!");
  }

  return stopnonliniter;
}


/*----------------------------------------------------------------------*
 | calculate velocities by a FD approximation               wirtz 12/15 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> EHL::Partitioned::CalcVelocity(
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
 | constructor                                              wirtz 12/15 |
 *----------------------------------------------------------------------*/
EHL::Partitioned_StrToLub_Relax::Partitioned_StrToLub_Relax(const Epetra_Comm& comm,
    const Teuchos::ParameterList& globaltimeparams,
    const Teuchos::ParameterList& lubricationparams,
    const Teuchos::ParameterList& structparams,
    const std::string struct_disname,
    const std::string lubrication_disname)
  : Partitioned(comm, globaltimeparams, lubricationparams, structparams, struct_disname, lubrication_disname)
{
  const Teuchos::ParameterList& ehlparamspart = DRT::Problem::Instance()->ElastoHydroDynamicParams().sublist("PARTITIONED");

  //Get minimal relaxation parameter from input file
  omega_ = ehlparamspart.get<double>("STARTOMEGA");
}

/*----------------------------------------------------------------------*
 | outer timeloop for EHL with relaxed displacements        wirtz 12/15 |
 *----------------------------------------------------------------------*/
void EHL::Partitioned_StrToLub_Relax::OuterLoop()
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
      dispnp->Update(1.0,*(structure_->Dispnp()),0.0);
      velnp->Update(1.0,*(structure_->Velnp()),0.0);
    }

    // store pressures and displacements for the convergence check later
    preincnp_->Update(1.0,*lubrication_->LubricationField()->Prenp(),0.0);
    dispincnp_->Update(1.0,*structure_->Dispnp(),0.0);


    // begin nonlinear solver / outer iteration ***************************

    // set relaxed mesh displacements and velocity field
    SetStructSolution( dispnp , velnp);

    // solve lubrication equation
    DoLubricationStep();

    // set pressure fields
    SetLubricationSolution(lubrication_->LubricationField()->Prenp());

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
 | calculate relaxation parameter                           wirtz 12/15 |
 *----------------------------------------------------------------------*/
void EHL::Partitioned_StrToLub_Relax::CalcOmega(double& omega, const int itnum)
{
  //nothing to do in here since we have a constant relaxation parameter: omega != startomega_;
  if (Comm().MyPID()==0 )
    std::cout<<"Fixed relaxation parameter omega is: " <<omega<<std::endl;
}

/*----------------------------------------------------------------------*
 | constructor                                              wirtz 12/15 |
 *----------------------------------------------------------------------*/
EHL::Partitioned_StrToLub_Aitken::Partitioned_StrToLub_Aitken(const Epetra_Comm& comm,
    const Teuchos::ParameterList& globaltimeparams,
    const Teuchos::ParameterList& lubricationparams,
    const Teuchos::ParameterList& structparams,
    const std::string struct_disname,
    const std::string lubrication_disname)
  : Partitioned_StrToLub_Relax(comm, globaltimeparams, lubricationparams, structparams, struct_disname, lubrication_disname),
    del_(LINALG::CreateVector(*structure_->DofRowMap(0),true)),
    delhist_(LINALG::CreateVector(*structure_->DofRowMap(0),true))
{
  const Teuchos::ParameterList& ehlparamspart = DRT::Problem::Instance()->ElastoHydroDynamicParams().sublist("PARTITIONED");;

  //Get maximal relaxation parameter from input file
  maxomega_ = ehlparamspart.get<double>("MAXOMEGA");

  //Get minimal relaxation parameter from input file
  minomega_ = ehlparamspart.get<double>("MINOMEGA");

  del_->PutScalar(1.0e05);
  delhist_->PutScalar(0.0);
}

/*----------------------------------------------------------------------*
 | calculate relaxation parameter via Aitken                wirtz 12/15 |
 *----------------------------------------------------------------------*/
void EHL::Partitioned_StrToLub_Aitken::CalcOmega(double& omega, const int itnum)
{
  // calculate difference of current (i+1) and old (i) residual vector
  // delhist = ( r^{i+1}_{n+1} - r^i_{n+1} )
  // update history vector old increment r^i_{n+1}

  delhist_->Update(1.0,*dispincnp_,(-1.0),*del_,0.0);  // update r^{i+1}_{n+1} - r^i_{n+1}

  // del_ = r^{i+1}_{n+1} = Disp^{i+1}_{n+1} - Disp^{i}_{n+1}
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
 | constructor                                              wirtz 12/15 |
 *----------------------------------------------------------------------*/
EHL::Partitioned_LubToStr_Relax::Partitioned_LubToStr_Relax(const Epetra_Comm& comm,
    const Teuchos::ParameterList& globaltimeparams,
    const Teuchos::ParameterList& lubricationparams,
    const Teuchos::ParameterList& structparams,
    const std::string struct_disname,
    const std::string lubrication_disname)
  : Partitioned(comm, globaltimeparams, lubricationparams, structparams, struct_disname, lubrication_disname)
{
  const Teuchos::ParameterList& ehlparamspart = DRT::Problem::Instance()->ElastoHydroDynamicParams().sublist("PARTITIONED");

  //Get start relaxation parameter from input file
  omega_ = ehlparamspart.get<double>("STARTOMEGA");

  if (Comm().MyPID()==0)
  {
    std::cout<<"\n#########################################################################\n  "<<std::endl;
    std::cout<<"The LubToStr relaxations are not well tested . Keep your eyes open!  "<<std::endl;
    std::cout<<"\n#########################################################################\n  "<<std::endl;
  }
}

/*----------------------------------------------------------------------*
 | outer timeloop for EHL with relaxed pressure             wirtz 12/15 |
 *----------------------------------------------------------------------*/
void EHL::Partitioned_LubToStr_Relax::OuterLoop()
{
  int  itnum = 0;
  bool stopnonliniter = false;

  if (Comm().MyPID()==0)
  {
    std::cout<<"\n****************************************\n          OUTER ITERATION LOOP\n****************************************\n";
  }

  //this is the relaxed input
  Teuchos::RCP<Epetra_Vector> prenp
    = LINALG::CreateVector(*lubrication_->LubricationField()->Discretization()->DofRowMap(0), true);

  while (stopnonliniter==false)
  {
    itnum++;

    if (itnum == 1)
    {
      prenp->Update(1.0,*(lubrication_->LubricationField()->Prenp()),0.0);
    }

    // store pressures and displacements for the convergence check later
    preincnp_->Update(1.0,*lubrication_->LubricationField()->Prenp(),0.0);
    dispincnp_->Update(1.0,*structure_->Dispnp(),0.0);


    // begin nonlinear solver / outer iteration ***************************

    // set relaxed pressures
    SetLubricationSolution(prenp);

    //prepare partioned tructure step
    if(itnum!=1)
      structure_->PreparePartitionStep();

    // solve structural system
    DoStructStep();

    // set mesh displacement and velocity fields
    SetStructSolution( structure_->Dispnp() , structure_->Velnp() );

    //update variable whic is going to be relaxed later
    prenp->Update(1.0,*lubrication_->LubricationField()->Prenp(),0.0);

    // solve lubrication equation
    DoLubricationStep();

    // check convergence for all fields and stop iteration loop if
    // convergence is achieved overall
    stopnonliniter = ConvergenceCheck(itnum);

    // get relaxation parameter
    CalcOmega(omega_,itnum);

    // do the relaxation
    // d^{i+1} = omega^{i+1} . d^{i+1} + (1- omega^{i+1}) d^i
    //         = d^i + omega^{i+1} * ( d^{i+1} - d^i )
    prenp->Update(omega_,*preincnp_,1.0);
    }
}

/*----------------------------------------------------------------------*
 | calculate relaxation parameter                           wirtz 12/15 |
 *----------------------------------------------------------------------*/
void EHL::Partitioned_LubToStr_Relax::CalcOmega(double& omega, const int itnum)
{
  //nothing to do in here since we have a constant relaxation parameter: omega != startomega_;
  if (Comm().MyPID()==0 )
    std::cout<<"Fixed relaxation parameter omega is: " <<omega<<std::endl;
}

/*----------------------------------------------------------------------*
 | constructor                                              wirtz 12/15 |
 *----------------------------------------------------------------------*/
EHL::Partitioned_LubToStr_Aitken::Partitioned_LubToStr_Aitken(const Epetra_Comm& comm,
    const Teuchos::ParameterList& globaltimeparams,
    const Teuchos::ParameterList& lubricationparams,
    const Teuchos::ParameterList& structparams,
    const std::string struct_disname,
    const std::string lubrication_disname)
  : Partitioned_LubToStr_Relax(comm, globaltimeparams, lubricationparams, structparams, struct_disname, lubrication_disname),
    del_(LINALG::CreateVector(*lubrication_->LubricationField()->Discretization()->DofRowMap(0),true)),
    delhist_(LINALG::CreateVector(*lubrication_->LubricationField()->Discretization()->DofRowMap(0),true))
{
  const Teuchos::ParameterList& ehlparamspart = DRT::Problem::Instance()->ElastoHydroDynamicParams().sublist("PARTITIONED");;

  //Get maximal relaxation parameter from input file
  maxomega_ = ehlparamspart.get<double>("MAXOMEGA");

  //Get minimal relaxation parameter from input file
  minomega_ = ehlparamspart.get<double>("MINOMEGA");

  del_->PutScalar(1.0e05);
  delhist_->PutScalar(0.0);
}

/*----------------------------------------------------------------------*
 | calculate relaxation parameter via Aitken                wirtz 12/15 |
 *----------------------------------------------------------------------*/
void EHL::Partitioned_LubToStr_Aitken::CalcOmega(double& omega, const int itnum)
{
  //delhist =  r^{i+1}_{n+1} - r^i_{n+1}
  delhist_->Update(1.0,*preincnp_,(-1.0),*del_,0.0);

  // del_ = r^{i+1}_{n+1} = Pre^{i+1}_{n+1} - Pre^{i}_{n+1}
  del_->Update(1.0,*preincnp_,0.0);

  double delhistnorm = 0.0;
  delhist_->Norm2(&delhistnorm);

  if (delhistnorm <=1e-05 and Comm().MyPID()==0 )
    std::cout<<"Warning: The pressure increment is to small in order to use it for Aitken relaxation. Using the previous omega instead!"<<std::endl;

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


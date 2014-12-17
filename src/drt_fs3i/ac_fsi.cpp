/*!----------------------------------------------------------------------
\file ac_fsi.H
\brief H-file associated with algorithmic routines for two-way coupled partitioned
       solution approaches to fluid-structure-scalar-scalar interaction
       (FS3I). Specifically related version for multiscale approches

<pre>
Maintainers: Moritz Thon
             thon@mhpc.mw.tum.de
             089/289-10364
</pre>

*----------------------------------------------------------------------*/


#include "ac_fsi.H"

#include "../drt_fsi/fsi_monolithic.H"
#include "../drt_scatra/scatra_algorithm.H"
#include "../drt_inpar/inpar_scatra.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_adapter/ad_str_fsiwrapper.H"
#include "../drt_adapter/ad_fld_fluid_fsi.H"
#include "../drt_adapter/ad_ale_fsi.H"
#include "../linalg/linalg_utils.H"

/*----------------------------------------------------------------------*
 | constructor                                               Thon 12/14 |
 *----------------------------------------------------------------------*/
FS3I::ACFSI::ACFSI(const Epetra_Comm& comm)
  :PartFS3I(comm)
{
  const Teuchos::ParameterList& fs3idyn = DRT::Problem::Instance()->FS3IDynamicParams();
  const Teuchos::ParameterList& fs3idynpart = fs3idyn.sublist("PARTITIONED");
  const Teuchos::ParameterList& fs3idynac = fs3idyn.sublist("AC");

  // get input parameters for partioned FS3I problems
  ittol_ = fs3idynpart.get<double>("CONVTOL");
  itmax_ = fs3idynpart.get<int>("ITEMAX");
  couplingalgo_ = DRT::INPUT::IntegralValue<INPAR::FS3I::SolutionSchemeOverFields>(fs3idynpart,"COUPALGO");

  // get input parameters for AC FS3I problems
  fsiperssisteps_ = fs3idynac.get<int>("FSISTEPSPERSCATRASTEP");
  //  if (fsiperssisteps_ != 1) //at this point this case is already been dealed by Manipulating the FSI parameter list

  /// initialize increment vectors
  structureincrement_ = LINALG::CreateVector(*fsi_->StructureField()->DofRowMap(0),true);
  fluidincrement_ = LINALG::CreateVector(*fsi_->FluidField()->DofRowMap(0),true);
  aleincrement_ = LINALG::CreateVector(*fsi_->AleField()->DofRowMap(),true);

  // build a proxy of the scatra discretization for the structure field
  Teuchos::RCP<DRT::DofSet> scatradofset
    = scatravec_[1]->ScaTraField()->Discretization()->GetDofSetProxy();

  // check if scatra field has 2 discretizations, so that coupling is possible
  if (fsi_->StructureField()->Discretization()->AddDofSet(scatradofset)!=1)
    dserror("unexpected dof sets in structure field");
}

/*----------------------------------------------------------------------*
 | Timeloop                                                  Thon 12/14 |
 *----------------------------------------------------------------------*/
void FS3I::ACFSI::Timeloop()
{
  // output of initial state
  ScatraOutput();
  fsi_->PrepareTimeloop();

  while (NotFinished())
  {
    PrepareTimeStep();
    OuterLoop();
    UpdateAndOutput();
  }

}

/*----------------------------------------------------------------------*
 | Prepare time step                                         Thon 12/14 |
 *----------------------------------------------------------------------*/
void FS3I::ACFSI::PrepareTimeStep()
{
  IncrementTimeAndStep();

  fsi_->PrepareTimeStep();

  // prepare time step for both fluid- and structure-based scatra field
  for (unsigned i=0; i<scatravec_.size(); ++i)
  {
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra = scatravec_[i];
    scatra->ScaTraField()->PrepareTimeStep();
  }
}

/*----------------------------------------------------------------------*
 | OuterLoop                                                 Thon 12/14 |
 *----------------------------------------------------------------------*/
void FS3I::ACFSI::OuterLoop()
{
  switch (couplingalgo_)
  {
  case INPAR::FS3I::fs3i_SequStagg:
    OuterLoopSequStagg();
    break;
  case INPAR::FS3I::fs3i_IterStagg:
    OuterLoopIterStagg();
    break;
  default:
    dserror("partitioned FS3I coupling scheme not implemented!");
    break;
  }
}

/*----------------------------------------------------------------------*
 | OuterLoop for sequentially staggered FS3I scheme          Thon 12/14 |
 *----------------------------------------------------------------------*/
void FS3I::ACFSI::OuterLoopSequStagg()
{
  SetStructScatraSolution();

  DoFSIStep();

  SetFSISolution();

  DoScatraStep();
}

/*----------------------------------------------------------------------*
 | OuterLoop for iterative staggered FS3I scheme             Thon 12/14 |
 *----------------------------------------------------------------------*/
void FS3I::ACFSI::OuterLoopIterStagg()
{
  int itnum = 0;

  bool stopnonliniter = false;

  while ( stopnonliniter==false )
  {
    itnum++;

    structureincrement_->Update(1.0,*fsi_->StructureField()->Dispnp(),0.0);
    fluidincrement_->Update(1.0,*fsi_->FluidField()->Velnp(),0.0);
    aleincrement_->Update(1.0,*fsi_->AleField()->Dispnp(),0.0);

    if (Comm().MyPID()==0)
    {
      std::cout<<"\n***********************************************************************************"
                 "\n                   OUTER FS3I ITERATION LOOP "<<itnum <<"/"<<itmax_ <<" START..."
                 "\n***********************************************************************************"<<std::endl;
    }

    SetStructScatraSolution();

    DoFSIStep();

    SetFSISolution();

    DoScatraStep();

    stopnonliniter = PartFs3iConvergenceCkeck(itnum);
  }

}

/*----------------------------------------------------------------------*
 | Convergence check for iterative staggered FS3I scheme     Thon 12/14 |
 *----------------------------------------------------------------------*/
bool FS3I::ACFSI::PartFs3iConvergenceCkeck( const int itnum)
{
  // convergence check based on the scalar increment
  bool stopnonliniter = false;


  //calculate fsi increments. scatra increment is already done in the scatra fields convergence check
  structureincrement_->Update(1.0,*fsi_->StructureField()->Dispnp(),-1.0);
  fluidincrement_->Update(1.0,*fsi_->FluidField()->Velnp(),-1.0);
  aleincrement_->Update(1.0,*fsi_->AleField()->Dispnp(),-1.0);

  //L2-norms of increment vectors
  double scatraincconnorm_L2(0.0);
  scatraincrement_->Norm2(&scatraincconnorm_L2);
  double structureincconnorm_L2(0.0);
  structureincrement_->Norm2(&structureincconnorm_L2);
  double fluidincconnorm_L2(0.0);
  fluidincrement_->Norm2(&fluidincconnorm_L2);
  double aleincconnorm_L2(0.0);
  aleincrement_->Norm2(&aleincconnorm_L2);

  // set up vector of absolute concentrations
  Teuchos::RCP<Epetra_Vector> scatra = Teuchos::rcp(new Epetra_Vector(scatraincrement_->Map()));
  Teuchos::RCP<const Epetra_Vector> scatra1 = scatravec_[0]->ScaTraField()->Phinp(); //fluidscatra
  Teuchos::RCP<const Epetra_Vector> scatra2 = scatravec_[1]->ScaTraField()->Phinp(); //structurescatra
  SetupCoupledScatraVector(scatra,scatra1,scatra2);

  //norms of solution vectors
  double scatranorm_L2(0.0);
  scatra->Norm2(&scatranorm_L2);
  if (scatranorm_L2 < 1e-05) scatranorm_L2=1.0;
  double structurenorm_L2(0.0);
  fsi_->StructureField()->Dispnp()->Norm2(&structurenorm_L2);
  if (structurenorm_L2 < 1e-05) structurenorm_L2=1.0;
  double fluidnorm_L2(0.0);
  fsi_->FluidField()->Velnp()->Norm2(&fluidnorm_L2);
  if (fluidnorm_L2 < 1e-05) fluidnorm_L2=1.0;
  double alenorm_L2(0.0);
  fsi_->AleField()->Dispnp()->Norm2(&alenorm_L2);
  if (alenorm_L2 < 1e-05) alenorm_L2=1.0;

  // print the incremental based convergence check to the screen
  if (Comm().MyPID()==0 )
  {
    std::cout<<"\n***********************************************************************************\n"
               "               OUTER ITERATION STEP  "<<itnum<<"/"<<itmax_<<" ...DONE\n"<<std::endl;

    printf("+--------------+---------------------+----------------+---------------+---------------+--------------+\n");
    printf("|   step/max   |   tol      [norm]   |   scalar-inc   |   disp-inc    |   vel-inc     |   ale-inc    |\n");
    printf("|   %3d/%3d    |  %10.3E[L_2 ]   |  %10.3E    |  %10.3E   |  %10.3E   | %10.3E   |\n",
         itnum,itmax_,ittol_,scatraincconnorm_L2/scatranorm_L2,structureincconnorm_L2/structurenorm_L2,fluidincconnorm_L2/fluidnorm_L2,aleincconnorm_L2/alenorm_L2);
    printf("+--------------+---------------------+----------------+---------------+---------------+--------------+\n");
    std::cout<<"***********************************************************************************\n";
  }

  if ( scatraincconnorm_L2/scatranorm_L2 <= ittol_ and structureincconnorm_L2/structurenorm_L2 <= ittol_
     and fluidincconnorm_L2/fluidnorm_L2 <= ittol_ and aleincconnorm_L2/alenorm_L2 <= ittol_ ) //converged!
   {
     stopnonliniter = true;
   }
  else if (itnum == itmax_)
  {
    stopnonliniter = true;
    if ((Comm().MyPID()==0) )
    {
      std::cout<<"\n***********************************************************************************\n"
                 "               OUTER ITERATION STEP NOT CONVERGED IN ITEMAX"
                 "\n***********************************************************************************\n"<<std::endl;
    }
    dserror("The partitioned FS3I solver did not converged in ITEMAX steps!");
  }

  return stopnonliniter;
}

/*----------------------------------------------------------------------*
 | Read restart                                              Thon 12/14 |
 *----------------------------------------------------------------------*/
void FS3I::ACFSI::ReadRestart()
{
  // read restart information, set vectors and variables
  // (Note that dofmaps might have changed in a redistribution call!)
  const int restart = DRT::Problem::Instance()->Restart();
  if (restart)
  {
    fsi_->ReadRestart(restart*fsiperssisteps_);

    for (unsigned i=0; i<scatravec_.size(); ++i)
    {
      Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> currscatra = scatravec_[i];
      currscatra->ScaTraField()->ReadRestart(restart);
    }

    time_ = fsi_->FluidField()->Time();
    step_ = fsi_->FluidField()->Step()/fsiperssisteps_;
  }
}

/*----------------------------------------------------------------------*
 | Update and output everything                              Thon 12/14 |
 *----------------------------------------------------------------------*/
void FS3I::ACFSI::UpdateAndOutput()
{

  fsi_->PrepareOutput();
  fsi_->Update();
  fsi_->Output();

  UpdateScatraFields();
  ScatraOutput();
}


/*----------------------------------------------------------------------*
 | Do a single fsi step (including all subcycles)            Thon 12/14 |
 *----------------------------------------------------------------------*/
void FS3I::ACFSI::DoFSIStep()
{
  if (fsiperssisteps_ == 1) //no subcycling
  {
    if (Comm().MyPID()==0)
    {
      std::cout<<"\n************************************************************************"
                 "\n                               FSI SOLVER "
                 "\n************************************************************************"<<std::endl;
    }
    fsi_->TimeStep(fsi_);
  }
  else //subcycling
  {
    double subcyclingiter = 0;
    while (subcyclingiter < fsiperssisteps_)
    {
      subcyclingiter++;

      if (subcyclingiter != 1) //in the first step we have already done this
      {
        fsi_->PrepareOutput();
        fsi_->Update();
        fsi_->PrepareTimeStep();
      }

      if (Comm().MyPID()==0)
      {
        std::cout<<"\n************************************************************************"
                   "\n                     FSI SUBCYCLING SOLVER " << subcyclingiter <<"/"<<fsiperssisteps_
                 <<"\n************************************************************************"<<std::endl;
      }

      fsi_->TimeStep(fsi_); //dt_ is already adapted!
    }
  }
}


/*----------------------------------------------------------------------*
 | Do a single scatra step                                   Thon 12/14 |
 *----------------------------------------------------------------------*/
void FS3I::ACFSI::DoScatraStep()
{
  if (Comm().MyPID()==0)
  {
    std::cout<<"\n************************************************************************"
               "\n                           AC TRANSPORT SOLVER"
               "\n************************************************************************\n"<<std::endl;

    std::cout<<"+- step/max -+- tol ---- [norm] -+-- scal-res --+-- scal-inc --+"<<std::endl;
  }

  bool stopnonliniter=false;
  int itnum = 0;

  while (stopnonliniter==false)
  {
    ScatraEvaluateSolveIterUpdate();
    itnum++;
    if (ScatraConvergenceCheck(itnum))
      break;
  }

}

/*----------------------------------------------------------------------*
 | check convergence of scatra fields                        Thon 12/14 |
 *----------------------------------------------------------------------*/
bool FS3I::ACFSI::ScatraConvergenceCheck(const int itnum)
{
  // some input parameters for the scatra fields
  const Teuchos::ParameterList& scatradyn = DRT::Problem::Instance()->ScalarTransportDynamicParams();
  const int scatraitemax = scatradyn.sublist("NONLINEAR").get<int>("ITEMAX");
  const double scatraittol = scatradyn.sublist("NONLINEAR").get<double>("CONVTOL");
  const double scatraabstolres = scatradyn.sublist("NONLINEAR").get<double>("ABSTOLRES");


  double conresnorm(0.0);
  scatrarhs_->Norm2(&conresnorm);
  double incconnorm(0.0);
  scatraincrement_->Norm2(&incconnorm);

  // set up vector of absolute concentrations
  Teuchos::RCP<Epetra_Vector> con = Teuchos::rcp(new Epetra_Vector(scatraincrement_->Map()));
  Teuchos::RCP<const Epetra_Vector> scatra1 = scatravec_[0]->ScaTraField()->Phinp(); //fluidscatra
  Teuchos::RCP<const Epetra_Vector> scatra2 = scatravec_[1]->ScaTraField()->Phinp(); //structurescatra
  SetupCoupledScatraVector(con,scatra1,scatra2);

  double connorm(0.0);
  con->Norm2(&connorm);

  // care for the case that nothing really happens in the concentration field
  if (connorm < 1e-5)
    connorm = 1.0;

  // print the screen info
  if (Comm().MyPID()==0)
  {
    printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E   | %10.3E   |\n",
           itnum,scatraitemax,scatraittol,conresnorm,incconnorm/connorm);
  }

  // this is the convergence check
  // We always require at least one solve. We test the L_2-norm of the
  // current residual. Norm of residual is just printed for information
  if (conresnorm <= scatraittol and incconnorm/connorm <= scatraittol)
  {
    if (Comm().MyPID()==0)
    {
      // print 'finish line'
      printf("+------------+-------------------+--------------+--------------+\n");
    }
    return true;
  }

  // abort iteration, when there's nothing more to do! -> more robustness
  else if (conresnorm < scatraabstolres)
  {
    // print 'finish line'
    if (Comm().MyPID()==0)
    {
      printf("+------------+-------------------+--------------+--------------+\n");
    }
    return true;
  }

  // if itemax is reached without convergence stop the simulation
  else if (itnum == scatraitemax)
  {
    if (Comm().MyPID()==0)
    {
      printf("+---------------------------------------------------------------+\n");
      printf("|            >>>>>> not converged in itemax steps!              |\n");
      printf("+---------------------------------------------------------------+\n");
    }
    // yes, we stop!
    dserror("Scatra not converged in itemax steps!");
    return true;
  }
  else
    return false;

}


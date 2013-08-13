/*----------------------------------------------------------------------*/
/*!
 \file poro_scatra_monolithic.cpp

 \brief

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
 </pre>
 *----------------------------------------------------------------------*/

#include "poro_scatra_monolithic.H"

#include "poro_base.H"
#include "../drt_adapter/adapter_scatra_base_algorithm.H"

#include "../drt_adapter/ad_fld_poro.H"
#include "../drt_adapter/ad_str_fsiwrapper.H"

#include "../drt_fluid_ele/fluid_ele_action.H"

#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_mapextractor.H"
#include "../linalg/linalg_blocksparsematrix.H"

#include "../drt_io/io_control.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_assemblestrategy.H"

#include <Teuchos_TimeMonitor.hpp>

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
POROELAST::PORO_SCATRA_Mono::PORO_SCATRA_Mono(
    const Epetra_Comm&            comm,
    const Teuchos::ParameterList& timeparams):
    PORO_SCATRA_Base(comm,timeparams),
    printscreen_(true), // ADD INPUT PARAMETER
    printiter_(true), // ADD INPUT PARAMETER
    printerrfile_(true), // ADD INPUT PARAMETER FOR 'true'
    errfile_(DRT::Problem::Instance()->ErrorFile()->Handle()),
    timer_(comm),
    iterinc_(Teuchos::null),
    zeros_(Teuchos::null),
    systemmatrix_(Teuchos::null),
    rhs_(Teuchos::null),
    directsolve_(true)
{
  // the problem is two way coupled, thus each discretization must know the other discretization
  AddDofSets();

  const Teuchos::ParameterList& sdynparams
  = DRT::Problem::Instance()->StructuralDynamicParams();

  //some solver paramaters are red form the structure dynamic list (this is not the best way to do it ...)
  solveradapttol_= (DRT::INPUT::IntegralValue<int>(sdynparams, "ADAPTCONV") == 1);
  solveradaptolbetter_ = (sdynparams.get<double> ("ADAPTCONV_BETTER"));
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Mono::Timeloop()
{
  while (NotFinished())
  {
    DoTimeStep();
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Mono::DoTimeStep()
{
  // counter and print header
  // predict solution of both field (call the adapter)
  PrepareTimeStep();

  // Newton-Raphson iteration
  Solve();

  // calculate stresses, strains, energies
  PrepareOutput();

  // update all single field solvers
  Update();

  // write output to screen and files
  Output();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Mono::ReadRestart(int restart)
{
  // read restart information, set vectors and variables
  // (Note that dofmaps might have changed in a redistribution call!)
  if (restart)
  {
    SetScatraSolution();
    SetPoroSolution();

    PoroField()->ReadRestart(restart);
    ScatraField().ReadRestart(restart);

    //in case of submeshes, we need to rebuild the subproxies, also (they are reset during restart)
    if(PoroField()->HasSubmeshes())
      AddDofSets(true);

    // the variables need to be set on other field
    SetScatraSolution();
    SetPoroSolution();

    //second restart needed due to two way coupling.
    ScatraField().ReadRestart(restart);
    PoroField()->ReadRestart(restart);

    //in case of submeshes, we need to rebuild the subproxies, also (they are reset during restart)
    if(PoroField()->HasSubmeshes())
      AddDofSets(true);

    SetTimeStep(PoroField()->Time(), restart);
  }
}

/*----------------------------------------------------------------------*/
//prepare time step
/*----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Mono::PrepareTimeStep()
{
  // the global control routine has its own time_ and step_ variables, as well as the single fields
  // keep them in sinc!
  IncrementTimeAndStep();

  ScatraField().PrepareTimeStep();
  // set structure-based scalar transport values
  SetScatraSolution();

  PoroField()-> PrepareTimeStep();
  SetPoroSolution();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Mono::PrepareOutput()
{
  PoroField()->PrepareOutput();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Mono::Output()
{
  PoroField()->Output();
  ScatraField().Output();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Mono::Update()
{
  PoroField()->Update();

  ScatraField().Update();
  ScatraField().EvaluateErrorComparedToAnalyticalSol();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Mono::Solve()
{

  // we do a Newton-Raphson iteration here.
  // the specific time integration has set the following
  // --> On #rhs_ is the positive force residuum
  // --> On #systemmatrix_ is the effective dynamic tangent matrix

  //initialize norms and increment
  SetupNewton();

  //---------------------------------------------- iteration loop

  // equilibrium iteration loop (loop over k)
  while (((not Converged()) and (iter_ <= itermax_)) or (iter_ <= itermin_))
  {
    timer_.ResetStartTime();
    Epetra_Time timer(Comm());

    // compute residual forces #rhs_ and tangent #tang_
    // whose components are globally oriented
    // build linear system stiffness matrix and rhs/force residual for each
    // field, here e.g. for structure field: field want the iteration increment
    // 1.) Update(iterinc_),
    // 2.) EvaluateForceStiffResidual(),
    // 3.) PrepareSystemForNewtonSolve()
    Evaluate(iterinc_);
    //cout << "  time for Evaluate diagonal blocks: " << timer.ElapsedTime() << "\n";
    //timer.ResetStartTime();

    // check whether we have a sanely filled tangent matrix
    if (not systemmatrix_->Filled())
    {
      dserror("Effective tangent matrix must be filled here");
    }

    // create full monolithic rhs vector
    SetupRHS(iter_==1);
    //cout << "  time for Evaluate SetupRHS: " << timer.ElapsedTime() << "\n";
    //timer.ResetStartTime();

    // (Newton-ready) residual with blanked Dirichlet DOFs (see adapter_timint!)
    // is done in PrepareSystemForNewtonSolve() within Evaluate(iterinc_)
    LinearSolve();
    //cout << "  time for Evaluate LinearSolve: " << timer.ElapsedTime() << "\n";
    //timer.ResetStartTime();

    // reset solver tolerance
    solver_->ResetTolerance();

    //build norms
    BuildCovergenceNorms();

    // print stuff
    PrintNewtonIter();

    // increment equilibrium loop index
    iter_ += 1;
  } // end equilibrium loop

  //---------------------------------------------- iteration loop

  // correct iteration counter
  iter_ -= 1;

  // test whether max iterations was hit
  if ( (Converged()) and (Comm().MyPID()==0) )
  {
    PrintNewtonConv();
  }
  else if (iter_ >= itermax_)
  {
    dserror("Newton unconverged in %d iterations", iter_);
  }

}    // Solve()

/*----------------------------------------------------------------------*
 | evaluate the single fields                              vuong 01/12   |
 *----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Mono::Evaluate(Teuchos::RCP<const Epetra_Vector> x)
{
  TEUCHOS_FUNC_TIME_MONITOR("PORO_SCATRA_Mono::Monolithic::Evaluate");

  // displacement and fluid velocity & pressure incremental vector
  Teuchos::RCP<const Epetra_Vector> poroinc;
  Teuchos::RCP<const Epetra_Vector> scatrainc;

  // if an increment vector exists
  if (x != Teuchos::null)
  {
    // process structure unknowns of the first field
    poroinc = Extractor().ExtractVector(x, 0);

    // process fluid unknowns of the second field
    scatrainc = Extractor().ExtractVector(x, 1);
  }

  // Newton update of the fluid field
  // update velocities and pressures before passed to the structural field
  //  UpdateIterIncrementally(fx),
  ScatraField().UpdateIter(scatrainc);

  // call all elements and assemble rhs and matrices
  /// structural field

  // structure Evaluate (builds tangent, residual and applies DBC)
  //Epetra_Time timerstructure(Comm());

  // apply current velocity and pressures to structure
  SetScatraSolution();

  // Monolithic Poroelasticity accesses the linearised structure problem:
  PoroField()->Evaluate(poroinc);
  //cout << "  structure time for calling Evaluate: " << timerstructure.ElapsedTime() << "\n";
  /// fluid field

  // fluid Evaluate
  // (builds tangent, residual and applies DBC and recent coupling values)
  //Epetra_Time timerfluid(Comm());

  SetScatraSolution();

  // monolithic Poroelasticity accesses the linearised fluid problem
  ScatraField().PrepareLinearSolve();
  //cout << "  fluid time for calling Evaluate: " << timerfluid.ElapsedTime() << "\n";

  // fill off diagonal blocks and build monolithic system matrix
  SetupSystemMatrix();

}

/*----------------------------------------------------------------------*
 | setup system (called in poro_dyn.cpp)                 vuong 01/12    |
 *----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Mono::SetupSystem()
{

  //setup the poro subsystem first
  PoroField()->SetupSystem();

  // create combined map
  std::vector<Teuchos::RCP<const Epetra_Map> > vecSpaces;

  {
    if(PoroField()->DofRowMap() == Teuchos::null)
      dserror("could not access DofRowMap of poro field!");
    // use its own DofRowMap, that is the 0th map of the discretization
    vecSpaces.push_back(PoroField()->DofRowMap());
    const Epetra_Map* dofrowmapscatra = (ScatraField().Discretization())->DofRowMap(0);
    vecSpaces.push_back(Teuchos::rcp(dofrowmapscatra,false));
  }

  if (vecSpaces[0]->NumGlobalElements() == 0)
    dserror("No structure equation. Panic.");
  if (vecSpaces[1]->NumGlobalElements()==0)
    dserror("No fluid equation. Panic.");

  // build dof row map of monolithic system
  SetDofRowMaps(vecSpaces);

  // build dbc map of monolithic system
  BuildCombinedDBCMap();

  // initialize Poroelasticity-systemmatrix_
  systemmatrix_ = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
                                      Extractor(),
                                      Extractor(),
                                      81,
                                      false,
                                      true));

  std::vector<Teuchos::RCP<const Epetra_Map> > scatravecSpaces;
  {
    const Epetra_Map* dofrowmapscatra = (ScatraField().Discretization())->DofRowMap(0);
    scatravecSpaces.push_back(Teuchos::rcp(dofrowmapscatra,false));
    scatrarowdofmap_.Setup(*dofrowmapscatra,scatravecSpaces);
  }

  k_sp_ = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
            scatrarowdofmap_,
            Extractor(),
            81,
            false,
            true));
  k_ps_ = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
                              scatrarowdofmap_,
                              Extractor(),
                              81,
                              false,
                              true));

  k_pss_ = Teuchos::rcp(new LINALG::SparseMatrix(
                      *(PoroField()->StructureField()->DofRowMap()), 81, true, true));
  k_pfs_ = Teuchos::rcp(new LINALG::SparseMatrix(
                        *(PoroField()->FluidField()->Discretization()->DofRowMap(0)),
                        //*(FluidField()->DofRowMap()),
                        81, true, true));
}

/*----------------------------------------------------------------------*
 | setup RHS (like fsimon)                                 vuong 01/12  |
 *----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Mono::SetupRHS(bool firstcall)
{
  // create full monolithic rhs vector
  if(rhs_==Teuchos::null)
    rhs_ = Teuchos::rcp(new Epetra_Vector(*DofRowMap(), true));

  PoroField()->SetupRHS(firstcall);

  // fill the Poroelasticity rhs vector rhs_ with the single field rhss
  SetupVector(*rhs_, PoroField()->RHS(), ScatraField().Residual());
}

/*----------------------------------------------------------------------*
 | setup vector of the structure and fluid field            vuong 01/12|
 *----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Mono::SetupVector(Epetra_Vector &f,
                                        Teuchos::RCP<const Epetra_Vector> sv,
                                        Teuchos::RCP<const Epetra_Vector> fv)
{
  // extract dofs of the two fields
  // and put the poro/scatra field vector into the global vector f
  // noticing the block number

  Extractor().InsertVector(*sv, 0, f);
  Extractor().InsertVector(*fv, 1, f);
}

/*----------------------------------------------------------------------*
 | setup system matrix of poroelasticity                   vuong 01/12  |
 *----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Mono::SetupSystemMatrix()
{
  // set loma block matrix to zero
  systemmatrix_->Zero();

  //----------------------------------------------------------------------
  // 1st diagonal block (upper left): poro weighting - poro solution
  //----------------------------------------------------------------------
  // get matrix block
  Teuchos::RCP<LINALG::SparseMatrix> mat_pp = PoroField()->SystemMatrix();

  // uncomplete matrix block (appears to be required in certain cases)
  mat_pp->UnComplete();

  // assign matrix block
  systemmatrix_->Assign(0,0,View,*mat_pp);

  //----------------------------------------------------------------------
  // 2nd diagonal block (lower right): scatra weighting - scatra solution
  //----------------------------------------------------------------------
  // get matrix block
  Teuchos::RCP<LINALG::SparseMatrix> mat_ss = ScatraField().SystemMatrix();

  // uncomplete matrix block (appears to be required in certain cases)
  mat_ss->UnComplete();

  // assign matrix block
  systemmatrix_->Assign(1,1,View,*mat_ss);

  // complete scatra block matrix
  systemmatrix_->Complete();

  //----------------------------------------------------------------------
  // 1st off-diagonal block (upper right): poro weighting - scatra solution
  //----------------------------------------------------------------------

  // evaluate off-diagonal matrix block in fluid
  EvaluateODBlockMatPoro();

  //k_ps_->Complete(mat_pp->DomainMap(),mat_ss->RangeMap());
  k_ps_->Complete();

  Teuchos::RCP<LINALG::SparseMatrix> k_ps_sparse = k_ps_->Merge();

  // uncomplete matrix block (appears to be required in certain cases)
  k_ps_sparse->UnComplete();

  // assign matrix block
  systemmatrix_->Assign(0,1,View,*(k_ps_sparse));

  //----------------------------------------------------------------------
  // 2nd off-diagonal block (lower left): scatra weighting - poro solution
  //----------------------------------------------------------------------

  // evaluate off-diagonal matrix block in scatra
  // (for present fixed-point-like iteration: no entries)
  EvaluateODBlockMatScatra();

  // uncomplete matrix block (appears to be required in certain cases)
 // k_sp_->UnComplete();

  // assign matrix block
 // systemmatrix_->Assign(1,0,View,*(k_sp_->Merge()));

  // complete block matrix
  systemmatrix_->Complete();
}

/*----------------------------------------------------------------------*
 | Solve linear Poroelasticity system                      vuong 01/12   |
 *----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Mono::LinearSolve()
{
  // Solve for inc_ = [disi_,tempi_]
  // Solve K_Teffdyn . IncX = -R  ===>  IncX_{n+1} with X=[d,T]
  // \f$x_{i+1} = x_i + \Delta x_i\f$
  if (solveradapttol_ and (iter_ > 1))
  {
    double worst = normrhs_;
    double wanted = tolfres_;
    solver_->AdaptTolerance(wanted, worst, solveradaptolbetter_);
  }
  // apply Dirichlet BCs to system of equations
  iterinc_->PutScalar(0.0);  // Useful? depends on solver and more

  if(directsolve_)
  {
    // merge blockmatrix to SparseMatrix and solve
    Teuchos::RCP<LINALG::SparseMatrix> sparse = systemmatrix_->Merge();

    LINALG::ApplyDirichlettoSystem(
        sparse,
        iterinc_,
        rhs_,
        Teuchos::null,
        zeros_,
        *CombinedDBCMap()
        );
    //  if ( Comm().MyPID()==0 ) { cout << " DBC applied to system" << endl; }

    // standard solver call
    solver_->Solve( sparse->EpetraOperator(),
                    iterinc_, rhs_,
                    true,
                    iter_ == 1
                    );
    //  if ( Comm().MyPID()==0 ) { cout << " Solved" << endl; }
  }
  else // use bgs2x2_operator
  {
    // in case of inclined boundary conditions
    // rotate systemmatrix_ using GetLocSysTrafo()!=Teuchos::null
    LINALG::ApplyDirichlettoSystem(
      systemmatrix_,
      iterinc_,
      rhs_,
      Teuchos::null,
      zeros_,
      *CombinedDBCMap()
      );

    solver_->Solve(
               systemmatrix_->EpetraOperator(),
               iterinc_,
               rhs_,
               true,
               iter_==1
               );
  }
}

/*----------------------------------------------------------------------*
 | setup solver for monolithic system                    vuong 01/12     |
 *----------------------------------------------------------------------*/
bool POROELAST::PORO_SCATRA_Mono::SetupSolver()
{
  //  solver
  // create a linear solver
  // get dynamic section of poroelasticity
  const Teuchos::ParameterList& poroscatradyn = DRT::Problem::Instance()->PoroScatraControlParams();
  // get the solver number used for linear poroelasticity solver
  const int linsolvernumber = poroscatradyn.get<int>("LINEAR_SOLVER");
  // check if the poroelasticity solver has a valid solver number
  if (linsolvernumber == (-1))
    dserror("no linear solver defined for scalar transport in porous media. Please set LINEAR_SOLVER in POROSCATRA CONTROL to a valid number!");
  const Teuchos::ParameterList& solverparams =
    DRT::Problem::Instance()->SolverParams(linsolvernumber);
  const int solvertype = DRT::INPUT::IntegralValue<INPAR::SOLVER::SolverType>(
    solverparams, "SOLVER");

  directsolve_ = (solvertype == INPAR::SOLVER::umfpack);

  if (directsolve_)
  {
    solver_ = Teuchos::rcp(new LINALG::Solver( solverparams,
                                      Comm(),
                                      DRT::Problem::Instance()->ErrorFile()->Handle())
                 );
  }
  else
    // create a linear solver
    //CreateLinearSolver();
    dserror("no implicit solver supported yet!");

  // Get the parameters for the Newton iteration
  itermax_ = poroscatradyn.get<int> ("ITEMAX");
  itermin_ = poroscatradyn.get<int> ("ITEMIN");
  normtypeinc_ = DRT::INPUT::IntegralValue<INPAR::PORO_SCATRA::ConvNorm>(
      poroscatradyn, "NORM_INC");
  normtypefres_ = DRT::INPUT::IntegralValue<INPAR::PORO_SCATRA::ConvNorm>(
      poroscatradyn, "NORM_RESF");
  combincfres_ = DRT::INPUT::IntegralValue<INPAR::PORO_SCATRA::BinaryOp>(
      poroscatradyn, "NORMCOMBI_RESFINC");

  tolinc_ = poroscatradyn.get<double> ("INCTOL");
  tolfres_ = poroscatradyn.get<double> ("RESTOL");

  return true;
}

/*----------------------------------------------------------------------*
 | Setup Newton-Raphson iteration            vuong 01/12   |
 *----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Mono::SetupNewton()
{

  // initialize variables needed in newton loop
  iter_ = 1;
  normrhs_ = 0.0;
  norminc_ = 0.0;

  // incremental solution vector with length of all dofs
  iterinc_ = LINALG::CreateVector(*DofRowMap(), true);
  iterinc_->PutScalar(0.0);

  // a zero vector of full length
  zeros_ = LINALG::CreateVector(*DofRowMap(), true);
  zeros_->PutScalar(0.0);

  return;
}

/*----------------------------------------------------------------------*
 | check convergence of Newton iteration (public)         vuong 01/12   |
 *----------------------------------------------------------------------*/
bool POROELAST::PORO_SCATRA_Mono::Converged()
{
  // check for single norms
  bool convinc = false;
  bool convfres = false;

  // residual increments
  switch (normtypeinc_)
  {
    case INPAR::PORO_SCATRA::convnorm_abs:
      convinc = norminc_ < tolinc_;
      break;
    default:
      dserror("Cannot check for convergence of residual values!");
      break;
  }

  // residual forces
  switch (normtypefres_)
  {
    case INPAR::PORO_SCATRA::convnorm_abs:
    convfres = normrhs_ < tolfres_;
    break;
    default:
    dserror("Cannot check for convergence of residual forces!");
    break;
  }

  // combine increments and forces
  bool conv = false;
  if (combincfres_==INPAR::PORO_SCATRA::bop_and)
    conv = convinc and convfres;
  else if (combincfres_==INPAR::PORO_SCATRA::bop_or)
    conv = convinc or convfres;
  else
    dserror("Something went terribly wrong with binary operator!");

  // return things
  return conv;
}  // Converged()

/*----------------------------------------------------------------------*
 | print Newton-Raphson iteration to screen and error file              |
 | originally by lw 12/07, tk 01/08                                     |
 *----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Mono::PrintNewtonIter()
{
  // print to standard out
  // replace myrank_ here general by Comm().MyPID()
  if ((Comm().MyPID() == 0) and printscreen_ and (Step()%printscreen_==0) and printiter_)
  {
    if (iter_ == 1)
      PrintNewtonIterHeader(stdout);
    PrintNewtonIterText(stdout);
  }

  // print to error file
  if (printerrfile_ and printiter_)
  {
    if (iter_ == 1)
      PrintNewtonIterHeader(errfile_);
    PrintNewtonIterText(errfile_);
  }

  return;
} // PrintNewtonIter()


/*----------------------------------------------------------------------*
 | print Newton-Raphson iteration to screen and error file              |
 | originally by lw 12/07, tk 01/08                                     |
 *----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Mono::PrintNewtonIterHeader(FILE* ofile)
{
  // open outstringstream
  std::ostringstream oss;

  // enter converged state etc
  oss << std::setw(6) << "numiter";

  // different style due relative or absolute error checking
  // displacement
  switch (normtypefres_)
  {
    case INPAR::PORO_SCATRA::convnorm_abs:
      oss << std::setw(18) << "abs-res";
      break;
    default:
      dserror("You should not turn up here.");
      break;
  }

  switch ( normtypeinc_ )
  {
    case INPAR::PORO_SCATRA::convnorm_abs :
    oss <<std::setw(18)<< "abs-inc";
    break;
    default:
    dserror("You should not turn up here.");
    break;
  }

  switch ( normtypefres_ )
  {
    case INPAR::PORO_SCATRA::convnorm_abs :
      oss <<std::setw(18)<< "abs-poro-res";
    //oss <<std::setw(18)<< "abs-s-res";
    //oss <<std::setw(18)<< "abs-fvel-res";
    //oss <<std::setw(18)<< "abs-fpres-res";
      oss <<std::setw(18)<< "abs-scalar-res";
    break;
    default:
    dserror("You should not turn up here.");
    break;
  }

  switch ( normtypeinc_ )
  {
    case INPAR::PORO_SCATRA::convnorm_abs :
      oss <<std::setw(18)<< "abs-poro-inc";
    //oss <<std::setw(18)<< "abs-s-inc";
    //oss <<std::setw(18)<< "abs-fvel-inc";
    //oss <<std::setw(18)<< "abs-fpres-inc";
      oss <<std::setw(18)<< "abs-scalar-inc";
    break;
    default:
    dserror("You should not turn up here.");
    break;
  }

  // add solution time
  oss << std::setw(14)<< "wct";

  // finish oss
  oss << std::ends;

  // print to screen (could be done differently...)
  if (ofile==NULL)
  dserror("no ofile available");
  fprintf(ofile, "%s\n", oss.str().c_str());

  // print it, now
  fflush(ofile);

  // nice to have met you
  return;
}  // PrintNewtonIterHeader()


/*----------------------------------------------------------------------*
 | print Newton-Raphson iteration to screen                       |
 | originally by lw 12/07, tk 01/08                                     |
 *----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Mono::PrintNewtonIterText(FILE* ofile)
{

  // open outstringstream
  std::ostringstream oss;

  // enter converged state etc
  oss << std::setw(7) << iter_;

  // different style due relative or absolute error checking
  // displacement
  switch (normtypefres_)
  {
    case INPAR::POROELAST::convnorm_abs:
      oss << std::setw(18) << std::setprecision(5) << std::scientific
          << normrhs_;
      break;
    default:
      dserror("You should not turn up here.");
      break;
  }

  switch ( normtypeinc_ )
  {
    case INPAR::POROELAST::convnorm_abs :
      oss << std::setw(18) << std::setprecision(5) << std::scientific << norminc_;
    break;
    default:
    dserror("You should not turn up here.");
    break;
  }

  switch ( normtypefres_ )
  {
    case INPAR::POROELAST::convnorm_abs :
      oss << std::setw(18) << std::setprecision(5) << std::scientific << normrhsporo_;
      oss << std::setw(18) << std::setprecision(5) << std::scientific << normrhsscatra_;
    break;
    default:
    dserror("You should not turn up here.");
    break;
  }

  switch ( normtypeinc_ )
  {
    case INPAR::POROELAST::convnorm_abs :
      oss << std::setw(18) << std::setprecision(5) << std::scientific << normincporo_;
      oss << std::setw(18) << std::setprecision(5) << std::scientific << normincscatra_;
    break;
    default:
    dserror("You should not turn up here.");
    break;
  }

  // add solution time
  oss << std::setw(14) << std::setprecision(2) << std::scientific << timer_.ElapsedTime();

  // finish oss
  oss << std::ends;

  // print to screen (could be done differently...)
  if (ofile==NULL)
  dserror("no ofile available");
  fprintf(ofile, "%s\n", oss.str().c_str());

  // print it, now
  fflush(ofile);

  // nice to have met you
  return;

}  // PrintNewtonIterText

/*----------------------------------------------------------------------*
 | print statistics of converged NRI                                     |
 | orignially by bborn 08/09                                            |
 *----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Mono::PrintNewtonConv()
{
  // somebody did the door
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Mono::BuildCovergenceNorms()
{
  // build residual force norm
  // for now use for simplicity only L2/Euclidian norm
  rhs_->Norm2(&normrhs_);
  Teuchos::RCP<const Epetra_Vector> rhs_poro;
  Teuchos::RCP<const Epetra_Vector> rhs_scatra;


  // process structure unknowns of the first field
  rhs_poro = Extractor().ExtractVector(rhs_, 0);
  // process fluid unknowns of the second field
  rhs_scatra = Extractor().ExtractVector(rhs_, 1);

  rhs_poro->Norm2(&normrhsporo_);
  rhs_scatra->Norm2(&normrhsscatra_);

  // build residual increment norm
  iterinc_->Norm2(&norminc_);

  // displacement and fluid velocity & pressure incremental vector
  Teuchos::RCP<const Epetra_Vector> interincporo;
  Teuchos::RCP<const Epetra_Vector> interincscatra;

  // process structure unknowns of the first field
  interincporo = Extractor().ExtractVector(iterinc_, 0);
  // process fluid unknowns of the second field
  interincscatra = Extractor().ExtractVector(iterinc_, 1);

  interincporo->Norm2(&normincporo_);
  interincscatra->Norm2(&normincscatra_);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> POROELAST::PORO_SCATRA_Mono::DofRowMap() const
{
  return blockrowdofmap_.FullMap();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> POROELAST::PORO_SCATRA_Mono::CombinedDBCMap() const
{
  return dbcmaps_->FullMap();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> POROELAST::PORO_SCATRA_Mono::SystemMatrix() const
{
  return systemmatrix_->Merge();
}

/*----------------------------------------------------------------------*
 | put the single maps to one full                                      |
 | Poroelasticity map together                              vuong 01/12 |
 *----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Mono::SetDofRowMaps(const std::vector<Teuchos::RCP<
    const Epetra_Map> >& maps)
{
  Teuchos::RCP<Epetra_Map> fullmap = LINALG::MultiMapExtractor::MergeMaps(maps);

  // full monolithic-blockmap
  blockrowdofmap_.Setup(*fullmap, maps);
}

/*----------------------------------------------------------------------*
 |  map containing the dofs with Dirichlet BC               vuong 01/12 |
 *----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Mono::BuildCombinedDBCMap()
{
  const Teuchos::RCP<const Epetra_Map> porocondmap =
      PoroField()->CombinedDBCMap();
  const Teuchos::RCP<const Epetra_Map> scatracondmap =
      ScatraField().DirichMaps()->CondMap();
  Teuchos::RCP<const Epetra_Map> dbcmap = LINALG::MergeMap(porocondmap, scatracondmap, false);

  // Finally, create the global FSI Dirichlet map extractor
  dbcmaps_ = Teuchos::rcp(new LINALG::MapExtractor(*DofRowMap(),dbcmap,true));
  if (dbcmaps_ == Teuchos::null) { dserror("Creation of Dirichlet map extractor failed."); }

} // BuildCombinedDBCMap()

/*----------------------------------------------------------------------*
 |  Evaluate off diagonal matrix in poro row                     |
 *----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Mono::EvaluateODBlockMatPoro()
{
  k_pfs_->Zero();

  // create the parameters for the discretization
  Teuchos::ParameterList fparams;
  // action for elements
  fparams.set<int>("action", FLD::calc_poroscatra_mono_odblock);
  // physical type
  fparams.set<int>("physical type", PoroField()->FluidField()->PhysicalType());
  // other parameters that might be needed by the elements
  fparams.set("delta time", Dt());
  fparams.set("total time", Time());

  const Teuchos::RCP<DRT::Discretization>& porofluiddis = PoroField()->FluidField()->Discretization();
  porofluiddis->ClearState();

  // set general vector values needed by elements
  porofluiddis->SetState(0,"dispnp",PoroField()->FluidField()->Dispnp());
  porofluiddis->SetState(0,"gridv",PoroField()->FluidField()->GridVel());
  porofluiddis->SetState(0,"veln",PoroField()->FluidField()->Veln());
  porofluiddis->SetState(0,"accnp",PoroField()->FluidField()->Accnp());

  PoroField()->FluidField()->Discretization()->SetState(0,"scaaf",PoroField()->FluidField()->Scaaf());

  // set scheme-specific element parameters and vector values
  //TODO
  //if (is_genalpha_)
  //    discret_->SetState("velaf",velaf_);
  //else

  porofluiddis->SetState(0,"velaf",PoroField()->FluidField()->Velnp());
  porofluiddis->SetState(0,"velnp",PoroField()->FluidField()->Velnp());

  // build specific assemble strategy for the fluid-mechanical system matrix
  // from the point of view of FluidField:
  // fluiddofset = 0, structdofset = 1
  DRT::AssembleStrategy fluidstrategy(
      0,              // fluiddofset for row
      2,              // structdofset for column
      k_pfs_,   // fluid-mechanical matrix
      Teuchos::null,  // no other matrix or vectors
      Teuchos::null ,
      Teuchos::null,
      Teuchos::null
  );

  // evaluate the fluid-mechanical system matrix on the fluid element
  //PoroField()->FluidField()->Discretization()->EvaluateCondition( fparams, fluidstrategy,"PoroCoupling" );
  PoroField()->FluidField()->Discretization()->Evaluate( fparams, fluidstrategy );

  porofluiddis->ClearState();

  //************************************************************************************
  //************************************************************************************

  k_pss_->Zero();


  // create the parameters for the discretization
  Teuchos::ParameterList sparams;

  const Teuchos::ParameterList& sdynparams
  = DRT::Problem::Instance()->StructuralDynamicParams();

  enum INPAR::STR::DynamicType strmethodname = DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(sdynparams,"DYNAMICTYP");

  switch (strmethodname)
  {
   case  INPAR::STR::dyna_statics :
   {
     sparams.set("theta", 1.0);
   // continue
   break;
   }
  case INPAR::STR::dyna_onesteptheta:
  {
    double theta = sdynparams.sublist("ONESTEPTHETA").get<double> ("THETA");
    sparams.set("theta", theta);
    break;
  }
    // TODO: time factor for genalpha
    /*
     case INPAR::STR::dyna_genalpha :
     {
     double alphaf_ = sdynparams.sublist("GENALPHA").get<double>("ALPHA_F");
     // K_Teffdyn(T_n+1) = (1-alphaf_) . kst
     // Lin(dT_n+1-alphaf_/ dT_n+1) = (1-alphaf_)
     k_st->Scale(1.0 - alphaf_);
     }*/
  default:
  {
    dserror("Don't know what to do... only one-step theta time integration implemented");
    break;
  }
  } // end of switch(strmethodname_)

  const std::string action = "calc_struct_poroscatracoupling";
  sparams.set("action", action);
  // other parameters that might be needed by the elements
  sparams.set("delta time", Dt());
  sparams.set("total time", Time());

  PoroField()->StructureField()->Discretization()->ClearState();
  PoroField()->StructureField()->Discretization()->SetState(0,"displacement",PoroField()->StructureField()->Dispnp());
  PoroField()->StructureField()->Discretization()->SetState(0,"velocity",PoroField()->StructureField()->ExtractVelnp());

  PoroField()->StructureField()->SetCouplingState();

  // build specific assemble strategy for mechanical-fluid system matrix
  // from the point of view of StructureField:
  // structdofset = 0, fluiddofset = 1
  DRT::AssembleStrategy structuralstrategy(
      0,               // structdofset for row
      2,               // fluiddofset for column
      k_pss_,            // mechanical-fluid coupling matrix
      Teuchos::null ,
      Teuchos::null ,
      Teuchos::null,
      Teuchos::null
  );

  // evaluate the mechancial-fluid system matrix on the structural element
  PoroField()->StructureField()->Discretization()->EvaluateCondition( sparams, structuralstrategy,"PoroCoupling" );
  //StructureField()->Discretization()->Evaluate( sparams, structuralstrategy);
  PoroField()->StructureField()->Discretization()->ClearState();


  //************************************************************************************
  //************************************************************************************

  //Teuchos::RCP<LINALG::SparseMatrix> k_pss = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(k_pss_);
  //Teuchos::RCP<LINALG::SparseMatrix> k_pfs = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(k_pfs_);

  // assign matrix blocks
  k_ps_->Assign(0,0,View,*k_pss_);
  k_ps_->Assign(1,0,View,*k_pfs_);

  return;
}

/*----------------------------------------------------------------------*
 |  Evaluate off diagonal matrix in scatra row                    |
 *----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Mono::EvaluateODBlockMatScatra()
{

}

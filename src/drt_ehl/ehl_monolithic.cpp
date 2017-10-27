/*--------------------------------------------------------------------------*/
/*!
\file ehl_monolithic.cpp

\brief basis of all monolithic EHL algorithms that perform a coupling between
       the structure field equation and lubrication field equations

\level 3
\maintainer Alexander Seitz
*/
/*--------------------------------------------------------------------------*/



/*----------------------------------------------------------------------*
 | headers                                                  wirtz 01/16 |
 *----------------------------------------------------------------------*/
#include "ehl_monolithic.H"

#include "../drt_lubrication/lubrication_timint_implicit.H"

#include "../drt_lubrication_ele/lubrication_ele_action.H"

#include <Teuchos_TimeMonitor.hpp>

#include "../drt_adapter/adapter_lubrication.H"
#include "../drt_adapter/ad_str_structure.H"
#include "../drt_lib/drt_assemblestrategy.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_locsys.H"

#include "../linalg/linalg_blocksparsematrix.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_utils.H"

#include "../drt_inpar/inpar_solver.H"

#include "../drt_io/io_control.H"


//! Note: The order of calling the two BaseAlgorithm-constructors is
//! important here! In here control file entries are written. And these entries
//! define the order in which the filters handle the Discretizations, which in
//! turn defines the dof number ordering of the Discretizations.


/*----------------------------------------------------------------------*
 | destructor (public)                                      wirtz 01/16 |
 *----------------------------------------------------------------------*/
EHL::Monolithic::~Monolithic()
{
}


/*----------------------------------------------------------------------*
 | monolithic                                               wirtz 01/16 |
 *----------------------------------------------------------------------*/
EHL::Monolithic::Monolithic(
  const Epetra_Comm& comm,
  const Teuchos::ParameterList& globaltimeparams,
  const Teuchos::ParameterList& lubricationparams,
  const Teuchos::ParameterList& structparams,
  const std::string struct_disname,
  const std::string lubrication_disname
  )
: Base(comm, globaltimeparams,lubricationparams,structparams,struct_disname,lubrication_disname),
//: Alglrithm(comm),
  solveradapttol_(DRT::INPUT::IntegralValue<int>(((DRT::Problem::Instance()->ElastoHydroDynamicParams()).sublist("MONOLITHIC")),"ADAPTCONV") == 1),
  solveradaptolbetter_(((DRT::Problem::Instance()->ElastoHydroDynamicParams()).sublist("MONOLITHIC")).get<double>("ADAPTCONV_BETTER")),
  printiter_(true),  // ADD INPUT PARAMETER
  printerrfile_(false),  // ADD INPUT PARAMETER FOR 'true'
  errfile_(NULL),
  zeros_(Teuchos::null),
  strmethodname_(DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(structparams,"DYNAMICTYP")),
  ehldyn_(DRT::Problem::Instance()->ElastoHydroDynamicParams()),
  ehldynmono_( (DRT::Problem::Instance()->ElastoHydroDynamicParams()).sublist("MONOLITHIC") ),
  blockrowdofmap_(Teuchos::null),
  systemmatrix_(Teuchos::null),
  k_sl_(Teuchos::null),
  k_ls_(Teuchos::null),
  merge_ehl_blockmatrix_(DRT::INPUT::IntegralValue<bool>(ehldynmono_,"MERGE_EHL_BLOCK_MATRIX")),
  soltech_(DRT::INPUT::IntegralValue<INPAR::EHL::NlnSolTech>(ehldynmono_,"NLNSOL")),
  iternorm_(DRT::INPUT::IntegralValue<INPAR::EHL::VectorNorm>(ehldynmono_,"ITERNORM")),
  iter_(0),
  sdyn_(structparams),
  timernewton_(comm),
  ptcdt_(ehldynmono_.get<double>("PTCDT")),
  dti_(1.0/ptcdt_),
  vel_(Teuchos::null)
{

  errfile_ = DRT::Problem::Instance()->ErrorFile()->Handle();
  if (errfile_)
    printerrfile_ = true;

  blockrowdofmap_ = Teuchos::rcp(new LINALG::MultiMapExtractor);

  // initialise internal varible with new velocities V_{n+1} at t_{n+1}
  vel_ = LINALG::CreateVector(*(StructureField()->DofRowMap(0)), true);

  // --------------------------------- EHL solver: create a linear solver

  // get iterative solver
  if (merge_ehl_blockmatrix_ == false)
    CreateLinearSolver();
  // get direct solver, e.g. UMFPACK
  else  // (merge_ehl_blockmatrix_ == true)
  {
    if (Comm().MyPID() == 0)
      std::cout << "Merged EHL block matrix is used!\n" << std::endl;

    Teuchos::RCP<Teuchos::ParameterList> solverparams = Teuchos::rcp(new Teuchos::ParameterList);
    solverparams->set("solver","umfpack");

    solver_ = Teuchos::rcp(new LINALG::Solver(
                                 solverparams,
                                 Comm(),
                                 DRT::Problem::Instance()->ErrorFile()->Handle()
                                 )
                );
  }  // end BlockMatrixMerge

  // StructureField: check whether we have locsys BCs, i.e. inclined structural
  //  Dirichlet BC
  {
    std::vector<DRT::Condition*> locsysconditions(0);
    (StructureField()->Discretization())->GetCondition("Locsys", locsysconditions);

    // if there are inclined structural Dirichlet BC, get the structural LocSysManager
    if (locsysconditions.size())
    {
      locsysman_ = StructureField()->LocsysManager();
    }
  }

}  // Monolithic()


/*----------------------------------------------------------------------*
 | read restart information for given time step (public)    wirtz 01/16 |
 *----------------------------------------------------------------------*/
void EHL::Monolithic::ReadRestart(int step)
{
  lubrication_->LubricationField()->ReadRestart(step);
  StructureField()->ReadRestart(step);

  // pass the current coupling variables to the respective field
  ApplyStructCouplingState(StructureField()->Dispnp(),StructureField()->Velnp());
  ApplyLubricationCouplingState(lubrication_->LubricationField()->Prenp());

  // second ReadRestart needed due to the coupling variables
  lubrication_->LubricationField()->ReadRestart(step);
  StructureField()->ReadRestart(step);

  SetTimeStep(StructureField()->TimeOld(),step);

  return;
}  // ReadRestart()


/*----------------------------------------------------------------------*
 | prepare time step (public)                               wirtz 01/16 |
 *----------------------------------------------------------------------*/
void EHL::Monolithic::PrepareTimeStep()
{
  // counter and print header
  // increment time and step counter
  IncrementTimeAndStep();
  PrintHeader();

  // pass the current coupling variables to the respective fields
  ApplyStructCouplingState(StructureField()->Dispnp(),StructureField()->Velnp());
  ApplyLubricationCouplingState(lubrication_->LubricationField()->Prenp());

  // call the predictor
  StructureField()->PrepareTimeStep();
  lubrication_->LubricationField()->PrepareTimeStep();

}  // PrepareTimeStep()


/*----------------------------------------------------------------------*
 | create linear solver                                     wirtz 01/16 |
 *----------------------------------------------------------------------*/
void EHL::Monolithic::CreateLinearSolver()
{
  // get the solver number used for linear EHL solver
  const int linsolvernumber = ehldynmono_.get<int>("LINEAR_SOLVER");
  // check if the EHL solver has a valid solver number
  if (linsolvernumber == (-1))
    dserror("no linear solver defined for monolithic EHL. Please set LINEAR_SOLVER in ELASTO HYDRO DYNAMIC to a valid number!");

  // get parameter list of structural dynamics
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
  // use solver blocks for structure
  // get the solver number used for structural solver
  const int slinsolvernumber = sdyn.get<int>("LINEAR_SOLVER");
  // check if the structural solver has a valid solver number
  if (slinsolvernumber == (-1))
    dserror("no linear solver defined for structural field. Please set LINEAR_SOLVER in STRUCTURAL DYNAMIC to a valid number!");

  // get parameter list of lubrication dynamics
  const Teuchos::ParameterList& ldyn = DRT::Problem::Instance()->LubricationDynamicParams();
  // use solver blocks for pressure (lubrication field)
  // get the solver number used for lubrication solver
  const int tlinsolvernumber = ldyn.get<int>("LINEAR_SOLVER");
  // check if the EHL solver has a valid solver number
  if (tlinsolvernumber == (-1))
    dserror("no linear solver defined for lubrication field. Please set LINEAR_SOLVER in THERMAL DYNAMIC to a valid number!");

  // get solver parameter list of linear EHL solver
  const Teuchos::ParameterList& ehlsolverparams
    = DRT::Problem::Instance()->SolverParams(linsolvernumber);

  const int solvertype
    = DRT::INPUT::IntegralValue<INPAR::SOLVER::SolverType>(
        ehlsolverparams,
        "SOLVER"
        );

  if ( (solvertype != INPAR::SOLVER::aztec_msr) and (solvertype != INPAR::SOLVER::belos) )
  {
    std::cout << "!!!!!!!!!!!!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!!!!!!!!!" << std::endl;
    std::cout << " Note: the BGS2x2 preconditioner now "                  << std::endl;
    std::cout << " uses the structural solver and lubrication solver blocks"  << std::endl;
    std::cout << " for building the internal inverses"                    << std::endl;
    std::cout << " Remove the old BGS PRECONDITIONER BLOCK entries "      << std::endl;
    std::cout << " in the dat files!"                                     << std::endl;
    std::cout << "!!!!!!!!!!!!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!!!!!!!!!" << std::endl;
    dserror("aztec solver expected");
  }
  const int azprectype
    = DRT::INPUT::IntegralValue<INPAR::SOLVER::AzPrecType>(
        ehlsolverparams,
        "AZPREC"
        );

  // plausibility check
  switch (azprectype)
  {
  case INPAR::SOLVER::azprec_BGS2x2 :
    break;
  case INPAR::SOLVER::azprec_BGSnxn :
  case INPAR::SOLVER::azprec_TekoSIMPLE :
  {
    dserror("Teko preconditioners only available with HAVE_TEKO flag for TRILINOS_DEV (>Q1/2011)");
    break;
  }
  case INPAR::SOLVER::azprec_MueLuAMG_sym :
  case INPAR::SOLVER::azprec_AMGnxn :
  case INPAR::SOLVER::azprec_CheapSIMPLE :
  {
    // no plausibility checks here
    // if you forget to declare an xml file you will get an error message anyway
  }
  break;
  default:
    dserror("Block Gauss-Seidel BGS2x2 preconditioner expected. Alternatively you can define your own AMG block preconditioner (using an xml file). This is experimental.");
    break;
  }


  // prepare linear solvers and preconditioners
  switch (azprectype)
  {
  case INPAR::SOLVER::azprec_BGS2x2 :
  case INPAR::SOLVER::azprec_BGSnxn :
  case INPAR::SOLVER::azprec_TekoSIMPLE :
  case INPAR::SOLVER::azprec_AMGnxn :
  case INPAR::SOLVER::azprec_CheapSIMPLE :
  {
    // This should be the default case (well-tested and used)
    solver_ = Teuchos::rcp(new LINALG::Solver(
                                 ehlsolverparams,
                                 // ggfs. explizit Comm von STR wie lungscatra
                                 Comm(),
                                 DRT::Problem::Instance()->ErrorFile()->Handle()
                                 )
                );

    // use solver blocks for structure and pressure (lubrication field)
    const Teuchos::ParameterList& ssolverparams = DRT::Problem::Instance()->SolverParams(slinsolvernumber);
    const Teuchos::ParameterList& tsolverparams = DRT::Problem::Instance()->SolverParams(tlinsolvernumber);

    solver_->PutSolverParamsToSubParams("Inverse1", ssolverparams);
    solver_->PutSolverParamsToSubParams("Inverse2", tsolverparams);

    // prescribe rigid body modes
    StructureField()->Discretization()->ComputeNullSpaceIfNecessary(
      solver_->Params().sublist("Inverse1")
      );
    lubrication_->LubricationField()->Discretization()->ComputeNullSpaceIfNecessary(
      solver_->Params().sublist("Inverse2")
      );


    if(azprectype==INPAR::SOLVER::azprec_CheapSIMPLE)
    {
    // Tell to the LINALG::SOLVER::SimplePreconditioner that we use the general implementation
      solver_->Params().set<bool>("GENERAL",true);
    }

    break;
  }
  case INPAR::SOLVER::azprec_MueLuAMG_sym:
  {
    solver_ = Teuchos::rcp(new LINALG::Solver(
                                 ehlsolverparams,
                                 // ggfs. explizit Comm von STR wie lungscatra
                                 Comm(),
                                 DRT::Problem::Instance()->ErrorFile()->Handle()
                                 )
                );

    // use solver blocks for structure and pressure (lubrication field)
    const Teuchos::ParameterList& ssolverparams = DRT::Problem::Instance()->SolverParams(slinsolvernumber);
    const Teuchos::ParameterList& tsolverparams = DRT::Problem::Instance()->SolverParams(tlinsolvernumber);

    // This is not very elegant:
    // first read in solver parameters. These have to contain ML parameters such that...
    solver_->PutSolverParamsToSubParams("Inverse1", ssolverparams);
    solver_->PutSolverParamsToSubParams("Inverse2", tsolverparams);

    // ... BACI calculates the null space vectors. These are then stored in the sublists
    //     Inverse1 and Inverse2 from where they...
    StructureField()->Discretization()->ComputeNullSpaceIfNecessary(
      solver_->Params().sublist("Inverse1")
      );
    lubrication_->LubricationField()->Discretization()->ComputeNullSpaceIfNecessary(
      solver_->Params().sublist("Inverse2")
      );

    // ... are copied from here to ...
    const Teuchos::ParameterList& inv1source = solver_->Params().sublist("Inverse1").sublist("ML Parameters");
    const Teuchos::ParameterList& inv2source = solver_->Params().sublist("Inverse2").sublist("ML Parameters");

    // ... here. The "MueLu Parameters" sublists "Inverse1" and "Inverse2" only contain the basic
    //     information about the corresponding null space vectors, which are actually copied ...
    Teuchos::ParameterList& inv1 = solver_->Params().sublist("MueLu Parameters").sublist("Inverse1");
    Teuchos::ParameterList& inv2 = solver_->Params().sublist("MueLu Parameters").sublist("Inverse2");

    // ... here.
    inv1.set<int>("PDE equations", inv1source.get<int>("PDE equations"));
    inv2.set<int>("PDE equations", inv2source.get<int>("PDE equations"));
    inv1.set<int>("null space: dimension", inv1source.get<int>("null space: dimension"));
    inv2.set<int>("null space: dimension", inv2source.get<int>("null space: dimension"));
    inv1.set<double*>("null space: vectors", inv1source.get<double*>("null space: vectors"));
    inv2.set<double*>("null space: vectors", inv2source.get<double*>("null space: vectors"));
    inv1.set<Teuchos::RCP<std::vector<double> > >("nullspace", inv1source.get<Teuchos::RCP<std::vector<double> > >("nullspace"));
    inv2.set<Teuchos::RCP<std::vector<double> > >("nullspace", inv2source.get<Teuchos::RCP<std::vector<double> > >("nullspace"));

    solver_->Params().sublist("MueLu Parameters").set("EHL",true);
    break;
  }
  default:
    dserror("Block Gauss-Seidel BGS2x2 preconditioner expected");
    break;
  }

}  // CreateLinearSolver()


/*----------------------------------------------------------------------*
 | non-linear solve, i.e. (multiple) corrector (public)     wirtz 01/16 |
 *----------------------------------------------------------------------*/
void EHL::Monolithic::Solve()
{
  // choose solution technique according to input file
  switch (soltech_)
  {
  // Newton-Raphson iteration
  case INPAR::EHL::soltech_newtonfull :
    NewtonFull();
    break;
  // Pseudo-transient continuation
  case INPAR::EHL::soltech_ptc :
    PTC();
    break;
  // catch problems
  default :
    dserror(
      "Solution technique \"%s\" is not implemented",
      INPAR::EHL::NlnSolTechString(soltech_).c_str()
      );
    break;
  }  // end switch (soltechnique_)

  return;
}  // Solve()


/*----------------------------------------------------------------------*
 | time loop of the monolithic system                       wirtz 01/16 |
 *----------------------------------------------------------------------*/
void EHL::Monolithic::Timeloop()
{
  // time loop
  while (NotFinished())
  {
    // counter and print header
    // predict solution of both field (call the adapter)
    PrepareTimeStep();

    // integrate time step
    Solve();

    // calculate stresses, strains, energies
    PrepareOutput();

    // update all single field solvers
    Update();

    // write output to screen and files
    Output();

  }  // NotFinished
}  // TimeLoop()


/*----------------------------------------------------------------------*
 | solution with full Newton-Raphson iteration              wirtz 01/16 |
 | in ehl_algorithm: NewtonFull()                                       |
 *----------------------------------------------------------------------*/
void EHL::Monolithic::NewtonFull()
{

  // we do a Newton-Raphson iteration here.
  // the specific time integration has set the following
  // --> On #rhs_ is the positive force residuum
  // --> On #systemmatrix_ is the effective dynamic EHL tangent matrix

  // initialise equilibrium loop
  iter_ = 1;

  // incremental solution vector with length of all EHL dofs
  iterinc_ = LINALG::CreateVector(*DofRowMap(), true);
  iterinc_->PutScalar(0.0);
  // a zero vector of full length
  zeros_ = LINALG::CreateVector(*DofRowMap(), true);
  zeros_->PutScalar(0.0);

  //------------------------------------------------------ iteration loop

  // equilibrium iteration loop (loop over k)
  while ( ( (not Converged()) and (iter_ <= itermax_) ) or (iter_ <= itermin_) )
  {
    // reset timer
    timernewton_.ResetStartTime();

    // compute residual forces #rhs_ and tangent #systemmatrix_
    // whose components are globally oriented
    // build linear system stiffness matrix and rhs/force residual for each
    // field, here e.g. for structure field: field want the iteration increment
    // 1.) Update(iterinc_),
    // 2.) EvaluateForceStiffResidual(),
    // 3.) PrepareSystemForNewtonSolve() --> if (locsysman_!=null) k_ss is rotated
    Evaluate(iterinc_);

    // create the linear system
    // \f$J(x_i) \Delta x_i = - R(x_i)\f$
    // create the systemmatrix
    SetupSystemMatrix();

    // check whether we have a sanely filled tangent matrix
    if (not systemmatrix_->Filled())
      dserror("Effective tangent matrix must be filled here");

    // create full monolithic rhs vector
    // make negative residual not necessary: rhs_ is already negative
    // (STR/LUBRICATION)-RHS is put negative in PrepareSystemForNewtonSolve()
    SetupRHS();

    // *********** time measurement ***********
    double dtcpu = timernewton_.WallTime();
    // *********** time measurement ***********
    // (Newton-ready) residual with blanked Dirichlet DOFs (see adapter_timint!)
    // is done in PrepareSystemForNewtonSolve() within Evaluate(iterinc_)
    LinearSolve();
    // *********** time measurement ***********
    dtsolve_ = timernewton_.WallTime() - dtcpu;
    // *********** time measurement ***********

    // reset solver tolerance
    solver_->ResetTolerance();

    // --------------------------------------------- build residual norms
    // include all stuff here related with convergence test
    normrhs_ = CalculateVectorNorm(iternorm_, rhs_);
    // vector of displacement and pressure residual
    Teuchos::RCP<Epetra_Vector> strrhs;
    Teuchos::RCP<Epetra_Vector> lubricationrhs;
    // extract field vectors
    ExtractFieldVectors(rhs_, strrhs, lubricationrhs);
    normstrrhs_ = CalculateVectorNorm(iternormstr_, strrhs);
    normlubricationrhs_ = CalculateVectorNorm(iternormlubrication_, lubricationrhs);

    // --------------------------------- build residual incremental norms
    // vector of displacement and pressure increments
    Teuchos::RCP<Epetra_Vector> sx;
    Teuchos::RCP<Epetra_Vector> lx;
    // extract field vectors
    ExtractFieldVectors(iterinc_,sx,lx);
    norminc_ = CalculateVectorNorm(iternorm_, iterinc_);
    normdisi_ = CalculateVectorNorm(iternormstr_, sx);
    normprei_ = CalculateVectorNorm(iternormlubrication_, lx);

    // in case of 'Mix'-convergence criterion: save the norm of the 1st
    // iteration in (norm . iter0_)
    if (iter_ == 1)
    {
      // save initial residual norms
      normrhsiter0_ = normrhs_;
      normstrrhsiter0_ = normstrrhs_;
      normlubricationrhsiter0_ = normlubricationrhs_;
      // save initial incremental norms
      norminciter0_ = norminc_;
      normdisiiter0_ = normdisi_;
      normpreiiter0_ = normprei_;

      // set the minimum of iter0_ and tolrhs_, because we want to prevent the
      // case of a zero characteristic initial norm
      if (normrhsiter0_ == 0.0) normrhsiter0_ = tolrhs_;
      if (normstrrhsiter0_ == 0.0) normstrrhsiter0_ = tolstrrhs_;
      if (normlubricationrhsiter0_ == 0.0) normlubricationrhsiter0_ = tollubricationrhs_;
      if (norminciter0_ == 0.0) norminciter0_ = tolinc_;
      if (normdisiiter0_ == 0.0) normdisiiter0_ = toldisi_;
      if (normpreiiter0_ == 0.0) normpreiiter0_ = tolprei_;
    }

    // print stuff
    PrintNewtonIter();

    // increment equilibrium loop index
    iter_ += 1;

  }  // end equilibrium loop

  // ----------------------------------------------------- iteration loop

  // correct iteration counter
  iter_ -= 1;

  // test whether max iterations was hit
  if ( (Converged()) and (Comm().MyPID() == 0) )
  {
    PrintNewtonConv();
  }
  else if (iter_ >= itermax_)
    dserror("Newton unconverged in %d iterations", iter_);

}  // NewtonFull()


/*----------------------------------------------------------------------*
 | solution with pseudo-transient continuation              wirtz 01/16 |
 *----------------------------------------------------------------------*/
void EHL::Monolithic::PTC()
{
  // do a PTC iteration here
  // implementation is based on the work of Gee, Kelley, Lehouq (2009):
  // "Pseudo-transient continuation for nonlinear transient elasticity"

  // the specific time integration has set the following
  // --> On #rhs_ is the positive force residuum
  // --> On #systemmatrix_ is the effective dynamic EHL tangent matrix

  // initialise equilibrium loop
  iter_ = 1;

  // incremental solution vector with length of all EHL DOFs
  iterinc_ = LINALG::CreateVector(*DofRowMap(), true);
  iterinc_->PutScalar(0.0);
  // a zero vector of full length
  zeros_ = LINALG::CreateVector(*DofRowMap(), true);
  zeros_->PutScalar(0.0);

  // create full monolithic rhs vector
  // make negative residual not necessary: rhs_ is already negative
  // (STR/LUBRICATION)-RHS is put negative in PrepareSystemForNewtonSolve()
  SetupRHS();

  // ----------------------------------------------- special stuff of PTC

  // compute the PTC parameters
  double ptcdt = ptcdt_;
  // norm of residual of old iteration step
  double nc = 0.0;
  // as recommended by Michael, apply PTC to the whole EHL system, i.e. use EHL
  // residual here
  // here the new convergence test stuff has to be included

  normrhs_ = CalculateVectorNorm(iternorm_, rhs_);
  rhs_->NormInf(&nc);
  // define the pseudo time step delta^{-1}
  double dti = 1/ptcdt;

  //------------------------------------------------------ iteration loop

  // equilibrium iteration loop (loop over k)
  while ( ( (not Converged()) and (iter_ <= itermax_) ) or (iter_ <= itermin_) )
  {
    // reset timer
    timernewton_.ResetStartTime();

    // compute residual forces #rhs_ and tangent #systemmatrix_
    // whose components are globally oriented
    // build linear EHL tangent matrix and rhs/force residual for each field,
    // here e.g. for structure field: STR field wants the iteration increment
    // 1.) Update(iterinc_),
    // 2.) EvaluateForceStiffResidual(),
    // 3.) PrepareSystemForNewtonSolve() --> if (locsysman_ != null) k_ss is rotated
    Evaluate(iterinc_);

    // ---------- modify diagonal blocks of systemmatrix according to PTC

    // modify structural diagonal block k_ss
    {
      Teuchos::RCP<Epetra_Vector> tmp_SS
        = LINALG::CreateVector(StructureField()->SystemMatrix()->RowMap(),false);
      tmp_SS->PutScalar(dti);
      Teuchos::RCP<Epetra_Vector> diag_SS
        = LINALG::CreateVector(StructureField()->SystemMatrix()->RowMap(),false);
      StructureField()->SystemMatrix()->ExtractDiagonalCopy(*diag_SS);

      diag_SS->Update(1.0, *tmp_SS, 1.0);

      StructureField()->SystemMatrix()->ReplaceDiagonalValues(*diag_SS);
    }
    // modify lubrication diagonal block k_ll
    {
      Teuchos::RCP<Epetra_Vector> tmp_ll
        = LINALG::CreateVector(lubrication_->LubricationField()->SystemMatrix()->RowMap(),false);
      tmp_ll->PutScalar(dti);
      Teuchos::RCP<Epetra_Vector> diag_ll
        = LINALG::CreateVector(lubrication_->LubricationField()->SystemMatrix()->RowMap(),false);
      lubrication_->LubricationField()->SystemMatrix()->ExtractDiagonalCopy(*diag_ll);
      diag_ll->Update(1.0, *tmp_ll, 1.0);
      lubrication_->LubricationField()->SystemMatrix()->ReplaceDiagonalValues(*diag_ll);
    }

    // create the linear system including PTC-modified systemmatrices k_ss and k_ll
    // \f$J(x_i) \Delta x_i = - R(x_i)\f$
    // create the systemmatrix
    SetupSystemMatrix();

    // check whether we have a sanely filled tangent matrix
    if (not systemmatrix_->Filled())
      dserror("Effective tangent matrix must be filled here");

    // create full monolithic rhs vector
    // make negative residual not necessary: rhs_ is already negative
    // (STR/LUBRICATION)-RHS is put negative in PrepareSystemForNewtonSolve()
    SetupRHS();

    // (Newton-ready) residual with blanked Dirichlet DOFs (see adapter_timint!)
    // is done in PrepareSystemForNewtonSolve() within Evaluate(iterinc_)
    LinearSolve();

    // reset solver tolerance
    solver_->ResetTolerance();

    // --------------------------------------------- build residual norms
    // include all stuff here related with convergence test
    normrhs_ = CalculateVectorNorm(iternorm_, rhs_);
    // vector of displacement and pressure residual
    Teuchos::RCP<Epetra_Vector> strrhs;
    Teuchos::RCP<Epetra_Vector> lubricationrhs;
    // extract field vectors
    ExtractFieldVectors(rhs_, strrhs, lubricationrhs);
    normstrrhs_ = CalculateVectorNorm(iternormstr_, strrhs);
    normlubricationrhs_ = CalculateVectorNorm(iternormlubrication_, lubricationrhs);

    // --------------------------------- build residual incremental norms
    // vector of displacement and pressure increments
    Teuchos::RCP<Epetra_Vector> sx;
    Teuchos::RCP<Epetra_Vector> lx;
    // extract field vectors
    ExtractFieldVectors(iterinc_,sx,lx);
    norminc_ = CalculateVectorNorm(iternorm_, iterinc_);
    normdisi_ = CalculateVectorNorm(iternormstr_, sx);
    normprei_ = CalculateVectorNorm(iternormlubrication_, lx);

    // in case of 'Mix'-convergence criterion: save the norm of the 1st
    // iteration in norm*iter0_
    if (iter_ == 1)
    {
      // set residuals
      normrhsiter0_ = normrhs_;
      normstrrhsiter0_ = normstrrhs_;
      normlubricationrhsiter0_ = normlubricationrhs_;
      // set increments
      norminciter0_ = norminc_;
      normdisiiter0_ = normdisi_;
      normpreiiter0_ = normprei_;

      // we set the minimum of iter0_ and tolrhs_, because
      // we want to prevent the case of a zero characteristic initial norm
      if (normrhsiter0_ == 0.0) normrhsiter0_ = tolrhs_;
      if (normstrrhsiter0_ == 0.0) normstrrhsiter0_ = tolstrrhs_;
      if (normlubricationrhsiter0_ == 0.0) normlubricationrhsiter0_ = tollubricationrhs_;
      if (norminciter0_ == 0.0) norminciter0_ = tolinc_;
      if (normdisiiter0_ == 0.0) normdisiiter0_ = toldisi_;
      if (normpreiiter0_ == 0.0) normpreiiter0_ = tolprei_;
    }

    // print stuff
    PrintNewtonIter();

    // increment equilibrium loop index
    iter_ += 1;

    // save old pseudo-time step in dti_
    dti_ = dti;

    // update ptc
    {
      double np = 0.0;
      rhs_->NormInf(&np);
      dti *= (np/nc);
      dti = std::max(dti,0.0);
      nc = np;
    }

  }  // end equilibrium loop

  // ----------------------------------------------------- iteration loop

  // correct iteration counter
  iter_ -= 1;

  // test whether max iterations was hit
  if ( (Converged()) and (Comm().MyPID() == 0) )
  {
    PrintNewtonConv();
  }
  else if (iter_ >= itermax_)
    dserror("PTC unconverged in %d iterations", iter_);

}  // PTC()


/*----------------------------------------------------------------------*
 | evaluate the single fields                               wirtz 01/16 |
 *----------------------------------------------------------------------*/
void EHL::Monolithic::Evaluate(Teuchos::RCP<Epetra_Vector> x)
{

  TEUCHOS_FUNC_TIME_MONITOR("EHL::Monolithic::Evaluate");

  // displacement and pressure incremental vector
  Teuchos::RCP<Epetra_Vector> sx;
  Teuchos::RCP<Epetra_Vector> lx;

  // if an increment vector exists
  if (x != Teuchos::null)
  {
    // extract displacement sx and pressure lx incremental vector of global
    // unknown incremental vector x
    ExtractFieldVectors(x,sx,lx);
  }

  // else (x == Teuchos::null): initialise the system

  // Newton update of the lubrication field
  // update pressure before passed to the structural field
  //   UpdateIterIncrementally(lx),
  lubrication_->LubricationField()->UpdateNewton(lx);

  // call all elements and assemble rhs and matrices
  /// structural field

  // structure Evaluate (builds tangent, residual and applies DBC)
  Epetra_Time timerstructure(Comm());

  // apply current pressure to structure
  ApplyLubricationCouplingState(lubrication_->LubricationField()->Prenp(),lx);

  // Monolithic EHL accesses the linearised structure problem:
  //   UpdaterIterIncrementally(sx),
  //   EvaluateForceStiffResidual()
  //   PrepareSystemForNewtonSolve()
  //     blank residual DOFs that are on Dirichlet BC
  //     in case of local coordinate systems rotate the residual forth and back
  //     Be AWARE: ApplyDirichlettoSystem has to be called with rotated stiff_!
  StructureField()->Evaluate(sx);
  StructureField()->Discretization()->ClearState(true);

  /// lubrication field

  // lubrication Evaluate
  // (builds tangent, residual and applies DBC and recent coupling values)
  Epetra_Time timerlubrication(Comm());

  // apply current displacements and velocities to the lubrication field
//  if (strmethodname_ == INPAR::STR::dyna_statics)
//  {
//    // calculate velocity V_n+1^k = (D_n+1^k-D_n)/Dt()
//    vel_ = CalcVelocity(StructureField()->Dispnp());
//  }
  // else: use velnp
//  else
    vel_ = StructureField()->Velnp();

  // pass the structural values to the lubrication field
  ApplyStructCouplingState(StructureField()->Dispnp(),vel_);

  // monolithic EHL accesses the linearised lubrication problem
  //   EvaluateRhsTangResidual() and
  //   PrepareSystemForNewtonSolve()
  lubrication_->LubricationField()->Evaluate();
  lubrication_->LubricationField()->Discretization()->ClearState(true);

}  // Evaluate()


/*----------------------------------------------------------------------*
 | extract field vectors for calling Evaluate() of the      wirtz 01/16 |
 | single fields                                                        |
 *----------------------------------------------------------------------*/
void EHL::Monolithic::ExtractFieldVectors(
  Teuchos::RCP<Epetra_Vector> x,
  Teuchos::RCP<Epetra_Vector>& sx,
  Teuchos::RCP<Epetra_Vector>& lx
  )
{
  TEUCHOS_FUNC_TIME_MONITOR("EHL::Monolithic::ExtractFieldVectors");

  // process structure unknowns of the first field
  sx = Extractor()->ExtractVector(x,0);

  // process lubrication unknowns of the second field
  lx = Extractor()->ExtractVector(x,1);
}  // ExtractFieldVectors()


/*----------------------------------------------------------------------*
 | full monolithic dof row map                              wirtz 01/16 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> EHL::Monolithic::DofRowMap() const
{
  return Extractor()->FullMap();
}  // DofRowMap()


/*----------------------------------------------------------------------*
 | setup system (called in ehl_dyn)                         wirtz 01/16 |
 *----------------------------------------------------------------------*/
void EHL::Monolithic::SetupSystem()
{

  // set parameters that remain the same in the whole calculation
  SetDefaultParameters();

  // create combined map
  std::vector<Teuchos::RCP<const Epetra_Map> > vecSpaces;

  // use its own DofRowMap, that is the 0th map of the discretization
  vecSpaces.push_back(StructureField()->DofRowMap(0));
  vecSpaces.push_back(lubrication_->LubricationField()->DofRowMap(0));

  if (vecSpaces[0]->NumGlobalElements() == 0)
    dserror("No structure equation. Panic.");
  if (vecSpaces[1]->NumGlobalElements() == 0)
    dserror("No pressure equation. Panic.");

  SetDofRowMaps(vecSpaces);

  /*----------------------------------------------------------------------*/
  // initialise EHL-systemmatrix_
  systemmatrix_
    = Teuchos::rcp(
        new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
              *Extractor(),
              *Extractor(),
              81,
              false,
              true
              )
        );

  // create empty matrix
  k_sl_ = Teuchos::rcp(
           new LINALG::SparseMatrix(
                 *(StructureField()->Discretization()->DofRowMap(0)),
                 81,
                 true,
                 true
                 )
           );

  // create empty matrix
  k_ls_ = Teuchos::rcp(
           new LINALG::SparseMatrix(
                 *(lubrication_->LubricationField()->Discretization()->DofRowMap(0)),
                 81,
                 true,
                 true
                 )
           );

}  // SetupSystem()


/*----------------------------------------------------------------------*
 | put the single maps to one full EHL map together         wirtz 01/16 |
 *----------------------------------------------------------------------*/
void EHL::Monolithic::SetDofRowMaps(
  const std::vector<Teuchos::RCP<const Epetra_Map> >& maps
  )
{
  Teuchos::RCP<Epetra_Map> fullmap
    = LINALG::MultiMapExtractor::MergeMaps(maps);

  // full EHL-blockmap
  Extractor()->Setup(*fullmap,maps);
}  // SetDofRowMaps()


/*----------------------------------------------------------------------*
 | setup system matrix of EHL                               wirtz 01/16 |
 *----------------------------------------------------------------------*/
void EHL::Monolithic::SetupSystemMatrix()
{

  TEUCHOS_FUNC_TIME_MONITOR("EHL::Monolithic::SetupSystemMatrix");

  /*----------------------------------------------------------------------*/
  // pure structural part k_ss (3nx3n)

  // build pure structural block k_ss
  // build block matrix
  // The maps of the block matrix have to match the maps of the blocks we
  // insert here. Extract Jacobian matrices and put them into composite system
  // matrix W
  // k_ss: already applied ApplyDirichletWithTrafo() in PrepareSystemToNewtonSolve
  Teuchos::RCP<LINALG::SparseMatrix> k_ss = StructureField()->SystemMatrix();

  // build block matrix
  // The maps of the block matrix have to match the maps of the blocks we
  // insert here.
  k_ss->UnComplete();

  // assign structure part to the EHL matrix
  systemmatrix_->Assign(0,0,LINALG::View,*k_ss);

  /*----------------------------------------------------------------------*/
  // structural block k_sl (3nxn)
  // build mechanical-lubrication block
#if(0)
  k_sl_->Reset();
  // call the element and calculate the matrix block
  ApplyStrCouplMatrix(k_sl_);

  // apply dirichlet boundary conditions properly on matrix k_sl, i.e. blank row
  // if dof is a structural DBC
  // Normally, DBC should be applied on complete systemmatrix k_EHL, but for
  // diagonal blocks (here k_ss, k_ll) DBC are ALREADY applied in
  // PrepareSystemForNewtonSolve() included in Evaluate(sx)
  //
  // to avoid double work, we only call ApplyDirichlet for the off-diagonal blocks,
  // here k_sl
  // k_sl is an off-diagonal block --> pass the bool diagonal==false
  // ApplyDirichlet*() expect filled matrix
  //
  // in case of inclined STR-DBC
  //   1.) transform the off-diagonal block k_sl to the local system --> k_sl^{~}
  //   2.) apply ApplyDirichletWithTrafo() on rotated block k_sl^{~}
  //              --> blank the row, which has a DBC
  //
  // to apply Multiply in LocSys, k_sl has to be FillCompleted

  k_sl_->Complete(
          *(StructureField()->Discretization()->DofRowMap(1)),
          *(StructureField()->Discretization()->DofRowMap(0))
          );

  if (locsysman_ != Teuchos::null)
  {
    // rotate k_sl to local coordinate system --> k_sl^{~}
    locsysman_->RotateGlobalToLocal(k_sl_);
    // apply ApplyDirichletWithTrafo() on rotated block k_sl^{~}
    // --> if dof has an inclined DBC: blank the complete row, the '1.0' is set
    //     on diagonal of row, i.e. on diagonal of k_ss
    k_sl_->ApplyDirichletWithTrafo(
            locsysman_->Trafo(),
            *StructureField()->GetDBCMapExtractor()->CondMap(),
            false
            );
  }  // end locsys
  // default: (locsysman_ == Teuchos::null), i.e. NO inclined Dirichlet BC
  else
    k_sl_->ApplyDirichlet(*StructureField()->GetDBCMapExtractor()->CondMap(),false);

  k_sl_->UnComplete();

  // assign lubrication part to the EHL matrix
  systemmatrix_->Assign(0,1,LINALG::View,*(k_sl_));
#endif
  /*----------------------------------------------------------------------*/
  // pure lubrication part k_ll (nxn)

  // build pure lubrication block k_ll
  // build block matrix
  // The maps of the block matrix have to match the maps of the blocks we
  // insert here. Extract Jacobian matrices and put them into composite system
  // matrix systemmatrix_
  Teuchos::RCP<LINALG::SparseMatrix> k_ll = lubrication_->LubricationField()->SystemMatrix();

  // Uncomplete lubrication matrix to be able to deal with slightly defective
  // interface meshes.
  k_ll->UnComplete();

  // assign lubrication part to the EHL matrix
  systemmatrix_->Assign(1,1,LINALG::View,*(k_ll));

  /*----------------------------------------------------------------------*/
  // lubrication part k_ls (nx3n)
  // build lubrication-mechanical block
#if(0)
  k_ls_->Reset();

  // call the element and calculate the matrix block
  ApplyLubricationCouplMatrix(k_ls_);

  // apply dirichlet boundary conditions properly on matrix k_ls, i.e. blank row
  // if dof is a lubrication DBC
  // Normally, DBC should be applied on full systemmatrix k_EHL, but on diagonal
  // blocks (here k_ss, k_ll) DBC are already applied in PrepareSystemForNewtonSolve()
  // to avoid double work, we only call ApplyDirichlet at the off-diagonal blokcs,
  // here k_ls
  // k_ls is an off-diagonal block --> pass the bool diagonal==false
  // ApplyDirichlet() expect filled matrix

  k_ls_->Complete(
          *(lubrication_->LubricationField()->Discretization()->DofRowMap(1)),
          *(lubrication_->LubricationField()->Discretization()->DofRowMap(0))
          );

  // assign lubrication part to the EHL matrix
  k_ls_->ApplyDirichlet(*lubrication_->LubricationField()->GetDBCMapExtractor()->CondMap(),false);
  k_ls_->UnComplete();
  systemmatrix_->Assign(1,0,LINALG::View,*k_ls_);
#endif
  /*----------------------------------------------------------------------*/
  // done. make sure all blocks are filled.
  systemmatrix_->Complete();

}  // SetupSystemMatrix()


/*----------------------------------------------------------------------*
 | setup RHS (like fsimon)                                  wirtz 01/16 |
 *----------------------------------------------------------------------*/
void EHL::Monolithic::SetupRHS()
{

  TEUCHOS_FUNC_TIME_MONITOR("EHL::Monolithic::SetupRHS");

  // create full monolithic rhs vector
  rhs_ = Teuchos::rcp(new Epetra_Vector(*DofRowMap(), true));

  // fill the EHL rhs vector rhs_ with the single field rhs
  SetupVector(
    *rhs_,
    StructureField()->RHS(),
    lubrication_->LubricationField()->RHS()
    );

}  // SetupRHS()


/*----------------------------------------------------------------------*
 | solve linear EHL system                                  wirtz 01/16 |
 *----------------------------------------------------------------------*/
void EHL::Monolithic::LinearSolve()
{

  // Solve for inc_ = [disi_,prei_]
  // Solve K_Teffdyn . IncX = -R  ===>  IncX_{n+1} with X=[d,T]
  // \f$x_{i+1} = x_i + \Delta x_i\f$
  if (solveradapttol_ and (iter_ > 1))
  {
    double worst = normrhs_;
    double wanted = tolrhs_;
    solver_->AdaptTolerance(wanted, worst, solveradaptolbetter_);
  }

  // Dirichlet boundary conditions are already applied to EHL system, i.e. EHL
  // system is prepared for solve, i.e. EHL systemmatrix, EHL rhs, EHL inc
  // --> in PrepareSystemForNewtonSolve(): done for rhs and diagonal blocks
  // --> in SetupSystemMatrix() done for off-diagonal blocks k_sl, k_ls

  // apply Dirichlet BCs to system of equations
  iterinc_->PutScalar(0.0);  // Useful? depends on solver and more

  // default: use block matrix
  if (merge_ehl_blockmatrix_ == false)
  {
  // Infnormscaling: scale system before solving
  ScaleSystem(*systemmatrix_,*rhs_);

  // solve the problem, work is done here!
  solver_->Solve(
             systemmatrix_->EpetraOperator(),
             iterinc_,
             rhs_,
             true,
             iter_==1
             );

  // Infnormscaling: unscale system after solving
  UnscaleSolution(*systemmatrix_, *iterinc_, *rhs_);

  }  // use block matrix

  else // (merge_ehl_blockmatrix_ == true)
  {
    // merge blockmatrix to SparseMatrix and solve
    Teuchos::RCP<LINALG::SparseMatrix> sparse = systemmatrix_->Merge();

    // standard solver call
    solver_->Solve(
               sparse->EpetraOperator(),
               iterinc_,
               rhs_,
               true,
               iter_==1
               );
  }  // MergeBlockMatrix

}  // LinearSolve()


/*----------------------------------------------------------------------*
 | setup vector of the structure and lubrication field      wirtz 01/16 |
 *----------------------------------------------------------------------*/
void EHL::Monolithic::SetupVector(
  Epetra_Vector &f,
  Teuchos::RCP<const Epetra_Vector> sv,
  Teuchos::RCP<const Epetra_Vector> lv
  )
{
  // extract dofs of the two fields
  // and put the structural/lubrication field vector into the global vector f
  // noticing the block number
  Extractor()->InsertVector(*sv,0,f);
  Extractor()->InsertVector(*lv,1,f);

}  // SetupVector()


/*----------------------------------------------------------------------*
 | check convergence of Newton iteration (public)           wirtz 01/16 |
 *----------------------------------------------------------------------*/
bool EHL::Monolithic::Converged()
{
  // check for single norms
  bool convrhs = false;
  bool convinc = false;
  bool convstrrhs = false;
  bool convdisp = false;
  bool convlubricationrhs = false;
  bool convpre = false;

  // ----------------------------------------------------------- EHL test
  // residual EHL forces
  switch (normtyperhs_)
  {
  case INPAR::EHL::convnorm_abs :
    convrhs = normrhs_ < tolrhs_;
    break;
  case INPAR::EHL::convnorm_rel :
    convrhs = normrhs_ < std::max(tolrhs_*normrhsiter0_, 1.0e-15);
    break;
  case INPAR::EHL::convnorm_mix :
    convrhs = ( (normrhs_ < tolrhs_) and (normrhs_ < std::max(normrhsiter0_*tolrhs_, 1.0e-15)) );
    break;
  default:
    dserror("Cannot check for convergence of residual forces!");
    break;
  }

  // residual EHL increments
  switch (normtypeinc_)
  {
  case INPAR::EHL::convnorm_abs :
    convinc = norminc_ < tolinc_;
    break;
  case INPAR::EHL::convnorm_rel :
    convinc = norminc_ < std::max(norminciter0_*tolinc_,1e-15);
    break;
  case INPAR::EHL::convnorm_mix :
    convinc = norminc_ < std::max(norminciter0_*tolinc_,1e-15);
    break;
  default:
    dserror("Cannot check for convergence of increments!");
    break;
  }  // switch (normtypeinc_)

  // -------------------------------------------------- single field test
  // ---------------------------------------------------------- structure
  // structural residual forces
  switch (normtypestrrhs_)
  {
  case INPAR::STR::convnorm_abs :
    convstrrhs = normstrrhs_ < tolstrrhs_;
    break;
  case INPAR::STR::convnorm_rel :
    convstrrhs = normstrrhs_ < std::max(normstrrhsiter0_*tolstrrhs_,1e-15);
    break;
  case INPAR::STR::convnorm_mix :
    convstrrhs = ( (normstrrhs_ < tolstrrhs_) or
                   (normstrrhs_ < std::max(normstrrhsiter0_*tolstrrhs_,1e-15))
                 );
    break;
  default :
    dserror("Cannot check for convergence of residual forces!");
    break;
  }  // switch (normtypestrrhs_)

  // residual displacements
  switch (normtypedisi_)
  {
  case INPAR::STR::convnorm_abs :
    convdisp = normdisi_ < toldisi_;
    break;
  case INPAR::STR::convnorm_rel :
    convdisp = normdisi_ < std::max(normdisiiter0_*toldisi_,1e-15);
    break;
  case INPAR::STR::convnorm_mix :
    convdisp = ( (normdisi_ < toldisi_) or
                 (normdisi_ < std::max(normdisiiter0_*toldisi_,1e-15))
               );
    break;
  default :
    dserror("Cannot check for convergence of displacements!");
    break;
  }  // switch (normtypedisi_)

  // ------------------------------------------------------------- lubrication
  // lubrication residual forces
  switch (normtypelubricationrhs_)
  {
  case INPAR::LUBRICATION::convnorm_abs :
    convlubricationrhs = normlubricationrhs_ < tollubricationrhs_;
    break;
  case INPAR::LUBRICATION::convnorm_rel :
    convlubricationrhs = normlubricationrhs_ < normlubricationrhsiter0_*tollubricationrhs_;
    break;
  case INPAR::LUBRICATION::convnorm_mix :
    convlubricationrhs = ( (normlubricationrhs_ < tollubricationrhs_) or (normlubricationrhs_ < normlubricationrhsiter0_*tollubricationrhs_) );
    break;
  default :
    dserror("Cannot check for convergence of residual forces!");
    break;
  }  // switch (normtypelubricationrhs_)

  // residual pressures
  switch (normtypeprei_)
  {
  case INPAR::LUBRICATION::convnorm_abs :
    convpre = normprei_ < tolprei_;
    break;
  case INPAR::LUBRICATION::convnorm_rel :
    convpre = normprei_ < normpreiiter0_*tolprei_;
    break;
  case INPAR::LUBRICATION::convnorm_mix :
    convpre = ( (normprei_ < tolprei_) or (normprei_ < normpreiiter0_*tolprei_) );
    break;
  default :
    dserror("Cannot check for convergence of pressures!");
    break;
  }  // switch (normtypeprei_)

  // -------------------------------------------------------- convergence
  // combine increment-like and force-like residuals, combine EHL and single
  // field values
  bool conv = false;
  if (combincrhs_ == INPAR::EHL::bop_and)
    conv = convinc and convrhs;
  else if (combincrhs_ == INPAR::EHL::bop_or)
    conv = convinc or convrhs;
  else if (combincrhs_ == INPAR::EHL::bop_coupl_and_singl)
    conv = convinc and convrhs and convdisp and convstrrhs and convpre and convlubricationrhs;
  else if (combincrhs_ == INPAR::EHL::bop_coupl_or_singl)
    conv = (convinc and convrhs) or (convdisp and convstrrhs and convpre and convlubricationrhs);
  else if (combincrhs_ == INPAR::EHL::bop_and_singl)
    conv = convdisp and convstrrhs and convpre and convlubricationrhs;
  else if (combincrhs_ == INPAR::EHL::bop_or_singl)
    conv = (convdisp or convstrrhs or convpre or convlubricationrhs);
  else
    dserror("Something went terribly wrong with binary operator!");

  // return things
  return conv;

}  // Converged()


/*----------------------------------------------------------------------*
 | print Newton-Raphson iteration to screen and error file  wirtz 01/16 |
 | originally by lw 12/07, tk 01/08                                     |
 *----------------------------------------------------------------------*/
void EHL::Monolithic::PrintNewtonIter()
{
  // print to standard out
  // replace myrank_ here general by Comm().MyPID()
  if ( (Comm().MyPID() == 0) and PrintScreenEvry() and
       (Step()%PrintScreenEvry() == 0) and printiter_
     )
  {
    if (iter_== 1)
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

  // see you
  return;
}  // PrintNewtonIter()


/*----------------------------------------------------------------------*
 | print Newton-Raphson iteration to screen and error file  wirtz 01/16 |
 | originally by lw 12/07, tk 01/08                                     |
 *----------------------------------------------------------------------*/
void EHL::Monolithic::PrintNewtonIterHeader(FILE* ofile)
{
  // open outstringstream
  std::ostringstream oss;

  // enter converged state etc
  oss << std::setw(6)<< "numiter";

  // ---------------------------------------------------------------- EHL
  // different style due relative or absolute error checking
  // displacement
  switch (normtyperhs_)
  {
  case INPAR::EHL::convnorm_abs :
    oss <<std::setw(15)<< "abs-res-norm";
    break;
  case INPAR::EHL::convnorm_rel :
    oss <<std::setw(15)<< "rel-res-norm";
    break;
  case INPAR::EHL::convnorm_mix :
    oss << std::setw(15)<< "mix-res-norm";
    break;
  default:
    dserror("You should not turn up here.");
    break;
  }

  switch (normtypeinc_)
  {
  case INPAR::EHL::convnorm_abs :
    oss << std::setw(15)<< "abs-inc-norm";
    break;
  case INPAR::EHL::convnorm_rel :
    oss << std::setw(15)<< "rel-inc-norm";
    break;
  default:
    dserror("You should not turn up here.");
    break;
  }

  // -------------------------------------------------- single field test
  // ---------------------------------------------------------- structure
  switch (normtypestrrhs_)
  {
  case INPAR::STR::convnorm_rel :
    oss << std::setw(18)<< "rel-str-res-norm";
    break;
  case INPAR::STR::convnorm_abs :
    oss << std::setw(18)<< "abs-str-res-norm";
    break;
  case INPAR::STR::convnorm_mix :
    oss << std::setw(18)<< "mix-str-res-norm";
    break;
  default:
    dserror("You should not turn up here.");
    break;
  }  // switch (normtypestrrhs_)

  switch (normtypedisi_)
  {
  case INPAR::STR::convnorm_rel :
    oss << std::setw(16)<< "rel-dis-norm";
    break;
  case INPAR::STR::convnorm_abs :
    oss << std::setw(16)<< "abs-dis-norm";
    break;
  case INPAR::STR::convnorm_mix :
    oss << std::setw(16)<< "mix-dis-norm";
    break;
  default :
    dserror("You should not turn up here.");
    break;
  }  // switch (normtypedisi_)

  // ------------------------------------------------------------- lubrication
  switch (normtypelubricationrhs_)
  {
  case INPAR::LUBRICATION::convnorm_rel :
    oss << std::setw(18)<< "rel-lub-res-norm";
    break;
  case INPAR::LUBRICATION::convnorm_abs :
    oss << std::setw(18)<< "abs-lub-res-norm";
    break;
  case INPAR::LUBRICATION::convnorm_mix :
    oss << std::setw(18)<< "mix-lub-res-norm";
    break;
  default :
    dserror("You should not turn up here.");
    break;
  }  // switch (normtypelubricationrhs_)

  switch (normtypeprei_)
  {
  case INPAR::LUBRICATION::convnorm_rel :
    oss << std::setw(16)<< "rel-pre-norm";
    break;
  case INPAR::LUBRICATION::convnorm_abs :
    oss << std::setw(16)<< "abs-pre-norm";
    break;
  case INPAR::LUBRICATION::convnorm_mix :
    oss << std::setw(16)<< "mix-pre-norm";
    break;
  default :
    dserror("You should not turn up here.");
    break;
  }  // switch (normtypeprei_)

  if (soltech_ == INPAR::EHL::soltech_ptc)
  {
    oss << std::setw(16)<< "        PTC-dti";
  }

  // add solution time
  oss << std::setw(12)<< "ts";
  oss << std::setw(12)<< "wct";

  // finish oss
  oss << std::ends;

  // print to screen (could be done differently...)
  if (ofile == NULL)
    dserror("no ofile available");
  fprintf(ofile, "%s\n", oss.str().c_str());

  // print it, now
  fflush(ofile);

  // nice to have met you
  return;

}  // PrintNewtonIterHeader()


/*----------------------------------------------------------------------*
 | print Newton-Raphson iteration to screen                 wirtz 01/16 |
 | originally by lw 12/07, tk 01/08                                     |
 *----------------------------------------------------------------------*/
void EHL::Monolithic::PrintNewtonIterText(FILE* ofile)
{
  // open outstringstream
  std::ostringstream oss;

  // enter converged state etc
  oss << std::setw(7)<< iter_;

  // different style due relative or absolute error checking

  // ----------------------------------------------- test coupled problem
  switch (normtyperhs_)
  {
  case INPAR::EHL::convnorm_abs :
    oss << std::setw(15) << std::setprecision(5) << std::scientific << normrhs_;
    break;
  case INPAR::EHL::convnorm_rel :
    oss << std::setw(15) << std::setprecision(5) << std::scientific << normrhs_/normrhsiter0_;
    break;
  case INPAR::EHL::convnorm_mix :
    oss << std::setw(15) << std::setprecision(5) << std::scientific << std::min(normrhs_, normrhs_/normrhsiter0_);
    break;
  default :
    dserror("You should not turn up here.");
    break;
  }

  switch (normtypeinc_)
  {
  case INPAR::EHL::convnorm_abs :
    oss << std::setw(15) << std::setprecision(5) << std::scientific << norminc_;
    break;
  case INPAR::EHL::convnorm_rel :
    oss << std::setw(15) << std::setprecision(5) << std::scientific << norminc_/norminciter0_;
    break;
  case INPAR::EHL::convnorm_mix :
    oss << std::setw(15) << std::setprecision(5) << std::scientific << std::min(norminc_, norminc_/norminciter0_);
    break;
  default :
    dserror("You should not turn up here.");
    break;
  }  // switch (normtypeinc_)

  // ------------------------------------------------- test single fields
  // ---------------------------------------------------------- structure
  // different style due relative or absolute error checking
  // displacement
  switch (normtypestrrhs_)
  {
  case INPAR::STR::convnorm_abs :
    oss << std::setw(18) << std::setprecision(5) << std::scientific << normstrrhs_;
    break;
  case INPAR::STR::convnorm_rel :
    oss << std::setw(18) << std::setprecision(5) << std::scientific << normstrrhs_/normstrrhsiter0_;
    break;
  case INPAR::STR::convnorm_mix :
    oss << std::setw(18) << std::setprecision(5) << std::scientific << std::min(normstrrhs_, normstrrhs_/normstrrhsiter0_);
    break;
  default :
    dserror("You should not turn up here.");
    break;
  }  // switch (normtypestrrhs_)

  switch (normtypedisi_)
  {
  case INPAR::STR::convnorm_abs :
    oss << std::setw(16) << std::setprecision(5) << std::scientific << normdisi_;
    break;
  case INPAR::STR::convnorm_rel :
    oss << std::setw(16) << std::setprecision(5) << std::scientific << normdisi_/normdisiiter0_;
    break;
  case INPAR::STR::convnorm_mix :
   oss << std::setw(16) << std::setprecision(5) << std::scientific << std::min(normdisi_, normdisi_/normdisiiter0_);
    break;
  default :
    dserror("You should not turn up here.");
    break;
  }  // switch (normtypedisi_)

  // ------------------------------------------------------------- lubrication
  switch (normtypelubricationrhs_)
  {
  case INPAR::LUBRICATION::convnorm_abs :
    oss << std::setw(18) << std::setprecision(5) << std::scientific << normlubricationrhs_;
    break;
  case INPAR::LUBRICATION::convnorm_rel :
    oss << std::setw(18) << std::setprecision(5) << std::scientific << normlubricationrhs_/normlubricationrhsiter0_;
    break;
  case INPAR::LUBRICATION::convnorm_mix :
    oss << std::setw(18) << std::setprecision(5) << std::scientific << std::min(normlubricationrhs_, normlubricationrhs_/normlubricationrhsiter0_);
    break;
  default :
    dserror("You should not turn up here.");
    break;
  }  // switch (normtypelubricationrhs_)

  switch (normtypeprei_)
  {
  case INPAR::LUBRICATION::convnorm_abs :
    oss << std::setw(16) << std::setprecision(5) << std::scientific << normprei_;
    break;
  case INPAR::LUBRICATION::convnorm_rel :
    oss << std::setw(16) << std::setprecision(5) << std::scientific << normprei_/normpreiiter0_;
    break;
  case INPAR::LUBRICATION::convnorm_mix :
    oss << std::setw(16) << std::setprecision(5) << std::scientific << std::min(normprei_, normprei_/normpreiiter0_);
    break;
  default :
    dserror("You should not turn up here.");
    break;
  }  // switch (normtypeprei_)

  if (soltech_ == INPAR::EHL::soltech_ptc)
  {
    oss << std::setw(16) << std::setprecision(5) << std::scientific << dti_;
  }

  // add solution time of to print to screen
  oss << std::setw(12) << std::setprecision(2) << std::scientific << dtsolve_;
  oss << std::setw(12) << std::setprecision(2) << std::scientific << timernewton_.ElapsedTime();

  // finish oss
  oss << std::ends;

  // print to screen (could be done differently...)
  if (ofile == NULL)
    dserror("no ofile available");
  fprintf(ofile, "%s\n", oss.str().c_str());

  // print it, now
  fflush(ofile);

  // nice to have met you
  return;

}  // PrintNewtonIterText


/*----------------------------------------------------------------------*
 | print statistics of converged NRI                        wirtz 01/16 |
 | orignially by bborn 08/09                                            |
 *----------------------------------------------------------------------*/
void EHL::Monolithic::PrintNewtonConv()
{
  // somebody did the door
  return;
}  // PrintNewtonConv()


/*----------------------------------------------------------------------*
 | evaluate mechanical-lubrication system matrix at state   wirtz 01/16 |
 *----------------------------------------------------------------------*/
void EHL::Monolithic::ApplyStrCouplMatrix(
  Teuchos::RCP<LINALG::SparseMatrix> k_sl  //!< off-diagonal tangent matrix term
  )
{

  // create the parameters for the discretization
  Teuchos::ParameterList sparams;
  const std::string action = "calc_struct_stiffpre";
  sparams.set("action", action);
  // other parameters that might be needed by the elements
  sparams.set("delta time", Dt());
  sparams.set("total time", Time());

  StructureField()->Discretization()->ClearState(true);
  StructureField()->Discretization()->SetState(0,"displacement",StructureField()->Dispnp());

  ApplyLubricationCouplingState(lubrication_->LubricationField()->Prenp());

  // build specific assemble strategy for mechanical-lubrication system matrix
  // from the point of view of StructureField:
  // structdofset = 0, thermdofset = 1
  DRT::AssembleStrategy structuralstrategy(
                          0,  // structdofset for row
                          1,  // thermdofset for column
                          k_sl,  // build mechanical-lubrication matrix
                          Teuchos::null,  // no other matrix or vectors
                          Teuchos::null,
                          Teuchos::null,
                          Teuchos::null
                          );

  // evaluate the mechanical-lubrication system matrix on the structural element
  StructureField()->Discretization()->Evaluate(sparams,structuralstrategy);
  StructureField()->Discretization()->ClearState(true);

  // TODO 2013-11-11 move scaling to the so3_lubrication element
  // --> consistent with lubrication element and clearer, more consistent

  // for consistent linearisation scale k_sl with time factor
  // major switch to different time integrators
  switch (strmethodname_)
  {
  case INPAR::STR::dyna_statics :
  {
    // continue
    break;
  }
  case INPAR::STR::dyna_onesteptheta :
  {
    double theta = sdyn_.sublist("ONESTEPTHETA").get<double>("THETA");
    // K_Teffdyn(T_n+1^i) = theta * k_sl
    k_sl->Scale(theta);
    break;
  }
  case INPAR::STR::dyna_genalpha :
  {
    double alphaf = sdyn_.sublist("GENALPHA").get<double>("ALPHA_F");
    // K_Teffdyn(T_n+1) = (1-alphaf_) . kst
    // Lin(dT_n+1-alphaf_/ dT_n+1) = (1-alphaf_)
    k_sl->Scale(1.0 - alphaf);
    break;
  }
  default :
    dserror("Don't know what to do...");
    break;
  }  // end of switch(strmethodname_)

}  // ApplyStrCouplMatrix()


/*----------------------------------------------------------------------*
 | evaluate lubrication-mechanical system matrix at state   wirtz 01/16 |
 *----------------------------------------------------------------------*/
void EHL::Monolithic::ApplyLubricationCouplMatrix(
  Teuchos::RCP<LINALG::SparseMatrix> k_ls  //!< off-diagonal tangent matrix term
  )
{

  // create the parameters for the discretization
  Teuchos::ParameterList lparams;
  // action for elements
  const LUBRICATION::Action action = LUBRICATION::calc_lubrication_coupltang;
  lparams.set<int>("action", action);
  // other parameters that might be needed by the elements
  lparams.set("delta time", Dt());
  lparams.set("total time", Time());

  lubrication_->LubricationField()->Discretization()->ClearState(true);
  // set the variables that are needed by the elements
  lubrication_->LubricationField()->Discretization()->SetState(0,"pressure",lubrication_->LubricationField()->Prenp());

  ApplyStructCouplingState(StructureField()->Dispnp(),vel_);

  // build specific assemble strategy for the lubrication-mechanical system matrix
  // from the point of view of lubrication_->LubricationField:
  // thermdofset = 0, structdofset = 1
  DRT::AssembleStrategy lubricationstrategy(
                          0,  // thermdofset for row
                          1,  // structdofset for column
                          k_ls,  // lubrication-mechanical matrix
                          Teuchos::null,  // no other matrix or vectors
                          Teuchos::null,
                          Teuchos::null,
                          Teuchos::null
                          );

  // evaluate the lubrication-mechanical system matrix on the lubrication element
  lubrication_->LubricationField()->Discretization()->Evaluate(lparams,lubricationstrategy);
  lubrication_->LubricationField()->Discretization()->ClearState(true);

}  // ApplyLubricationCouplMatrix()

/*----------------------------------------------------------------------*
 | map containing the dofs with Dirichlet BC                wirtz 01/16 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> EHL::Monolithic::CombinedDBCMap()
{
  const Teuchos::RCP<const Epetra_Map > scondmap
    = StructureField()->GetDBCMapExtractor()->CondMap();
  const Teuchos::RCP<const Epetra_Map > lcondmap
    = lubrication_->LubricationField()->GetDBCMapExtractor()->CondMap();
  Teuchos::RCP<Epetra_Map> condmap = LINALG::MergeMap(scondmap, lcondmap, false);
  return condmap;

}  // CombinedDBCMap()

/*----------------------------------------------------------------------*
 | scale system, i.e. apply infnorm scaling to linear       wirtz 01/16 |
 | block system before solving system                                   |
 *----------------------------------------------------------------------*/
void EHL::Monolithic::ScaleSystem(
  LINALG::BlockSparseMatrixBase& mat,
  Epetra_Vector& b
  )
{
  //should we scale the system?
  const bool scaling_infnorm
    = (bool)DRT::INPUT::IntegralValue<int>(ehldynmono_,"INFNORMSCALING");

  if (scaling_infnorm)
  {
    // The matrices are modified here. Do we have to change them back later on?

    Teuchos::RCP<Epetra_CrsMatrix> A = mat.Matrix(0,0).EpetraMatrix();
    srowsum_ = Teuchos::rcp(new Epetra_Vector(A->RowMap(),false));
    scolsum_ = Teuchos::rcp(new Epetra_Vector(A->RowMap(),false));
    A->InvRowSums(*srowsum_);
    A->InvColSums(*scolsum_);
    if ( (A->LeftScale(*srowsum_)) or (A->RightScale(*scolsum_)) or
         (mat.Matrix(0,1).EpetraMatrix()->LeftScale(*srowsum_)) or
         (mat.Matrix(1,0).EpetraMatrix()->RightScale(*scolsum_))
       )
      dserror("structure scaling failed");

    A = mat.Matrix(1,1).EpetraMatrix();
    lrowsum_ = Teuchos::rcp(new Epetra_Vector(A->RowMap(),false));
    lcolsum_ = Teuchos::rcp(new Epetra_Vector(A->RowMap(),false));
    A->InvRowSums(*lrowsum_);
    A->InvColSums(*lcolsum_);
    if ( (A->LeftScale(*lrowsum_)) or (A->RightScale(*lcolsum_)) or
         (mat.Matrix(1,0).EpetraMatrix()->LeftScale(*lrowsum_)) or
         (mat.Matrix(0,1).EpetraMatrix()->RightScale(*lcolsum_))
       )
      dserror("lubrication scaling failed");

    Teuchos::RCP<Epetra_Vector> sx = Extractor()->ExtractVector(b,0);
    Teuchos::RCP<Epetra_Vector> lx = Extractor()->ExtractVector(b,1);

    if (sx->Multiply(1.0, *srowsum_, *sx, 0.0))
      dserror("structure scaling failed");
    if (lx->Multiply(1.0, *lrowsum_, *lx, 0.0))
      dserror("lubrication scaling failed");

    Extractor()->InsertVector(*sx,0,b);
    Extractor()->InsertVector(*lx,1,b);
  }
}  // ScaleSystem


/*----------------------------------------------------------------------*
 | unscale solution after solving the linear system         wirtz 01/16 |
 *----------------------------------------------------------------------*/
void EHL::Monolithic::UnscaleSolution(
  LINALG::BlockSparseMatrixBase& mat,
  Epetra_Vector& x,
  Epetra_Vector& b
  )
{
  const bool scaling_infnorm
    = (bool)DRT::INPUT::IntegralValue<int>(ehldynmono_,"INFNORMSCALING");

  if (scaling_infnorm)
  {
    Teuchos::RCP<Epetra_Vector> sy = Extractor()->ExtractVector(x,0);
    Teuchos::RCP<Epetra_Vector> ly = Extractor()->ExtractVector(x,1);

    if (sy->Multiply(1.0, *scolsum_, *sy, 0.0))
      dserror("structure scaling failed");
    if (ly->Multiply(1.0, *lcolsum_, *ly, 0.0))
      dserror("lubrication scaling failed");

    Extractor()->InsertVector(*sy,0,x);
    Extractor()->InsertVector(*ly,1,x);

    Teuchos::RCP<Epetra_Vector> sx = Extractor()->ExtractVector(b,0);
    Teuchos::RCP<Epetra_Vector> lx = Extractor()->ExtractVector(b,1);

    if (sx->ReciprocalMultiply(1.0, *srowsum_, *sx, 0.0))
      dserror("structure scaling failed");
    if (lx->ReciprocalMultiply(1.0, *lrowsum_, *lx, 0.0))
      dserror("lubrication scaling failed");

    Extractor()->InsertVector(*sx,0,b);
    Extractor()->InsertVector(*lx,1,b);

    Teuchos::RCP<Epetra_CrsMatrix> A = mat.Matrix(0,0).EpetraMatrix();
    srowsum_->Reciprocal(*srowsum_);
    scolsum_->Reciprocal(*scolsum_);
    if ( (A->LeftScale(*srowsum_)) or (A->RightScale(*scolsum_)) or
         (mat.Matrix(0,1).EpetraMatrix()->LeftScale(*srowsum_)) or
         (mat.Matrix(1,0).EpetraMatrix()->RightScale(*scolsum_))
       )
      dserror("structure scaling failed");

    A = mat.Matrix(1,1).EpetraMatrix();
    lrowsum_->Reciprocal(*lrowsum_);
    lcolsum_->Reciprocal(*lcolsum_);
    if ( (A->LeftScale(*lrowsum_)) or (A->RightScale(*lcolsum_)) or
         (mat.Matrix(1,0).EpetraMatrix()->LeftScale(*lrowsum_)) or
         (mat.Matrix(0,1).EpetraMatrix()->RightScale(*lcolsum_))
       )
      dserror("lubrication scaling failed");

  }  // if (scaling_infnorm)

}  // UnscaleSolution()


/*----------------------------------------------------------------------*
 | calculate vector norm                                    wirtz 01/16 |
 *----------------------------------------------------------------------*/
double EHL::Monolithic::CalculateVectorNorm(
  const enum INPAR::EHL::VectorNorm norm,
  const Teuchos::RCP<Epetra_Vector> vect
  )
{
  // L1 norm
  // norm = sum_0^i vect[i]
  if (norm == INPAR::EHL::norm_l1)
  {
    double vectnorm;
    vect->Norm1(&vectnorm);
    return vectnorm;
  }
  // L2/Euclidian norm
  // norm = sqrt{sum_0^i vect[i]^2 }
  else if (norm == INPAR::EHL::norm_l2)
  {
    double vectnorm;
    vect->Norm2(&vectnorm);
    return vectnorm;
  }
  // RMS norm
  // norm = sqrt{sum_0^i vect[i]^2 }/ sqrt{length_vect}
  else if (norm == INPAR::EHL::norm_rms)
  {
    double vectnorm;
    vect->Norm2(&vectnorm);
    return vectnorm/sqrt((double) vect->GlobalLength());
  }
  // infinity/maximum norm
  // norm = max( vect[i] )
  else if (norm == INPAR::EHL::norm_inf)
  {
    double vectnorm;
    vect->NormInf(&vectnorm);
    return vectnorm;
  }
  // norm = sum_0^i vect[i]/length_vect
  else if (norm == INPAR::EHL::norm_l1_scaled)
  {
    double vectnorm;
    vect->Norm1(&vectnorm);
    return vectnorm/((double) vect->GlobalLength());
  }
  else
  {
    dserror("Cannot handle vector norm");
    return 0;
  }
}  // CalculateVectorNorm()


/*----------------------------------------------------------------------*
 | set parameters for EHL remaining constant over whole     wirtz 01/16 |
 | simulation                                                           |
 *----------------------------------------------------------------------*/
void EHL::Monolithic::SetDefaultParameters()
{
  // time parameters
  // call the EHL parameter list
  const Teuchos::ParameterList& ldyn
    = DRT::Problem::Instance()->LubricationDynamicParams();

  // get the parameters for the Newton iteration
  itermax_ = ehldyn_.get<int>("ITEMAX");
  itermin_ = ehldyn_.get<int>("ITEMIN");

  // what kind of norm do we wanna test for coupled EHL problem
  normtypeinc_
    = DRT::INPUT::IntegralValue<INPAR::EHL::ConvNorm>(ehldynmono_,"NORM_INC");
  normtyperhs_
    = DRT::INPUT::IntegralValue<INPAR::EHL::ConvNorm>(ehldynmono_,"NORM_RESF");
  // what kind of norm do we wanna test for the single fields
  normtypedisi_
    = DRT::INPUT::IntegralValue<INPAR::STR::ConvNorm>(sdyn_,"NORM_DISP");
  normtypestrrhs_
    = DRT::INPUT::IntegralValue<INPAR::STR::ConvNorm>(sdyn_,"NORM_RESF");
  enum INPAR::STR::VectorNorm striternorm
    = DRT::INPUT::IntegralValue<INPAR::STR::VectorNorm>(sdyn_,"ITERNORM");
  normtypeprei_
    = DRT::INPUT::IntegralValue<INPAR::LUBRICATION::ConvNorm>(ldyn,"NORM_PRE");
  normtypelubricationrhs_
    = DRT::INPUT::IntegralValue<INPAR::LUBRICATION::ConvNorm>(ldyn,"NORM_RESF");
  enum INPAR::LUBRICATION::VectorNorm lubricationiternorm
    = DRT::INPUT::IntegralValue<INPAR::LUBRICATION::VectorNorm>(ldyn,"ITERNORM");
  // in total when do we reach a converged state for complete problem
  combincrhs_
    = DRT::INPUT::IntegralValue<INPAR::EHL::BinaryOp>(ehldynmono_,"NORMCOMBI_RESFINC");

  switch (combincrhs_)
  {
  case INPAR::EHL::bop_and :
  {
    if (Comm().MyPID() == 0)
      std::cout << "Convergence test of EHL:\n res, inc with 'AND'." << std::endl;
    break;
  }
  case INPAR::EHL::bop_or :
  {
    if (Comm().MyPID() == 0)
      std::cout << "Convergence test of EHL:\n res, inc with 'OR'." << std::endl;
    break;
  }
  case INPAR::EHL::bop_coupl_and_singl :
  {
    if (Comm().MyPID() == 0)
      std::cout << "Convergence test of EHL:\n res, inc, str-res, lub-res, dis, pre with 'AND'." << std::endl;
    break;
  }
  case INPAR::EHL::bop_coupl_or_singl :
  {
    if (Comm().MyPID() == 0)
      std::cout << "Convergence test of EHL:\n (res, inc) or (str-res, lub-res, dis, pre)." << std::endl;
    break;
  }
  case INPAR::EHL::bop_and_singl :
  {
    if (Comm().MyPID() == 0)
      std::cout << "Convergence test of EHL:\n str-res, lub-res, dis, pre with 'AND'." << std::endl;
    break;
  }
  case INPAR::EHL::bop_or_singl :
  {
    if (Comm().MyPID() == 0)
      std::cout << "Convergence test of EHL:\n str-res, lub-res, dis, pre with 'OR'." << std::endl;
    break;
  }
  default :
  {
    dserror("Something went terribly wrong with binary operator!");
    break;
  }
  }  // switch (combincrhs_)

  // convert the single field norms to be used within EHL
  // what norm is used for structure
  switch (striternorm)
  {
  case INPAR::STR::norm_l1 :
    iternormstr_ = INPAR::EHL::norm_l1;
    break;
  case INPAR::STR::norm_l2 :
    iternormstr_ = INPAR::EHL::norm_l2;
    break;
  case INPAR::STR::norm_rms :
    iternormstr_ = INPAR::EHL::norm_rms;
    break;
  case INPAR::STR::norm_inf :
    iternormstr_ = INPAR::EHL::norm_inf;
    break;
  case INPAR::STR::norm_vague :
  default :
    dserror("STR norm is not determined");
    break;
  }  // switch (striternorm)

  // what norm is used for lubrication
  switch (lubricationiternorm)
  {
  case INPAR::LUBRICATION::norm_l1 :
    iternormlubrication_ = INPAR::EHL::norm_l1;
    break;
  case INPAR::LUBRICATION::norm_l2 :
    iternormlubrication_ = INPAR::EHL::norm_l2;
    break;
  case INPAR::LUBRICATION::norm_rms :
    iternormlubrication_ = INPAR::EHL::norm_rms;
    break;
  case INPAR::LUBRICATION::norm_inf :
    iternormlubrication_ = INPAR::EHL::norm_inf;
    break;
  case INPAR::LUBRICATION::norm_vague :
  default :
  {
    dserror("LUBRICATION norm is not determined.");
    break;
  }
  }  // switch (lubricationiternorm)

  // if scaled L1-norm is wished to be used
  if ( (iternorm_ == INPAR::EHL::norm_l1_scaled) and
       ( (combincrhs_ == INPAR::EHL::bop_coupl_and_singl) or
         (combincrhs_ == INPAR::EHL::bop_coupl_or_singl)
       )
     )
  {
    iternormstr_ = INPAR::EHL::norm_l1_scaled;
    iternormlubrication_ = INPAR::EHL::norm_l1_scaled;
  }

  // test the EHL-residual and the EHL-increment
  tolinc_ = ehldynmono_.get<double>("TOLINC");
  tolrhs_ = ehldynmono_.get<double>("CONVTOL");

  // get the single field tolerances from this field itselves
  toldisi_ = sdyn_.get<double>("TOLDISP");
  tolstrrhs_ = sdyn_.get<double>("TOLRES");
  tolprei_ = ldyn.get<double>("TOLPRE");
  tollubricationrhs_ = ldyn.get<double>("TOLRES");

  // initialise norms for coupled EHL problem
  normrhs_ = 0.0;
  normrhsiter0_ = 0.0;
  norminc_ = 0.0;
  norminciter0_ = 0.0;

  // initialise norms for single field tests
  normdisi_ = 0.0;
  normstrrhs_ = 0.0;
  normstrrhsiter0_ = 0.0;
  normprei_ = 0.0;
  normlubricationrhs_ = 0.0;
  normlubricationrhsiter0_ = 0.0;

  return;

}  // SetDefaultParameter()

/*----------------------------------------------------------------------*
 | calculate stresses, strains, energies                    wirtz 01/16 |
 *----------------------------------------------------------------------*/
void EHL::Monolithic::PrepareOutput()
{
  //set pressures on structure field for evaluating stresses
  ApplyLubricationCouplingState(lubrication_->LubricationField()->Prenp());
  // prepare output (i.e. calculate stresses, strains, energies)
  StructureField()->PrepareOutput();

  //reset states
  StructureField()->Discretization()->ClearState(true);
}
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | apply pressure state on structure discretization         wirtz 01/16 |
 *----------------------------------------------------------------------*/
void EHL::Monolithic::ApplyLubricationCouplingState(Teuchos::RCP<const Epetra_Vector> pre,
                                              Teuchos::RCP<const Epetra_Vector> pre_res)
{
//  EHL::Algorithm::ApplyLubricationCouplingState(pre,pre_res);
}  // ApplyLubricationCouplingState()


/*----------------------------------------------------------------------*
 | apply structural displacements and velocities on         wirtz 01/16 |
 | lubrication discretization                                           |
 *----------------------------------------------------------------------*/
void EHL::Monolithic::ApplyStructCouplingState(Teuchos::RCP<const Epetra_Vector> disp,
                                              Teuchos::RCP<const Epetra_Vector> vel)
{

  Teuchos::RCP<Epetra_Vector> vec_lub = Teuchos::rcp(new Epetra_Vector(
      *lubrication_->LubricationField()->Discretization()->DofRowMap(1)));
  if (disp != Teuchos::null)
  {
    LINALG::Export(*disp,*vec_lub);
    lubrication_->LubricationField()->Discretization()->SetState(1, "displacement", vec_lub);
  }
  if (vel != Teuchos::null)
  {
    LINALG::Export(*vel,*vec_lub);
    lubrication_->LubricationField()->Discretization()->SetState(1, "velocity", vec_lub);
  }

}  // ApplyStructCouplingState()

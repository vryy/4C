/*----------------------------------------------------------------------*/
/*!
\file tsi_monolithic.cpp

\brief  Basis of all monolithic TSI algorithms that perform a coupling between
        the linear momentum equation and the heat conduction equation

<pre>
Maintainer: Caroline Danowski
            danowski@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>
*/

/*----------------------------------------------------------------------*
 | headers                                                   dano 11/10 |
 *----------------------------------------------------------------------*/
#include "tsi_monolithic.H"
#include "tsi_defines.H"

#include <Teuchos_TimeMonitor.hpp>

// include this header for coupling stiffness terms
#include "../drt_lib/drt_assemblestrategy.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_blocksparsematrix.H"
#include "../linalg/linalg_utils.H"
#include "../drt_inpar/inpar_solver.H"
#include "../linalg/linalg_solver.H"

// contact
#include "../drt_contact/contact_abstract_strategy.H"
#include "../drt_contact/contact_interface.H"
#include "../drt_contact/contact_node.H"
#include "../drt_thermo/thr_contact.H"
#include "../drt_mortar/mortar_manager_base.H"

//! Note: The order of calling the two BaseAlgorithm-constructors is
//! important here! In here control file entries are written. And these entries
//! define the order in which the filters handle the Discretizations, which in
//! turn defines the dof number ordering of the Discretizations.

/*----------------------------------------------------------------------*
 | destructor (public)                                       dano 11/10 |
 *----------------------------------------------------------------------*/
TSI::Monolithic::~Monolithic()
{
}


/*----------------------------------------------------------------------*
 | monolithic                                                dano 11/10 |
 *----------------------------------------------------------------------*/
TSI::Monolithic::Monolithic(
  const Epetra_Comm& comm,
  const Teuchos::ParameterList& sdynparams
  )
: Algorithm(comm),
  solveradapttol_(DRT::INPUT::IntegralValue<int>(sdynparams,"ADAPTCONV")==1),
  solveradaptolbetter_(sdynparams.get<double>("ADAPTCONV_BETTER")),
  printiter_(true),  // ADD INPUT PARAMETER
  printerrfile_(true and errfile_),  // ADD INPUT PARAMETER FOR 'true'
  errfile_(NULL),
  zeros_(Teuchos::null),
  strmethodname_(DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(sdynparams,"DYNAMICTYP")),
  blockrowdofmap_(Teuchos::null),
  systemmatrix_(Teuchos::null),
  sdyn_(sdynparams),
  veln_(Teuchos::null)
{
  // monolithic TSI must know the other discretization
  // build a proxy of the structure discretization for the temperature field
  Teuchos::RCP<DRT::DofSet> structdofset
    = StructureField()->Discretization()->GetDofSetProxy();
  // build a proxy of the temperature discretization for the structure field
  Teuchos::RCP<DRT::DofSet> thermodofset
    = ThermoField()->Discretization()->GetDofSetProxy();

  // check if ThermoField has 2 discretizations, so that coupling is possible
  if (ThermoField()->Discretization()->AddDofSet(structdofset)!=1)
    dserror("unexpected dof sets in thermo field");
  if (StructureField()->Discretization()->AddDofSet(thermodofset)!=1)
    dserror("unexpected dof sets in structure field");

  // access the thermal parameter lists
  const Teuchos::ParameterList& tdyn
    = DRT::Problem::Instance()->ThermalDynamicParams();

  // check time integration algo -> currently only one-step-theta scheme supported
  INPAR::STR::DynamicType structtimealgo
    = DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(sdyn_,"DYNAMICTYP");
  INPAR::THR::DynamicType thermotimealgo
    = DRT::INPUT::IntegralValue<INPAR::THR::DynamicType>(tdyn,"DYNAMICTYP");

  // use for both fields the same time integrator
  if ( ( structtimealgo!=INPAR::STR::dyna_onesteptheta ||
         thermotimealgo!=INPAR::THR::dyna_onesteptheta )
         &&
       ( structtimealgo!=INPAR::STR::dyna_statics ||
         thermotimealgo!=INPAR::THR::dyna_statics )
      )
    dserror("same time integration scheme for STR and THR required for monolithic.");

  // add extra parameters (a kind of work-around)
  Teuchos::RCP<Teuchos::ParameterList> xparams
    = Teuchos::rcp(new Teuchos::ParameterList());
  xparams->set<FILE*>("err file", DRT::Problem::Instance()->ErrorFile()->Handle());
  errfile_ = xparams->get<FILE*>("err file");

  blockrowdofmap_ = Teuchos::rcp(new LINALG::MultiMapExtractor);

  // velocities V_{n+1} at t_{n+1}
  veln_ = LINALG::CreateVector(*(StructureField()->DofRowMap(0)), true);
  veln_->PutScalar(0.0);

  // tsi solver
#ifdef TSIBLOCKMATRIXMERGE
  // create a linear solver
  // get UMFPACK...
  Teuchos::RCP<Teuchos::ParameterList> solverparams = rcp(new Teuchos::ParameterList);
  solverparams->set("solver","umfpack");

  solver_ = rcp(new LINALG::Solver(
                      solverparams,
                      Comm(),
                      DRT::Problem::Instance()->ErrorFile()->Handle()
                      )
                );
#else
  // create a linear solver
  CreateLinearSolver();

#endif

  // structural and thermal contact
  if(StructureField()->ContactManager()!=null)
  {
    cmtman_ = StructureField()->ContactManager();

    // initialize thermal contact manager
    ThermoField()->PrepareThermoContact(StructureField()->ContactManager(),StructureField()->Discretization());

    // get thermal contact manager
    thermcontman_ = ThermoField()->ThermoContactManager();

    // check input
    if (cmtman_->GetStrategy().Friction())
      dserror ("TSI with contact only for frictionless contact so far!");
  }

  // StructureField: check whether we have locsys BCs and create LocSysManager
  // if so after checking
  {
    std::vector<DRT::Condition*> locsysconditions(0);
    (StructureField()->Discretization())->GetCondition("Locsys", locsysconditions);
    if (locsysconditions.size())
    {
      locsysman_ = StructureField()->LocsysManager();
    }
  }

}  // Monolithic()


/*----------------------------------------------------------------------*
 | read restart information for given time step (public)     dano 11/10 |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::ReadRestart(int step)
{
  ThermoField()->ReadRestart(step);
  StructureField()->ReadRestart(step);

  // pass the current coupling varibles to the respective field
  ThermoField()->ApplyStructVariables(StructureField()->Dispnp(),
                                     StructureField()->ExtractVelnp());
  StructureField()->ApplyTemperatures(ThermoField()->Tempnp());

  // second ReadRestart needed due to the coupling variables
  ThermoField()->ReadRestart(step);
  StructureField()->ReadRestart(step);

  SetTimeStep(ThermoField()->GetTime(),step);

  return;
}


/*----------------------------------------------------------------------*
 | prepare time step (public)                                dano 11/10 |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::PrepareTimeStep()
{
  // counter and print header
  // increment time and step counter
  IncrementTimeAndStep();
  PrintHeader();

  // call the predictor
  StructureField()->PrepareTimeStep();
  ThermoField()->PrepareTimeStep();
}


/*----------------------------------------------------------------------*
 | create linear solver                                   wiesner 07/11 |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::CreateLinearSolver()
{
  // get dynamic section of TSI
  const Teuchos::ParameterList& tsidyn = DRT::Problem::Instance()->TSIDynamicParams();
  // get the solver number used for linear TSI solver
  const int linsolvernumber = tsidyn.get<int>("LINEAR_SOLVER");
  // check if the TSI solver has a valid solver number
  if (linsolvernumber == (-1))
    dserror("no linear solver defined for monolithic TSI. Please set LINEAR_SOLVER in TSI DYNAMIC to a valid number!");

  // get parameter list of structural dynamics
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
  // use solver blocks for structure
  // get the solver number used for structural solver
  const int slinsolvernumber = sdyn.get<int>("LINEAR_SOLVER");
  // check if the structural solver has a valid solver number
  if (slinsolvernumber == (-1))
    dserror("no linear solver defined for structural field. Please set LINEAR_SOLVER in STRUCTURAL DYNAMIC to a valid number!");

  // get parameter list of thermal dynamics
  const Teuchos::ParameterList& tdyn = DRT::Problem::Instance()->ThermalDynamicParams();
  // use solver blocks for temperature (thermal field)
  // get the solver number used for thermal solver
  const int tlinsolvernumber = tdyn.get<int>("LINEAR_SOLVER");
  // check if the TSI solver has a valid solver number
  if (tlinsolvernumber == (-1))
    dserror("no linear solver defined for thermal field. Please set LINEAR_SOLVER in THERMAL DYNAMIC to a valid number!");

  // get solver parameter list of linear TSI solver
  const Teuchos::ParameterList& tsisolverparams
    = DRT::Problem::Instance()->SolverParams(linsolvernumber);

  const int solvertype
    = DRT::INPUT::IntegralValue<INPAR::SOLVER::SolverType>(
        tsisolverparams,
        "SOLVER"
        );

  if (solvertype != INPAR::SOLVER::aztec_msr &&
      solvertype != INPAR::SOLVER::belos)
  {
    cout << "!!!!!!!!!!!!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!!!!!!!!!" << endl;
    cout << " Note: the BGS2x2 preconditioner now "                  << endl;
    cout << " uses the structural solver and thermal solver blocks"  << endl;
    cout << " for building the internal inverses"                    << endl;
    cout << " Remove the old BGS PRECONDITIONER BLOCK entries "      << endl;
    cout << " in the dat files!"                                     << endl;
    cout << "!!!!!!!!!!!!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!!!!!!!!!" << endl;
    dserror("aztec solver expected");
  }
  const int azprectype
    = DRT::INPUT::IntegralValue<INPAR::SOLVER::AzPrecType>(
        tsisolverparams,
        "AZPREC"
        );

  // plausibility check
  switch (azprectype)
  {
    case INPAR::SOLVER::azprec_BGS2x2:
      break;
    case INPAR::SOLVER::azprec_BGSnxn:
    case INPAR::SOLVER::azprec_TekoSIMPLE:
    {
#ifdef HAVE_TEKO
      // check if structural solver and thermal solver are Stratimikos based (Teko expects stratimikos)
      int solvertype = DRT::INPUT::IntegralValue<INPAR::SOLVER::SolverType>(DRT::Problem::Instance()->SolverParams(slinsolvernumber), "SOLVER");
      if (solvertype != INPAR::SOLVER::stratimikos_amesos &&
          solvertype != INPAR::SOLVER::stratimikos_aztec  &&
          solvertype != INPAR::SOLVER::stratimikos_belos)
      dserror("Teko expects a STRATIMIKOS solver object in STRUCTURE SOLVER");

      solvertype = DRT::INPUT::IntegralValue<INPAR::SOLVER::SolverType>(DRT::Problem::Instance()->SolverParams(tlinsolvernumber), "SOLVER");
      if (solvertype != INPAR::SOLVER::stratimikos_amesos &&
          solvertype != INPAR::SOLVER::stratimikos_aztec  &&
          solvertype != INPAR::SOLVER::stratimikos_belos)
        dserror("Teko expects a STRATIMIKOS solver object in thermal solver %3d",tlinsolvernumber);
#else
      dserror("Teko preconditioners only available with HAVE_TEKO flag for TRILINOS_DEV (>Q1/2011)");
#endif
    }
    break;
    default:
          dserror("Block Gauss-Seidel BGS2x2 preconditioner expected");
          break;
  }

  solver_ = rcp(new LINALG::Solver(
                         tsisolverparams,
                         // ggfs. explizit Comm von STR wie lungscatra
                         Comm(),
                         DRT::Problem::Instance()->ErrorFile()->Handle()
                         )
                     );

  // use solver blocks for structure and temperature (thermal field)
  const Teuchos::ParameterList& ssolverparams = DRT::Problem::Instance()->SolverParams(slinsolvernumber);
  const Teuchos::ParameterList& tsolverparams = DRT::Problem::Instance()->SolverParams(tlinsolvernumber);

  solver_->PutSolverParamsToSubParams("Inverse1", ssolverparams);
  solver_->PutSolverParamsToSubParams("Inverse2", tsolverparams);

  // prescribe rigid body modes
  StructureField()->Discretization()->ComputeNullSpaceIfNecessary(
                                       solver_->Params().sublist("Inverse1")
                                       );
  ThermoField()->Discretization()->ComputeNullSpaceIfNecessary(
                                    solver_->Params().sublist("Inverse2")
                                    );
}


/*----------------------------------------------------------------------*
 | time loop of the monolithic system                        dano 11/10 |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::TimeLoop()
{
  // time loop
  while (NotFinished())
  {
    // counter and print header
    // predict solution of both field (call the adapter)
    PrepareTimeStep();

    // Newton-Raphson iteration
    NewtonFull();

    // calculate stresses, strains, energies
    PrepareOutput();

    // update all single field solvers
    Update();

    // write output to screen and files
    Output();

#ifdef TSIMONOLITHASOUTPUT
    printf("Ende Timeloop ThermoField()->ExtractTempnp[0] %12.8f\n",(*ThermoField()->ExtractTempnp())[0]);
    printf("Ende Timeloop ThermoField()->ExtractTempn[0] %12.8f\n",(*ThermoField()->ExtractTempn())[0]);

    printf("Ende Timeloop disp %12.8f\n",(*StructureField()->Dispn())[0]);
    cout << "dispn\n" << *(StructureField()->Dispn()) << endl;
#endif // TSIMONOLITHASOUTPUT

  }  // NotFinished
}  // TimeLoop


/*----------------------------------------------------------------------*
 | solution with full Newton-Raphson iteration               dano 10/10 |
 | in tsi_algorithm: NewtonFull()                                       |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::NewtonFull()
{
#ifndef TFSI
  if ( Comm().MyPID()==0 )
    cout << "TSI::Monolithic::NewtonFull()" << endl;
#endif

  // we do a Newton-Raphson iteration here.
  // the specific time integration has set the following
  // --> On #rhs_ is the positive force residuum
  // --> On #systemmatrix_ is the effective dynamic tangent matrix

  // time parameters
  // call the TSI parameter list
  const Teuchos::ParameterList& tsidyn =
    DRT::Problem::Instance()->TSIDynamicParams();
  // Get the parameters for the Newton iteration
  itermax_ = tsidyn.get<int>("ITEMAX");
  itermin_ = tsidyn.get<int>("ITEMIN");
  normtypeinc_
    = DRT::INPUT::IntegralValue<INPAR::TSI::ConvNorm>(tsidyn,"NORM_INC");
  normtypefres_
    = DRT::INPUT::IntegralValue<INPAR::TSI::ConvNorm>(tsidyn,"NORM_RESF");
  combincfres_
    = DRT::INPUT::IntegralValue<INPAR::TSI::BinaryOp>(tsidyn,"NORMCOMBI_RESFINC");

  // FIRST STEP: to test the residual and the increments use the same tolerance
  tolinc_ =  tsidyn.get<double>("CONVTOL");
  tolfres_ = tsidyn.get<double>("CONVTOL");

  // initialise equilibrium loop
  iter_ = 1;
  normrhs_ = 0.0;
  normrhsiter0_ = 0.0;
  normstrrhs_ = 0.0;
  normthrrhs_ = 0.0;
  norminc_ = 0.0;
  // get length of the structural and thermal vector
  ns_ = (*(StructureField()->RHS())).GlobalLength();
  nt_ = (*(ThermoField()->RHS())).GlobalLength();
  // get length of all TSI dofs
  ntsi_ = ns_ + nt_;

  Epetra_Time timerthermo(Comm());
  timerthermo.ResetStartTime();

  // incremental solution vector with length of all TSI dofs
  iterinc_ = LINALG::CreateVector(*DofRowMap(), true);
  iterinc_->PutScalar(0.0);
  // a zero vector of full length
  zeros_ = LINALG::CreateVector(*DofRowMap(), true);
  zeros_->PutScalar(0.0);

  //---------------------------------------------- iteration loop

  // equilibrium iteration loop (loop over k)
  while ( ( (not Converged()) and (iter_ <= itermax_) ) or (iter_ <= itermin_) )
  {
    // compute residual forces #rhs_ and tangent #tang_
    // whose components are globally oriented
    // build linear system stiffness matrix and rhs/force residual for each
    // field, here e.g. for structure field: field want the iteration increment
    // 1.) Update(iterinc_),
    // 2.) EvaluateForceStiffResidual(),
    // 3.) PrepareSystemForNewtonSolve()
    Evaluate(iterinc_);

    // create the linear system
    // \f$J(x_i) \Delta x_i = - R(x_i)\f$
    // create the systemmatrix
    SetupSystemMatrix();

    // check whether we have a sanely filled tangent matrix
    if (not systemmatrix_->Filled())
    {
      dserror("Effective tangent matrix must be filled here");
    }

    // create full monolithic rhs vector
    // make negative residual not necessary: rhs_ is yet TSI negative
    SetupRHS();

    // (Newton-ready) residual with blanked Dirichlet DOFs (see adapter_timint!)
    // is done in PrepareSystemForNewtonSolve() within Evaluate(iterinc_)
    LinearSolve();

    // recover LM in the case of contact
    RecoverStructThermLM();

    // reset solver tolerance
    solver_->ResetTolerance();

    // build residual force norm
    // for now use for simplicity only L2/Euclidian norm
    rhs_->Norm2(&normrhs_);
    if (iter_ == 1)
      normrhsiter0_ = normrhs_;
    StructureField()->RHS()->Norm2(&normstrrhs_);
    ThermoField()->RHS()->Norm2(&normthrrhs_);

    // build residual increment norm
    iterinc_->Norm2(&norminc_);

    // print stuff
    PrintNewtonIter();

    // increment equilibrium loop index
    iter_ += 1;

  }  // end equilibrium loop

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

}  // NewtonFull()


/*----------------------------------------------------------------------*
 | evaluate the single fields                                dano 11/10 |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::Evaluate(Teuchos::RCP<Epetra_Vector> x)
{
#ifndef TFSI
  if (Comm().MyPID() == 0)
    cout << "\n TSI::Monolithic::Evaluate()" << endl;
#endif
  TEUCHOS_FUNC_TIME_MONITOR("TSI::Monolithic::Evaluate");

  // displacement and temperature incremental vector
  Teuchos::RCP<Epetra_Vector> sx;
  Teuchos::RCP<Epetra_Vector> tx;

  // if an increment vector exists
  if (x!=Teuchos::null)
  {
    // extract displacement sx and temperature tx incremental vector of global
    // unknown incremental vector x
    ExtractFieldVectors(x,sx,tx);

#ifdef TSIMONOLITHASOUTPUT
    cout << "Recent thermal increment DT_n+1^i\n" << *(tx) << endl;
    cout << "Recent structural increment Dd_n+1^i\n" << *(sx) << endl;

    cout << "Until here only old solution of Newton step. No update applied\n" << *(ThermoField()->Tempnp()) << endl;
#endif // TSIMONOLITHASOUTPUT
  }
  // else(x=Teuchos::null): initialize the system

#ifdef TSIMONOLITHASOUTPUT
  cout << "Tempnp vor UpdateNewton\n" << *(ThermoField()->Tempnp()) << endl;
  printf("Tempnp vor UpdateNewton ThermoField()->ExtractTempnp[0] %12.8f\n",(*ThermoField()->ExtractTempnp())[0]);
#endif // TSIMONOLITHASOUTPUT

  // Newton update of the thermo field
  // update temperature before passed to the structural field
  //   UpdateIterIncrementally(tx),
  ThermoField()->UpdateNewton(tx);

#ifdef TSIMONOLITHASOUTPUT
  cout << "Tempnp nach UpdateNewton\n" << *(ThermoField()->Tempnp()) << endl;
  printf("Tempnp nach UpdateNewton ThermoField()->ExtractTempnp[0] %12.8f\n",(*ThermoField()->ExtractTempnp())[0]);
#endif // TSIMONOLITHASOUTPUT

  // call all elements and assemble rhs and matrices
  /// structural field

  // structure Evaluate (builds tangent, residual and applies DBC)
  Epetra_Time timerstructure(Comm());

#ifndef MonTSIwithoutTHR
  // apply current temperature to structure
  StructureField()->ApplyTemperatures(ThermoField()->Tempnp());
#endif

#ifdef TSIPARALLEL
  cout << Comm().MyPID() << " nach ApplyTemp!!" << endl;
#endif // TSIPARALLEL

#ifdef TSIMONOLITHASOUTPUT
//    Teuchos::RCP<Epetra_Vector> tempera = rcp(new Epetra_Vector(ThermoField()->Tempn()->Map(),true));
//    if (ThermoField()->Tempnp() != Teuchos::null)
//      tempera->Update(1.0, *ThermoField()->Tempnp(), 0.0);
//    StructureField()->ApplyTemperatures(tempera);
//    StructureField()->ApplyTemperatures(ThermoField()->Tempn());
#endif // TSIMONOLITHASOUTPUT

  // Monolithic TSI accesses the linearised structure problem:
  //   UpdaterIterIncrementally(sx),
  //   EvaluateForceStiffResidual()
  //   PrepareSystemForNewtonSolve()
  //     blank residual DOFs that are on Dirichlet BC
  //     in the case of local systems we have to rotate forth and back
  //     ApplyDirichlettoSystem is called as well with roated stiff_!
  StructureField()->Evaluate(sx);
#ifndef TFSI
  cout << "  structure time for calling Evaluate: " << timerstructure.ElapsedTime() << "\n";
#endif

#ifdef TSIMONOLITHASOUTPUT
  cout << "STR fres_" << *StructureField()->RHS() << endl;
#endif // TSIMONOLITHASOUTPUT

  /// thermal field

  // thermo Evaluate
  // (builds tangent, residual and applies DBC and recent coupling values)
  Epetra_Time timerthermo(Comm());

  // apply current displacements and velocities to the thermo field
  if (strmethodname_==INPAR::STR::dyna_statics)
  {
    // calculate velocity V_n+1^k = (D_n+1^k-D_n)/Dt()
    veln_ = CalcVelocity(StructureField()->Dispnp());
  }
  else
  {
    veln_ = StructureField()->ExtractVelnp();
  }
#ifndef MonTSIwithoutSTR
  // pass the structural values to the thermo field
  ThermoField()->ApplyStructVariables(StructureField()->Dispnp(),veln_);
#endif

#ifdef TSIMONOLITHASOUTPUT
  cout << "d_n+1 inserted in THR field\n" << *(StructureField()->Dispnp()) << endl;
  cout << "v_n+1\n" << *veln_ << endl;
#endif // TSIMONOLITHASOUTPUT

  // monolithic TSI accesses the linearised thermo problem
  //   EvaluateRhsTangResidual() and
  //   PrepareSystemForNewtonSolve()
  ThermoField()->Evaluate();
#ifndef TFSI
  cout << "  thermo time for calling Evaluate: " << timerthermo.ElapsedTime() << "\n";
#endif

}  // Evaluate()


/*----------------------------------------------------------------------*
 | extract field vectors for calling Evaluate() of the       dano 11/10 |
 | single fields                                                        |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::ExtractFieldVectors(
  Teuchos::RCP<Epetra_Vector> x,
  Teuchos::RCP<Epetra_Vector>& sx,
  Teuchos::RCP<Epetra_Vector>& tx
  )
{
  TEUCHOS_FUNC_TIME_MONITOR("TSI::Monolithic::ExtractFieldVectors");

  // process structure unknowns of the first field
  sx = Extractor()->ExtractVector(x,0);

  // process thermo unknowns of the second field
  tx = Extractor()->ExtractVector(x,1);
}


/*----------------------------------------------------------------------*
 | full monolithic dof row map                               dano 05/12 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> TSI::Monolithic::DofRowMap() const
{
  return Extractor()->FullMap();
}


/*----------------------------------------------------------------------*
 | setup system (called in tsi_dyn)                          dano 11/10 |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::SetupSystem()
{
#ifndef TFSI
  if ( Comm().MyPID()==0 )
    cout << " TSI::Monolithic::SetupSystem()" << endl;
#endif

  // create combined map
  std::vector<Teuchos::RCP<const Epetra_Map> > vecSpaces;

#ifdef TSIPARALLEL
  cout << Comm().MyPID() << " :PID" << endl;
  cout << "structure dofmap" << endl;
  cout << *StructureField()->DofRowMap(0) << endl;
  cout << "thermo dofmap" << endl;
  cout << *StructureField()->DofRowMap(1) << endl;
#endif // TSIPARALLEL

  // use its own DofRowMap, that is the 0th map of the discretization
  vecSpaces.push_back(StructureField()->DofRowMap(0));
  vecSpaces.push_back(ThermoField()->DofRowMap(0));

  if (vecSpaces[0]->NumGlobalElements()==0)
    dserror("No structure equation. Panic.");
  if (vecSpaces[1]->NumGlobalElements()==0)
    dserror("No temperature equation. Panic.");

  SetDofRowMaps(vecSpaces);

}  // SetupSystem()


/*----------------------------------------------------------------------*
 | put the single maps to one full TSI map together          dano 11/10 |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::SetDofRowMaps(
  const std::vector<Teuchos::RCP<const Epetra_Map> >& maps
  )
{
  Teuchos::RCP<Epetra_Map> fullmap
    = LINALG::MultiMapExtractor::MergeMaps(maps);

  // full TSI-blockmap
  Extractor()->Setup(*fullmap,maps);
}


/*----------------------------------------------------------------------*
 | setup system matrix of TSI                                dano 11/10 |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::SetupSystemMatrix()
{
#ifndef TFSI
  if ( Comm().MyPID()==0 )
    cout << " TSI::Monolithic::SetupSystemMatrix()" << endl;
#endif
  TEUCHOS_FUNC_TIME_MONITOR("TSI::Monolithic::SetupSystemMatrix");

  /*----------------------------------------------------------------------*/
  // initialize TSI-systemmatrix_
  systemmatrix_
    = rcp(
        new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
              *Extractor(),
              *Extractor(),
              81,
              false,
              true
              )
        );

  /*----------------------------------------------------------------------*/
  // pure structural part k_ss (3nx3n)

  // build pure structural block k_ss
  // build block matrix
  // The maps of the block matrix have to match the maps of the blocks we
  // insert here. Extract Jacobian matrices and put them into composite system
  // matrix W
  Teuchos::RCP<LINALG::SparseMatrix> k_ss = StructureField()->SystemMatrix();

  // build block matrix
  // The maps of the block matrix have to match the maps of the blocks we
  // insert here.

  // uncomplete because the fluid interface can have more connections than the
  // structural one. (Tet elements in fluid can cause this.) We should do
  // this just once...
  k_ss->UnComplete();

  // assign structure part to the TSI matrix
  systemmatrix_->Assign(0,0,View,*k_ss);

  /*----------------------------------------------------------------------*/
  // structural part k_st (3nxn)
  // build mechanical-thermal block

  // create empty matrix
  Teuchos::RCP<LINALG::SparseMatrix> k_st = Teuchos::null;
  k_st = Teuchos::rcp(
           new LINALG::SparseMatrix(
                 *(StructureField()->Discretization()->DofRowMap(0)),
                 81,
                 true,
                 true
                 )
           );

  // call the element and calculate the matrix block
#ifndef MonTSIwithoutTHR
  ApplyStrCouplMatrix(k_st);

#endif // MonTSIwithoutTHR

  // modify towards contact
  ApplyStructContact(k_st);

  // Uncomplete mechanical-thermal matrix to be able to deal with slightly
  // defective interface meshes.
  k_st->UnComplete();

  // assign thermo part to the TSI matrix
  systemmatrix_->Assign(0,1,View,*(k_st));

  /*----------------------------------------------------------------------*/
  // pure thermo part k_tt (nxn)

  // build pure thermal block k_tt
  // build block matrix
  // The maps of the block matrix have to match the maps of the blocks we
  // insert here. Extract Jacobian matrices and put them into composite system
  // matrix systemmatrix_
  Teuchos::RCP<LINALG::SparseMatrix> k_tt = ThermoField()->SystemMatrix();

  // Uncomplete thermo matrix to be able to deal with slightly defective
  // interface meshes.
  k_tt->UnComplete();

  // assign thermo part to the TSI matrix
  systemmatrix_->Assign(1,1,View,*(k_tt));

  /*----------------------------------------------------------------------*/
  // thermo part k_ts (nx3n)
  // build thermal-mechanical block

  // create empty matrix
  Teuchos::RCP<LINALG::SparseMatrix> k_ts = Teuchos::null;
  k_ts = Teuchos::rcp(
           new LINALG::SparseMatrix(
                 *(ThermoField()->Discretization()->DofRowMap(0)),
                 81,
                 true,
                 true
                 )
           );

  // call the element and calculate the matrix block
#if !defined(MonTSIwithoutSTR) && !defined(COUPLEINITTEMPERATURE)
  ApplyThrCouplMatrix(k_ts);
#endif

  // modify towards contact
  ApplyThermContact(k_ts);

  // Uncomplete thermo matrix to be able to deal with slightly defective
  // interface meshes.
  k_ts->UnComplete();

  // assign thermo part to the TSI matrix
  systemmatrix_->Assign(1,0,View,*(k_ts));

  /*----------------------------------------------------------------------*/
  // done. make sure all blocks are filled.
  systemmatrix_->Complete();
}  // SetupSystemMatrix


/*----------------------------------------------------------------------*
 | setup RHS (like fsimon)                                   dano 11/10 |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::SetupRHS()
{
#ifndef TFSI
  if ( Comm().MyPID()==0 )
    cout << " TSI::Monolithic::SetupRHS()" << endl;
#endif
  TEUCHOS_FUNC_TIME_MONITOR("TSI::Monolithic::SetupRHS");

  // create full monolithic rhs vector
  rhs_ = rcp(new Epetra_Vector(*DofRowMap(), true));

  // fill the TSI rhs vector rhs_ with the single field rhss
  SetupVector(
    *rhs_,
    StructureField()->RHS(),
    ThermoField()->RHS()
    );

}  // SetupRHS()


/*----------------------------------------------------------------------*
 | Solve linear TSI system                                   dano 04/11 |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::LinearSolve()
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

#ifdef TSIBLOCKMATRIXMERGE
  // merge blockmatrix to SparseMatrix and solve
  Teuchos::RCP<LINALG::SparseMatrix> sparse = systemmatrix_->Merge();

  // apply Dirichlet BCs to system of equations
  iterinc_->PutScalar(0.0);  // Useful? depends on solver and more
  LINALG::ApplyDirichlettoSystem(
    sparse,
    iterinc_,
    rhs_,
    Teuchos::null,
    zeros_,
    *CombinedDBCMap()
    );
#ifndef TFSI
  if ( Comm().MyPID()==0 ) { cout << " DBC applied to TSI system" << endl; }
#endif

  // standard solver call
  solver_->Solve(
             sparse->EpetraOperator(),
             iterinc_,
             rhs_,
             true,
             iter_==1
             );

#ifndef TFSI
  if ( Comm().MyPID()==0 ) { cout << " Solved" << endl; }
#endif

#else // use bgs2x2_operator

  // apply Dirichlet BCs to system of equations
  iterinc_->PutScalar(0.0);  // Useful? depends on solver and more

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

#ifndef TFSI
  if ( Comm().MyPID()==0 )
  { cout << " DBC applied to TSI system on proc" << Comm().MyPID() << endl; }
#endif

  solver_->Solve(
             systemmatrix_->EpetraOperator(),
             iterinc_,
             rhs_,
             true,
             iter_==1
             );

#ifndef TFSI
  if ( Comm().MyPID()==0 ) { cout << " Solved" << endl; }
#endif

#endif  // TSIBLOCKMATRIXMERGE

}  // LinearSolve()


/*----------------------------------------------------------------------*
 | initial guess of the displacements/temperatures           dano 11/10 |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::InitialGuess(Teuchos::RCP<Epetra_Vector> ig)
{
  TEUCHOS_FUNC_TIME_MONITOR("TSI::Monolithic::InitialGuess");

  // InitalGuess() is called of the single fields and results are put in TSI
  // increment vector ig
  SetupVector(
    *ig,
    // returns residual displacements \f$\Delta D_{n+1}^{<k>}\f$ - disi_
    StructureField()->InitialGuess(),
    // returns residual temperatures or iterative thermal increment - tempi_
    ThermoField()->InitialGuess()
    );
} // InitialGuess()


/*----------------------------------------------------------------------*
 | setup vector of the structure and thermo field            dano 11/10 |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::SetupVector(
  Epetra_Vector &f,
  Teuchos::RCP<const Epetra_Vector> sv,
  Teuchos::RCP<const Epetra_Vector> tv
  )
{
  // extract dofs of the two fields
  // and put the structural/thermal field vector into the global vector f
  // noticing the block number
  Extractor()->InsertVector(*sv,0,f);
  Extractor()->InsertVector(*tv,1,f);
}  // SetupVector()


/*----------------------------------------------------------------------*
 | check convergence of Newton iteration (public)            dano 11/10 |
 *----------------------------------------------------------------------*/
bool TSI::Monolithic::Converged()
{
  // check for single norms
  bool convinc = false;
  bool convfres = false;

  // residual increments
  switch (normtypeinc_)
  {
    case INPAR::TSI::convnorm_abs:
      convinc = norminc_ < tolinc_;
      break;
    default:
      dserror("Cannot check for convergence of increments!");
  }

  // structural and thermal residual forces
  switch (normtypefres_)
  {
  case INPAR::TSI::convnorm_abs:
    convfres = (normrhs_/ntsi_) < tolfres_;
    break;
  case INPAR::TSI::convnorm_rel:
    convfres = ( ((normstrrhs_/ns_) < tolfres_) and ((normthrrhs_/nt_) < tolfres_) );
    break;
  case INPAR::TSI::convnorm_reliter0:
    // when Converged() is called the first time, normrhsiter0_ is not filled
    // with the norm includes its initial value, i.e. 0.0 used in the denominator
    // in convfres leading to infinity
    // convfres := normrhs_/normrhsiter0_ < tolfres_
    convfres = ( normrhs_ < (tolfres_ * normrhsiter0_ ));
    break;
  case INPAR::TSI::convnorm_mix:
    convfres = ( (normstrrhs_ < tolfres_) and (normthrrhs_ < tolfres_) );
    break;
  default:
    dserror("Cannot check for convergence of residual forces!");
  }

  // combine temperature-like and force-like residuals
  bool conv = false;
  if (combincfres_==INPAR::TSI::bop_and)
    conv = convinc and convfres;
  else
    dserror("Something went terribly wrong with binary operator!");

  // return things
  return conv;

}  // Converged()


/*----------------------------------------------------------------------*
 | print Newton-Raphson iteration to screen and error file   dano 11/10 |
 | originally by lw 12/07, tk 01/08                                     |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::PrintNewtonIter()
{
  // print to standard out
  // replace myrank_ here general by Comm().MyPID()
  if ( (Comm().MyPID()==0) and PrintScreenEvry() and (Step()%PrintScreenEvry()==0) and printiter_ )
  {
    if (iter_== 1)
      PrintNewtonIterHeader(stdout);
    PrintNewtonIterText(stdout);
  }

  // print to error file
  if ( printerrfile_ and printiter_ )
  {
    if (iter_== 1)
      PrintNewtonIterHeader(errfile_);
    PrintNewtonIterText(errfile_);
  }

  // see you
  return;
}  // PrintNewtonIter()


/*----------------------------------------------------------------------*
 | print Newton-Raphson iteration to screen and error file   dano 11/10 |
 | originally by lw 12/07, tk 01/08                                     |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::PrintNewtonIterHeader(FILE* ofile)
{
  // open outstringstream
  std::ostringstream oss;

  // enter converged state etc
  oss << std::setw(6)<< "numiter";

  // different style due relative or absolute error checking
  // displacement
  switch ( normtypefres_ )
  {
  case INPAR::TSI::convnorm_abs :
    oss <<std::setw(18)<< "abs-res-norm";
    break;
  case INPAR::TSI::convnorm_rel :
    oss <<std::setw(18)<< "rel-res-norm";
    break;
  case INPAR::TSI::convnorm_reliter0 :
    oss <<std::setw(18)<< "reliter0-res-norm";
    break;
  case INPAR::TSI::convnorm_mix :
    oss <<std::setw(18)<< "mix-res-norm";
    break;
  default:
    dserror("You should not turn up here.");
  }

  switch ( normtypeinc_ )
  {
  case INPAR::TSI::convnorm_abs :
    oss <<std::setw(18)<< "abs-inc-norm";
    break;
  default:
    dserror("You should not turn up here.");
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
 | print Newton-Raphson iteration to screen                  dano 11/10 |
 | originally by lw 12/07, tk 01/08                                     |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::PrintNewtonIterText(FILE* ofile)
{
  // open outstringstream
  std::ostringstream oss;

  // enter converged state etc
  oss << std::setw(7)<< iter_;

  // different style due relative or absolute error checking
  // displacement
  switch ( normtypefres_ )
  {
  case INPAR::TSI::convnorm_abs :
    oss << std::setw(18) << std::setprecision(5) << std::scientific << normrhs_/ntsi_;
    break;
  case INPAR::TSI::convnorm_rel :
    oss << std::setw(18) << std::setprecision(5) << std::scientific << max( (normstrrhs_/ns_), (normthrrhs_/nt_) );
    break;
  case INPAR::TSI::convnorm_reliter0 :
    oss << std::setw(18) << std::setprecision(5) << std::scientific << (normrhs_/normrhsiter0_);
    break;
  case INPAR::TSI::convnorm_mix :
    oss << std::setw(18) << std::setprecision(5) << std::scientific << max( (normstrrhs_), (normthrrhs_) );
    break;
  default:
    dserror("You should not turn up here.");
  }

  switch ( normtypeinc_ )
  {
  case INPAR::TSI::convnorm_abs :
    oss << std::setw(18) << std::setprecision(5) << std::scientific << norminc_;
    break;
  default:
    dserror("You should not turn up here.");
  }

  // TODO: 26.11.10 double declaration for the timer (first in Evaluate for thermal field)
  Epetra_Time timerthermo(Comm());
  // add solution time
  oss << std::setw(14) << std::setprecision(2) << std::scientific << timerthermo.ElapsedTime();

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
 | print statistics of converged NRI                         dano 11/10 |
 | orignially by bborn 08/09                                            |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::PrintNewtonConv()
{
  // somebody did the door
  return;
}


/*----------------------------------------------------------------------*
 |  evaluate mechanical-thermal system matrix at state       dano 03/11 |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::ApplyStrCouplMatrix(
  Teuchos::RCP<LINALG::SparseMatrix> k_st  //!< off-diagonal tangent matrix term
  )
{
#ifndef TFSI
  if ( Comm().MyPID()==0 )
    cout << " TSI::Monolithic::ApplyStrCouplMatrix()" << endl;
#endif

  // create the parameters for the discretization
  Teuchos::ParameterList sparams;
  const std::string action = "calc_struct_stifftemp";
  sparams.set("action", action);
  // other parameters that might be needed by the elements
  sparams.set("delta time", Dt());
  sparams.set("total time", Time());
//  cout << "STR Parameterliste\n " <<  sparams << endl;
  StructureField()->Discretization()->ClearState();
  StructureField()->Discretization()->SetState(0,"displacement",StructureField()->Dispnp());

  // build specific assemble strategy for mechanical-thermal system matrix
  // from the point of view of StructureField:
  // structdofset = 0, thermdofset = 1
  DRT::AssembleStrategy structuralstrategy(
                          0,  // structdofset for row
                          1,  // thermdofset for column
                          k_st,  // build mechanical-thermal matrix
                          Teuchos::null,  // no other matrix or vectors
                          Teuchos::null,
                          Teuchos::null,
                          Teuchos::null
                          );

  // evaluate the mechancial-thermal system matrix on the structural element
  StructureField()->Discretization()->Evaluate( sparams, structuralstrategy );
  StructureField()->Discretization()->ClearState();

  // for consistent linearisation scale k_st with time factor
  // major switch to different time integrators
  switch (strmethodname_)
  {
   case  INPAR::STR::dyna_statics :
   {
     // continue
     break;
   }
   case  INPAR::STR::dyna_onesteptheta :
   {
     double theta_ = sdyn_.sublist("ONESTEPTHETA").get<double>("THETA");
     // K_Teffdyn(T_n+1^i) = theta * k_st
     k_st->Scale(theta_);
     break;
   }
   // TODO: time factor for genalpha
   case  INPAR::STR::dyna_genalpha :
   {
     double alphaf_ = sdyn_.sublist("GENALPHA").get<double>("ALPHA_F");
     // K_Teffdyn(T_n+1) = (1-alphaf_) . kst
     // Lin(dT_n+1-alphaf_/ dT_n+1) = (1-alphaf_)
     k_st->Scale(1.0 - alphaf_);
   }
   default :
   {
     dserror("Don't know what to do...");
     break;
   }
  }  // end of switch(strmethodname_)

}  // ApplyStrCouplMatrix()


/*----------------------------------------------------------------------*
 |  evaluate thermal-mechanical system matrix at state       dano 03/11 |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::ApplyThrCouplMatrix(
  Teuchos::RCP<LINALG::SparseMatrix> k_ts  //!< off-diagonal tangent matrix term
  )
{
#ifndef TFSI
  if ( Comm().MyPID()==0 )
    cout << " TSI::Monolithic::ApplyThrCouplMatrix()" << endl;
#endif

  // create the parameters for the discretization
  Teuchos::ParameterList tparams;
  // action for elements
  const std::string action = "calc_thermo_coupltang";
  tparams.set("action", action);
  // other parameters that might be needed by the elements
  tparams.set("delta time", Dt());
  tparams.set("total time", Time());
  // create specific time integrator
  const Teuchos::ParameterList& tdyn
    = DRT::Problem::Instance()->ThermalDynamicParams();
  tparams.set<int>(
    "time integrator",
    DRT::INPUT::IntegralValue<INPAR::THR::DynamicType>(tdyn,"DYNAMICTYP")
    );
  switch (DRT::INPUT::IntegralValue<INPAR::THR::DynamicType>(tdyn, "DYNAMICTYP"))
  {
    // static analysis
    case INPAR::THR::dyna_statics :
    {
      break;
    }
    // dynamic analysis
    case INPAR::THR::dyna_onesteptheta :
    {
      double theta_ = tdyn.sublist("ONESTEPTHETA").get<double>("THETA");
      tparams.set("theta",theta_);
      break;
    }
    case INPAR::THR::dyna_genalpha :
    {
      dserror("Genalpha not yet implemented");
      break;
    }
    case INPAR::THR::dyna_undefined :
    default :
    {
      dserror("Don't know what to do...");
      break;
    }
  }

  ThermoField()->Discretization()->ClearState();
  // set the variables that are needed by the elements
  ThermoField()->Discretization()->SetState(0,"temperature",ThermoField()->Tempnp());
  ThermoField()->Discretization()->SetState(1,"displacement",StructureField()->Dispnp());
  ThermoField()->Discretization()->SetState(1,"velocity",veln_);

  // build specific assemble strategy for the thermal-mechanical system matrix
  // from the point of view of ThermoField:
  // thermdofset = 0, structdofset = 1
  DRT::AssembleStrategy thermostrategy(
                          0,  // thermdofset for row
                          1,  // structdofset for column
                          k_ts,  // thermal-mechancial matrix
                          Teuchos::null,  // no other matrix or vectors
                          Teuchos::null,
                          Teuchos::null,
                          Teuchos::null
                          );
  // evaluate the thermal-mechancial system matrix on the thermal element
  ThermoField()->Discretization()->Evaluate(tparams,thermostrategy);
  ThermoField()->Discretization()->ClearState();

  // consider linearisation of velocities due to displacements
  // major switch to different time integrators
  switch (strmethodname_)
  {
    case  INPAR::STR::dyna_statics :
    {
      // Lin (v_n+1) . \Delta d_n+1 = 1/dt, cf. Diss N. Karajan (2009) for quasistatic approach
      const double fac = 1.0/Dt();
      k_ts->Scale(fac);
      break;
    }
    case  INPAR::STR::dyna_onesteptheta :
    {
      double theta = sdyn_.sublist("ONESTEPTHETA").get<double>("THETA");
      const double fac = 1.0 / ( theta * Dt() );
      k_ts->Scale(fac);
      break;
    }
    case  INPAR::STR::dyna_genalpha :
    {
      double beta = sdyn_.sublist("GENALPHA").get<double>("BETA");
      double gamma = sdyn_.sublist("GENALPHA").get<double>("GAMMA");
      // Lin (v_n+1) . \Delta d_n+1 = (gamma) / (beta . dt)
      const double fac =  gamma / ( beta * Dt() );
      k_ts->Scale(fac);
    }
    default :
    {
      dserror("Don't know what to do...");
      break;
    }
  }  // end of switch(strmethodname_)

}  // ApplyThrCouplMatrix()


/*----------------------------------------------------------------------*
 |  map containing the dofs with Dirichlet BC                dano 03/11 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> TSI::Monolithic::CombinedDBCMap()
{
  const Teuchos::RCP<const Epetra_Map > scondmap = StructureField()->GetDBCMapExtractor()->CondMap();
  const Teuchos::RCP<const Epetra_Map > tcondmap = ThermoField()->GetDBCMapExtractor()->CondMap();
  Teuchos::RCP<Epetra_Map> condmap = LINALG::MergeMap(scondmap, tcondmap, false);
  return condmap;
} // CombinedDBCMap()


/*----------------------------------------------------------------------*
 |  apply contact to off diagonal block (k_st)                mgit 05/11 |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::ApplyStructContact(Teuchos::RCP<LINALG::SparseMatrix>& k_st)
{
  // only in the case of contact
  if (cmtman_==Teuchos::null)
    return;

  // contact strategy
  MORTAR::StrategyBase& strategy = cmtman_->GetStrategy();
  CONTACT::CoAbstractStrategy& cstrategy = static_cast<CONTACT::CoAbstractStrategy&>(strategy);

  // check if contact contributions are present,
  // if not we can skip this routine to speed things up
  if (!cstrategy.IsInContact() && !cstrategy.WasInContact() && !cstrategy.WasInContactLastTimeStep())
    return;

  //**********************************************************************
  // maps/matrices form structural and thermal field
  //**********************************************************************
  // necessary maps from structural problem
  RCP<Epetra_Map> sdofs,adofs,idofs,mdofs,amdofs,ndofs,smdofs;
  const Epetra_Map* structprobrowmap = StructureField()->Discretization()->DofRowMap(0);
  sdofs = cstrategy.SlaveRowDofs();
  adofs = cstrategy.ActiveRowDofs();
  mdofs = cstrategy.MasterRowDofs();
  smdofs = LINALG::MergeMap(sdofs,mdofs,false);
  ndofs = LINALG::SplitMap(*structprobrowmap,*smdofs);
  idofs =  LINALG::SplitMap(*sdofs,*adofs);

  // necessary matrices from structural problem
  RCP<LINALG::SparseMatrix> dmatrix = cstrategy.DMatrix();
  RCP<LINALG::SparseMatrix> mmatrix = cstrategy.MMatrix();

  // necessary maps from thermal problem
  RCP<Epetra_Map> thermoprobrowmap = rcp(new Epetra_Map(*(ThermoField()->Discretization()->DofRowMap(0))));

  // abbreviations for active set
  int aset = adofs->NumGlobalElements();

  //**********************************************************************
  // splitting of matrix k_st
  //**********************************************************************
  // complete (needed for splitmap)
  k_st->Complete(*thermoprobrowmap,*structprobrowmap);

  // matrix to split
  RCP<LINALG::SparseMatrix> k_struct_temp = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(k_st);

  RCP<Epetra_Map> tmp;
  RCP<LINALG::SparseMatrix> ksmt,knt,kst,kmt,kat,kit,tmp1,tmp2;

  // first split: k_struct_temp -> ksmt, knt
  LINALG::SplitMatrix2x2(k_struct_temp,smdofs,ndofs,thermoprobrowmap,tmp,ksmt,tmp1,knt,tmp2);

  // second split: ksmt -> kst, kmt
  LINALG::SplitMatrix2x2(ksmt,sdofs,mdofs,thermoprobrowmap,tmp,kst,tmp1,kmt,tmp2);

  // third split: kst -> kat,kit
  LINALG::SplitMatrix2x2(kst,adofs,idofs,thermoprobrowmap,tmp,kat,tmp1,kit,tmp2);

  /**********************************************************************/
  /* evaluation of the inverse of D, active part of M                   */
  /**********************************************************************/
   RCP<LINALG::SparseMatrix> invd = rcp(new LINALG::SparseMatrix(*dmatrix));
   RCP<Epetra_Vector> diag = LINALG::CreateVector(*sdofs,true);
   int err = 0;

   // extract diagonal of invd into diag
   invd->ExtractDiagonalCopy(*diag);

   // set zero diagonal values to dummy 1.0
   for (int i=0;i<diag->MyLength();++i)
     if ((*diag)[i]==0.0) (*diag)[i]=1.0;

   // scalar inversion of diagonal values
   err = diag->Reciprocal(*diag);
   if (err>0) dserror("ERROR: Reciprocal: Zero diagonal entry!");

   // re-insert inverted diagonal into invd
   err = invd->ReplaceDiagonalValues(*diag);

   // store some stuff for condensation of LM
   kst_ = kst;
   invd_ = invd;

   // active part of invd
   RCP<LINALG::SparseMatrix> invda,tempmtx1,tempmtx2,tempmtx3;
   LINALG::SplitMatrix2x2(invd,adofs,idofs,adofs,idofs,invda,tempmtx1,tempmtx2,tempmtx3);

   // active part of mmatrix
   RCP<Epetra_Map> tempmap;
   RCP<LINALG::SparseMatrix> mmatrixa;
   LINALG::SplitMatrix2x2(mmatrix,adofs,idofs,mdofs,tempmap,mmatrixa,tempmtx1,tempmtx2,tempmtx3);

  /**********************************************************************/
  /* additional entries in master row                                   */
  /**********************************************************************/
  // do the multiplication mhataam = invda * mmatrixa
  RCP<LINALG::SparseMatrix> mhataam = rcp(new LINALG::SparseMatrix(*adofs,10));
  mhataam = LINALG::MLMultiply(*invda,false,*mmatrixa,false,false,false,true);
  mhataam->Complete(*mdofs,*adofs);

  // kmn: add T(mhataam)*kat
  RCP<LINALG::SparseMatrix> kmtadd = LINALG::MLMultiply(*mhataam,true,*kat,false,false,false,true);

  /**********************************************************************/
  /* additional entries in active tangential row                        */
  /**********************************************************************/
  // matrix T
  RCP<LINALG::SparseMatrix> tmatrix = cstrategy.TMatrix();

  // kaa: multiply tmatrix with invda and kaa
  RCP<LINALG::SparseMatrix> katadd;
  if (aset)
  {
    katadd = LINALG::MLMultiply(*tmatrix,false,*invda,true,false,false,true);
    katadd = LINALG::MLMultiply(*katadd,false,*kat,false,false,false,true);
  }

  /**********************************************************************/
  /* global setup of k_st_new                                           */
  /**********************************************************************/
  RCP<LINALG::SparseMatrix> k_st_new = rcp(new LINALG::SparseMatrix(*(StructureField()->Discretization()->DofRowMap(0)),81,true,false,k_st->GetMatrixtype()));
  k_st_new->Add(*knt,false,1.0,0.0);
  k_st_new->Add(*kmt,false,1.0,0.0);
  k_st_new->Add(*kmtadd,false,1.0,1.0);
  k_st_new->Add(*kit,false,1.0,1.0);
  if(aset) k_st_new->Add(*katadd,false,1.0,1.0);

  // FillComplete k_st_new
  k_st_new->Complete(*thermoprobrowmap,*structprobrowmap);

  // finally, do the replacement
  k_st = k_st_new;

  return;

} // ApplyStructContact()


/*----------------------------------------------------------------------*
 |  apply contact to diagonal block k_ts                      mgit 05/11 |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::ApplyThermContact(Teuchos::RCP<LINALG::SparseMatrix>& k_ts)
{
  // only in the case of contact
  if (cmtman_==Teuchos::null)
    return;

  // contact strategy
  MORTAR::StrategyBase& strategy = cmtman_->GetStrategy();
  CONTACT::CoAbstractStrategy& cstrategy = static_cast<CONTACT::CoAbstractStrategy&>(strategy);

  // check if contact contributions are present,
  // if not we can skip this routine to speed things up
  if (!cstrategy.IsInContact() && !cstrategy.WasInContact() && !cstrategy.WasInContactLastTimeStep())
    return;

  //**********************************************************************
  // maps/matrices form structural and thermal field
  //**********************************************************************
  // FIXGIT: This should be obtained from thermal field (and not build again)
  // convert maps (from structure discretization to thermo discretization)
  RCP<Epetra_Map> sdofs,adofs,idofs,mdofs,amdofs,ndofs,smdofs;
  RCP<Epetra_Map> thermoprobrowmap = rcp(new Epetra_Map(*(ThermoField()->Discretization()->DofRowMap(0))));
  thermcontman_->ConvertMaps(sdofs,adofs,mdofs);
  smdofs = LINALG::MergeMap(sdofs,mdofs,false);
  ndofs = LINALG::SplitMap(*(ThermoField()->Discretization()->DofRowMap(0)),*smdofs);

  // FIXGIT: This should be obtained form thermal field (and not build again)
  // structural mortar matrices, converted to thermal dofs
  RCP<LINALG::SparseMatrix> dmatrix = rcp(new LINALG::SparseMatrix(*sdofs,10));
  RCP<LINALG::SparseMatrix> mmatrix = rcp(new LINALG::SparseMatrix(*sdofs,100));
  thermcontman_->TransformDM(*dmatrix,*mmatrix,sdofs,mdofs);

  // FillComplete() global Mortar matrices
  dmatrix->Complete();
  mmatrix->Complete(*mdofs,*sdofs);

  // necessary map from structural problem
  RCP<Epetra_Map> structprobrowmap = rcp(new Epetra_Map(*(StructureField()->Discretization()->DofRowMap(0))));

  // abbreviations for active and inactive set
  int aset = adofs->NumGlobalElements();

  //**********************************************************************
  // linearization entries from mortar additional terms in balance equation
  // (lindmatrix, linmmatrix) with respect to displacements
  // and linearization entries from thermal contact condition (lindismatrix)
  // with respect to displacements
  //**********************************************************************

  // respective matrices
  RCP<LINALG::SparseMatrix> lindmatrix = rcp(new LINALG::SparseMatrix(*sdofs,100,true,false,LINALG::SparseMatrix::FE_MATRIX));
  RCP<LINALG::SparseMatrix> linmmatrix = rcp(new LINALG::SparseMatrix(*mdofs,100,true,false,LINALG::SparseMatrix::FE_MATRIX));
  RCP<LINALG::SparseMatrix> lindismatrix = rcp(new LINALG::SparseMatrix(*adofs,100,true,false,LINALG::SparseMatrix::FE_MATRIX));

  // assemble them
  AssembleLinDM(*lindmatrix,*linmmatrix);
  AssembleThermContCondition(*lindismatrix);

  // complete
  lindmatrix->Complete(*(cstrategy.SlaveMasterRowDofs()),*sdofs);
  linmmatrix->Complete(*(cstrategy.SlaveMasterRowDofs()),*mdofs);
  lindismatrix->Complete(*(cstrategy.SlaveMasterRowDofs()),*adofs);

  // add them to the off-diagonal block
  k_ts->Add(*lindmatrix,false,1.0,1.0);
  k_ts->Add(*linmmatrix,false,1.0,1.0);

  //**********************************************************************
  // splitting of matrix k_ts
  //**********************************************************************
  // complete
  k_ts->Complete(*structprobrowmap,*thermoprobrowmap);

  // matrix to split
  RCP<LINALG::SparseMatrix> k_temp_struct = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(k_ts);

  RCP<Epetra_Map> tmp;
  RCP<LINALG::SparseMatrix> ksmstruct,knstruct,ksstruct,kmstruct,kastruct,kistruct,tmp1,tmp2;

  // first split: k_temp_struct -> ksmstruct, knstruct
  LINALG::SplitMatrix2x2(k_temp_struct,smdofs,ndofs,structprobrowmap,tmp,ksmstruct,tmp1,knstruct,tmp2);

  // second split: ksmstruct -> ksstruct, kmstruct
  LINALG::SplitMatrix2x2(ksmstruct,sdofs,mdofs,structprobrowmap,tmp,ksstruct,tmp1,kmstruct,tmp2);

  // third split: ksstruct -> kastruct,kistruct
  LINALG::SplitMatrix2x2(ksstruct,adofs,idofs,structprobrowmap,tmp,kastruct,tmp1,kistruct,tmp2);

  /**********************************************************************/
  /* evaluation of the inverse of D                                     */
  /**********************************************************************/
  RCP<LINALG::SparseMatrix> invd = rcp(new LINALG::SparseMatrix(*dmatrix));
  RCP<Epetra_Vector> diag = LINALG::CreateVector(*sdofs,true);
  int err = 0;

  // extract diagonal of invd into diag
  invd->ExtractDiagonalCopy(*diag);

  // set zero diagonal values to dummy 1.0
  for (int i=0;i<diag->MyLength();++i)
    if ((*diag)[i]==0.0) (*diag)[i]=1.0;

  // scalar inversion of diagonal values
  err = diag->Reciprocal(*diag);
  if (err>0) dserror("ERROR: Reciprocal: Zero diagonal entry!");

  // re-insert inverted diagonal into invd
  err = invd->ReplaceDiagonalValues(*diag);

  // store some stuff for condensation of LM
  kts_ = ksstruct;
  invdtherm_ = invd;

  /**********************************************************************/
  /* evaluation of mhatmatrix, active parts                             */
  /**********************************************************************/
 // do the multiplication M^ = inv(D) * M
  RCP<LINALG::SparseMatrix> mhatmatrix = rcp(new LINALG::SparseMatrix(*sdofs,10));
  mhatmatrix = LINALG::MLMultiply(*invd,false,*mmatrix,false,false,false,true);
  mhatmatrix->Complete(*mdofs,*sdofs);

  // maps
  RCP<Epetra_Map> tempmap1,tmpmap;

  // active part of mhatmatrix and invd
  RCP<LINALG::SparseMatrix> mhata,invda,tempmtx1,tempmtx2,tempmtx3,tmp3;
  LINALG::SplitMatrix2x2(mhatmatrix,adofs,idofs,mdofs,tmpmap,mhata,tmp1,tmp2,tmp3);
  LINALG::SplitMatrix2x2(invd,sdofs,tempmap1,adofs,idofs,invda,tempmtx1,tempmtx2,tempmtx3);

  /**********************************************************************/
  /* additional entries in master row                                   */
  /**********************************************************************/
  // kmstructadd: add T(mhataam)*kan
  RCP<LINALG::SparseMatrix> kmstructadd = LINALG::MLMultiply(*mhata,true,*kastruct,false,false,false,true);

  /**********************************************************************/
  /* additional entries in active tangential row                        */
  /**********************************************************************/
  // thermcondLMmatrix
  RCP<LINALG::SparseMatrix> thermcondLMMatrix = thermcontman_->ThermCondLMMatrix();

  // kastructadd: multiply thermcontLMmatrix with invda and kastruct
  RCP<LINALG::SparseMatrix> kastructadd;
  if (aset)
  {
    kastructadd = LINALG::MLMultiply(*thermcondLMMatrix,false,*invda,false,false,false,true);
    kastructadd = LINALG::MLMultiply(*kastructadd,false,*kastruct,false,false,false,true);
  }

  /**********************************************************************/
  /* Global setup of k_ts_new                                           */
  /**********************************************************************/
  RCP<LINALG::SparseMatrix> k_ts_new = rcp(new LINALG::SparseMatrix(*(ThermoField()->Discretization()->DofRowMap(0)),81,true,false,k_ts->GetMatrixtype()));
  k_ts_new->Add(*knstruct,false,1.0,0.0);
  k_ts_new->Add(*kmstruct,false,1.0,0.0);
  k_ts_new->Add(*kmstructadd,false,1.0,1.0);
  k_ts_new->Add(*kistruct,false,1.0,1.0);
  if(aset) k_ts_new->Add(*kastructadd,false,1.0,1.0);
  if(aset) k_ts_new->Add(*lindismatrix,false,-1.0,1.0);

  // FillComplete k_ts_newteffnew (square)
  k_ts_new->Complete(*structprobrowmap,*thermoprobrowmap);

  // finally, do the replacement
  k_ts = k_ts_new;

  return;

} // ApplyThermContact()


/*----------------------------------------------------------------------*
 | recover structural and thermal Lagrange multipliers from   mgit 04/10|
 | displacements and temperature                                        |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::RecoverStructThermLM()
{
  // only in the case of contact
  if (cmtman_ == Teuchos::null)
    return;

  // initialize thermal Lagrange multiplier
  // FIXGIT: this should be done before
  // for structural LM, this is done in within the structural field
  RCP<Epetra_Map> sthermdofs,athermdofs,mthermdofs;
  thermcontman_->ConvertMaps (sthermdofs,athermdofs,mthermdofs);

  thermcontman_->InitializeThermLM(sthermdofs);

  // check if contact contributions are present,
  // if not we can skip this routine to speed things up
  // static cast of mortar strategy to contact strategy
  MORTAR::StrategyBase& strategy = cmtman_->GetStrategy();
  CONTACT::CoAbstractStrategy& cstrategy = static_cast<CONTACT::CoAbstractStrategy&>(strategy);

  if (!cstrategy.IsInContact() && !cstrategy.WasInContact() && !cstrategy.WasInContactLastTimeStep())
    return;

  // vector of displacement and temperature increments
  Teuchos::RCP<Epetra_Vector> sx;
  Teuchos::RCP<Epetra_Vector> tx;

  // extract field vectors
  ExtractFieldVectors(iterinc_,sx,tx);
  Teuchos::RCP<Epetra_Vector> siterinc = rcp(new Epetra_Vector((sx->Map())));
  siterinc->Update(1.0,*sx,0.0);
  Teuchos::RCP<Epetra_Vector> titerinc = rcp(new Epetra_Vector((tx->Map())));
  titerinc->Update(1.0,*tx,0.0);

  /**********************************************************************/
  /* recover of structural LM                                           */
  /**********************************************************************/
  // this requires two step
  // 1. recover structural LM from displacement dofs
  // 2. additionally evaluate part from thermal dofs

  // 1. recover structural LM form displacement dofs
  cmtman_->GetStrategy().Recover(siterinc);

  // 2. additionally evaluate part from thermal dofs
  // evaluate part from thermal dofs

  // matrices and maps
  RCP<LINALG::SparseMatrix> invda;
  RCP<Epetra_Map> tempmap;
  RCP<LINALG::SparseMatrix> tempmtx1, tempmtx2, tempmtx3;

  // necessary maps
  RCP<Epetra_Map> sdofs,adofs,idofs,mdofs,amdofs,ndofs,smdofs;
  sdofs = cstrategy.SlaveRowDofs();
  adofs = cstrategy.ActiveRowDofs();
  mdofs = cstrategy.MasterRowDofs();
  smdofs = LINALG::MergeMap(sdofs,mdofs,false);
  idofs =  LINALG::SplitMap(*sdofs,*adofs);
  ndofs = LINALG::SplitMap(*(StructureField()->Discretization()->DofRowMap(0)),*smdofs);

  // multiplication
  RCP<Epetra_Vector> mod = rcp(new Epetra_Vector(*sdofs));
  kst_->Multiply(false,*tx,*mod);

  // active part of invd
  LINALG::SplitMatrix2x2(invd_,adofs,tempmap,adofs,tempmap,invda,tempmtx1,tempmtx2,tempmtx3);
  RCP<LINALG::SparseMatrix> invdmod = rcp(new LINALG::SparseMatrix(*sdofs,10));
  invdmod->Add(*invda,false,1.0,1.0);
  invdmod->Complete();

  // vector to add
  RCP<Epetra_Vector> zadd = rcp(new Epetra_Vector(*sdofs));
  invdmod->Multiply(true,*mod,*zadd);

  // lagrange multipliers from structural field to be modified
  RCP<Epetra_Vector> lagrmult = cmtman_->GetStrategy().LagrMult();

  // modify structural Lagrange multipliers and store them to nodes
  lagrmult->Update(-1.0,*zadd,1.0);
  cmtman_->GetStrategy().StoreNodalQuantities(MORTAR::StrategyBase::lmupdate);


  /**********************************************************************/
  /* recover of thermal LM                                              */
  /**********************************************************************/
  // this requires two step
  // 1. recover thermal LM from temperature dofs
  // 2. additionally evaluate part from structural dofs

  // 1. recover thermal LM from temperature dofs
  thermcontman_->RecoverThermLM(titerinc);

  // 2. additionally evaluate part from structural dofs

  // matrices and maps
  RCP<LINALG::SparseMatrix> invdatherm;
  RCP<Epetra_Map> tempmaptherm;
  RCP<LINALG::SparseMatrix> tempmtx4, tempmtx5, tempmtx6;

  // necessary maps
  RCP<Epetra_Map> sdofstherm,adofstherm,idofstherm,mdofstherm;
  thermcontman_->ConvertMaps(sdofstherm,adofstherm,mdofstherm);
  idofstherm = LINALG::SplitMap(*sdofstherm,*adofstherm);

  // multiplication
  RCP<Epetra_Vector> modtherm = rcp(new Epetra_Vector(*sdofs));
  kts_->Multiply(false,*sx,*modtherm);

  // active part of invdtherm
  LINALG::SplitMatrix2x2(invdtherm_,adofstherm,tempmaptherm,adofstherm,tempmaptherm,invdatherm,tempmtx4,tempmtx5,tempmtx6);
  RCP<LINALG::SparseMatrix> invdmodtherm = rcp(new LINALG::SparseMatrix(*sdofstherm,10));
  invdmodtherm->Add(*invdatherm,false,1.0,1.0);
  invdmodtherm->Complete();

  // vector to add
  RCP<Epetra_Vector> zaddtherm = rcp(new Epetra_Vector(*sdofstherm));
  invdmodtherm->Multiply(true,*modtherm,*zaddtherm);

  // lagrange multipliers from thermal field to be modified
  RCP<Epetra_Vector> thermlagrmult = thermcontman_->ThermLM();

  // modify thermal Lagrange multiplier
  thermlagrmult->Update(+1.0,*zaddtherm,-1.0);

  return;

}  // RecoverStructThermLM()


/*----------------------------------------------------------------------*
 | linearization of D and M with respect to displacements     mgit 06/11 |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::AssembleLinDM(LINALG::SparseMatrix& lindglobal,
                                    LINALG::SparseMatrix& linmglobal)
{
  // stactic cast of mortar strategy to contact strategy
  MORTAR::StrategyBase& strategy = cmtman_->GetStrategy();
  CONTACT::CoAbstractStrategy& cstrategy = static_cast<CONTACT::CoAbstractStrategy&>(strategy);

  // get vector of contact interfaces
  vector<RCP<CONTACT::CoInterface> > interface = cstrategy.ContactInterfaces();

  // this currently works only for one interface yet
  if (interface.size()>1)
    dserror("Error in TSI::Algorithm::AssembleLinDM: Only for one interface yet.");

  // slave nodes
 const RCP<Epetra_Map> slavenodes = interface[0]->SlaveRowNodes();

  // loop over all slave nodes (row map)
  for (int j=0;j<slavenodes->NumMyElements();++j)
  {
    int gid = slavenodes->GID(j);
    DRT::Node* node = (interface[0]->Discret()).gNode(gid);
    DRT::Node* nodeges = ThermoField()->Discretization()->gNode(gid);

    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CONTACT::CoNode* cnode = static_cast<CONTACT::CoNode*>(node);

    int rowtemp = StructureField()->Discretization()->Dof(1,nodeges)[0];
    int locid = (thermcontman_->ThermLM()->Map()).LID(rowtemp);
    double lm = (*thermcontman_->ThermLM())[locid];

    // Mortar matrix D and M derivatives
    map<int,map<int,double> >& dderiv = cnode->CoData().GetDerivD();
    map<int,map<int,double> >& mderiv = cnode->CoData().GetDerivM();

    // get sizes and iterator start
    int slavesize = (int)dderiv.size();
    int mastersize = (int)mderiv.size();
    map<int,map<int,double> >::iterator scurr = dderiv.begin();
    map<int,map<int,double> >::iterator mcurr = mderiv.begin();

    /********************************************** LinDMatrix **********/
    // loop over all DISP slave nodes in the DerivD-map of the current LM slave node
    for (int k=0;k<slavesize;++k)
    {
      int sgid = scurr->first;
      ++scurr;

      DRT::Node* snode = interface[0]->Discret().gNode(sgid);
      DRT::Node* snodeges = ThermoField()->Discretization()->gNode(sgid);
      if (!snode) dserror("ERROR: Cannot find node with gid %",sgid);

      // Mortar matrix D derivatives
      map<int,double>& thisdderiv = cnode->CoData().GetDerivD()[sgid];
      int mapsize = (int)(thisdderiv.size());

      int row = StructureField()->Discretization()->Dof(1,snodeges)[0];

      map<int,double>::iterator scolcurr = thisdderiv.begin();

      // loop over all directional derivative entries
      for (int c=0;c<mapsize;++c)
      {
        int col = scolcurr->first;
        double val = lm*(scolcurr->second);
        ++scolcurr;

        if (abs(val)>1.0e-12) lindglobal.FEAssemble(-val,row,col);
      }

        // check for completeness of DerivD-Derivatives-iteration
        if (scolcurr!=thisdderiv.end())
          dserror("ERROR: AssembleLinDM: Not all derivative entries of DerivD considered!");
    }

    // check for completeness of DerivD-Slave-iteration
    if (scurr!=dderiv.end())
      dserror("ERROR: AssembleLinDM: Not all DISP slave entries of DerivD considered!");
    /******************************** Finished with LinDMatrix **********/


    /********************************************** LinMMatrix **********/
    // loop over all master nodes in the DerivM-map of the current LM slave node
    for (int l=0;l<mastersize;++l)
    {
      int mgid = mcurr->first;
      ++mcurr;

      DRT::Node* mnode = interface[0]->Discret().gNode(mgid);
      DRT::Node* mnodeges = ThermoField()->Discretization()->gNode(mgid);
      if (!mnode) dserror("ERROR: Cannot find node with gid %",mgid);

      // Mortar matrix M derivatives
      map<int,double>&thismderiv = cnode->CoData().GetDerivM()[mgid];
      int mapsize = (int)(thismderiv.size());

      int row = StructureField()->Discretization()->Dof(1,mnodeges)[0];
      map<int,double>::iterator mcolcurr = thismderiv.begin();

      // loop over all directional derivative entries
      for (int c=0;c<mapsize;++c)
      {
        int col = mcolcurr->first;
        double val =  lm * (mcolcurr->second);
        ++mcolcurr;

        // owner of LM slave node can do the assembly, although it actually
        // might not own the corresponding rows in lindglobal (DISP slave node)
        // (FE_MATRIX automatically takes care of non-local assembly inside!!!)
        //cout << "Assemble LinM: " << row << " " << col << " " << val << endl;
        if (abs(val)>1.0e-12) linmglobal.FEAssemble(val,row,col);
      }

      // check for completeness of DerivM-Derivatives-iteration
      if (mcolcurr!=thismderiv.end())
        dserror("ERROR: AssembleLinDM: Not all derivative entries of DerivM considered!");
    }

    // check for completeness of DerivM-Master-iteration
    if (mcurr!=mderiv.end())
      dserror("ERROR: AssembleLinDM: Not all master entries of DerivM considered!");
    /******************************** Finished with LinMMatrix **********/
  }
  return;

}  // AssembleLinDM()


/*----------------------------------------------------------------------*
 | linearization of thermal contact condition with respect to           |
 | displacements                                             mgit 06/11 |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::AssembleThermContCondition(LINALG::SparseMatrix& lindisglobal)
{
  // stactic cast of mortar strategy to contact strategy
  MORTAR::StrategyBase& strategy = cmtman_->GetStrategy();
  CONTACT::CoAbstractStrategy& cstrategy = static_cast<CONTACT::CoAbstractStrategy&>(strategy);

  // get vector of contact interfaces
  vector<RCP<CONTACT::CoInterface> > interface = cstrategy.ContactInterfaces();

  // this currently works only for one interface yet
  if (interface.size()>1)
    dserror("Error in TSI::Algorithm::AssembleThermContCondition: Only for one interface yet.");

  // heat transfer coefficient for slave and master surface
  double heattranss = interface[0]->IParams().get<double>("HEATTRANSSLAVE");
  double heattransm = interface[0]->IParams().get<double>("HEATTRANSMASTER");
  double beta = heattranss*heattransm/(heattranss+heattransm);

  // slave nodes
  const RCP<Epetra_Map> slavenodes = interface[0]->SlaveRowNodes();

  // loop over all LM slave nodes (row map)
  for (int j=0;j<slavenodes->NumMyElements();++j)
  {
    int gid = slavenodes->GID(j);

    DRT::Node* node = (interface[0]->Discret()).gNode(gid);
    DRT::Node* nodeges = ThermoField()->Discretization()->gNode(gid);

    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CONTACT::CoNode* cnode = static_cast<CONTACT::CoNode*>(node);

    if(cnode->Active()!=true)
      break;

    int row = StructureField()->Discretization()->Dof(1,nodeges)[0];

    // Mortar matrix D and M derivatives
    map<int,map<int,double> >& dderiv = cnode->CoData().GetDerivD();
    map<int,map<int,double> >& mderiv = cnode->CoData().GetDerivM();

    // get sizes and iterator start
    int slavesize = (int)dderiv.size();
    int mastersize = (int)mderiv.size();
    map<int,map<int,double> >::iterator scurr = dderiv.begin();
    map<int,map<int,double> >::iterator mcurr = mderiv.begin();

    /********************************************** LinDMatrix **********/
    // loop over all DISP slave nodes in the DerivD-map of the current LM slave node
    for (int k=0;k<slavesize;++k)
    {
      int sgid = scurr->first;
      ++scurr;

      DRT::Node* snode = interface[0]->Discret().gNode(sgid);
      DRT::Node* snodeges = ThermoField()->Discretization()->gNode(gid);
      if (!snode) dserror("ERROR: Cannot find node with gid %",sgid);

      int rowtemp = StructureField()->Discretization()->Dof(1,snodeges)[0];

      int locid = (thermcontman_->ThermLM()->Map()).LID(rowtemp);
      int locid1 = (ThermoField()->Tempnp()->Map()).LID(rowtemp);

      double lm = (*thermcontman_->ThermLM())[locid];
      double Ts = (*ThermoField()->Tempnp())[locid1];

      // Mortar matrix D derivatives
      map<int,double>& thisdderiv = cnode->CoData().GetDerivD()[sgid];
      int mapsize = (int)(thisdderiv.size());

      map<int,double>::iterator scolcurr = thisdderiv.begin();

      // loop over all directional derivative entries
      for (int c=0;c<mapsize;++c)
      {
        int col = scolcurr->first;
        double val = lm*(scolcurr->second);
        double val1 = -beta*Ts*(scolcurr->second);
        ++scolcurr;

        if (abs(val)>1.0e-12) lindisglobal.FEAssemble(val,row,col);
        if (abs(val1)>1.0e-12) lindisglobal.FEAssemble(val1,row,col);
      }

      // check for completeness of DerivD-Derivatives-iteration
      if (scolcurr!=thisdderiv.end())
        dserror("ERROR: AssembleThermContCondition: Not all derivative entries of DerivD considered!");
    }

    // check for completeness of DerivD-Slave-iteration
    if (scurr!=dderiv.end())
      dserror("ERROR: AssembleThermContCondition: Not all DISP slave entries of DerivD considered!");
    /******************************** Finished with LinDMatrix **********/


    /********************************************** LinMMatrix **********/
    // loop over all master nodes in the DerivM-map of the current LM slave node
    for (int l=0;l<mastersize;++l)
    {
      int mgid = mcurr->first;
      ++mcurr;

      DRT::Node* mnode = interface[0]->Discret().gNode(mgid);
      DRT::Node* mnodeges = ThermoField()->Discretization()->gNode(mgid);
      if (!mnode) dserror("ERROR: Cannot find node with gid %",mgid);

      int rowtemp = StructureField()->Discretization()->Dof(1,mnodeges)[0];

      int locid = (ThermoField()->Tempnp()->Map()).LID(rowtemp);
      double Tm = (*ThermoField()->Tempnp())[locid];

      // Mortar matrix M derivatives
      map<int,double>&thismderiv = cnode->CoData().GetDerivM()[mgid];
      int mapsize = (int)(thismderiv.size());

      map<int,double>::iterator mcolcurr = thismderiv.begin();

      // loop over all directional derivative entries
      for (int c=0;c<mapsize;++c)
      {
        int col = mcolcurr->first;
        double val = beta*Tm*(mcolcurr->second);
        ++mcolcurr;

        if (abs(val)>1.0e-12) lindisglobal.FEAssemble(val,row,col);
      }

      // check for completeness of DerivM-Derivatives-iteration
      if (mcolcurr!=thismderiv.end())
        dserror("ERROR: AssembleThermContCondition: Not all derivative entries of DerivM considered!");
    }

    // check for completeness of DerivM-Master-iteration
    if (mcurr!=mderiv.end())
      dserror("ERROR: AssembleThermContCondition: Not all master entries of DerivM considered!");
    /******************************** Finished with LinMMatrix **********/
  }
  return;

}  // AssembleThermContCondition()


/*----------------------------------------------------------------------*/

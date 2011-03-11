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
 | definitions                                               dano 12/09 |
 *----------------------------------------------------------------------*/
#ifdef CCADISCRET

#ifdef PARALLEL
#include <mpi.h>
#endif


/*----------------------------------------------------------------------*
 | headers                                                   dano 11/10 |
 *----------------------------------------------------------------------*/
#include "tsi_monolithic.H"
#include "tsi_defines.H"

#include <Teuchos_TimeMonitor.hpp>
// needed for PrintNewton
#include <sstream>

// include this header for coupling stiffness terms
#include "../drt_lib/drt_assemblestrategy.H"


//! Note: The order of calling the two BaseAlgorithm-constructors is
//! important here! In here control file entries are written. And these entries
//! define the order in which the filters handle the Discretizations, which in
//! turn defines the dof number ordering of the Discretizations.
/*----------------------------------------------------------------------*
 | constructor (public)                                      dano 11/10 |
 *----------------------------------------------------------------------*/
TSI::MonolithicBase::MonolithicBase(Epetra_Comm& comm)
  : AlgorithmBase(comm,DRT::Problem::Instance()->TSIDynamicParams()),
    StructureBaseAlgorithm(DRT::Problem::Instance()->TSIDynamicParams()),
    ThermoBaseAlgorithm(DRT::Problem::Instance()->TSIDynamicParams())
{
  // monolithic TSI must know the other discretization
  // build a proxy of the structure discretization for the temperature field
  Teuchos::RCP<DRT::DofSet> structdofset
    = StructureField().Discretization()->GetDofSetProxy();
  // build a proxy of the temperature discretization for the structure field
  Teuchos::RCP<DRT::DofSet> thermodofset
    = ThermoField().Discretization()->GetDofSetProxy();

  // check if ThermoField has 2 discretizations, so that coupling is possible
  if (ThermoField().Discretization()->AddDofSet(structdofset)!=1)
    dserror("unexpected dof sets in thermo field");
  if (StructureField().Discretization()->AddDofSet(thermodofset)!=1)
    dserror("unexpected dof sets in structure field");
}


/*----------------------------------------------------------------------*
 | destructor (public)                                       dano 11/10 |
 *----------------------------------------------------------------------*/
TSI::MonolithicBase::~MonolithicBase()
{
}


/*----------------------------------------------------------------------*
 | read restart information for given time step (public)     dano 11/10 |
 *----------------------------------------------------------------------*/
void TSI::MonolithicBase::ReadRestart(int step)
{
  ThermoField().ReadRestart(step);
  StructureField().ReadRestart(step);
  SetTimeStep(ThermoField().GetTime(),step);

  return;
}


/*----------------------------------------------------------------------*
 | prepare time step (public)                                dano 11/10 |
 *----------------------------------------------------------------------*/
void TSI::MonolithicBase::PrepareTimeStep()
{
  // counter and print header
  IncrementTimeAndStep();
  PrintHeader();

  // call the predictor
  StructureField().PrepareTimeStep();
  ThermoField().PrepareTimeStep();
}

/*----------------------------------------------------------------------*
 | update (protected)                                        dano 11/10 |
 *----------------------------------------------------------------------*/
void TSI::MonolithicBase::Update()
{
  StructureField().Update();
  ThermoField().Update();
}


/*----------------------------------------------------------------------*
 | output (protected)                                        dano 11/10 |
 *----------------------------------------------------------------------*/
void TSI::MonolithicBase::Output()
{
  // Note: The order is important here! In here control file entries are
  // written. And these entries define the order in which the filters handle
  // the Discretizations, which in turn defines the dof number ordering of the
  // Discretizations.
  StructureField().Output();

  // write the thermo output (temperatures at the moment) to the the structure output
  // get disc writer from structure field
  Teuchos::RCP<IO::DiscretizationWriter> output = StructureField().DiscWriter();

  // get the temperature and the noderowmap of thermo discretization
  Epetra_Vector temperature = *(ThermoField().Tempn());
  const Epetra_Map* temprowmap = ThermoField().Discretization()->NodeRowMap();

  // replace map and write it to output
  temperature.ReplaceMap(*temprowmap);
  RCP<Epetra_Vector> temp = rcp(new Epetra_Vector(temperature));
  output->WriteVector("temperature",temp);

  ThermoField().Output();
} // MonolithicBase::Output()





/*----------------------------------------------------------------------*
 | monolithic                                                dano 11/10 |
 *----------------------------------------------------------------------*/
TSI::Monolithic::Monolithic(
  Epetra_Comm& comm,
  const Teuchos::ParameterList& sdynparams
  )
  : MonolithicBase(comm),
    solveradapttol_(DRT::INPUT::IntegralValue<int>(sdynparams,"ADAPTCONV")==1),
    solveradaptolbetter_(sdynparams.get<double>("ADAPTCONV_BETTER")),
    printscreen_(true),  // ADD INPUT PARAMETER
    printiter_(true),  // ADD INPUT PARAMETER
    printerrfile_(true and errfile_),  // ADD INPUT PARAMETER FOR 'true'
    errfile_(NULL),
    zeros_(Teuchos::null),
    strmethodname_(DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(sdynparams,"DYNAMICTYP")),
    veln_(StructureField().ExtractVeln())
{
  // add extra parameters (a kind of work-around)
  Teuchos::RCP<Teuchos::ParameterList> xparams
    = Teuchos::rcp(new Teuchos::ParameterList());
  xparams->set<FILE*>("err file", DRT::Problem::Instance()->ErrorFile()->Handle());
  errfile_ = xparams->get<FILE*>("err file");

  veln_->PutScalar(0.0);
}

/*----------------------------------------------------------------------*
 | time loop of the monolithic system                        dano 11/10 |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::TimeLoop(
  const Teuchos::ParameterList& sdynparams
  )
{
  // time loop
  while (NotFinished())
  {
    // calculate initial linear system at current position
    // This initializes the field algorithms and creates the first linear
    // systems.
    // coupling parameters must be passed to the fields
    Evaluate(Teuchos::null); // pass null vector as unknown increments

    // counter and print header
    // predict solution of both field (call the adapter)
    PrepareTimeStep();

    // Newton-Raphson iteration
    NewtonFull(sdynparams);

    // update all single field solvers
    Update();

    // write output to screen and files
    Output();

#ifdef TSIMONOLITHASOUTPUT
    printf("Ende Timeloop ThermoField().ExtractTempnp[0] %12.8f\n",(*ThermoField().ExtractTempnp())[0]);
    printf("Ende Timeloop ThermoField().ExtractTempn[0] %12.8f\n",(*ThermoField().ExtractTempn())[0]);
#endif // TSIMONOLITHASOUTPUT

  }  // timeloop
}


/*----------------------------------------------------------------------*
 | solution with full Newton-Raphson iteration               dano 10/10 |
 | in tsi_algorithm: NewtonFull()                                       |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::NewtonFull(
  const Teuchos::ParameterList& sdynparams
  )
{
  cout << "TSI::Monolithic::NewtonFull()" << endl;

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
  norminc_ = 0.0;
  Epetra_Time timerthermo(Comm());
  timerthermo.ResetStartTime();

  // incremental solution vector with length of all TSI dofs
  iterinc_ = LINALG::CreateVector(*DofRowMap(), true);
  iterinc_->PutScalar(0.0);
  // a zero vector of full length
  zeros_ = LINALG::CreateVector(*DofRowMap(), true);
  zeros_->PutScalar(0.0);

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
    SetupSystemMatrix(sdynparams);

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

    // apply Dirichlet BCs to system of equations
    iterinc_->PutScalar(0.0);  // Useful? depends on solver and more

    // create the solver
    // get UMFPACK...
    Teuchos::ParameterList solverparams
      = DRT::Problem::Instance()->ThermalSolverParams();

    solver_ = rcp(new LINALG::Solver(
                        solverparams,
                        Comm(),
                        DRT::Problem::Instance()->ErrorFile()->Handle()
                        )
                  );

    //---------------------------------------------- solve linear system

    // Solve for inc_ = [disi_,tempi_]
    // Solve K_Teffdyn . IncX = -R  ===>  IncX_{n+1} with X=[d,T]
    // \f$x_{i+1} = x_i + \Delta x_i\f$
    if (solveradapttol_ and (iter_ > 1))
    {
      double worst = normrhs_;
      double wanted = tolfres_;
      solver_->AdaptTolerance(wanted, worst, solveradaptolbetter_);
    }
    // mere blockmatrix to SparseMatrix and solve
    // --> very expensive, but ok for development of new code 29.11.10
    Teuchos::RCP<LINALG::SparseMatrix> m = systemmatrix_->Merge();

    // standard solver call
    solver_->Solve(m->EpetraOperator(), iterinc_, rhs_, true, iter_==1);
    cout << " Solved" << endl;

    // blank all increments that have Dirichlet BC, because they must be zero
    Teuchos::RCP<Epetra_Vector> tmp = LINALG::CreateVector(*DofRowMap(), false);
    // blank the increments of the TSI system that have Dirichlet BC
    LINALG::ApplyDirichlettoSystem(iterinc_, tmp, zeros_, *CombinedDBCMap());
    tmp = Teuchos::null;
    cout << " DBC applied to TSI system" << endl;

    // reset solver tolerance
    solver_->ResetTolerance();

    // build residual force norm
    // for now use for simplicity only L2/Euclidian norm
    rhs_->Norm2(&normrhs_);
    // build residual increment norm
    iterinc_->Norm2(&norminc_);

    // print stuff
    PrintNewtonIter();

    // increment equilibrium loop index
    iter_ += 1;

  }  // end equilibrium loop

  // correct iteration counter
  iter_ -= 1;

  // test whether max iterations was hit
  if ( (iter_ >= itermax_) and (not iterdivercont_) )
  {
    dserror("Newton unconverged in %d iterations", iter_);
  }
  else if ( (iter_ >= itermax_) and (iterdivercont_) and (Comm().MyPID() == 0) )
  {
    printf("Newton unconverged in %d iterations ... continuing\n", iter_);
  }
  else if ( (Converged()) and (Comm().MyPID() == 0) )
  {
    PrintNewtonConv();
  }

}  // NewtonFull()


/*----------------------------------------------------------------------*
 | evaluate the single fields                                dano 11/10 |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::Evaluate(Teuchos::RCP<const Epetra_Vector> x)
{
  TEUCHOS_FUNC_TIME_MONITOR("TSI::Monolithic::Evaluate");

  // displacement and temperature incremental vector
  Teuchos::RCP<const Epetra_Vector> sx;
  Teuchos::RCP<const Epetra_Vector> tx;

  // if an increment vector exists
  if (x!=Teuchos::null)
  {
    // extract displacement sx and temperature tx incremental vector of global
    // unknown incremental vector x
    ExtractFieldVectors(x,sx,tx);

#ifdef TSIASOUTPUT
    cout << "Recent thermal increment DT_n+1^i\n" << *(tx) << endl;
    cout << "Recent structural increment Dd_n+1^i\n" << *(sx) << endl;

    cout << "Until here only old solution of Newton step. No update applied\n" << *(ThermoField().Tempnp()) << endl;
#endif // TSIASOUTPUT
  }
  // else(x=Teuchos::null): initialize the system

#ifdef TSIASOUTPUT
  cout << "Tempnp vor UpdateNewton\n" << *(ThermoField().Tempnp()) << endl;
  printf("Tempnp vor UpdateNewton ThermoField().ExtractTempnp[0] %12.8f\n",(*ThermoField().ExtractTempnp())[0]);
#endif // TSIASOUTPUT

  // Newton update of the thermo field
  // update temperature before passed to the structural field
  //   UpdateIterIncrementally(tx),
  ThermoField().UpdateNewton(tx);

#ifdef TSIASOUTPUT
  cout << "Tempnp nach UpdateNewton\n" << *(ThermoField().Tempnp()) << endl;
  printf("Tempnp nach UpdateNewton ThermoField().ExtractTempnp[0] %12.8f\n",(*ThermoField().ExtractTempnp())[0]);
#endif // TSIASOUTPUT

  // call all elements and assemble rhs and matrices
  cout << " \nEvaluate elements\n" << endl;


  /// structural field

  // structure Evaluate (builds tangent, residual and applies DBC)
  Epetra_Time timerstructure(Comm());

  // apply current temperature to structure
  StructureField().ApplyTemperatures(ThermoField().Tempnp());

#ifdef TSIPARALLEL
  cout << Comm().MyPID() << " nach ApplyTemp!!" << endl;
#endif // TSIPARALLEL

#ifdef TSIASOUTPUT
//    Teuchos::RCP<Epetra_Vector> tempera = rcp(new Epetra_Vector(ThermoField().Tempn()->Map(),true));
//    if (ThermoField().Tempnp() != Teuchos::null)
//      tempera->Update(1.0, *ThermoField().Tempnp(), 0.0);
//    StructureField().ApplyTemperatures(tempera);
//    StructureField().ApplyTemperatures(ThermoField().Tempn());
#endif // TSIASOUTPUT

  // Monolithic TSI accesses the linearised structure problem:
  //   UpdaterIterIncrementally(sx),
  //   EvaluateForceStiffResidual()
  //   PrepareSystemForNewtonSolve()
  StructureField().Evaluate(sx);
  cout << " structure time for calling Evaluate: " << timerstructure.ElapsedTime() << "\n";

#ifdef TSIASOUTPUT
  cout << "STR fres_" << *StructureField().RHS() << endl;
#endif // TSIASOUTPUT

  /// thermal field

  // thermo Evaluate
  // (builds tangent, residual and applies DBC and recent coupling values)
  Epetra_Time timerthermo(Comm());

  // apply current displacements and velocities to the thermo field
  if (strmethodname_ == INPAR::STR::dyna_statics)
  {
    // calculate velocity V_n+1^k = (D_n+1^k-D_n)/Dt()
    veln_ = CalcVelocity(StructureField().Dispnp());
  }
  else
  {
    veln_ = StructureField().ExtractVelnp();
  }
  // pass the structural values to the thermo field
  ThermoField().ApplyStructVariables(StructureField().Dispnp(),veln_);

#ifdef TSIASOUTPUT
  cout << "d_n+1 inserted in THR field\n" << *(StructureField().Dispnp()) << endl;
  cout << "v_n+1\n" << *veln_ << endl;
#endif // TSIASOUTPUT

  // monolithic TSI accesses the linearised thermo problem
  //   EvaluateRhsTangResidual() and
  //   PrepareSystemForNewtonSolve()
  ThermoField().Evaluate();
  cout << " thermo time for calling Evaluate: " << timerthermo.ElapsedTime() << "\n";

}  // Evaluate()


/*----------------------------------------------------------------------*
 | extract field vectors for calling Evaluate() of the       dano 11/10 |
 | single fields                                                        |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::ExtractFieldVectors(
  Teuchos::RCP<const Epetra_Vector> x,
  Teuchos::RCP<const Epetra_Vector>& sx,
  Teuchos::RCP<const Epetra_Vector>& tx
  )
{
  TEUCHOS_FUNC_TIME_MONITOR("TSI::Monolithic::ExtractFieldVectors");

  // process structure unknowns of the first field
  sx = Extractor().ExtractVector(x,0);

  // process thermo unknowns of the second field
  tx = Extractor().ExtractVector(x,1);
}


/*----------------------------------------------------------------------*
 | calculate velocities                                      dano 12/10 |
 | like InterfaceVelocity(disp) in FSI::DirichletNeumann                |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> TSI::Monolithic::CalcVelocity(
  Teuchos::RCP<const Epetra_Vector> sx
  )
{
  Teuchos::RCP<Epetra_Vector> vel = Teuchos::null;
  // copy D_n onto V_n+1
  vel = rcp(new Epetra_Vector( *(StructureField().ExtractDispn()) ) );
  // calculate velocity with timestep Dt()
  //  V_n+1^k = (D_n+1^k - D_n) / Dt
  vel->Update(1./Dt(), *sx, -1./Dt());

  return vel;
}  // CalcVelocity()


/*----------------------------------------------------------------------*
 | setup system (called in tsi_dyn)                          dano 11/10 |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::SetupSystem()
{
  cout << " TSI::Monolithic::SetupSystem()" << endl;

  // create combined map
  std::vector<Teuchos::RCP<const Epetra_Map> > vecSpaces;

#ifdef TSIPARALLEL
  cout << Comm().MyPID() << " :PID" << endl;
  cout << "structure dofmap" << endl;
  cout << *StructureField().DofRowMap(0) << endl;
  cout << "thermo dofmap" << endl;
  cout << *StructureField().DofRowMap(1) << endl;
#endif // TSIPARALLEL

  // use its own DofRowMap, that is the 0th map of the discretization
  vecSpaces.push_back(StructureField().DofRowMap(0));
  vecSpaces.push_back(ThermoField().DofRowMap(0));

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
  blockrowdofmap_.Setup(*fullmap,maps);
}


/*----------------------------------------------------------------------*
 | setup system matrix of TSI                                dano 11/10 |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::SetupSystemMatrix(
  const Teuchos::ParameterList& sdynparams
  )
{
  cout << " TSI::Monolithic::SetupSystemMatrix()" << endl;
  TEUCHOS_FUNC_TIME_MONITOR("TSI::Monolithic::SetupSystemMatrix");

  /*----------------------------------------------------------------------*/
  // initialize TSI-systemmatrix_
  systemmatrix_
    = rcp(
        new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
          Extractor(),
          Extractor(),
          81,
          false
          )
        );

  /*----------------------------------------------------------------------*/
  // pure structural part k_ss (3nx3n)

  // build pure structural block k_ss
  // build block matrix
  // The maps of the block matrix have to match the maps of the blocks we
  // insert here. Extract Jacobian matrices and put them into composite system
  // matrix W
  Teuchos::RCP<LINALG::SparseMatrix> k_ss = StructureField().SystemMatrix();

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
                     *(StructureField().Discretization()->DofRowMap(0)),
                     81,
                     true,
                     true
                     )
                   );

  // call the element and calculate the matrix block
  ApplyStrCouplMatrix(k_st);

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
  // matrix W
  Teuchos::RCP<LINALG::SparseMatrix> k_tt = ThermoField().SystemMatrix();

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
                     *(ThermoField().Discretization()->DofRowMap(0)),
                     81,
                     true,
                     true
                     )
                   );

  // call the element and calculate the matrix block
  ApplyThrCouplMatrix(k_ts,sdynparams);

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
  cout << " TSI::Monolithic::SetupRHS()" << endl;
  TEUCHOS_FUNC_TIME_MONITOR("TSI::Monolithic::SetupRHS");

  // create full monolithic rhs vector
  rhs_ = rcp(new Epetra_Vector(*DofRowMap(), true));

  // fill the TSI rhs vector rhs_ with the single field rhss
  SetupVector(
    *rhs_,
    StructureField().RHS(),
    ThermoField().RHS()
    );

}  // SetupRHS()


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
    StructureField().InitialGuess(),
    // returns residual temperatures or iterative thermal increment - tempi_
    ThermoField().InitialGuess()
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
  Extractor().InsertVector(*sv,0,f);
  Extractor().InsertVector(*tv,1,f);
}


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
      dserror("Cannot check for convergence of residual values!");
  }

  // residual forces
  switch (normtypefres_)
  {
    case INPAR::TSI::convnorm_abs:
      convfres = normrhs_ < tolfres_;
      break;
    default:
      dserror("Cannot check for convergence of residual forces!");
  }

  // combine temperature-like and force-like residuals
  bool conv = false;
  if (combincfres_ == INPAR::TSI::bop_and)
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
  if ( (Comm().MyPID() == 0) and printscreen_ and printiter_ )
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
  default:
    dserror("You should not turn up here.");
  }

  switch ( normtypeinc_ )
  {
  case INPAR::TSI::convnorm_abs :
    oss <<std::setw(18)<< "abs-temp-norm";
    break;
  default:
    dserror("You should not turn up here.");
  }

  // add solution time
  oss << std::setw(14)<< "wct";

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
    oss << std::setw(18) << std::setprecision(5) << std::scientific << normrhs_;
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
  if (ofile == NULL)
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
  // create the parameters for the discretization
  Teuchos::ParameterList sparams;
  const std::string action = "calc_struct_stifftemp";
  sparams.set("action", action);
//  cout << "STR Parameterliste\n " <<  sparams << endl;
  StructureField().Discretization()->SetState(0,"displacement",StructureField().Dispnp());

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
  StructureField().Discretization()->Evaluate( sparams, structuralstrategy );

}  // ApplyStrCouplMatrix()


/*----------------------------------------------------------------------*
 |  evaluate thermal-mechanical system matrix at state       dano 03/11 |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::ApplyThrCouplMatrix(
  Teuchos::RCP<LINALG::SparseMatrix> k_ts,  //!< off-diagonal tangent matrix term
  const Teuchos::ParameterList& sdynparams
  )
{
  // create the parameters for the discretization
  Teuchos::ParameterList tparams;
  // type of calling structural time integrator
  tparams.set<int>("time integrator", strmethodname_);
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
      double theta_ = sdynparams.sublist("ONESTEPTHETA").get<double>("THETA");
      tparams.set<double>("theta", theta_);
      break;
    }
    case  INPAR::STR::dyna_genalpha :
    {
      double alphaf_ = sdynparams.sublist("GENALPHA").get<double>("ALPHA_F");
      double beta_ = sdynparams.sublist("GENALPHA").get<double>("BETA");
      double gamma_ = sdynparams.sublist("GENALPHA").get<double>("GAMMA");
      tparams.set<double>("ALPHA_F", alphaf_);
      tparams.set<double>("BETA", beta_);
      tparams.set<double>("GAMMA", gamma_);
    }
    default :
    {
      dserror("Don't know what to do...");
      break;
    }
  }  // end of switch(strmethodname_)
  // action for elements
  const std::string action = "calc_thermo_coupltang";
  tparams.set("action", action);
  // other parameters that might be needed by the elements
  tparams.set("delta time", Dt());
  tparams.set("total time", Time());

  // set the variables that are needed by the elements
  ThermoField().Discretization()->SetState(0,"temperature",ThermoField().Tempnp());
  ThermoField().Discretization()->SetState(1,"displacement",StructureField().Dispnp());
  ThermoField().Discretization()->SetState(1,"velocity",veln_);

//  cout << "veln_" << *veln_ << endl;

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
  ThermoField().Discretization()->Evaluate(tparams,thermostrategy);

}  // ApplyThrCouplMatrix()


/*----------------------------------------------------------------------*
 |  map containing the dofs with Dirichlet BC                dano 03/11 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> TSI::Monolithic::CombinedDBCMap()
{
  const Teuchos::RCP<const Epetra_Map > scondmap = StructureField().GetDBCMapExtractor()->CondMap();
  const Teuchos::RCP<const Epetra_Map > tcondmap = ThermoField().GetDBCMapExtractor()->CondMap();
  Teuchos::RCP<Epetra_Map> condmap = LINALG::MergeMap(scondmap, tcondmap, false);
  return condmap;
} // CombinedDBCMap()


/*----------------------------------------------------------------------*/
#endif  // CCADISCRET

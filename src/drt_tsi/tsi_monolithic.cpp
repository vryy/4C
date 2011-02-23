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
#include <Teuchos_TimeMonitor.hpp>
// needed for PrintNewton
#include <sstream>


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
TSI::Monolithic::Monolithic(Epetra_Comm& comm)
  : MonolithicBase(comm),
    printscreen_(true),  // ADD INPUT PARAMETER
    printiter_(true),  // ADD INPUT PARAMETER
    printerrfile_(true and errfile_),  // ADD INPUT PARAMETER FOR 'true'
    errfile_(NULL)
{
  // add extra parameters (a kind of work-around)
  Teuchos::RCP<Teuchos::ParameterList> xparams
    = Teuchos::rcp(new Teuchos::ParameterList());
  xparams->set<FILE*>("err file", DRT::Problem::Instance()->ErrorFile()->Handle());
  errfile_ = xparams->get<FILE*>("err file");

  // if structure field is quasi-static --> CalcVelocity
  const Teuchos::ParameterList& sdyn
    = DRT::Problem::Instance()->StructuralDynamicParams();
  // major switch to different time integrators
  quasistatic_
    = DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(sdyn,"DYNAMICTYP")
        ==INPAR::STR::dyna_statics;
}

/*----------------------------------------------------------------------*
 | time loop of the monolithic system                        dano 11/10 |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::TimeLoop()
{
  zeros_ = Teuchos::null;
  // a zero vector of full length of all TSI GID
  zeros_ = LINALG::CreateVector(*DofRowMap(), true);

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
    NewtonFull();

    // update all single field solvers
    Update();

    // write output to screen and files
    Output();

  }  // timeloop
}


/*----------------------------------------------------------------------*
 | solution with full Newton-Raphson iteration               dano 10/10 |
 | in tsi_algorithm: NewtonFull()                                       |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::NewtonFull()
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
  Epetra_Time timerthermo(Comm());
  timerthermo.ResetStartTime();

  // incremental solution vector with length of all TSI dofs
  iterinc_ = LINALG::CreateVector(*DofRowMap(), true);

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
    SetupRHS();

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

    // make negative residual not necessary: rhs_ is yet TSI negative
    // (Newton-ready) residual with blanked Dirichlet DOFs (see adapter_timint!)

    // apply Dirichlet BCs to system of equations
    iterinc_->PutScalar(0.0);  // Useful? depends on solver and more

    //---------------------------------------------- solve linear system

    // Solve for inc_ = [disi_,tempi_]
    // Solve K_Teffdyn . IncX = -R  ===>  IncX_{n+1} with X=[d,T]
    // \f$x_{i+1} = x_i + \Delta x_i\f$

    // call the thermo parameter list
    const Teuchos::ParameterList& tdynparams
      = DRT::Problem::Instance()->ThermalDynamicParams();
    solveradapttol_ = DRT::INPUT::IntegralValue<int>(tdynparams,"ADAPTCONV")==1;
    solveradaptolbetter_ = tdynparams.get<double>("ADAPTCONV_BETTER");
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
    cout << " solved" << endl;

    // reset solver tolerance
    solver_->ResetTolerance();

    //-------------------------------------- update configuration values
    // update end-point temperatures etc.
//    NewtonUpdate();

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
  }
  // else(x=Teuchos::null): initialize the system

  // call all elements and assemble rhs and matrices
  cout << "\nEvaluate elements\n" << endl;

  // structure Evaluate (builds tangent, residual and applies DBC)
  {
    Epetra_Time timerstructure(Comm());

    // apply current temperature increments to structure
    StructureField().ApplyTemperatures(ThermoField().Tempnp());

#ifdef TSIASOUTPUT
    cout << "T_n+1 inserted in STR field\n" << *(ThermoField().Tempnp()) << endl;
#endif // TSIASOUTPUT

    // Monolithic TSI accesses the linearised structure problem:
    //   UpdaterIterIncrementally(sx),
    //   EvaluateForceStiffResidual()
    //   PrepareSystemForNewtonSolve()
    StructureField().Evaluate(sx);
    cout << "structure time for calling Evaluate: " << timerstructure.ElapsedTime() << "\n";
  }

  // thermo Evaluate (builds tangent, residual and applies DBC)
  {
    Epetra_Time timerthermo(Comm());

    // apply current displacements and velocities to the thermo field
    Teuchos::RCP<const Epetra_Vector> velnp;
    if (quasistatic_)
    {
      // calculate velocity V_n+1^k = (D_n+1^k-D_n)/Dt()
      velnp = CalcVelocity(StructureField().Dispnp());
    }
    else
    {
      velnp = StructureField().ExtractVelnp();
    }
    // pass the structural values to the thermo field
    ThermoField().ApplyStructVariables(StructureField().Dispnp(),velnp);

#ifdef TSIASOUTPUT
    cout << "d_n+1 inserted in THR field\n" << *(StructureField().Dispnp()) << endl;
    cout << "v_n+1\n" << *velnp << endl;
#endif // TSIASOUTPUT

    // monolithic TSI accesses the linearised thermo problem
    //   UpdateIterIncrementally(tx),
    //   EvaluateRhsTangResidual() and
    //   PrepareSystemForNewtonSolve()
    ThermoField().Evaluate(tx);
    cout << "thermo time for calling Evaluate: " << timerthermo.ElapsedTime() << "\n";
  }
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
  cout << "TSI::Monolithic::SetupSystem()" << endl;

  // create combined map
  std::vector<Teuchos::RCP<const Epetra_Map> > vecSpaces;
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
  Teuchos::RCP<Epetra_Map> fullmap = LINALG::MultiMapExtractor::MergeMaps(maps);
  blockrowdofmap_.Setup(*fullmap,maps);
}


/*----------------------------------------------------------------------*
 | setup system matrix of TSI                                dano 11/10 |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::SetupSystemMatrix()
{
  cout << "TSI::Monolithic::SetupSystemMatrix()" << endl;
  TEUCHOS_FUNC_TIME_MONITOR("TSI::Monolithic::SetupSystemMatrix");

  /*----------------------------------------------------------------------*/
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
  // structure part

  Teuchos::RCP<LINALG::SparseMatrix> s = StructureField().SystemMatrix();

  // build block matrix
  // The maps of the block matrix have to match the maps of the blocks we
  // insert here.

  // uncomplete because the fluid interface can have more connections than the
  // structural one. (Tet elements in fluid can cause this.) We should do
  // this just once...
  s->UnComplete();

  // assign structure part to the TSI matrix
  systemmatrix_->Assign(0,0,View,*s);

  /*----------------------------------------------------------------------*/
  // thermo part

  // build block matrix
  // The maps of the block matrix have to match the maps of the blocks we
  // insert here. Extract Jacobian matrices and put them into composite system
  // matrix W
  Teuchos::RCP<LINALG::SparseMatrix> t = ThermoField().SystemMatrix();

  // Uncomplete fluid matrix to be able to deal with slightly defective
  // interface meshes.
  t->UnComplete();

  // assign thermo part to the TSI matrix
  systemmatrix_->Assign(1,1,View,*t);

  /*----------------------------------------------------------------------*/
  // done. make sure all blocks are filled.
  systemmatrix_->Complete();

}  // SetupSystemMatrix


/*----------------------------------------------------------------------*
 | setup RHS (like fsimon)                                   dano 11/10 |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::SetupRHS()
{
  cout << "TSI::Monolithic::SetupRHS()" << endl;
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
 | iteration update of state                                 dano 11/10 |
 | originally by bborn 08/09                                            |                   |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::NewtonUpdate()
{
//  //! new end-point temperatures
//  //! T_{n+1}^{<k+1>} := T_{n+1}^{<k>} + IncT_{n+1}^{<k>}
//  timestepincn_->Update(1.0, *iterinc_, 1.0);
//
//  //! bye
//  return;
}


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


/*----------------------------------------------------------------------*/
#endif  // CCADISCRET

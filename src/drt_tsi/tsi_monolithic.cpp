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
}





/*----------------------------------------------------------------------*
 | monolithic                                                dano 11/10 |
 *----------------------------------------------------------------------*/
TSI::Monolithic::Monolithic(Epetra_Comm& comm)
  : MonolithicBase(comm)
{
}

/*----------------------------------------------------------------------*
 | time loop of the monolithic system                        dano 11/10 |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::TimeLoop()
{
  // time loop
  while (NotFinished())
  {
    PrepareTimeStep();

    // calculate initial linear system at current position
    // This initializes the field algorithms and creates the first linear
    // systems.
    Evaluate(Teuchos::null); // pass null vector as unknown increments

    // get initial guess
    // The initial system is there, so we can happily extract the
    // initial guess. (The Dirichlet conditions are already build in!)
    Teuchos::RCP<Epetra_Vector> initial_guess
      = Teuchos::rcp(new Epetra_Vector(*DofRowMap()));
    InitialGuess(initial_guess);

    // predictor step
    PredictConstValueRate();

    // create the linear system
    // \f$J(x_i) \Delta x_i = - R(x_i)\f$

    // create the solver

    // create the systemmatrix/Jacobian (FSI monolithic: computeJacobian() )
//    TEUCHOS_FUNC_TIME_MONITOR("TSI::Monolithic::ComputeSystemmatrix");
//    Evaluate(Teuchos::rcp(&x,false));
//    LINALG::BlockSparseMatrixBase& mat = Teuchos::dyn_cast<LINALG::BlockSparseMatrixBase>(Jac);
//    SetupSystemMatrix(mat);

    // Newton-Raphson iteration
    NewtonFull();

    // update all single field solvers
    Update();

    // write output to screen and files
    Output();

  } // time loop
}

/*----------------------------------------------------------------------*
 | evaluate the single fields                                dano 11/10 |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::Evaluate(Teuchos::RCP<const Epetra_Vector> x)
{
  TEUCHOS_FUNC_TIME_MONITOR("TSI::Monolithic::Evaluate");

  // displacement and temperature incremental vector
  Teuchos::RCP<const Epetra_Vector> sx;
  Teuchos::RCP<const Epetra_Vector> tx;

  // 24.11.10
  cout << "sx\n" << sx << endl;
  cout << "tx\n" << tx << endl;

  // extract displacement sx and temperature tx incremental vector of global
  // unknown incremental vector x
  ExtractFieldVectors(x,sx,tx);

  // call all elements and assemble rhs and matrices
  cout << "\nEvaluate elements\n" << endl;

  {
    Epetra_Time ts(Comm());
    // builds tangent, residual and applies DBC
    // Monolithic TSI accesses the linearised structure problem
    // UpdaterIterIncrementally(sx), EvaluateForceStiffResidual() and
    // PrepareSystemForNewtonSolve()
    StructureField().Evaluate(sx);
    cout << "structure: " << ts.ElapsedTime() << "\n";
  }

  {
    Epetra_Time tt(Comm());
    // builds tangent, residual and applies DBC
    // monolithic TSI accesses the linearised thermo problem
    // UpdateIterIncrementally, EvaluateRhsTangResidual() and
    // PrepareSystemForNewtonSolve()
    ThermoField().Evaluate(tx);
    cout << "thermo: " << tt.ElapsedTime() << "\n";
  }
}


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

  // process structure unknowns
  sx = Extractor().ExtractVector(x,0);

  // process thermo unknowns
  tx = Extractor().ExtractVector(x,1);

}


/*----------------------------------------------------------------------*
 | initial guess of the displacements/temperatures           dano 11/10 |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::InitialGuess(Teuchos::RCP<Epetra_Vector> ig)
{
  TEUCHOS_FUNC_TIME_MONITOR("TSI::Monolithic::InitialGuess");

  SetupVector(*ig,
              StructureField().InitialGuess(),
              ThermoField().InitialGuess()
              );
}


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
 | solution with full Newton-Raphson iteration               dano 10/10 |
 | in tsi_algorithm: NewtonFull()                               |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::NewtonFull()
{
//  // we do a Newton-Raphson iteration here.
//  // the specific time integration has set the following
//  // --> On #fres_ is the positive force residuum
//  // --> On #systemmatrix_ is the effective dynamic tangent matrix
//
//  // check whether we have a sanely filled tangent matrix
//  if (not systemmatrix_->Filled())
//  {
//    dserror("Effective tangent matrix must be filled here");
//  }
//
//  // initialise equilibrium loop
//  iter_ = 1;
//  normfres_ = CalcRefNormForce();
//  // normtempi_ was already set in predictor; this is strictly >0
//  timer_.ResetStartTime();
//
//  // equilibrium iteration loop
//  while ( ( (not Converged()) and (iter_ <= itermax_) ) or (iter_ <= itermin_) )
//  {
//    // make negative residual
//    fres_->Scale(-1.0);
//
//    // apply Dirichlet BCs to system of equations
//    tempi_->PutScalar(0.0);  // Useful? depends on solver and more
//    LINALG::ApplyDirichlettoSystem(tang_, tempi_, fres_,
//                                   Teuchos::null, zeros_, *(dbcmaps_->CondMap()));
//
//    //---------------------------------------------- solve linear system
//
//    // Solve for inc_ = [disi_,tempi_]
//    // Solve K_Teffdyn . IncX = -R  ===>  IncX_{n+1} with X=[d,T]
//    // \f$x_{i+1} = x_i + \Delta x_i\f$
//    // Solve K_Teffdyn . IncT = -R  ===>  IncT_{n+1}
//    if (solveradapttol_ and (iter_ > 1))
//    {
//      double worst = normfres_;
//      double wanted = tolfres_;
//      solver_->AdaptTolerance(wanted, worst, solveradaptolbetter_);
//    }
//    // standard solver call
//    solver_->Solve(systemmatrix_->EpetraMatrix(), tempi_, fres_, true, iter_==1);
//    // reset solver tolerance
//    solver_->ResetTolerance();
//
//    //-------------------------------------- update configuration values
//    // update end-point temperatures etc
//    UpdateIter(iter_);
//
//    // compute residual forces #fres_ and tangent #tang_
//    // whose components are globally oriented
//    EvaluateRhsTangResidual();
//
//    // extract reaction forces
//    // reactions are negative to balance residual on DBC
//    freact_->Update(-1.0, *fres_, 0.0);
//    dbcmaps_->InsertOtherVector(dbcmaps_->ExtractOtherVector(zeros_), freact_);
//
//    // blank residual at DOFs on Dirichlet BC
//    dbcmaps_->InsertCondVector(dbcmaps_->ExtractCondVector(zeros_), fres_);
//
//    // build residual force norm
//    normfres_ = THR::AUX::CalculateVectorNorm(iternorm_, fres_);
//    // build residual temperature norm
//    normtempi_ = THR::AUX::CalculateVectorNorm(iternorm_, tempi_);
//
//    // print stuff
//    PrintNewtonIter();
//
//    // increment equilibrium loop index
//    iter_ += 1;
//
//  }  // end equilibrium loop
//
}  // NewtonFull


/*----------------------------------------------------------------------*
 | setup system (called in tsi_dyn)                          dano 11/10 |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::SetupSystem()
{

  // create combined map
  std::vector<Teuchos::RCP<const Epetra_Map> > vecSpaces;
  vecSpaces.push_back(StructureField().DofRowMap());
  vecSpaces.push_back(ThermoField().DofRowMap());

  if (vecSpaces[0]->NumGlobalElements()==0)
    dserror("No structure equation. Panic.");
  if (vecSpaces[1]->NumGlobalElements()==0)
    dserror("No temperature equation. Panic.");

  SetDofRowMaps(vecSpaces);

  // create block system matrix
  systemmatrix_ = Teuchos::rcp(new TSIBlockMatrix(
                                     Extractor(),
                                     StructureField(),
                                     ThermoField()
                                     )
                                );
} // SetupSystem()


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
 | setup RHS                                                 dano 11/10 |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::SetupRHS(Epetra_Vector& f, bool firstcall)
{
  TEUCHOS_FUNC_TIME_MONITOR("TSI::Monolithic::SetupRHS");

  SetupVector(f,
              StructureField().RHS(),
              ThermoField().RHS()
              );
}


/*----------------------------------------------------------------------*
 | setup system matrix of TSI                                dano 11/10 |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::SetupSystemMatrix(LINALG::BlockSparseMatrixBase& mat)
{
  TEUCHOS_FUNC_TIME_MONITOR("TSI::MonolithicStructureSplit::SetupSystemMatrix");

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
  mat.Assign(0,0,View,*s);

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
  mat.Assign(1,1,View,*t);

  /*----------------------------------------------------------------------*/
  // done. make sure all blocks are filled.
  mat.Complete();

}  // SetupSystemMatrix


/*----------------------------------------------------------------------*
 | tsi block matrix                                          dano 11/10 |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::PredictConstValueRate()
{
  // constant predictor
  // d_n+1^0 = d_n
//  disn_->Update(1.0, *(*dis_)(0), 0.0);
//  veln_->Update(1.0, *(*vel_)(0), 0.0);
//  accn_->Update(1.0, *(*acc_)(0), 0.0);

  // T_n+1^0 = T_n
//  tempn_->Update(1.0, *(*temp_)(0), 0.0);
//  raten_->Update(1.0, *(*rate_)(0), 0.0);

  // see you next time step
}


/*----------------------------------------------------------------------*
 | convergence check for both fields (thermo & structure)    dano 11/10 |
 | originally for partitioned TSI (by vg 01/09)                         |
 *----------------------------------------------------------------------*/
bool TSI::Monolithic::ConvergenceCheck(
  int itnum,
  const int itmax,
  const double ittol
  )
{
  // convergence check based on the temperature increment
  bool stopnonliniter = false;

  //    | temperature increment |_2
  //  -------------------------------- < Tolerance
  //     | temperature_n+1 |_2

  // Variables to save different L2 - Norms
  // define L2-norm of incremental temperature and temperature
  // here: only the temperature field is checked for convergence!!!
  double tempincnorm_L2(0.0);
  double tempnorm_L2(0.0);
  double dispincnorm_L2(0.0);
  double dispnorm_L2(0.0);

  // build the current temperature increment Inc T^{i+1}
  // \f Delta T^{k+1} = Inc T^{k+1} = T^{k+1} - T^{k}  \f
  tempincnp_->Update(1.0,*(ThermoField().Tempnp()),-1.0);
  dispincnp_->Update(1.0,*(StructureField().Dispnp()),-1.0);

  // build the L2-norm of the temperature increment and the temperature
  tempincnp_->Norm2(&tempincnorm_L2);
  ThermoField().Tempnp()->Norm2(&tempnorm_L2);
  dispincnp_->Norm2(&dispincnorm_L2);
  StructureField().Dispnp()->Norm2(&dispnorm_L2);

  // care for the case that there is (almost) zero temperature
  // (usually not required for temperature)
  if (tempnorm_L2 < 1e-6) tempnorm_L2 = 1.0;
  if (dispnorm_L2 < 1e-6) dispnorm_L2 = 1.0;

  // Print the incremental based convergence check to the screen
  if (Comm().MyPID() == 0)
  {
    cout<<"\n";
    cout<<"***********************************************************************************\n";
    cout<<"    OUTER ITERATION STEP    \n";
    cout<<"***********************************************************************************\n";
    printf("+--------------+------------------------+--------------------+--------------------+\n");
    printf("|-  step/max  -|-  tol      [norm]     -|--  temp-inc      --|--  disp-inc      --|\n");
    printf("|   %3d/%3d    |  %10.3E[L_2 ]      | %10.3E         | %10.3E         |",
         itnum,itmax,ittol,tempincnorm_L2/tempnorm_L2,dispincnorm_L2/dispnorm_L2);
    printf("\n");
    printf("+--------------+------------------------+--------------------+--------------------+\n");
  }

  // Converged
  if ((tempincnorm_L2/tempnorm_L2 <= ittol) &&
      (dispincnorm_L2/dispnorm_L2 <= ittol))
  {
    stopnonliniter = true;
    if (Comm().MyPID() == 0)
    {
      printf("\n");
      printf("|  Outer Iteration loop converged after iteration %3d/%3d !                       |\n", itnum,itmax);
      printf("+--------------+------------------------+--------------------+--------------------+\n");
    }
  }

  // warn if itemax is reached without convergence, but proceed to next
  // timestep
  if ((itnum == itmax) and
       ((tempincnorm_L2/tempnorm_L2 > ittol) || (dispincnorm_L2/dispnorm_L2 > ittol))
     )
  {
    stopnonliniter = true;
    if ((Comm().MyPID() == 0))
    {
      printf("|     >>>>>> not converged in itemax steps!                                       |\n");
      printf("+--------------+------------------------+--------------------+--------------------+\n");
      printf("\n");
      printf("\n");
    }
  }

  return stopnonliniter;
}  // TSI::Monolithic::ConvergenceCheck




/*----------------------------------------------------------------------*
 | tsi block matrix                                          dano 11/10 |
 *----------------------------------------------------------------------*/
TSI::TSIBlockMatrix::TSIBlockMatrix(
  const LINALG::MultiMapExtractor& maps,
  ADAPTER::Structure& structure,
  ADAPTER::Thermo& thermo
  )
  : LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(maps,maps,81,false,true)
{
  // 25.11.10
  //structuresolver_ = Teuchos::rcp(new LINALG::Preconditioner(structure.LinearSolver()));
  //thermosolver_ = Teuchos::rcp(new LINALG::Preconditioner(thermo.LinearSolver()));

}  // TSIBlockMatrix


/*----------------------------------------------------------------------*
 | tsi block matrix                                          dano 11/10 |
 *----------------------------------------------------------------------*/
void TSI::TSIBlockMatrix::SetupTSIBlockMatrix()
{
//  const LINALG::SparseMatrix& structOp = Matrix(0,0);
//  const LINALG::SparseMatrix& thermoOp  = Matrix(1,1);
//
//  RCP<LINALG::MapExtractor> tsidofmapex = null;
//  RCP<Epetra_Map> irownodes = null;
//
//  structuresolver_->Setup(structOp.EpetraMatrix());
//  thermosolver_->Setup(thermoOp.EpetraMatrix());

}  // SetupTSIBlockMatrix


/*----------------------------------------------------------------------*
 | label the tsi block matrix                                dano 11/10 |
 *----------------------------------------------------------------------*/
const char* TSI::TSIBlockMatrix::Label() const
{
  return "TSI::TSIBlockMatrix";
}


/*----------------------------------------------------------------------*/
#endif  // CCADISCRET

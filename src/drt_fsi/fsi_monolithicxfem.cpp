#ifdef CCADISCRET

#include "fsi_monolithicxfem.H"
#include "fsi_overlapprec_fsiamg.H"
#include "fsi_statustest.H"
#include "fsi_nox_linearsystem_bgs.H"
#include "fsi_monolithic_linearsystem.H"
#include "../drt_lib/linalg_mapextractor.H"

#include "fsi_nox_group.H"
#include "fsi_nox_newton.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/inpar_fsi.H"

#include "../drt_io/io_control.H"


/*----------------------------------------------------------------------*/
// Note: The order of calling the three BaseAlgorithm-constructors is
// important here! In here control file entries are written. And these
// entries define the order in which the filters handle the
// Discretizations, which in turn defines the dof number ordering of the
// Discretizations.
/*----------------------------------------------------------------------*/
FSI::MonolithicBaseXFEM::MonolithicBaseXFEM(Epetra_Comm& comm)
  : AlgorithmBase(comm,DRT::Problem::Instance()->FSIDynamicParams()),
    StructureBaseAlgorithm(DRT::Problem::Instance()->FSIDynamicParams()),
    fluidfield_(DRT::Problem::Instance()->FSIDynamicParams(),"FSICoupling")
{
  ADAPTER::Coupling& coupsf = StructureFluidCoupling();

  // structure to fluid

  coupsf.SetupConditionCoupling(*StructureField().Discretization(),
                                 StructureField().Interface().CondMap(),
                                *FluidField().Discretization(),
                                 FluidField().Interface().FSICondMap(),
                                 "FSICoupling");

}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::MonolithicBaseXFEM::~MonolithicBaseXFEM()
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicBaseXFEM::ReadRestart(int step)
{
  StructureField().ReadRestart(step);
  FluidField()    .ReadRestart(step);

  //SetTimeStep(FluidField().Time(),FluidField().Step());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicBaseXFEM::PrepareTimeStep()
{
  cout << "FSI::MonolithicBaseXFEM::PrepareTimeStep()" << endl;
  IncrementTimeAndStep();

  PrintHeader();

  StructureField().PrepareTimeStep();
  FluidField().    PrepareTimeStep();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicBaseXFEM::Update()
{
  cout << "FSI::MonolithicBaseXFEM::Update()" << endl;
  StructureField().Update();
  FluidField().    Update();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicBaseXFEM::Output()
{
  // Note: The order is important here! In here control file entries are
  // written. And these entries define the order in which the filters handle
  // the Discretizations, which in turn defines the dof number ordering of the
  // Discretizations.
  cout << "FSI::MonolithicBaseXFEM::Output()" << endl;
  StructureField().Output();
  FluidField().    Output();
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::MonolithicXFEM::MonolithicXFEM(Epetra_Comm& comm)
  : MonolithicBaseXFEM(comm)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::Timeloop()
{
  cout << "FSI::MonolithicXFEM::Timeloop()" << endl;

  Teuchos::RCP<std::ofstream> log;
  if (Comm().MyPID()==0)
  {
    std::string s = DRT::Problem::Instance()->OutputControlFile()->FileName();
    s.append(".iteration");
    log = Teuchos::rcp(new std::ofstream(s.c_str()));
    (*log) << "# num procs      = " << Comm().NumProc() << "\n"
           << "#\n"
      ;
  }

  Teuchos::Time timer("time step timer");

  while (NotFinished())
  {
    PrepareTimeStep();

    // start time measurement
    Teuchos::RCP<Teuchos::TimeMonitor> timemonitor = rcp(new Teuchos::TimeMonitor(timer,true));

    Newton();

    // stop time measurement
    timemonitor = Teuchos::null;

    if (Comm().MyPID()==0)
    {
      (*log) << Step()
             << " " << timer.totalElapsedTime()
        ;
      (*log) << std::endl;
    }

    Update();
    Output();
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::Newton()
{
  cout << "FSI::MonolithicXFEM::Newton()" << endl;

  dispincsum_ = Teuchos::null;
  x_ = Teuchos::null;
  matrix_ = Teuchos::null;

  for (int i =0;i<2;++i)
  {
    Evaluate();

    SetupRHS();

    if (ConverganceTest())
    {
      break;
    }


    SetupSystemMatrix();

    // solve
    LinearSolve();
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::Evaluate()
{
  // take current increment and setup linear system
  cout << "FSI::MonolithicXFEM::Evaluate()" << endl;

  if (x_!=Teuchos::null)
  {
    cout << "later step" << endl;

    // extract structure displacement
    Teuchos::RCP<Epetra_Vector> sx = extractor_->ExtractStructureVector(x_);

    // residual displacements (or iteration increments or iteratively incremental displacements)
    Teuchos::RCP<Epetra_Vector> disp = Teuchos::rcp(new Epetra_Vector(*sx));

    // update incremental displacement member to provided step increments
    // shortly: disp^<i+1> := disinc_^<i>
    dispincsum_->Update(1.0, *disp, 1.0);

    // expects x^{n+1}_{i+1}
    StructureField().Evaluate(dispincsum_);



    Teuchos::RCP<Epetra_Vector> fx = extractor_->ExtractFluidVector(x_);

    // expects \Delta x^{n+1}_{i+1} ???
    Teuchos::RCP<Epetra_Vector> idisp = StructureField().Interface().ExtractCondVector(sx);

    FluidField().Evaluate(StructToFluid(idisp), fx);
  }
  else
  {
    cout << "first step" << endl;
    StructureField().Evaluate(Teuchos::null);
    FluidField().Evaluate(Teuchos::null, Teuchos::null);

    dispincsum_ = LINALG::CreateVector(*StructureField().DofRowMap(), true);
  }


  // setup global (full monolithic) map
  std::vector<Teuchos::RCP<const Epetra_Map> > maps;
  maps.push_back(StructureField().DofRowMap());
  maps.push_back(FluidField().DofRowMap());
  Teuchos::RCP<Epetra_Map> fullmapnew = LINALG::MultiMapExtractor::MergeMaps(maps);

  extractor_ = rcp(new MonolithicXFEMExtractor);
  extractor_->Setup(*fullmapnew,maps);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::SetupRHS()
{
  // get rhs from fields

  // create full monolithic rhs vector
  Teuchos::RCP<const Epetra_Vector> srhs = StructureField().RHS();
  Teuchos::RCP<const Epetra_Vector> frhs = FluidField().RHS();

  rhs_ = extractor_->InsertStructureVector(srhs);
  extractor_->InsertFluidVector(frhs,rhs_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool FSI::MonolithicXFEM::ConverganceTest()
{
  cout << "FSI::MonolithicXFEM::ConverganceTest()" << endl;


  const double ittol = 1.0e-8;

  const RCP<const Epetra_Map> dofrowmap = extractor_->FullMap();
  Epetra_Vector full(*dofrowmap);
  Epetra_Import importer(*dofrowmap,rhs_->Map());

//  double incfullnorm_L2;
//  double fullnorm_L2;
  double fullresnorm;

  int err = full.Import(*rhs_,importer,Insert);
  if (err) dserror("Import using importer returned err=%d",err);
  full.Norm2(&fullresnorm);

//  err = full.Import(*x_,importer,Insert);
//  if (err) dserror("Import using importer returned err=%d",err);
//  full.Norm2(&incfullnorm_L2);
//
//  err = full.Import(*state_.velnp_,importer,Insert);
//  if (err) dserror("Import using importer returned err=%d",err);
//  full.Norm2(&fullnorm_L2);


  // test if rhs is zero (equilibrium)
  bool converged = false;
  cout << "fullresnorm = " << fullresnorm << endl;
  if (fullresnorm <= ittol
//      and incfullnorm_L2/fullnorm_L2 <= ittol
      )
  {
    converged = true;
  }

  return converged;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::SetupSystemMatrix()
{
  cout << "FSI::MonolithicXFEM::SetupSystemMatrix()" << endl;
  // build global block matrix

  const ADAPTER::Coupling& coupsf = StructureFluidCoupling();

  Teuchos::RCP<LINALG::SparseMatrix> s = StructureField().SystemMatrix();
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> blocks =
    s->Split<LINALG::DefaultBlockMatrixStrategy>(StructureField().Interface(),
                                                 StructureField().Interface());
  blocks->Complete();

  Teuchos::RCP<LINALG::SparseMatrix> f = FluidField().SystemMatrix();

  /// direct access to system matrix
  const map<std::string, Teuchos::RCP<LINALG::SparseMatrix> > cmats = FluidField().CouplingMatrices();

  Teuchos::RCP<LINALG::SparseMatrix> Cuu = cmats.find("Cuu")->second;
  Teuchos::RCP<LINALG::SparseMatrix> Mud = cmats.find("Mud")->second;
  Teuchos::RCP<LINALG::SparseMatrix> Mdu = cmats.find("Mdu")->second;
  Teuchos::RCP<LINALG::SparseMatrix> Cdd = cmats.find("Cdd")->second;


  /*----------------------------------------------------------------------*/

  double scale     = FluidField().ResidualScaling();
  double timescale = FluidField().TimeScaling();

  matrix_ = rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(*extractor_,*extractor_));

  matrix_->Assign(0,0,View,*s);
  matrix_->Assign(1,1,View,*f);

//  matrix_->Matrix(1,0).Add(*Mud,false,1.0,0.0);
  sigtransform_(*blocks,
                *Mud,
                1./timescale,
                ADAPTER::Coupling::SlaveConverter(coupsf),
                matrix_->Matrix(1,0),
                true,
                false);

  //  matrix_->Matrix(0,0).Add(*Cdd,false,1.0,1.0);
  sggtransform_(*Cdd,
                1./(scale*timescale),
                ADAPTER::Coupling::SlaveConverter(coupsf),
                ADAPTER::Coupling::SlaveConverter(coupsf),
                matrix_->Matrix(0,0),
                true,
                true);

  //  matrix_->Matrix(0,1).Add(*Mdu,false,1.0,0.0);
  sgitransform_(*Mdu,
                1./scale,
                ADAPTER::Coupling::SlaveConverter(coupsf),
                matrix_->Matrix(0,1));


//  matrix_->Matrix(1,1).Add(*Cuu,false,1.0,1.0);
  matrix_->Complete();


}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::LinearSolve()
{
  cout << "FSI::MonolithicXFEM::LinearSolve()" << endl;

  // get initial guess from fields
  // solve linear system

  // x^{n+1}_{0} = x^{n}
  //
  // Newton increment \Delta x^{n+1}_{i+1} = x^{n+1}_{i+1} - x^{n+1}_{i}

  //Teuchos::RCP<const Epetra_Vector> sx = StructureField().InitialGuess();
  //Teuchos::RCP<const Epetra_Vector> fx = FluidField().InitialGuess();

  x_ = LINALG::CreateVector(*extractor_->FullMap(),true);

  x_ = extractor_->InsertStructureVector(StructureField().InitialGuess());

  Teuchos::RCP<LINALG::SparseMatrix> m = matrix_->Merge();

  LINALG::PrintMatrixInMatlabFormat("monomatrix.txt", *m->EpetraMatrix(), true);

  // get UMFPACK...
  Teuchos::ParameterList solverparams = DRT::Problem::Instance()->FluidSolverParams();


  Teuchos::RCP<LINALG::Solver> solver =
      rcp(new LINALG::Solver(solverparams,
                             Comm(),
                             DRT::Problem::Instance()->ErrorFile()->Handle()));

  solver->Solve(m->EpetraOperator(), x_, rhs_, true, true);

  cout << x_ << endl;

}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::MonolithicXFEM::StructToFluid(Teuchos::RCP<Epetra_Vector> iv)
{
  const ADAPTER::Coupling& coupsf = StructureFluidCoupling();
  return coupsf.MasterToSlave(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::MonolithicXFEM::FluidToStruct(Teuchos::RCP<Epetra_Vector> iv)
{
  const ADAPTER::Coupling& coupsf = StructureFluidCoupling();
  return coupsf.SlaveToMaster(iv);
}
#endif

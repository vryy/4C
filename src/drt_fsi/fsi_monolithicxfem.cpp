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

  x_ = Teuchos::null;
  systemmatrix_ = Teuchos::null;

  for (int i =0;i<3;++i)
  {
    Evaluate();

    const bool firstcall = i==0;
    SetupRHS(firstcall);

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

    StructureField().Evaluate(sx);


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
void FSI::MonolithicXFEM::SetupRHS(bool firstcall)
{
  cout << "FSI::MonolithicXFEM::SetupRHS" << endl;

  // get rhs from fields

  // create full monolithic rhs vector
  Teuchos::RCP<const Epetra_Vector> srhs = StructureField().RHS();
  Teuchos::RCP<const Epetra_Vector> frhs = FluidField().RHS();

  rhs_ = extractor_->InsertStructureVector(srhs);
  extractor_->InsertFluidVector(frhs,rhs_);

//  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicXFEM::SetupRHS");
//
//  SetupVector(*rhs_,
//              StructureField().RHS(),
//              FluidField().RHS(),
//              FluidField().ResidualScaling());
//
  const map<std::string, Teuchos::RCP<Epetra_Vector> > cvecs = FluidField().CouplingVectors();
  Teuchos::RCP<Epetra_Vector> rhsd = cvecs.find("rhsd")->second;
  cout << *rhsd << endl;

  if (firstcall)
  {
//    Teuchos::RCP<LINALG::BlockSparseMatrixBase> blockf = FluidField().BlockSystemMatrix();

//    LINALG::SparseMatrix& fig = blockf->Matrix(0,1);
//    LINALG::SparseMatrix& fgg = blockf->Matrix(1,1);
    /// direct access to system matrix
    const map<std::string, Teuchos::RCP<LINALG::SparseMatrix> > cmats = FluidField().CouplingMatrices();


  //  Teuchos::RCP<LINALG::SparseMatrix> Cuu = cmats.find("Cuu")->second; // schon erledigt im Fluid
    Teuchos::RCP<LINALG::SparseMatrix> Mud = cmats.find("Mud")->second;
    Teuchos::RCP<LINALG::SparseMatrix> Mdu = cmats.find("Mdu")->second;
    Teuchos::RCP<LINALG::SparseMatrix> Cdd = cmats.find("Cdd")->second;


    LINALG::SparseMatrix& fig = *cmats.find("Mdu")->second;
    LINALG::SparseMatrix& fgg = *FluidField().SystemMatrix();

    Teuchos::RCP<Epetra_Vector> fveln = FluidField().ExtractInterfaceVeln();
    double timescale = FluidField().TimeScaling();
    double scale     = FluidField().ResidualScaling();

//    Teuchos::RCP<Epetra_Vector> rhs_fig = Teuchos::rcp(new Epetra_Vector(fig.RowMap()));
//    fig.Apply(*fveln,*rhs_fig);
//    rhs_fig->Scale(timescale*Dt());
//    Extractor().AddVector(*rhs_fig,1,*rhs_);
//
//    Teuchos::RCP<Epetra_Vector> rhs = Teuchos::rcp(new Epetra_Vector(fgg.RowMap()));
//    fgg.Apply(*fveln,*rhs);
//    rhs->Scale(scale*timescale*Dt());
//    Extractor().AddVector(*StructureField().Interface().InsertCondVector(FluidToStruct(rhs)),0,*rhs_);

  }

  Extractor().AddVector(*StructureField().Interface().InsertCondVector(FluidToStruct(rhsd)),0,*rhs_);

  // NOX expects a different sign here.
//  f.Scale(-1.);
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

  // extract Jacobian matrices and put them into composite system
  // matrix W

  const ADAPTER::Coupling& coupsf = StructureFluidCoupling();
  //const ADAPTER::Coupling& coupsa = StructureAleCoupling();

  Teuchos::RCP<LINALG::SparseMatrix> s = StructureField().SystemMatrix();
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> blocks =
    s->Split<LINALG::DefaultBlockMatrixStrategy>(StructureField().Interface(),
                                                 StructureField().Interface());
  blocks->Complete();

  Teuchos::RCP<LINALG::SparseMatrix> f = FluidField().SystemMatrix();

  /// direct access to system matrix
  const map<std::string, Teuchos::RCP<LINALG::SparseMatrix> > cmats = FluidField().CouplingMatrices();

//  Teuchos::RCP<LINALG::SparseMatrix> Cuu = cmats.find("Cuu")->second; // schon erledigt im Fluid
  Teuchos::RCP<LINALG::SparseMatrix> Mud = cmats.find("Mud")->second;
  Teuchos::RCP<LINALG::SparseMatrix> Mdu = cmats.find("Mdu")->second;
  Teuchos::RCP<LINALG::SparseMatrix> Cdd = cmats.find("Cdd")->second;


//  // split fluid matrix
//
//  Teuchos::RCP<LINALG::BlockSparseMatrixBase> blockf = FluidField().BlockSystemMatrix();
//
//  LINALG::SparseMatrix& fii = blockf->Matrix(0,0);
//  LINALG::SparseMatrix& fig = blockf->Matrix(0,1);
//  LINALG::SparseMatrix& fgi = blockf->Matrix(1,0);
//  LINALG::SparseMatrix& fgg = blockf->Matrix(1,1);



  /*----------------------------------------------------------------------*/

  double scale     = FluidField().ResidualScaling();
  double timescale = FluidField().TimeScaling();

  systemmatrix_ = rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(*extractor_,*extractor_));

  systemmatrix_->Assign(0,0,View,*s);
  systemmatrix_->Assign(1,1,View,*f);

//  matrix_->Matrix(1,0).Add(*Mud,false,1.0,0.0);
  sigtransform_(*blocks,
                *Mud,
                1./timescale,
                ADAPTER::Coupling::SlaveConverter(coupsf),
                systemmatrix_->Matrix(1,0),
                true,
                false);

  //  matrix_->Matrix(0,0).Add(*Cdd,false,1.0,1.0);
  sggtransform_(*Cdd,
                1./(scale*timescale),
                ADAPTER::Coupling::SlaveConverter(coupsf),
                ADAPTER::Coupling::SlaveConverter(coupsf),
                systemmatrix_->Matrix(0,0),
                true,
                true);

  //  matrix_->Matrix(0,1).Add(*Mdu,false,1.0,0.0);
  sgitransform_(*Mdu,
                1./scale,
                ADAPTER::Coupling::SlaveConverter(coupsf),
                systemmatrix_->Matrix(0,1));


//  matrix_->Matrix(1,1).Add(*Cuu,false,1.0,1.0);

  // done. make sure all blocks are filled.
  systemmatrix_->Complete();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::InitialGuess(Teuchos::RCP<Epetra_Vector> ig)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicXFEM::InitialGuess");

  SetupVector(*ig,
              StructureField().InitialGuess(),
              FluidField().InitialGuess(),
              0.0);
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

  Teuchos::RCP<LINALG::SparseMatrix> m = systemmatrix_->Merge();

  LINALG::PrintMatrixInMatlabFormat("monomatrix.txt", *m->EpetraMatrix(), true);

  // get UMFPACK...
  Teuchos::ParameterList solverparams = DRT::Problem::Instance()->FluidSolverParams();


  Teuchos::RCP<LINALG::Solver> solver =
      rcp(new LINALG::Solver(solverparams,
                             Comm(),
                             DRT::Problem::Instance()->ErrorFile()->Handle()));

//  incvel_->PutScalar(0.0);
//  solver->Solve(m->EpetraOperator(), incvel_, rhs_, true, true);
  solver->Solve(m->EpetraOperator(), x_, rhs_, true, true);
//  state_.velnp_->Update(1.0,*incvel_,1.0);
  cout << *x_ << endl;

}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::SetupVector(Epetra_Vector &f,
                                      Teuchos::RCP<const Epetra_Vector> sv,
                                      Teuchos::RCP<const Epetra_Vector> fv,
                                      double fluidscale)
{
  Extractor().InsertVector(*sv,0,f);
  Extractor().InsertVector(*fv,1,f);
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

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
                                 StructureField().Interface().FSICondMap(),
                                *FluidField().Discretization(),
                                 FluidField().Interface().FSICondMap(),
                                 "FSICoupling");

  // Use splitted structure matrix
  StructureField().UseBlockMatrix();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicBaseXFEM::SetDofRowMaps(const std::vector<Teuchos::RCP<const Epetra_Map> >& maps)
{
  cout << "FSI::MonolithicBaseXFEM::SetDofRowMaps" << endl;
  Teuchos::RCP<Epetra_Map> fullmap = LINALG::MultiMapExtractor::MergeMaps(maps);
  extractor_.Setup(*fullmap,maps);
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

//  // Get the top level parameter list
//  Teuchos::ParameterList& nlParams = NOXParameterList();
//
//  // sublists
//
//  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
//  //Teuchos::ParameterList& solverOptions = nlParams.sublist("Solver Options");
//  Teuchos::ParameterList& newtonParams = dirParams.sublist("Newton");
//  Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");
//
//  //Teuchos::ParameterList& searchParams = nlParams.sublist("Line Search");
//  Teuchos::ParameterList& printParams = nlParams.sublist("Printing");
//  printParams.set("MyPID", Comm().MyPID());
//
//#if 0
//  // there is a strange NOX bug...
//  Teuchos::ParameterList& oo = printParams.sublist("Output Information");
//
//  oo.set<bool>("Error",false);
//  oo.set<bool>("Warning",false);
//  oo.set<bool>("Outer Iteration",false);
//  oo.set<bool>("Inner Iteration",false);
//  oo.set<bool>("Parameters",false);
//  oo.set<bool>("Details",false);
//  oo.set<bool>("Outer Iteration StatusTest",false);
//  oo.set<bool>("Linear Solver Details",false);
//  oo.set<bool>("Test Details",false);
//  oo.set<bool>("Stepper Iteration",false);
//  oo.set<bool>("Stepper Details",false);
//  oo.set<bool>("Stepper Parameters",false);
//  oo.set<bool>("Debug",false);
//#else
//  printParams.set("Output Information",
//                  NOX::Utils::Error |
//                  NOX::Utils::Warning |
//                  NOX::Utils::OuterIteration |
//                  NOX::Utils::InnerIteration |
//                  //NOX::Utils::Parameters |
//                  NOX::Utils::Details |
//                  NOX::Utils::OuterIterationStatusTest |
//                  NOX::Utils::LinearSolverDetails |
//                  NOX::Utils::TestDetails |
//                  NOX::Utils::StepperIteration |
//                  NOX::Utils::StepperDetails |
//                  NOX::Utils::StepperParameters |
//                  NOX::Utils::Debug |
//                  0);
//#endif
//
//  // Create printing utilities
//  utils_ = Teuchos::rcp(new NOX::Utils(printParams));

  Teuchos::RCP<std::ofstream> log;
  if (Comm().MyPID()==0)
  {
    std::string s = DRT::Problem::Instance()->OutputControlFile()->FileName();
    s.append(".iteration");
    log = Teuchos::rcp(new std::ofstream(s.c_str()));
    (*log) << "# num procs      = " << Comm().NumProc() << "\n"
//           << "# Method         = " << nlParams.sublist("Direction").get("Method","Newton") << "\n"
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
//             << " " << nlParams.sublist("Output").get("Nonlinear Iterations",0)
//             << " " << nlParams.sublist("Output").get("2-Norm of Residual", 0.)
//             << " " << lsParams.sublist("Output").get("Total Number of Linear Iterations",0)
        ;
      (*log) << std::endl;
//      lsParams.sublist("Output").set("Total Number of Linear Iterations",0);
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

  rhs_ = Teuchos::null;
  stepinc_ = Teuchos::null;
  systemmatrix_ = Teuchos::null;

  for (int i =0;i<100;++i)
  {
    cout << "Newton step: " << i << endl;
    Evaluate();

    const bool firstcall = i==0;
    SetupRHS(firstcall);

    if (ConverganceTest() or i > 3)
    {
      break;
    }


    SetupSystemMatrix();

    LinearSolve();
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::Evaluate()
{
  // take current increment and setup linear system
  cout << "FSI::MonolithicXFEM::Evaluate()" << endl;

  if (stepinc_!=Teuchos::null)
  {
    cout << "later step" << endl;

    // extract structure displacement
    Teuchos::RCP<Epetra_Vector> sxstepinc_i = Extractor().ExtractStructureInteriorVector(stepinc_);
    Teuchos::RCP<Epetra_Vector> sxstepinc_b = Extractor().ExtractStructureBoundaryVector(stepinc_);

    Teuchos::RCP<Epetra_Vector> sxstepinc = StructureField().Interface().InsertOtherVector(sxstepinc_i);
    StructureField().Interface().InsertFSICondVector(sxstepinc_b,sxstepinc);

    StructureField().Evaluate(sxstepinc);


//    if (boundarydis_->Comm().MyPID() == 0 && itruerescol->MyLength() >= 3)
    {
      std::ofstream f;
      const std::string fname = DRT::Problem::Instance()->OutputControlFile()->FileName()
                              + ".outifacedispstepinc.txt";
      f.open(fname.c_str(),std::fstream::ate | std::fstream::app);
      f << (*sxstepinc)[0] << "  " << "\n";
      f.close();
    }

    cout << "solid stepinc:" << endl;
    cout << *sxstepinc << endl;

    Teuchos::RCP<Epetra_Vector> fxstepinc = Extractor().ExtractFluidVector(stepinc_);
    cout << "fluid stepinc:" << endl;
    cout << *fxstepinc << endl;

    FluidField().Evaluate(StructToFluid(sxstepinc_b), fxstepinc);
  }
  else
  {
    cout << "first step" << endl;
    StructureField().Evaluate(Teuchos::null);
    FluidField().Evaluate(Teuchos::null, Teuchos::null);
  }


  // setup global (full monolithic) map
  std::vector<Teuchos::RCP<const Epetra_Map> > maps;
  maps.push_back(StructureField().Interface().OtherMap());
  maps.push_back(StructureField().Interface().FSICondMap());
  maps.push_back(FluidField()    .DofRowMap());

  if (maps[0]->NumGlobalElements()==0)
    dserror("No inner structural equations. Splitting not possible. Panic.");

  SetDofRowMaps(maps);

  Extractor().CheckForValidMapExtractor();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::SetupRHS(bool firstcall)
{
  cout << "FSI::MonolithicXFEM::SetupRHS" << endl;

  // get rhs from fields
  Teuchos::RCP<const Epetra_Vector> srhs = StructureField().RHS();
  Teuchos::RCP<const Epetra_Vector> frhs = FluidField().RHS();
  Teuchos::RCP<Epetra_Vector> sxi = StructureField().Interface().ExtractOtherVector(srhs);
  Teuchos::RCP<Epetra_Vector> sxb = StructureField().Interface().ExtractFSICondVector(srhs);

  // create full monolithic rhs vector
  rhs_ = Extractor().InsertStructureInteriorVector(sxi);  // structure interior
  Extractor().InsertStructureBoundaryVector(sxb, rhs_);   // structure boundary
//  rhs_->Scale(-1.);
  Extractor().InsertFluidVector(frhs,rhs_);               // fluid all

  // add condensed contribution from stress Lagrange multiplier
  const map<std::string, Teuchos::RCP<Epetra_Vector> > cvecs = FluidField().CouplingVectors();
  Teuchos::RCP<Epetra_Vector> rhsd = cvecs.find("rhsd")->second;

//  rhsd->Scale(-1.0);
  Extractor().AddVector(*FluidToStruct(rhsd),1,*rhs_);

  cout << *rhs_ << endl;

  // NOX expects a different sign here.
//  rhs_->Scale(-1.);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool FSI::MonolithicXFEM::ConverganceTest()
{
  cout << "FSI::MonolithicXFEM::ConverganceTest()" << endl;


  const double ittol = 1.0e-8;

  const RCP<const Epetra_Map> dofrowmap = Extractor().FullMap();
  Epetra_Vector full(*dofrowmap);
  Epetra_Import importer(*dofrowmap,rhs_->Map());

  double fullresnorm;

  const int err = full.Import(*rhs_,importer,Insert);
  if (err) dserror("Import using importer returned err=%d",err);
  full.Norm2(&fullresnorm);

  // test if rhs is zero (equilibrium)
  const bool converged = fullresnorm <= ittol;
  cout << "fullresnorm = " << fullresnorm << endl;

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

  Teuchos::RCP<LINALG::BlockSparseMatrixBase> s = StructureField().BlockSystemMatrix();
  Teuchos::RCP<LINALG::SparseMatrix>          f = FluidField().SystemMatrix();

  // direct access to coupling matrices
  const map<std::string, Teuchos::RCP<LINALG::SparseMatrix> > cmats = FluidField().CouplingMatrices();

//  Teuchos::RCP<LINALG::SparseMatrix> Cuu = cmats.find("Cuu")->second; // schon erledigt im Fluid
  Teuchos::RCP<LINALG::SparseMatrix> Mud = cmats.find("Mud")->second;
  Teuchos::RCP<LINALG::SparseMatrix> Mdu = cmats.find("Mdu")->second;
  Teuchos::RCP<LINALG::SparseMatrix> Cdd = cmats.find("Cdd")->second;


  /*----------------------------------------------------------------------*/

//  double scale     = FluidField().ResidualScaling();
//  double timescale = FluidField().TimeScaling();

  systemmatrix_ = rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(Extractor(),Extractor()));

  systemmatrix_->Assign(0,0,View,s->Matrix(0,0));
  systemmatrix_->Assign(1,0,View,s->Matrix(1,0));
  systemmatrix_->Assign(0,1,View,s->Matrix(0,1));
  systemmatrix_->Assign(1,1,View,s->Matrix(1,1));
  systemmatrix_->Assign(2,2,View,*f);

//  systemmatrix_->Assign(1,2,View,*Mdu);
//  systemmatrix_->Assign(2,1,View,*Mud);
//  cout << "alles assigned" << endl;
//  systemmatrix_->Matrix(1,1).Add(*Cdd,false,1.0,1.0);
//cout << "alles addiert" << endl;


//  matrix_->Matrix(1,0).Add(*Mud,false,1.0,0.0);
  sigtransform_(Mud->RowMap(),
                Mud->ColMap(),
                *Mud,
                1.0,
                ADAPTER::Coupling::SlaveConverter(coupsf),
                systemmatrix_->Matrix(2,1),
                true,
                false);

  //  matrix_->Matrix(0,0).Add(*Cdd,false,1.0,1.0);
  sggtransform_(*Cdd,
                1.0,
                ADAPTER::Coupling::SlaveConverter(coupsf),
                ADAPTER::Coupling::SlaveConverter(coupsf),
                systemmatrix_->Matrix(1,1),
                true,
                true);


  //  matrix_->Matrix(0,1).Add(*Mdu,false,1.0,0.0);
  sgitransform_(*Mdu,
                1.0,
                ADAPTER::Coupling::SlaveConverter(coupsf),
                systemmatrix_->Matrix(1,2));
//  cout << "alles addiert" << endl;


//  matrix_->Matrix(1,1).Add(*Cuu,false,1.0,1.0);

  // done. make sure all blocks are filled.
  systemmatrix_->Complete();
  cout << "systemmatrix_->Complete() ... Done!" << endl;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::InitialGuess(Teuchos::RCP<Epetra_Vector> ig)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicXFEM::InitialGuess");

//  SetupVector(*ig,
//              StructureField().InitialGuess(),
//              FluidField().InitialGuess(),
//              0.0);
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
//  Teuchos::RCP<const Epetra_Vector> fx = FluidField().InitialGuess();

  Teuchos::RCP<Epetra_Vector> iterinc = LINALG::CreateVector(*Extractor().FullMap(),true);

//  Teuchos::RCP<const Epetra_Vector> sxi = StructureField().Interface().ExtractOtherVector(StructureField().InitialGuess());
//  Teuchos::RCP<const Epetra_Vector> sxb = StructureField().Interface().ExtractFSICondVector(StructureField().InitialGuess());

//  Extractor().InsertStructureInteriorVector(sxi, iterinc);
//  Extractor().InsertStructureBoundaryVector(sxb, iterinc);
//  Extractor().InsertFluidVector(fx, iterinc);

  Teuchos::RCP<LINALG::SparseMatrix> m = systemmatrix_->Merge();

//  LINALG::PrintMatrixInMatlabFormat("monomatrix.txt", *m->EpetraMatrix(), true);

  // get UMFPACK...
  Teuchos::ParameterList solverparams = DRT::Problem::Instance()->FluidSolverParams();


  Teuchos::RCP<LINALG::Solver> solver =
      rcp(new LINALG::Solver(solverparams,
                             Comm(),
                             DRT::Problem::Instance()->ErrorFile()->Handle()));

//  cout << *rhs_ << endl;

  solver->Solve(m->EpetraOperator(), iterinc, rhs_, true, true);
//  state_.velnp_->Update(1.0,*incvel_,1.0);
//  cout << *x_ << endl;

  if (stepinc_ == Teuchos::null)
    stepinc_ = LINALG::CreateVector(*Extractor().FullMap(),true);

  stepinc_->Update(1.0,*iterinc, 1.0);

}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
//void FSI::MonolithicXFEM::SetupVector(Epetra_Vector &f,
//                                      Teuchos::RCP<const Epetra_Vector> sv,
//                                      Teuchos::RCP<const Epetra_Vector> fv,
//                                      double fluidscale)
//{
//  Extractor().InsertVector(*sv,0,f);
//  Extractor().InsertVector(*fv,1,f);
//}

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

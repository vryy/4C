#ifdef CCADISCRET

#include "fsi_monolithicxfem.H"
#include "fsi_overlapprec_fsiamg.H"
#include "fsi_statustest.H"
#include "fsi_nox_linearsystem_bgs.H"
#include "fsi_monolithic_linearsystem.H"
#include "../linalg/linalg_mapextractor.H"

#include "fsi_nox_group.H"
#include "fsi_nox_newton.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/inpar_fsi.H"

#include "../drt_io/io_control.H"
#include "../drt_fem_general/debug_nan.H"

#include <Teuchos_Time.hpp>
#include <Teuchos_TimeMonitor.hpp>

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
    fluidfield_(DRT::Problem::Instance()->FSIDynamicParams(),"FSICoupling"),
    cout0_(fluidfield_.Discretization()->Comm(), std::cout)
{
  // structure to fluid
  coupsf_.SetupConditionCoupling(*StructureField().Discretization(),
                                  StructureField().Interface().FSICondMap(),
                                 *FluidField().Discretization(),
                                  FluidField().Interface().FSICondMap(),
                                 "FSICoupling",
                                  genprob.ndim);

  // Use splitted structure matrix
  StructureField().UseBlockMatrix();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::SetDofRowMaps(const std::vector<Teuchos::RCP<const Epetra_Map> >& maps)
{
  Teuchos::RCP<Epetra_Map> fullmap = LINALG::MultiMapExtractor::MergeMapsKeepOrder(maps);
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

  SetTimeStep(FluidField().Time(),FluidField().Step());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicBaseXFEM::PrepareTimeStep()
{
  IncrementTimeAndStep();

  PrintHeader();

  StructureField().PrepareTimeStep();
  FluidField().    PrepareTimeStep();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicBaseXFEM::Update()
{
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
  StructureField().Output();
  FluidField().    Output();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::SetupExtractor()
{
  // setup global (full monolithic) map
  std::vector<Teuchos::RCP<const Epetra_Map> > maps;
  maps.push_back(StructureField().Interface().OtherMap());
  maps.push_back(StructureField().Interface().FSICondMap());
  maps.push_back(FluidField()    .DofRowMap());

  if (maps[0]->NumGlobalElements()==0)
  {
//    dserror("No inner structural equations. Splitting not possible. Panic.");
    cout0_ <<"No inner structural equations... All structure DOFs are surface DOFs, Axel?" << endl;
  }

  SetDofRowMaps(maps);

  Extractor().CheckForValidMapExtractor();
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
  cout0_ << "FSI::MonolithicXFEM::Timeloop()" << endl;

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
  cout0_ << "FSI::MonolithicXFEM::Newton()" << endl;

  const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
  const int itemax = fsidyn.get<int>("ITEMAX");

  SetupExtractor();

  stepinc_ = LINALG::CreateVector(*Extractor().FullMap(), true);
  iterinc_ = LINALG::CreateVector(*Extractor().FullMap(), true);

  for (int i =0;i<itemax;++i)
  {
    cout0_ << endl << YELLOW << "Newton step: " << i << "/" << itemax << END_COLOR << endl;

    Evaluate();

    // setup global (full monolithic) map and the individual extractors
    SetupExtractor();

    SetupRHS();

    if (Converged())
    {
      break;
    }

    SetupSystemMatrix();

    LinearSolve();
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::InitialGuess()
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicXFEM::InitialGuess");

  // get rhs from fields
//  const Teuchos::RCP<const Epetra_Vector> siterinc = StructureField().InitialGuess();
//  const Teuchos::RCP<const Epetra_Vector> fiterinc = FluidField().InitialGuess();
//
//  const Teuchos::RCP<const Epetra_Vector> siterincinterior = StructureField().Interface().ExtractOtherVector(siterinc);
//  const Teuchos::RCP<const Epetra_Vector> siterincboundary = StructureField().Interface().ExtractFSICondVector(siterinc);

  // create full monolithic stepinc vector
  iterinc_ = rcp(new Epetra_Vector(*Extractor().FullMap(), true));

//  Extractor().InsertStructureInteriorVector(siterincinterior, iterinc_);  // structure interior
//  Extractor().InsertStructureBoundaryVector(siterincboundary, iterinc_);  // structure boundary
//  Extractor().InsertFluidVector(fiterinc, iterinc_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::Evaluate()
{
  // take current increment and setup linear system
  cout0_ << "FSI::MonolithicXFEM::Evaluate()" << endl;

  // extract structure displacement
  const Teuchos::RCP<const Epetra_Vector> sxstepinc_interior = Extractor().ExtractStructureInteriorVector(stepinc_);
  const Teuchos::RCP<const Epetra_Vector> sxstepinc_boundary = Extractor().ExtractStructureBoundaryVector(stepinc_);

  const Teuchos::RCP<Epetra_Vector> sxstepinc = LINALG::CreateVector(*StructureField().Interface().FullMap(),true);
  StructureField().Interface().InsertOtherVector(sxstepinc_interior,sxstepinc);
  StructureField().Interface().InsertFSICondVector(sxstepinc_boundary,sxstepinc);

  const Teuchos::RCP<const Epetra_Vector> fxstepinc = Extractor().ExtractFluidVector(stepinc_);

  //    cout << "solid interior step inc" << endl;
  //    cout << *sxstepinc_interior << endl;
  //    cout << "solid boundary step inc" << endl;
  //    cout << *sxstepinc_boundary << endl;
  //    cout << "fluid stepinc:" << endl;
  //    cout << *fxstepinc << endl;

  cout0_ << BLUE2_LIGHT << "Solid evaluation" << END_COLOR << endl;
  StructureField().Evaluate(sxstepinc);

  // compute step inc including the structure Dirichlet BC
  sxstepinc->Update(1.0, *StructureField().Dispnp(), -1.0, *StructureField().Dispn(), 0.0);
  // get surface displacement step inc
  const Teuchos::RCP<const Epetra_Vector> stepinc_solid_boundary_2 =
      StructureField().Interface().ExtractFSICondVector(sxstepinc);
  // fluid requires step inc including Dirichlet conditions
  cout0_ << BLUE2_LIGHT << "Fluid evaluation" << END_COLOR << endl;
  FluidField().Evaluate(StructToFluid(stepinc_solid_boundary_2), fxstepinc);

  // put the FSI interface displacement into a text file for debugging
  if (stepinc_solid_boundary_2->Comm().MyPID() == 0 && stepinc_solid_boundary_2->MyLength() >= 3)
  {
    std::ofstream f;
    const std::string fname = DRT::Problem::Instance()->OutputControlFile()->FileName()
                                  + ".outifacedispstepinc.txt";
    f.open(fname.c_str(),std::fstream::ate | std::fstream::app);
    f << (*stepinc_solid_boundary_2)[0] << "  " << "\n";
    f.close();
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::SetupRHS()
{
  cout0_ << "FSI::MonolithicXFEM::SetupRHS" << endl;

  // get rhs from fields
  const Teuchos::RCP<const Epetra_Vector> srhs = StructureField().RHS();
  DRT::DEBUGGING::NaNChecker(*srhs);
  const Teuchos::RCP<const Epetra_Vector> frhs = FluidField().RHS();
  DRT::DEBUGGING::NaNChecker(*frhs);

  const Teuchos::RCP<const Epetra_Vector> srhsinterior = StructureField().Interface().ExtractOtherVector(srhs);
  const Teuchos::RCP<const Epetra_Vector> srhsboundary = StructureField().Interface().ExtractFSICondVector(srhs);

//  cout << "srhsinterior" << endl;
//  cout << *srhsinterior << endl;
//  cout << "srhsboundary" << endl;
//  cout << *srhsboundary << endl;
//  cout << "frhs" << endl;
//  cout << *frhs << endl;

  // create full monolithic rhs vector
  rhs_ = rcp(new Epetra_Vector(*Extractor().FullMap(), true));

  Extractor().InsertStructureInteriorVector(srhsinterior, rhs_);  // structure interior
  Extractor().InsertStructureBoundaryVector(srhsboundary, rhs_);  // structure boundary
  Extractor().InsertFluidVector(frhs, rhs_);                      // fluid all

  // add condensed contribution from stress Lagrange multiplier
  const Teuchos::RCP<const Epetra_Vector> rhsd = FluidField().CouplingVectors().find("rhsd")->second;
  const Teuchos::RCP<const Epetra_Vector> structboundarycoupleforce = FluidToStruct(rhsd);

//  cout << "rhsd^sigma + rhs_condense" << endl;
//  cout << *structboundarycoupleforce << endl;
  Extractor().AddStructureBoundaryVector(structboundarycoupleforce,rhs_);

  const Teuchos::RCP<const Epetra_Vector> zeros = LINALG::CreateVector(*Extractor().FullMap(), true);
  Teuchos::RCP<Epetra_Vector> tmp = LINALG::CreateVector(*Extractor().FullMap(), true);

  LINALG::ApplyDirichlettoSystem(tmp, rhs_, zeros, *CombinedDBCMap());

  DRT::DEBUGGING::NaNChecker(*rhs_);

}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool FSI::MonolithicXFEM::Converged()
{
  cout0_ << "FSI::MonolithicXFEM::ConverganceTest()" << endl;


  const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
  const double ittol = fsidyn.get<double>("CONVTOL");

  double fullresnorm;
  rhs_->Norm2(&fullresnorm);
  fullresnorm /= rhs_->Map().NumGlobalElements();

  // test if rhs is zero (equilibrium)
  const bool converged = fullresnorm <= ittol;
  if (converged)
    cout0_ << GREEN   << "fullresnorm = " << fullresnorm << " / " << ittol << END_COLOR << endl;
  else
    cout0_ << RED     << "fullresnorm = " << fullresnorm << " / " << ittol << END_COLOR << endl;

  // put the FSI interface displacement into a text file for debugging
  if (stepinc_->Comm().MyPID() == 0)
  {
    std::ofstream f;
    const std::string fname = DRT::Problem::Instance()->OutputControlFile()->FileName()
                                  + ".fullresnorm.txt";
    f.open(fname.c_str(),std::fstream::ate | std::fstream::app);
    f << fullresnorm << "  " << "\n";
    f.close();
  }

  return converged;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::SetupSystemMatrix()
{
  cout0_ << "FSI::MonolithicXFEM::SetupSystemMatrix()" << endl;
  // build global block matrix

  // extract Jacobian matrices and put them into composite system
  // matrix W

  const ADAPTER::Coupling& coupsf = StructureFluidCoupling();

  const Teuchos::RCP<LINALG::BlockSparseMatrixBase> s = StructureField().BlockSystemMatrix();
  const Teuchos::RCP<LINALG::SparseMatrix>          f = FluidField().SystemMatrix();

  // direct access to coupling matrices
  const map<std::string, Teuchos::RCP<LINALG::SparseMatrix> > cmats = FluidField().CouplingMatrices();

  // Cuu schon erledigt im Fluid
  const Teuchos::RCP<LINALG::SparseMatrix> Cud = cmats.find("Cud")->second;
  const Teuchos::RCP<LINALG::SparseMatrix> Cdu = cmats.find("Cdu")->second;
  const Teuchos::RCP<LINALG::SparseMatrix> Cdd = cmats.find("Cdd")->second;


  /*----------------------------------------------------------------------*/
  systemmatrix_ = rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(Extractor(),Extractor(),81,false));

  s->Matrix(1,1).UnComplete(); // da kommt noch ein coupling block dazu

  systemmatrix_->Assign(0,0,View,s->Matrix(0,0)); // interior
  systemmatrix_->Assign(1,0,View,s->Matrix(1,0));
  systemmatrix_->Assign(0,1,View,s->Matrix(0,1));
  systemmatrix_->Assign(1,1,View,s->Matrix(1,1)); // fsi boundary
  systemmatrix_->Assign(2,2,View,*f);

//  LINALG::PrintMatrixInMatlabFormat("S00",*s->Matrix(0,0).EpetraMatrix(),true);
//  LINALG::PrintMatrixInMatlabFormat("S10",*s->Matrix(1,0).EpetraMatrix(),true);
//  LINALG::PrintMatrixInMatlabFormat("S01",*s->Matrix(0,1).EpetraMatrix(),true);
//  LINALG::PrintMatrixInMatlabFormat("S11",*s->Matrix(1,1).EpetraMatrix(),true);

  sigtransform_(Cud->RowMap(),
                Cud->ColMap(),
                *Cud,
                1.0,
                ADAPTER::Coupling::SlaveConverter(coupsf),
                systemmatrix_->Matrix(2,1),
                false,
                false);

  sgitransform_(*Cdu,
                1.0,
                ADAPTER::Coupling::SlaveConverter(coupsf),
                systemmatrix_->Matrix(1,2),
                true);

  //  matrix_->Matrix(0,0).Add(*Cdd,false,1.0,1.0);
  sggtransform_(*Cdd,
                1.0,
                ADAPTER::Coupling::SlaveConverter(coupsf),
                ADAPTER::Coupling::SlaveConverter(coupsf),
                systemmatrix_->Matrix(1,1),
                false,
                true);

  // done. make sure all blocks are filled.
  systemmatrix_->Complete();
  // if DBC are on the interface, the DBC have to be applied again
  systemmatrix_->ApplyDirichlet(*CombinedDBCMap(), true);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::LinearSolve()
{
  cout0_ << "FSI::MonolithicXFEM::LinearSolve()" << endl;

  // get initial guess from fields
  // solve linear system

  // x^{n+1}_{0} = x^{n}
  //
  // Newton increment \Delta x^{n+1}_{i+1} = x^{n+1}_{i+1} - x^{n+1}_{i}

  //Teuchos::RCP<const Epetra_Vector> sx = StructureField().InitialGuess();
//  Teuchos::RCP<const Epetra_Vector> fx = FluidField().InitialGuess();

//  iterinc_ = LINALG::CreateVector(*Extractor().FullMap(),true);
  InitialGuess();

//  Teuchos::RCP<const Epetra_Vector> sxi = StructureField().Interface().ExtractOtherVector(StructureField().InitialGuess());
//  Teuchos::RCP<const Epetra_Vector> sxb = StructureField().Interface().ExtractFSICondVector(StructureField().InitialGuess());

//  Extractor().InsertStructureInteriorVector(sxi, iterinc);
//  Extractor().InsertStructureBoundaryVector(sxb, iterinc);
//  Extractor().InsertFluidVector(fx, iterinc);

  Teuchos::RCP<LINALG::SparseMatrix> m = systemmatrix_->Merge();
  cout0_ << "  merged" << endl;

  // get UMFPACK...
  Teuchos::ParameterList solverparams = DRT::Problem::Instance()->FluidSolverParams();

  Teuchos::RCP<LINALG::Solver> solver =
      rcp(new LINALG::Solver(solverparams,
                             Comm(),
                             DRT::Problem::Instance()->ErrorFile()->Handle()));

  solver->Solve(m->EpetraOperator(), iterinc_, rhs_, true, true);
  cout0_ << "  solved" << endl;

  if (stepinc_ == Teuchos::null)
    dserror("schimpf!");

  if (not stepinc_->Map().SameAs(iterinc_->Map()))
  {
    cout0_ << RED_LIGHT << "  Resetting global FSI stepinc_... " << END_COLOR << endl;
    VectorRescue(stepinc_);
    cout0_ << "  recovered" << endl;
  }

  const Teuchos::RCP<const Epetra_Vector> zeros = LINALG::CreateVector(*Extractor().FullMap(), true);
  Teuchos::RCP<Epetra_Vector> tmp = LINALG::CreateVector(*Extractor().FullMap(), false);
  LINALG::ApplyDirichlettoSystem(iterinc_, tmp, zeros, *CombinedDBCMap());
  tmp = Teuchos::null;
  cout0_ << "  DBC applied" << endl;

  stepinc_->Update(1.0,*iterinc_, 1.0);
  cout0_ << "  updated" << endl;

  double fulliterincnorm;
  iterinc_->Norm2(&fulliterincnorm);
  fulliterincnorm /= iterinc_->Map().NumGlobalElements();
  // test if rhs is zero (equilibrium)
  cout0_ << RED     << "fulliterincnorm = " << fulliterincnorm << END_COLOR << endl;

  // put the FSI interface displacement into a text file for debugging
  if (stepinc_->Comm().MyPID() == 0)
  {
    std::ofstream f;
    const std::string fname = DRT::Problem::Instance()->OutputControlFile()->FileName()
                                  + ".fulliterincnorm.txt";
    f.open(fname.c_str(),std::fstream::ate | std::fstream::app);
    f << fulliterincnorm << "  " << "\n";
    f.close();
  }

}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::VectorRescue(
    Teuchos::RCP<Epetra_Vector>& oldv
    ) const
{

  Teuchos::RCP<Epetra_Vector> newv = LINALG::CreateVector(*Extractor().FullMap(),true);
  // recover step vector as much as possible
  for (int newLID = 0; newLID < newv->Map().NumMyElements(); newLID++)
  {
    const int newGID = newv->Map().GID(newLID);
    const int oldLID = oldv->Map().LID(newGID);
    if (oldLID == -1)
    {
      // not found
    }
    else
    {
      (*newv)[newLID] = (*oldv)[oldLID];
    }
  }
  oldv = newv;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> FSI::MonolithicXFEM::CombinedDBCMap()
{
  const Teuchos::RCP<const Epetra_Map > scondmap = StructureField().GetDBCMapExtractor()->CondMap();
  const Teuchos::RCP<const Epetra_Map > fcondmap = FluidField().GetDBCMapExtractor()->CondMap();
  Teuchos::RCP<Epetra_Map> condmap = LINALG::MergeMap(scondmap, fcondmap, false);
  return condmap;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::MonolithicXFEM::StructToFluid(const Teuchos::RCP<const Epetra_Vector> iv) const
{
  const ADAPTER::Coupling& coupsf = StructureFluidCoupling();
  return coupsf.MasterToSlave(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::MonolithicXFEM::FluidToStruct(const Teuchos::RCP<const Epetra_Vector> iv) const
{
  const ADAPTER::Coupling& coupsf = StructureFluidCoupling();
  return coupsf.SlaveToMaster(iv);
}
#endif

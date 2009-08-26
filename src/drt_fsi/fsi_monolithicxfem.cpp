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
FSI::MonolithicXFEM::MonolithicXFEM(Epetra_Comm& comm)
  : MonolithicBaseXFEM(comm)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::Timeloop()
{
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

  Teuchos::RCP<const Epetra_Map> structuremap = StructureField().DofRowMap();
  vector<Teuchos::RCP<const Epetra_Map> > maps;
  maps.push_back(structuremap);


  Teuchos::RCP<const Epetra_Map> fullmap = LINALG::MultiMapExtractor::MergeMaps(maps);

  for (;;)
  {
    Evaluate(fullmap);

    SetupRHS();

    if (ConverganceTest())
    {
//      break;
    }


    SetupSystemMatrix();

    // solve
    LinearSolve();
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::Evaluate(Teuchos::RCP<const Epetra_Map> fullmapold)
{
  // take current increment and setup linear system
  cout << "FSI::MonolithicXFEM::Evaluate" << endl;

  // extract structure displacement
//  Teuchos::RCP<Epetra_Vector> sx = Extract...
//  StructureField().Evaluate(sx);

//  Teuchos::RCP<Epetra_Vector> fx = Extract...
  FluidField().Evaluate(null);



  std::vector<Teuchos::RCP<const Epetra_Map> > maps;
  maps.push_back(StructureField().DofRowMap());
  maps.push_back(FluidField().DofRowMap());
  Teuchos::RCP<Epetra_Map> fullmapnew = LINALG::MultiMapExtractor::MergeMaps(maps);


}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::SetupRHS()
{
  // get rhs from fields
  // setup global (full monolithic) map
  std::vector<Teuchos::RCP<const Epetra_Map> > maps;
  maps.push_back(StructureField().DofRowMap());
  maps.push_back(FluidField().DofRowMap());
  Teuchos::RCP<Epetra_Map> fullmap = LINALG::MultiMapExtractor::MergeMaps(maps);

  // create full monolithic rhs vector
  rhs_ = rcp(new Epetra_Vector(*fullmap, false));

  Teuchos::RCP<const Epetra_Vector> srhs = StructureField().RHS();
  Teuchos::RCP<const Epetra_Vector> frhs = FluidField().RHS();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool FSI::MonolithicXFEM::ConverganceTest()
{
  // test if rhs is zero (equilibrium)
  return false;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::SetupSystemMatrix()
{
  // build global block matrix

  Teuchos::RCP<LINALG::SparseMatrix> s = StructureField().SystemMatrix();
  Teuchos::RCP<LINALG::SparseMatrix> f = FluidField().SystemMatrix();
  cout << *s << endl;
  cout << *f << endl;

  /// direct access to system matrix
  map<std::string, Teuchos::RCP<LINALG::SparseMatrix> > cmats = FluidField().CouplingMatrices();

  Teuchos::RCP<LINALG::SparseMatrix> Mud = cmats.find("Mud")->second;
  Teuchos::RCP<LINALG::SparseMatrix> Mdu = cmats.find("Mdu")->second;
  Teuchos::RCP<LINALG::SparseMatrix> Cdd = cmats.find("Cdd")->second;

  cout << *Mud << endl;

  cout << *Cdd << endl;

  exit(0);

}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicXFEM::LinearSolve()
{
  // get initial guess from fields
  // solve linear system

  Teuchos::RCP<const Epetra_Vector> sx = StructureField().InitialGuess();
  //Teuchos::RCP<const Epetra_Vector> fx = FluidField().InitialGuess();
}


#endif

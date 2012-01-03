#ifdef CCADISCRET

#include <Teuchos_TimeMonitor.hpp>

#include "../drt_lib/drt_utils_createdis.H"
#include "../drt_io/io_control.H"
#include "../drt_fsi/fsi_partitionedmonolithic.H"
#include "../drt_fsi/fsi_monolithicfluidsplit.H"
#include "../drt_fsi/fsi_monolithiclagrange.H"
#include "../drt_fsi/fsi_monolithicstructuresplit.H"
#include "../drt_fsi/fsi_utils.H"
#include "../drt_lib/drt_condition_selector.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../linalg/linalg_utils.H"

#include "../drt_fsi/fs_monolithic.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/drt_validparameters.H"
#include "../drt_lib/drt_colors.H"

#include "../drt_scatra/scatra_utils.H"

#include "../drt_lib/drt_condition_utils.H"

#ifdef PARALLEL
#include <mpi.h>
#endif

#ifdef PARALLEL
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include "thermo_fsi.H"

extern struct _GENPROB     genprob;

#define SCATRABLOCKMATRIXMERGE


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FS3I::ThermoFSI::ThermoFSI(const Epetra_Comm& comm)
  :FS3I_Base(),
   comm_(comm)
{
  //---------------------------------------------------------------------
  // ensure correct order of three discretizations, with dof-numbering
  // such that structure dof < fluid dof < ale dofs
  // (ordering required at certain non-intuitive points)
  //---------------------------------------------------------------------
  RCP<DRT::Problem> problem = DRT::Problem::Instance();
  problem->Dis(genprob.numsf,0)->FillComplete();
  problem->Dis(genprob.numff,0)->FillComplete();
  problem->Dis(genprob.numaf,0)->FillComplete();
  problem->Dis(genprob.numscatra,0)->FillComplete();
  problem->Dis(genprob.numscatra,1)->FillComplete();

  //---------------------------------------------------------------------
  // create ale elements if not yet existing
  //---------------------------------------------------------------------
  RCP<DRT::Discretization> aledis = problem->Dis(genprob.numaf,0);
  if (aledis->NumGlobalNodes()==0)
  {
    RCP<DRT::Discretization> fluiddis = problem->Dis(genprob.numff,0);
    Teuchos::RCP<DRT::UTILS::DiscretizationCreator<FSI::UTILS::AleFluidCloneStrategy> > alecreator =
      Teuchos::rcp(new DRT::UTILS::DiscretizationCreator<FSI::UTILS::AleFluidCloneStrategy>() );
    alecreator->CreateMatchingDiscretization(fluiddis,aledis,-1);
  }
  //FSI::UTILS::CreateAleDiscretization();

  //---------------------------------------------------------------------
  // access discretizations for structure, fluid as well as fluid- and
  // structure-based scalar transport and get material map for scalar
  // transport elements
  //---------------------------------------------------------------------
  RefCountPtr<DRT::Discretization> structdis       = problem->Dis(0,0);
  RefCountPtr<DRT::Discretization> fluiddis        = problem->Dis(1,0);
  RefCountPtr<DRT::Discretization> fluidscatradis  = problem->Dis(3,0);
  RefCountPtr<DRT::Discretization> structscatradis = problem->Dis(3,1);

  std::map<std::pair<string,string>,std::map<int,int> > clonefieldmatmap = problem->ClonedMaterialMap();
  if (clonefieldmatmap.size() < 2)
    dserror("At least two material lists required for scalar transport in TFSI!");

  //---------------------------------------------------------------------
  // create discretization for fluid-based scalar transport from and
  // according to fluid discretization
  //---------------------------------------------------------------------
  if (fluiddis->NumGlobalNodes()==0) dserror("Fluid discretization is empty!");

  if (fluidscatradis->NumGlobalNodes()==0)
  {
    Epetra_Time time(comm);
    {
      Teuchos::RCP<DRT::UTILS::DiscretizationCreator<SCATRA::ScatraFluidCloneStrategy> > clonewizard =
        Teuchos::rcp(new DRT::UTILS::DiscretizationCreator<SCATRA::ScatraFluidCloneStrategy>() );
      std::pair<string,string> key("fluid","scatra1");
      std::map<int,int> fluidmatmap = clonefieldmatmap[key];
      clonewizard->CreateMatchingDiscretization(fluiddis,fluidscatradis,fluidmatmap);
    }
    if (comm.MyPID()==0) cout << "Created fluid-based scalar transport discretization from fluid discretization in...." << time.ElapsedTime() << " secs\n\n";
  }
  else dserror("Fluid AND ScaTra discretization present. This is not supported.");

  //---------------------------------------------------------------------
  // create discretization for structure-based scalar transport from and
  // according to structure discretization
  //---------------------------------------------------------------------
  if (structdis->NumGlobalNodes()==0) dserror("Structure discretization is empty!");

  if (structscatradis->NumGlobalNodes()==0)
  {
    Epetra_Time time(comm);
    {
      Teuchos::RCP<DRT::UTILS::DiscretizationCreator<SCATRA::ScatraFluidCloneStrategy> > clonewizard =
        Teuchos::rcp(new DRT::UTILS::DiscretizationCreator<SCATRA::ScatraFluidCloneStrategy>() );
      std::pair<string,string> key("structure","scatra2");
      std::map<int,int> structmatmap = clonefieldmatmap[key];
      clonewizard->CreateMatchingDiscretization(structdis,structscatradis,structmatmap);
    }
    if (comm.MyPID()==0) cout << "Created structure-based scalar transport discretization from structure discretization in...." << time.ElapsedTime() << " secs\n\n";
  }
  else dserror("Structure AND ScaTra discretization present. This is not supported.");

  //---------------------------------------------------------------------
  // get FSI coupling algorithm
  //---------------------------------------------------------------------
  const Teuchos::ParameterList& fsidyn = problem->FSIDynamicParams();
  int coupling = Teuchos::getIntegralValue<int>(fsidyn,"COUPALGO");
  switch (coupling)
  {
    case fsi_iter_monolithicfluidsplit:
    case fsi_iter_monolithicstructuresplit:
    {
      INPAR::FSI::LinearBlockSolver linearsolverstrategy = DRT::INPUT::IntegralValue<INPAR::FSI::LinearBlockSolver>(fsidyn,"LINEARBLOCKSOLVER");

      // call constructor for initialization of base class
      if (linearsolverstrategy==INPAR::FSI::PartitionedAitken or
          linearsolverstrategy==INPAR::FSI::PartitionedVectorExtrapolation or
          linearsolverstrategy==INPAR::FSI::PartitionedJacobianFreeNewtonKrylov)
      {
        fsi_ = Teuchos::rcp(new FSI::PartitionedMonolithic(comm));
      }
      else if (coupling==fsi_iter_monolithicfluidsplit)
      {
        fsi_ = Teuchos::rcp(new FSI::MonolithicFluidSplit(comm));
      }
      else if (coupling==fsi_iter_monolithicstructuresplit)
      {
        fsi_ = Teuchos::rcp(new FSI::MonolithicStructureSplit(comm));
      }
      else
      {
        dserror("Cannot find appropriate monolithic solver for coupling %d and linear strategy %d",coupling,linearsolverstrategy);
      }
      break;
    }
    default:
    {
      dserror("Unknown FSI coupling algorithm!");
    }
  }

  //---------------------------------------------------------------------
  // create instances for fluid- and structure-based scalar transport
  // solver and arrange them in combined vector
  //---------------------------------------------------------------------
  const Teuchos::ParameterList& scatradyn = problem->ScalarTransportDynamicParams();
  Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> fluidscatra =
    Teuchos::rcp(new ADAPTER::ScaTraBaseAlgorithm(scatradyn,true,0,problem->ScalarTransportFluidSolverParams()));
  Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> structscatra =
    Teuchos::rcp(new ADAPTER::ScaTraBaseAlgorithm(scatradyn,true,1,problem->ScalarTransportStructureSolverParams()));

  scatravec_.push_back(fluidscatra);
  scatravec_.push_back(structscatra);

  //---------------------------------------------------------------------
  // check various input parameters
  //---------------------------------------------------------------------
  // check time-integration scheme -> currently only one-step-theta scheme supported
  //const Teuchos::ParameterList& structdyn = problem->StructuralDynamicParams();
  //const Teuchos::ParameterList& fluiddyn  = problem->FluidDynamicParams();
  /*INPAR::SCATRA::TimeIntegrationScheme scatratimealgo =
    DRT::INPUT::IntegralValue<INPAR::SCATRA::TimeIntegrationScheme>(scatradyn,"TIMEINTEGR");
  INPAR::FLUID::TimeIntegrationScheme fluidtimealgo = fsi_->FluidAdapter().TimIntScheme();
  INPAR::STR::DynamicType structtimealgo = DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(structdyn,"DYNAMICTYP");

  if (scatratimealgo != INPAR::SCATRA::timeint_one_step_theta or
      fluidtimealgo != INPAR::FLUID::timeint_one_step_theta or
      structtimealgo != INPAR::STR::dyna_onesteptheta)
    dserror("For the time being, only one-step-theta time-integration scheme available for thermo-fluid-structure interaction!");

  // check whether time-step length, number of time steps and theta
  // are chosen identically
  if (scatradyn.get<double>("TIMESTEP") != fsidyn.get<double>("TIMESTEP") or
      scatradyn.get<int>("NUMSTEP")     != fsidyn.get<int>("NUMSTEP") or
      scatradyn.get<double>("THETA")    != fluiddyn.get<double>("THETA") or
      scatradyn.get<double>("THETA")    != structdyn.sublist("ONESTEPTHETA").get<double>("THETA"))
    dserror("One or more time-integration parameters do not match!");*/

  // check solver type -> needs to be incremental, otherwise residual and
  // stiffness matrix determined by scatra fields do not match formulation below
  if (scatravec_[0]->ScaTraField().Incremental() == false)
    dserror("Incremental formulation required for thermo-fluid-structure interaction!");

  //if (DRT::INPUT::IntegralValue<INPAR::SCATRA::ConvForm>(scatradyn,"CONVFORM") != INPAR::SCATRA::convform_conservative)
    //dserror("Conservative formulation needs to be chosen for structure-based scalar -> set CONVFORM to conservative!");

  //---------------------------------------------------------------------
  // check existence of scatra coupling condit. for both discretizations
  //---------------------------------------------------------------------
  std::vector<std::set<int> > condIDs;
  std::set<int> fluidIDs;
  std::set<int> structIDs;
  condIDs.push_back(fluidIDs);
  condIDs.push_back(structIDs);

  for (unsigned i=0; i<scatravec_.size(); ++i)
  {
    Teuchos::RCP<DRT::Discretization> dis = (scatravec_[i])->ScaTraField().Discretization();
    std::vector<DRT::Condition*> coupcond;
    dis->GetCondition("ScaTraCoupling",coupcond);

    for (unsigned iter=0; iter<coupcond.size(); ++iter)
    {
      int myID = (coupcond[iter])->GetInt("coupling id");
      condIDs[i].insert(myID);
    }
  }

  if (condIDs[0].size() != condIDs[1].size())
    dserror("ScaTra coupling conditions need to be defined for both discretizations!");
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::ThermoFSI::SetupSystem()
{
  //---------------------------------------------------------------------
  // set up FSI system
  //---------------------------------------------------------------------
  //fsi_->SetupSystem();

  //---------------------------------------------------------------------
  // create map extractors needed for scatra condition coupling
  //---------------------------------------------------------------------
  for (unsigned i=0; i<scatravec_.size(); ++i)
  {
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> currscatra = scatravec_[i];
    Teuchos::RCP<DRT::Discretization> currdis = currscatra->ScaTraField().Discretization();
    LINALG::MultiMapExtractor mapex;
    DRT::UTILS::MultiConditionSelector mcs;
    mcs.AddSelector(rcp(new DRT::UTILS::NDimConditionSelector(*currdis,"ScaTraCoupling",0,genprob.ndim)));
    mcs.SetupExtractor(*currdis,*currdis->DofRowMap(),mapex);
    scatrafieldexvec_.push_back(mapex);
  }

  scatracoup_.SetupConditionCoupling(*(scatravec_[0]->ScaTraField().Discretization()),
                                     scatrafieldexvec_[0].Map(1),
                                     *(scatravec_[1]->ScaTraField().Discretization()),
                                     scatrafieldexvec_[1].Map(1),
                                     "ScaTraCoupling",
                                     1);

  //---------------------------------------------------------------------
  // create map extractor for coupled scatra fields:
  // - Temperature values at interface from both sides are
  //   constrained to be equal.
  // - The fluid scatra dofs are kept, whereas the structure scatra dofs
  //   are condensed, which resembles the so-called "structure split" in
  //   a monolithic FSI system.
  //---------------------------------------------------------------------
  std::vector<Teuchos::RCP<const Epetra_Map> > maps;
  maps.push_back(scatrafieldexvec_[0].FullMap());
  maps.push_back(scatrafieldexvec_[1].Map(0));
  Teuchos::RCP<Epetra_Map> fullmap = LINALG::MultiMapExtractor::MergeMaps(maps);
  scatraglobalex_.Setup(*fullmap,maps);

  //---------------------------------------------------------------------
  // create scatra block matrix, rhs vector and incremental vector
  //---------------------------------------------------------------------
  scatrasystemmatrix_ =
    Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(scatraglobalex_,scatraglobalex_,54,false,true));
  scatrarhs_       = rcp(new Epetra_Vector(*scatraglobalex_.FullMap(),true));
  scatraincrement_ = rcp(new Epetra_Vector(*scatraglobalex_.FullMap(),true));

  //---------------------------------------------------------------------
  // check whether potential Dirichlet conditions at the interface are
  // defined for both discretizations
  //---------------------------------------------------------------------
  CheckInterfaceDirichletBC();

  //---------------------------------------------------------------------
  // create scatra solver
  //---------------------------------------------------------------------
  // get input parameters for coupled solver from respective section
  const Teuchos::ParameterList& coupledscatrasolvparams =
    DRT::Problem::Instance()->CoupledFluidAndScalarTransportSolverParams();

  // set and check solver type, which needs to be of Aztec type
  const int solvertype = DRT::INPUT::IntegralValue<INPAR::SOLVER::SolverType>(coupledscatrasolvparams,"SOLVER");
  if (solvertype != INPAR::SOLVER::aztec_msr) dserror("Aztec solver expected!");

  // set and check preconditioner type, which needs to be of block Gauss-Seidel type
  const int azprectype = DRT::INPUT::IntegralValue<INPAR::SOLVER::AzPrecType>(coupledscatrasolvparams,"AZPREC");
  if (azprectype != INPAR::SOLVER::azprec_BGS2x2)
    dserror("Block Gauss-Seidel preconditioner expected!");

  // create coupled scatra solver object
  // (first of two scatra discretizations required for this)
  Teuchos::RCP<DRT::Discretization> firstscatradis = (scatravec_[0])->ScaTraField().Discretization();
  scatrasolver_ = rcp(new LINALG::Solver(coupledscatrasolvparams,
                                         firstscatradis->Comm(),
                                         DRT::Problem::Instance()->ErrorFile()->Handle()));

  // get parameters for fluid- and structure-based scatra solvers
  scatrasolver_->PutSolverParamsToSubParams("Inverse1",DRT::Problem::Instance()->ScalarTransportFluidSolverParams());
  scatrasolver_->PutSolverParamsToSubParams("Inverse2",DRT::Problem::Instance()->ScalarTransportStructureSolverParams());

  // compute null space if necessary for fluid- and structure-based scatra solvers
  (scatravec_[0])->ScaTraField().Discretization()->ComputeNullSpaceIfNecessary(scatrasolver_->Params().sublist("Inverse1"));
  (scatravec_[1])->ScaTraField().Discretization()->ComputeNullSpaceIfNecessary(scatrasolver_->Params().sublist("Inverse2"));

}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::ThermoFSI::Timeloop()
{
  // output of initial state for scalars
  ScatraOutput();

  //fsi_->PrepareTimeloop();

  while (fsi_->NotFinished())
  {
    //DoFsiStep();
    DoScatraStep();
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::ThermoFSI::DoScatraStep()
{
  if (Comm().MyPID()==0)
    cout<<"\n*********************\n TEMPERATURE SOLVER \n*********************\n";

  // first scatra field: fluid-based, second scatra field: structure-based
  bool stopnonliniter=false;
  int itnum = 0;

  // prepare time step
  PrepareTimeStep();

  while (stopnonliniter==false)
  {
    itnum++;

    // set mesh displacement field for present iteration step
    SetMeshDisp();

    // set velocity fields from fluid and structure solution
    // for present iteration step
    SetVelocityFields();

    // prepare linear solve for both fluid- and structure-based scatra field
    for (unsigned i=0; i<scatravec_.size(); ++i)
    {
      Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra = scatravec_[i];
      scatra->ScaTraField().PrepareLinearSolve();
    }

    // set up matrix and right-hand-side for coupled scatra system
    SetupCoupledScatraMatrix();
    SetupCoupledScatraRHS();

    // set incremental vector to zero and solve coupled scatra system
    scatraincrement_->PutScalar(0.0);
    scatrasolver_->Solve(scatrasystemmatrix_->EpetraOperator(),
                         scatraincrement_,
                         scatrarhs_,
                         true,
                         true);

    // update for next iteration step
    IterUpdate();

    // check convergence and stop iteration loop if convergence is achieved
    stopnonliniter = ConvergenceCheck(itnum);
  }

  // update for next time step
  for (unsigned i=0; i<scatravec_.size(); ++i)
  {
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra = scatravec_[i];
    scatra->ScaTraField().Update();
  }

  // output of scalar transport fields
  ScatraOutput();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::ThermoFSI::SetupCoupledScatraMatrix()
{
  // define and check fluid- and structure-based scatra block matrices
  Teuchos::RCP<LINALG::SparseMatrix> scatra1 = scatravec_[0]->ScaTraField().SystemMatrix();
  Teuchos::RCP<LINALG::SparseMatrix> scatra2 = scatravec_[1]->ScaTraField().SystemMatrix();
  if (scatra1==Teuchos::null) dserror("Fluid-based scatra block matrix expected!");
  if (scatra2==Teuchos::null) dserror("Structure-based scatra block matrix expected!");

  // uncomplete system matrix for enabling slightly defective interface meshes
  scatra1->UnComplete();

  // split matrix into 2x2 blocks (boundary vs. inner dofs)
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> blockscatra2 = scatra2->Split<LINALG::DefaultBlockMatrixStrategy>(scatrafieldexvec_[1],scatrafieldexvec_[1]);
  blockscatra2->Complete();

  // assign structure-based scatra block matrix to lower right block
  scatrasystemmatrix_->Assign(1,1,View,blockscatra2->Matrix(0,0));

  // do transforms for upper right and lower left as well as lower right block
  sibtransform_(blockscatra2->FullRowMap(),
                blockscatra2->FullColMap(),
                blockscatra2->Matrix(0,1),
                1.0,
                ADAPTER::Coupling::SlaveConverter(scatracoup_),
                scatrasystemmatrix_->Matrix(1,0));
  sbitransform_(blockscatra2->Matrix(1,0),
                1.0,
                ADAPTER::Coupling::SlaveConverter(scatracoup_),
                scatrasystemmatrix_->Matrix(0,1));
  sbbtransform_(blockscatra2->Matrix(1,1),
                1.0,
                ADAPTER::Coupling::SlaveConverter(scatracoup_),
                ADAPTER::Coupling::SlaveConverter(scatracoup_),
                *scatra1,
                true,
                true);

  // assign fluid-based scatra block matrix to upper left block
  scatrasystemmatrix_->Assign(0,0,View,*scatra1);

  // complete matrix of coupled scatra system
  scatrasystemmatrix_->Complete();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::ThermoFSI::SetupCoupledScatraRHS()
{
  // define residual vectors for fluid- and structure-based scatra fields
  Teuchos::RCP<const Epetra_Vector> scatra1 = scatravec_[0]->ScaTraField().Residual();
  Teuchos::RCP<const Epetra_Vector> scatra2 = scatravec_[1]->ScaTraField().Residual();

  // extract inner (uncoupled) and boundary (coupled) dofs from
  // structure-based residual vector and put in separate vectors
  Teuchos::RCP<Epetra_Vector> scatra2_in = scatrafieldexvec_[1].ExtractVector(scatra2,0);
  Teuchos::RCP<Epetra_Vector> scatra2_bo = scatrafieldexvec_[1].ExtractVector(scatra2,1);

  // insert boundary dofs from structure-based residual vector in temporary
  // vector and add fluid-based residual vector
  Teuchos::RCP<Epetra_Vector> temp = scatrafieldexvec_[0].InsertVector(Scatra2ToScatra1(scatra2_bo),1);
  temp->Update(1.0,*scatra1,1.0);

  // insert temporary vector and boundary dofs from structure-based
  // residual vector into residual vector for coupled scatra system
  scatraglobalex_.InsertVector(*temp,0,*scatrarhs_);
  scatraglobalex_.InsertVector(*scatra2_in,1,*scatrarhs_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::ThermoFSI::IterUpdate()
{
  // define incremental vectors for fluid- and structure-based scatra
  // fields and extract respective vectors
  Teuchos::RCP<const Epetra_Vector> inc1;
  Teuchos::RCP<const Epetra_Vector> inc2;
  ExtractScatraFieldVectors(scatraincrement_,inc1,inc2);

  // update both fluid- and structure-based solution vectors
  scatravec_[0]->ScaTraField().UpdateIter(inc1);
  scatravec_[1]->ScaTraField().UpdateIter(inc2);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::ThermoFSI::ExtractScatraFieldVectors(
  Teuchos::RCP<const Epetra_Vector>  globalvec,
  Teuchos::RCP<const Epetra_Vector>& vec1,
  Teuchos::RCP<const Epetra_Vector>& vec2
)
{
  // extract fluid-based vector from global vector to obtain vector 1
  vec1 = scatraglobalex_.ExtractVector(globalvec,0);

  // extract boundary (coupled) dofs from fluid-based vector
  Teuchos::RCP<Epetra_Vector> vec1_bo = scatrafieldexvec_[0].ExtractVector(vec1,1);

  // extract inner (uncoupled) dofs of structure-based vector
  // from global vector and get boundary (coupled) dofs
  Teuchos::RCP<const Epetra_Vector> vec2_in = scatraglobalex_.ExtractVector(globalvec,1);
  Teuchos::RCP<Epetra_Vector> vec2_bo = Scatra1ToScatra2(vec1_bo);

  // insert inner and boundary dofs of structure-based vector
  // in temporary vector and generate vector 2
  Teuchos::RCP<Epetra_Vector> vec2_temp = scatrafieldexvec_[1].InsertVector(vec2_in,0);
  scatrafieldexvec_[1].InsertVector(vec2_bo,1,vec2_temp);
  vec2 = vec2_temp;
}


#endif

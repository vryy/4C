#ifdef CCADISCRET

#include <Teuchos_TimeMonitor.hpp>

#include "../drt_lib/drt_utils_createdis.H"

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

#include "gas_fsi.H"

extern struct _GENPROB     genprob;

#define SCATRABLOCKMATRIXMERGE


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FS3I::GasFSI::GasFSI(Epetra_Comm& comm)
  :FS3I_Base(),
   comm_(comm)
{
  RCP<DRT::Problem> problem = DRT::Problem::Instance();

  // make sure the three discretizations are filled in the right order
  // this creates dof numbers with
  //
  //       structure dof < fluid dof < ale dof
  //
  // We rely on this ordering in certain non-intuitive places!

  problem->Dis(genprob.numsf,0)->FillComplete();
  problem->Dis(genprob.numff,0)->FillComplete();
  problem->Dis(genprob.numaf,0)->FillComplete();
  problem->Dis(genprob.numscatra,0)->FillComplete();
  problem->Dis(genprob.numscatra,1)->FillComplete();

  // create ale elements if the ale discretization is empty
  RCP<DRT::Discretization> aledis = problem->Dis(genprob.numaf,0);
  if (aledis->NumGlobalNodes()==0)
  {
    RCP<DRT::Discretization> fluiddis = DRT::Problem::Instance()->Dis(genprob.numff,0);

    Teuchos::RCP<DRT::UTILS::DiscretizationCreator<FSI::UTILS::AleFluidCloneStrategy> > alecreator =
      Teuchos::rcp(new DRT::UTILS::DiscretizationCreator<FSI::UTILS::AleFluidCloneStrategy>() );

    alecreator->CreateMatchingDiscretization(fluiddis,aledis,-1);
  }
  //FSI::UTILS::CreateAleDiscretization();

  // access the fluid discretization
  RefCountPtr<DRT::Discretization> fluiddis = DRT::Problem::Instance()->Dis(1,0);
  // access the structure discretization
  RefCountPtr<DRT::Discretization> structdis = DRT::Problem::Instance()->Dis(0,0);
  // access the fluid scatra discretization
  RefCountPtr<DRT::Discretization> fluidscatradis = DRT::Problem::Instance()->Dis(3,0);
  // access the fluid structure discretization
  RefCountPtr<DRT::Discretization> structscatradis = DRT::Problem::Instance()->Dis(3,1);

  // get material map for the transport elements
  std::map<std::pair<string,string>,std::map<int,int> > clonefieldmatmap = DRT::Problem::Instance()->ClonedMaterialMap();
  if (clonefieldmatmap.size() < 2)
    dserror("at least 2 matlists needed for lung gas exchange");

  // FLUID SCATRA
  // we use the fluid discretization as layout for the scalar transport discretization
  if (fluiddis->NumGlobalNodes()==0) dserror("Fluid discretization is empty!");

  // create fluid scatra elements if the fluid scatra discretization is empty
  if (fluidscatradis->NumGlobalNodes()==0)
  {
    Epetra_Time time(comm);

    // create the fluid scatra discretization
    {
      Teuchos::RCP<DRT::UTILS::DiscretizationCreator<SCATRA::ScatraFluidCloneStrategy> > clonewizard =
        Teuchos::rcp(new DRT::UTILS::DiscretizationCreator<SCATRA::ScatraFluidCloneStrategy>() );

      std::pair<string,string> key("fluid","scatra1");
      std::map<int,int> fluidmatmap = clonefieldmatmap[key];

      clonewizard->CreateMatchingDiscretization(fluiddis,fluidscatradis,fluidmatmap);
    }
    if (comm.MyPID()==0)
      cout<<"Created scalar transport discretization from fluid field in...."
          <<time.ElapsedTime() << " secs\n\n";
  }
  else
    dserror("Fluid AND ScaTra discretization present. This is not supported.");

  // STRUCTURE SCATRA
  // we use the structure discretization as layout for the scalar transport discretization
  if (structdis->NumGlobalNodes()==0) dserror("Structure discretization is empty!");

  // create structure scatra elements if the structure scatra discretization is empty
  if (structscatradis->NumGlobalNodes()==0)
  {
    Epetra_Time time(comm);

    // create the structure scatra discretization
    {
      Teuchos::RCP<DRT::UTILS::DiscretizationCreator<SCATRA::ScatraFluidCloneStrategy> > clonewizard =
        Teuchos::rcp(new DRT::UTILS::DiscretizationCreator<SCATRA::ScatraFluidCloneStrategy>() );

      std::pair<string,string> key("structure","scatra2");
      std::map<int,int> structmatmap = clonefieldmatmap[key];

      clonewizard->CreateMatchingDiscretization(structdis,structscatradis,structmatmap);
    }
    if (comm.MyPID()==0)
      cout<<"Created scalar transport discretization from structure field in...."
          <<time.ElapsedTime() << " secs\n\n";
  }
  else
    dserror("Structure AND ScaTra discretization present. This is not supported.");

  const Teuchos::ParameterList& fsidyn   = problem->FSIDynamicParams();

  int coupling = Teuchos::getIntegralValue<int>(fsidyn,"COUPALGO");
  switch (coupling)
  {
  case fsi_iter_monolithicfluidsplit:
  case fsi_iter_monolithicstructuresplit:
  {
    INPAR::FSI::LinearBlockSolver linearsolverstrategy = DRT::INPUT::IntegralValue<INPAR::FSI::LinearBlockSolver>(fsidyn,"LINEARBLOCKSOLVER");

    // call constructor to initialise the base class
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
    dserror("Unknown coupling FSI algorithm");
  }
  }

  // access the problem-specific parameter lists
  const Teuchos::ParameterList& scatradyn = DRT::Problem::Instance()->ScalarTransportDynamicParams();
  const Teuchos::ParameterList& structdyn = DRT::Problem::Instance()->StructuralDynamicParams();
  const Teuchos::ParameterList& fluiddyn = DRT::Problem::Instance()->FluidDynamicParams();

  permeablesurf_ = DRT::INPUT::IntegralValue<int>(scatradyn,"PERMEABLESURF");

  // create one-way coupling algorithm instances
  Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> fluidscatra =
    Teuchos::rcp(new ADAPTER::ScaTraBaseAlgorithm(scatradyn,true,0,DRT::Problem::Instance()->ScalarTransportFluidSolverParams()));
  Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> structscatra =
    Teuchos::rcp(new ADAPTER::ScaTraBaseAlgorithm(scatradyn,true,1,DRT::Problem::Instance()->ScalarTransportStructureSolverParams()));

  scatravec_.push_back(fluidscatra);
  scatravec_.push_back(structscatra);

  /*----------------------------------------------------------------------*/
  /*                      Check of input parameters                       */
  /*----------------------------------------------------------------------*/

  // check time integration algo -> currently only one-step-theta scheme supported
  INPAR::SCATRA::TimeIntegrationScheme scatratimealgo =
    DRT::INPUT::IntegralValue<INPAR::SCATRA::TimeIntegrationScheme>(scatradyn,"TIMEINTEGR");
  INPAR::FLUID::TimeIntegrationScheme fluidtimealgo = fsi_->FluidAdapter().TimIntScheme();
  INPAR::STR::DynamicType structtimealgo =
    DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(structdyn,"DYNAMICTYP");

  if (scatratimealgo != INPAR::SCATRA::timeint_one_step_theta or
      fluidtimealgo != INPAR::FLUID::timeint_one_step_theta or
      structtimealgo != INPAR::STR::dyna_onesteptheta)
    dserror("lung gas exchange is limited in functionality (only one-step-theta scheme possible)");

  // check solver type -> it must be incremental, otherwise residual and
  // stiffness matrix determined by the scatra fields do not match the
  // formulation implemented below
  if (scatravec_[0]->ScaTraField().Incremental() == false)
    dserror("Incremental formulation needed for coupled lung scatra simulations");

  // make sure that initial time derivative of concentration is not calculated
  // automatically (i.e. field-wise)
  //if (DRT::INPUT::IntegralValue<int>(scatradyn,"SKIPINITDER")==false)
  //  dserror("Initial time derivative of phi must not be calculated automatically -> set SKIPINITDER to false");

  if (DRT::INPUT::IntegralValue<INPAR::SCATRA::ConvForm>(scatradyn,"CONVFORM") != INPAR::SCATRA::convform_conservative)
    dserror("Conservative formulation needs to be chosen for solids -> set CONVFORM to conservative!");

  // check if relevant parameters are chosen the same for FSI and ScaTra
  // dynamics
  if (scatradyn.get<double>("TIMESTEP") != fsidyn.get<double>("TIMESTEP") or
      scatradyn.get<int>("NUMSTEP") != fsidyn.get<int>("NUMSTEP") or
      scatradyn.get<double>("THETA") != fluiddyn.get<double>("THETA") or
      scatradyn.get<double>("THETA") != structdyn.sublist("ONESTEPTHETA").get<double>("THETA"))
    dserror("Fix your input file! Time integration parameters for FSI and ScaTra fields not matching!");

  // check if scatra coupling conditions are defined on both discretizations
  // and have the same permeability coefficient
  std::vector<std::set<int> > condIDs;
  std::set<int> fluidIDs;
  std::set<int> structIDs;
  condIDs.push_back(fluidIDs);
  condIDs.push_back(structIDs);
  std::vector<std::map<int,double> > PermCoeffs;
  std::map<int,double> fluidcoeff;
  std::map<int,double> structcoeff;
  PermCoeffs.push_back(fluidcoeff);
  PermCoeffs.push_back(structcoeff);

  for (unsigned i=0; i<scatravec_.size(); ++i)
  {
    Teuchos::RCP<DRT::Discretization> dis = (scatravec_[i])->ScaTraField().Discretization();
    std::vector<DRT::Condition*> coupcond;
    dis->GetCondition("ScaTraCoupling",coupcond);

    for (unsigned iter=0; iter<coupcond.size(); ++iter)
    {
      int myID = (coupcond[iter])->GetInt("coupling id");
      condIDs[i].insert(myID);

      if (permeablesurf_)
      {
        double myperm = (coupcond[iter])->GetDouble("permeability coefficient");
        PermCoeffs[i].insert(pair<int,double>(myID,myperm));
      }
    }
  }
  if (condIDs[0].size() != condIDs[1].size())
    dserror("ScaTra coupling conditions need to be defined on both discretizations");

  if (permeablesurf_)
  {
    std::map<int,double> fluid_PermCoeffs = PermCoeffs[0];
    std::map<int,double> struct_PermCoeffs = PermCoeffs[1];

    for (std::map<int,double>::iterator fit=fluid_PermCoeffs.begin(); fit!=fluid_PermCoeffs.end(); ++fit)
    {
      int ID = (*fit).first;
      double fluid_permcoef = (*fit).second;

      std::map<int,double>::iterator sit = struct_PermCoeffs.find(ID);
      if ((*sit).second != fluid_permcoef)
        dserror("Permeability coefficient of ScaTra interface needs to be the same in both conditions");
    }
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::GasFSI::ReadRestart()
{
  // read the restart information, set vectors and variables ---
  // be careful, dofmaps might be changed here in a Redistribute call
  if (genprob.restart)
  {
    fsi_->ReadRestart(genprob.restart);

    for (unsigned i=0; i<scatravec_.size(); ++i)
    {
      Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> currscatra = scatravec_[i];
      currscatra->ScaTraField().ReadRestart(genprob.restart);
    }
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::GasFSI::SetupSystem()
{
  // now do the coupling setup and create the combined dofmap
  fsi_->SetupSystem();

  /*----------------------------------------------------------------------*/
  /*                            General set up                            */
  /*----------------------------------------------------------------------*/

  // create map extractors needed for scatra condition coupling
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

  // create map extractor for coupled scatra fields
  // the second field (currently structure) is always split
  std::vector<Teuchos::RCP<const Epetra_Map> > maps;

  // If we do not consider the permeability of the interface between different
  // scatra fields, the concentrations on both sides of the interface are
  // constrained to be equal. In this case, we keep the fluid scatra dofs at the
  // interface as unknowns in the overall system, whereas the structure scatra
  // dofs are condensed (cf. "structuresplit" in a monolithic FSI
  // system). Otherwise, both concentrations are kept in the overall system
  // and the equality of fluxes is considered explicitly.
  if (!permeablesurf_)
  {
    maps.push_back(scatrafieldexvec_[0].FullMap());
    maps.push_back(scatrafieldexvec_[1].Map(0));
  }
  else
  {
    maps.push_back(scatrafieldexvec_[0].FullMap());
    maps.push_back(scatrafieldexvec_[1].FullMap());
  }
  Teuchos::RCP<Epetra_Map> fullmap = LINALG::MultiMapExtractor::MergeMaps(maps);
  scatraglobalex_.Setup(*fullmap,maps);

  // create coupling vectors and matrices (only needed when surface permeability is
  // considered)
  if (permeablesurf_)
  {
    for (unsigned i=0; i<scatravec_.size(); ++i)
    {
      Teuchos::RCP<Epetra_Vector> scatracoupforce =
      Teuchos::rcp(new Epetra_Vector(*(scatraglobalex_.Map(i)),true));
      scatracoupforce_.push_back(scatracoupforce);

      Teuchos::RCP<LINALG::SparseMatrix> scatracoupmat =
        Teuchos::rcp(new LINALG::SparseMatrix(*(scatraglobalex_.Map(i)),27,false,true));
      scatracoupmat_.push_back(scatracoupmat);

      Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra = scatravec_[i];
      const Epetra_Map* dofrowmap = scatra->ScaTraField().Discretization()->DofRowMap();
      Teuchos::RCP<Epetra_Vector> zeros = LINALG::CreateVector(*dofrowmap,true);
      scatrazeros_.push_back(zeros);
    }
  }

  // create scatra block matrix
  scatrasystemmatrix_ =
    Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(scatraglobalex_,
                                                                                   scatraglobalex_,
                                                                                   27,
                                                                                   false,
                                                                                   true));

  // create scatra rhs vector
  scatrarhs_ = rcp(new Epetra_Vector(*scatraglobalex_.FullMap(),true));

  // create scatra increment vector
  scatraincrement_ = rcp(new Epetra_Vector(*scatraglobalex_.FullMap(),true));

  // check whether potential Dirichlet conditions at the scatra interface are
  // defined on both discretizations
  CheckInterfaceDirichletBC();

  // scatra solver
  Teuchos::RCP<DRT::Discretization> firstscatradis = (scatravec_[0])->ScaTraField().Discretization();
#ifdef SCATRABLOCKMATRIXMERGE
  Teuchos::RCP<Teuchos::ParameterList> scatrasolvparams = rcp(new Teuchos::ParameterList);
  scatrasolvparams->set("solver","umfpack");
  scatrasolver_ = rcp(new LINALG::Solver(scatrasolvparams,
                                         firstscatradis->Comm(),
                                         DRT::Problem::Instance()->ErrorFile()->Handle()));
#else
  const Teuchos::ParameterList& coupledscatrasolvparams =
    DRT::Problem::Instance()->CoupledFluidAndScalarTransportSolverParams();
  const int solvertype =
    DRT::INPUT::IntegralValue<INPAR::SOLVER::SolverType>(coupledscatrasolvparams,"SOLVER");
  if (solvertype != INPAR::SOLVER::aztec_msr)
    dserror("aztec solver expected");
  const int azprectype =
    DRT::INPUT::IntegralValue<INPAR::SOLVER::AzPrecType>(coupledscatrasolvparams,"AZPREC");
  if (azprectype != INPAR::SOLVER::azprec_BGS2x2)
    dserror("Block Gauss-Seidel preconditioner expected");

  // use coupled SCATRA solver object
  scatrasolver_ = rcp(new LINALG::Solver(coupledscatrasolvparams,
                                         firstscatradis->Comm(),
                                         DRT::Problem::Instance()->ErrorFile()->Handle()));

  scatrasolver_->PutSolverParamsToSubParams("Inverse1",DRT::Problem::Instance()->ScalarTransportFluidSolverParams());
  scatrasolver_->PutSolverParamsToSubParams("Inverse2",DRT::Problem::Instance()->ScalarTransportStructureSolverParams());

  (scatravec_[0])->ScaTraField().Discretization()->ComputeNullSpaceIfNecessary(scatrasolver_->Params().sublist("Inverse1"));
  (scatravec_[1])->ScaTraField().Discretization()->ComputeNullSpaceIfNecessary(scatrasolver_->Params().sublist("Inverse2"));
#endif
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::GasFSI::Timeloop()
{
  // output of initial state
  ScatraOutput();

  fsi_->PrepareTimeloop();

  while (fsi_->NotFinished())
  {
    DoFsiStep();
    DoScatraStep();
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::GasFSI::DoFsiStep()
{
  fsi_->PrepareTimeStep();
  fsi_->TimeStep(fsi_);
  fsi_->PrepareOutput();
  fsi_->Update();
  fsi_->Output();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::GasFSI::DoScatraStep()
{
  if (Comm().MyPID()==0)
  {
    cout<<"\n***********************\n GAS TRANSPORT SOLVER \n***********************\n";
  }

  // first scatra field is associated with fluid, second scatra field is
  // associated with structure

  bool stopnonliniter=false;
  int itnum = 0;

  PrepareTimeStep();

  while (stopnonliniter==false)
  {
    SetMeshDisp();
    SetVelocityFields();

    EvaluateScatraFields();

    SetupCoupledScatraSystem();

    stopnonliniter = AbortScatraNonlinIter(itnum);
    if (stopnonliniter)
      break;

    LinearSolveScatra();
    FieldUpdateIter();

    itnum++;
  }

  UpdateScatraFields();

  ScatraOutput();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::GasFSI::PrepareTimeStep()
{
  SetMeshDisp();
  SetVelocityFields();

  for (unsigned i=0; i<scatravec_.size(); ++i)
  {
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra = scatravec_[i];
    scatra->ScaTraField().PrepareTimeStep();
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::GasFSI::EvaluateScatraFields()
{
  for (unsigned i=0; i<scatravec_.size(); ++i)
  {
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra_adap = scatravec_[i];
    SCATRA::ScaTraTimIntImpl& scatra = scatra_adap->ScaTraField();
    scatra.PrepareLinearSolve();

    // add contributions due to permeable surface/interface
    if (permeablesurf_)
    {
      Teuchos::RCP<Epetra_Vector> coupforce = scatracoupforce_[i];
      Teuchos::RCP<LINALG::SparseMatrix> coupmat = scatracoupmat_[i];

      coupforce->PutScalar(0.0);
      coupmat->Zero();

      // evaluate interface flux condition
      scatra.SurfacePermeability(coupmat,coupforce);

      // apply Dirichlet BC to coupling matrix and vector
      Teuchos::RCP<Epetra_Vector> zeros = scatrazeros_[i];
      const Teuchos::RCP<const LINALG::MapExtractor> dbcmapex = scatra.DirichMaps();
      const Teuchos::RCP< const Epetra_Map > dbcmap = dbcmapex->CondMap();
      coupmat->ApplyDirichlet(*dbcmap,false);
      LINALG::ApplyDirichlettoSystem(coupforce,zeros,*dbcmap);

    }
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::GasFSI::SetupCoupledScatraSystem()
{
  // set up scatra rhs
  SetupCoupledScatraRHS();

  // set up scatra system matrix
  SetupCoupledScatraMatrix();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::GasFSI::SetupCoupledScatraRHS()
{
  Teuchos::RCP<const Epetra_Vector> scatra1 = scatravec_[0]->ScaTraField().Residual();
  Teuchos::RCP<const Epetra_Vector> scatra2 = scatravec_[1]->ScaTraField().Residual();
  SetupCoupledScatraVector(scatrarhs_,scatra1,scatra2);

  // additional contributions in case of interface permeability
  if (permeablesurf_)
  {
    Teuchos::RCP<Epetra_Vector> coup1 = scatracoupforce_[0];
    Teuchos::RCP<Epetra_Vector> coup2 = scatracoupforce_[1];

    // contribution of the same field
    scatraglobalex_.AddVector(*coup1,0,*scatrarhs_,1.0);
    scatraglobalex_.AddVector(*coup2,1,*scatrarhs_,1.0);

    // contribution of the respective other field
    Teuchos::RCP<Epetra_Vector> coup1_boundary = scatrafieldexvec_[0].ExtractVector(coup1,1);
    Teuchos::RCP<Epetra_Vector> temp = scatrafieldexvec_[1].InsertVector(Scatra1ToScatra2(coup1_boundary),1);
    temp->Scale(-1.0);
    scatraglobalex_.AddVector(*temp,1,*scatrarhs_);

    Teuchos::RCP<Epetra_Vector> coup2_boundary = scatrafieldexvec_[1].ExtractVector(coup2,1);
    temp = scatrafieldexvec_[0].InsertVector(Scatra2ToScatra1(coup2_boundary),1);
    temp->Scale(-1.0);
    scatraglobalex_.AddVector(*temp,0,*scatrarhs_);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::GasFSI::SetupCoupledScatraVector(Teuchos::RCP<Epetra_Vector> globalvec,
                                               const Teuchos::RCP<const Epetra_Vector> vec1,
                                               const Teuchos::RCP<const Epetra_Vector> vec2)
{
  if (!permeablesurf_)
  {
    // concentrations are assumed to be equal at the interface

    // extract the inner (uncoupled) dofs from second field
    Teuchos::RCP<Epetra_Vector> vec2_other = scatrafieldexvec_[1].ExtractVector(vec2,0);

    Teuchos::RCP<Epetra_Vector> vec2_boundary = scatrafieldexvec_[1].ExtractVector(vec2,1);
    Teuchos::RCP<Epetra_Vector> temp = scatrafieldexvec_[0].InsertVector(Scatra2ToScatra1(vec2_boundary),1);
    temp->Update(1.0,*vec1,1.0);

    scatraglobalex_.InsertVector(*temp,0,*globalvec);
    scatraglobalex_.InsertVector(*vec2_other,1,*globalvec);
  }
  else
  {
    scatraglobalex_.InsertVector(*vec1,0,*globalvec);
    scatraglobalex_.InsertVector(*vec2,1,*globalvec);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::GasFSI::ExtractScatraFieldVectors(Teuchos::RCP<const Epetra_Vector> globalvec,
                                                Teuchos::RCP<const Epetra_Vector>& vec1,
                                                Teuchos::RCP<const Epetra_Vector>& vec2)
{
  if (!permeablesurf_)
  {
    // process fluid scatra unknowns
    vec1 = scatraglobalex_.ExtractVector(globalvec,0);

    // process structure scatra unknowns at the boundary
    Teuchos::RCP<Epetra_Vector> vec1_boundary = scatrafieldexvec_[0].ExtractVector(vec1,1);
    Teuchos::RCP<const Epetra_Vector> vec2_inner = scatraglobalex_.ExtractVector(globalvec,1);
    Teuchos::RCP<Epetra_Vector> vec2_boundary = Scatra1ToScatra2(vec1_boundary);

    Teuchos::RCP<Epetra_Vector> vec2_temp = scatrafieldexvec_[1].InsertVector(vec2_inner,0);
    scatrafieldexvec_[1].InsertVector(vec2_boundary,1,vec2_temp);
    vec2 = vec2_temp;
  }
  else
  {
    vec1 = scatraglobalex_.ExtractVector(globalvec,0);
    vec2 = scatraglobalex_.ExtractVector(globalvec,1);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::GasFSI::SetupCoupledScatraMatrix()
{
  Teuchos::RCP<LINALG::SparseMatrix> scatra1 = scatravec_[0]->ScaTraField().SystemMatrix();
  Teuchos::RCP<LINALG::SparseMatrix> scatra2 = scatravec_[1]->ScaTraField().SystemMatrix();

  if (scatra1==Teuchos::null)
    dserror("expect fluid scatra block matrix");
  if (scatra2==Teuchos::null)
    dserror("expect structure scatra block matrix");

  if (!permeablesurf_)
  {
    // Uncomplete system matrix to be able to deal with slightly defective
    // interface meshes.
    scatra1->UnComplete();

    // structure scatra
    // first split the matrix into 2x2 blocks (boundary vs. inner dofs)
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> blockscatra2 =
      scatra2->Split<LINALG::DefaultBlockMatrixStrategy>(scatrafieldexvec_[1],scatrafieldexvec_[1]);
    blockscatra2->Complete();

    scatrasystemmatrix_->Assign(1,1,View,blockscatra2->Matrix(0,0));

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

    // fluid scatra
    scatrasystemmatrix_->Assign(0,0,View,*scatra1);
  }
  else
  {
    // conventional contributions
    scatrasystemmatrix_->Assign(0,0,View,*scatra1);
    scatrasystemmatrix_->Assign(1,1,View,*scatra2);

    // additional contributions due to interface permeability (-> coupling terms)
    // contribution of the same field
    Teuchos::RCP<LINALG::SparseMatrix> coup1 = scatracoupmat_[0];
    Teuchos::RCP<LINALG::SparseMatrix> coup2 = scatracoupmat_[1];

    scatrasystemmatrix_->Matrix(0,0).Add(*coup1,false,1.0,1.0);
    scatrasystemmatrix_->Matrix(1,1).Add(*coup2,false,1.0,1.0);

    // contribution of the respective other field
    // first split the matrix into 2x2 blocks (boundary vs. inner dofs)
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> coupblock1
      = coup1->Split<LINALG::DefaultBlockMatrixStrategy>(scatrafieldexvec_[0],scatrafieldexvec_[0]);
    coupblock1->Complete();
    fbitransform_(coupblock1->Matrix(1,1),
                  -1.0,
                  ADAPTER::Coupling::MasterConverter(scatracoup_),
                  scatrasystemmatrix_->Matrix(1,0));

    Teuchos::RCP<LINALG::BlockSparseMatrixBase> coupblock2
      = coup2->Split<LINALG::DefaultBlockMatrixStrategy>(scatrafieldexvec_[1],scatrafieldexvec_[1]);
    coupblock2->Complete();
    sbitransform_(coupblock2->Matrix(1,1),
                  -1.0,
                  ADAPTER::Coupling::SlaveConverter(scatracoup_),
                  scatrasystemmatrix_->Matrix(0,1));
  }

  scatrasystemmatrix_->Complete();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::GasFSI::LinearSolveScatra()
{
  scatraincrement_->PutScalar(0.0);
  CoupledScatraSolve();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::GasFSI::FieldUpdateIter()
{
  Teuchos::RCP<const Epetra_Vector> inc1;
  Teuchos::RCP<const Epetra_Vector> inc2;

  ExtractScatraFieldVectors(scatraincrement_,inc1,inc2);

  scatravec_[0]->ScaTraField().UpdateIter(inc1);
  scatravec_[1]->ScaTraField().UpdateIter(inc2);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::GasFSI::UpdateScatraFields()
{
  for (unsigned i=0; i<scatravec_.size(); ++i)
  {
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra = scatravec_[i];
    scatra->ScaTraField().Update();
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::GasFSI::ScatraOutput()
{
  for (unsigned i=0; i<scatravec_.size(); ++i)
  {
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra = scatravec_[i];
    scatra->ScaTraField().Output();
    scatra->ScaTraField().OutputMeanScalars();
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::GasFSI::CoupledScatraSolve()
{
#ifdef SCATRABLOCKMATRIXMERGE
  Teuchos::RCP<LINALG::SparseMatrix> sparse = scatrasystemmatrix_->Merge();

  scatrasolver_->Solve(sparse->EpetraMatrix(),
                       scatraincrement_,
                       scatrarhs_,
                       true);
#else
  scatrasolver_->Solve(scatrasystemmatrix_->EpetraOperator(),
                       scatraincrement_,
                       scatrarhs_,
                       true,
                       true);
#endif
}

#endif

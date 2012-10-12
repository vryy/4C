/*!----------------------------------------------------------------------
\file fs3i_partitioned.cpp
\brief General algorithmic routines for partitioned solution approaches
       to fluid-structure-scalar-scalar interaction (FS3I), that is,
       algorithmic routines not specifically related to partitioned
       solution approaches to one -or two-way-coupled problem
       configurations, respectively

<pre>
Maintainers: Lena Yoshihara & Volker Gravemeier
             {yoshihara,vgravem}@lnm.mw.tum.de
             089/289-15303,-15245
</pre>

*----------------------------------------------------------------------*/


#include <Teuchos_TimeMonitor.hpp>

#include "../drt_lib/drt_utils_createdis.H"
#include "../drt_io/io_control.H"
#include "../drt_fsi/fsi_monolithic.H"
#include "../drt_fsi/fsi_monolithicfluidsplit.H"
#include "../drt_fsi/fsi_monolithicstructuresplit.H"
#include "../drt_fsi/fsi_utils.H"
#include "../drt_fsi/fsi_matrixtransform.H"
#include "../drt_lib/drt_condition_selector.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../linalg/linalg_utils.H"


#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/drt_validparameters.H"
#include "../drt_lib/drt_colors.H"

#include "../drt_adapter/adapter_coupling.H"

#include "../drt_scatra/passive_scatra_algorithm.H"
#include "../drt_scatra/scatra_utils_clonestrategy.H"

#include "../drt_ale/ale_utils_clonestrategy.H"

#include "../drt_lib/drt_condition_utils.H"

#include "fs3i_partitioned.H"


//#define SCATRABLOCKMATRIXMERGE

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FS3I::PartFS3I::PartFS3I(const Epetra_Comm& comm)
  : FS3I_Base(),
    comm_(comm)
{
  DRT::Problem* problem = DRT::Problem::Instance();

  //---------------------------------------------------------------------
  // read input parameters for FS3I problem
  //---------------------------------------------------------------------
  const Teuchos::ParameterList& fs3icontrol = problem->FS3IControlParams();
  dt_      = fs3icontrol.get<double>("TIMESTEP");
  numstep_ = fs3icontrol.get<int>("NUMSTEP");
  timemax_ = fs3icontrol.get<double>("MAXTIME");

  infperm_ = DRT::INPUT::IntegralValue<int>(fs3icontrol,"INF_PERM");

  //---------------------------------------------------------------------
  // set step and time
  //---------------------------------------------------------------------
  step_ = 0;
  time_ = 0.;

  //---------------------------------------------------------------------
  // ensure correct order of three discretizations, with dof-numbering
  // such that structure dof < fluid dof < ale dofs
  // (ordering required at certain non-intuitive points)
  //---------------------------------------------------------------------
  problem->GetDis("structure")->FillComplete();
  problem->GetDis("fluid")->FillComplete();
  problem->GetDis("ale")->FillComplete();
  problem->GetDis("scatra1")->FillComplete();
  problem->GetDis("scatra2")->FillComplete();

  //---------------------------------------------------------------------
  // access discretizations for structure, fluid, ale as well as fluid-
  // and structure-based scalar transport and get material map for fluid
  // and scalar transport elements
  //---------------------------------------------------------------------
  RefCountPtr<DRT::Discretization> fluiddis = problem->GetDis("fluid");
  RefCountPtr<DRT::Discretization> structdis = problem->GetDis("structure");
  RefCountPtr<DRT::Discretization> fluidscatradis = problem->GetDis("scatra1");
  RefCountPtr<DRT::Discretization> structscatradis = problem->GetDis("scatra2");
  RefCountPtr<DRT::Discretization> aledis = problem->GetDis("ale");

  //---------------------------------------------------------------------
  // create ale discretization as a clone from fluid discretization
  //---------------------------------------------------------------------
  if (aledis->NumGlobalNodes()==0)
  {
    DRT::UTILS::CloneDiscretization<ALE::UTILS::AleCloneStrategy>(fluiddis,aledis);
  }
  else
    dserror("Providing an ALE mesh is not supported for problemtype FS3I.");

  std::map<std::pair<string,string>,std::map<int,int> > clonefieldmatmap = problem->CloningMaterialMap();
  if (clonefieldmatmap.size() < 2)
    dserror("At least two material lists required for partitioned FS3I!");

  //---------------------------------------------------------------------
  // create discretization for fluid-based scalar transport from and
  // according to fluid discretization
  //---------------------------------------------------------------------
  if (fluiddis->NumGlobalNodes()==0) dserror("Fluid discretization is empty!");

  // create fluid-based scalar transport elements if fluid-based scalar
  // transport discretization is empty
  if (fluidscatradis->NumGlobalNodes()==0)
  {
    DRT::UTILS::CloneDiscretization<SCATRA::ScatraFluidCloneStrategy>(fluiddis,fluidscatradis);
  }
  else
    dserror("Fluid AND ScaTra discretization present. This is not supported.");

  //---------------------------------------------------------------------
  // create discretization for structure-based scalar transport from and
  // according to structure discretization
  //---------------------------------------------------------------------
  if (structdis->NumGlobalNodes()==0) dserror("Structure discretization is empty!");

  // create structure-based scalar transport elements if structure-based
  // scalar transport discretization is empty
  if (structscatradis->NumGlobalNodes()==0)
  {
    DRT::UTILS::CloneDiscretization<SCATRA::ScatraFluidCloneStrategy>(structdis,structscatradis);
  }
  else
    dserror("Structure AND ScaTra discretization present. This is not supported.");

  //---------------------------------------------------------------------
  // get FSI coupling algorithm
  //---------------------------------------------------------------------
  const Teuchos::ParameterList& fsidyn   = problem->FSIDynamicParams();
  int coupling = Teuchos::getIntegralValue<int>(fsidyn,"COUPALGO");
  switch (coupling)
  {
    case fsi_iter_monolithicfluidsplit:
    {
      fsi_ = Teuchos::rcp(new FSI::MonolithicFluidSplit(comm,fs3icontrol));
      break;
    }
    case fsi_iter_monolithicstructuresplit:
    {
      fsi_ = Teuchos::rcp(new FSI::MonolithicStructureSplit(comm,fs3icontrol));
      break;
    }
    default:
      dserror("Unknown coupling FSI algorithm");
  }

  //---------------------------------------------------------------------
  // create instances for fluid- and structure-based scalar transport
  // solver and arrange them in combined vector
  //---------------------------------------------------------------------
  // get the solver number used for structural ScalarTransport solver
  const int linsolver1number = fs3icontrol.get<int>("LINEAR_SOLVER1");
  // get the solver number used for structural ScalarTransport solver
    const int linsolver2number = fs3icontrol.get<int>("LINEAR_SOLVER2");
  // check if the TSI solver has a valid solver number
  if (linsolver1number == (-1))
    dserror("no linear solver defined for fluid ScalarTransport solver. Please set LINEAR_SOLVER2 in FS3I CONTROL to a valid number!");
  if (linsolver2number == (-1))
    dserror("no linear solver defined for structural ScalarTransport solver. Please set LINEAR_SOLVER2 in FS3I CONTROL to a valid number!");

  Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> fluidscatra =
    Teuchos::rcp(new ADAPTER::ScaTraBaseAlgorithm(fs3icontrol,true,"scatra1",problem->SolverParams(linsolver1number)));
  Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> structscatra =
    Teuchos::rcp(new ADAPTER::ScaTraBaseAlgorithm(fs3icontrol,true,"scatra2",DRT::Problem::Instance()->SolverParams(linsolver2number)));

  scatravec_.push_back(fluidscatra);
  scatravec_.push_back(structscatra);

  //---------------------------------------------------------------------
  // check various input parameters
  //---------------------------------------------------------------------
  const Teuchos::ParameterList& structdyn = problem->StructuralDynamicParams();
  const Teuchos::ParameterList& fluiddyn  = problem->FluidDynamicParams();
  const Teuchos::ParameterList& scatradyn = problem->ScalarTransportDynamicParams();

  // check consistency of time-integration schemes in input file
  // (including parameter theta itself in case of one-step-theta scheme)
  // and rule out unsupported versions of generalized-alpha time-integration
  // scheme (as well as other inappropriate schemes) for fluid subproblem
  INPAR::SCATRA::TimeIntegrationScheme scatratimealgo = DRT::INPUT::IntegralValue<INPAR::SCATRA::TimeIntegrationScheme>(scatradyn,"TIMEINTEGR");
  INPAR::FLUID::TimeIntegrationScheme fluidtimealgo = fsi_->FluidField().TimIntScheme();
  INPAR::STR::DynamicType structtimealgo = DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(structdyn,"DYNAMICTYP");
  if (fluidtimealgo  == INPAR::FLUID::timeint_one_step_theta)
  {
    if (scatratimealgo != INPAR::SCATRA::timeint_one_step_theta or
        structtimealgo != INPAR::STR::dyna_onesteptheta)
      dserror("Partitioned FS3I computations should feature consistent time-integration schemes for the subproblems; in this case, a one-step-theta scheme is intended to be used for the fluid subproblem, and different schemes are intended to be used for the structure and/or scalar transport subproblems!");

    if (scatradyn.get<double>("THETA") != fluiddyn.get<double>("THETA") or
        scatradyn.get<double>("THETA") != structdyn.sublist("ONESTEPTHETA").get<double>("THETA"))
    dserror("Parameter(s) theta for one-step-theta time-integration scheme defined in one or more of the individual fields do(es) not match for partitioned FS3I computation.");
  }
  else if (fluidtimealgo  == INPAR::FLUID::timeint_afgenalpha)
  {
    if (scatratimealgo != INPAR::SCATRA::timeint_gen_alpha or
        structtimealgo != INPAR::STR::dyna_genalpha)
      dserror("Partitioned FS3I computations should feature consistent time-integration schemes for the subproblems; in this case, a (alpha_f-based) generalized-alpha scheme is intended to be used for the fluid subproblem, and different schemes are intended to be used for the structure and/or scalar transport subproblems!");
  }
  else if (fluidtimealgo  == INPAR::FLUID::timeint_npgenalpha or
           fluidtimealgo  == INPAR::FLUID::timeint_gen_alpha)
  {
      dserror("Partitioned FS3I computations do not support n+1-based generalized-alpha time-integration schemes for the fluid subproblem!");
  }
  else if (fluidtimealgo  == INPAR::FLUID::timeint_bdf2 or
           fluidtimealgo  == INPAR::FLUID::timeint_stationary)
  {
      dserror("Partitioned FS3I computations do not support stationary of BDF2 time-integration schemes for the fluid subproblem!");
  }

  // check that incremental formulation is used for scalar transport field,
  // according to structure and fluid field
  if (scatravec_[0]->ScaTraField().IsIncremental() == false)
    dserror("Incremental formulation required for partitioned FS3I computations!");

  // ensure that initial time derivative of scalar is not calculated
  //if (DRT::INPUT::IntegralValue<int>(scatradyn,"SKIPINITDER")==false)
  //  dserror("Initial time derivative of phi must not be calculated automatically -> set SKIPINITDER to false");

  //---------------------------------------------------------------------
  // check existence of scatra coupling conditions for both
  // discretizations and definition of the permeability coefficient
  //---------------------------------------------------------------------
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

      if (!infperm_)
      {
        double myperm = (coupcond[iter])->GetDouble("permeability coefficient");
        PermCoeffs[i].insert(pair<int,double>(myID,myperm));
      }
    }
  }
  if (condIDs[0].size() != condIDs[1].size())
    dserror("ScaTra coupling conditions need to be defined on both discretizations");

  if (!infperm_)
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

  scatracoup_ = Teuchos::rcp(new ADAPTER::Coupling());
  scatraglobalex_ = Teuchos::rcp(new LINALG::MultiMapExtractor());
  sbbtransform_ = Teuchos::rcp(new FSI::UTILS::MatrixRowColTransform());
  sbitransform_ = Teuchos::rcp(new FSI::UTILS::MatrixRowTransform());
  sibtransform_ = Teuchos::rcp(new FSI::UTILS::MatrixColTransform());
  fbitransform_ = Teuchos::rcp(new FSI::UTILS::MatrixRowTransform());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I::ReadRestart()
{
  // read restart information, set vectors and variables
  // (Note that dofmaps might have changed in a redistribution call!)
  const int restart = DRT::Problem::Instance()->Restart();
  if (restart)
  {
    fsi_->ReadRestart(restart);

    for (unsigned i=0; i<scatravec_.size(); ++i)
    {
      Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> currscatra = scatravec_[i];
      currscatra->ScaTraField().ReadRestart(restart);
    }

    time_ = fsi_->FluidField().Time();
    step_ = fsi_->FluidField().Step();
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I::SetupSystem()
{
  // now do the coupling setup and create the combined dofmap
  fsi_->SetupSystem();

  /*----------------------------------------------------------------------*/
  /*                            General set up                            */
  /*----------------------------------------------------------------------*/

  // create map extractors needed for scatra condition coupling
  const int ndim = DRT::Problem::Instance()->NDim();
  for (unsigned i=0; i<scatravec_.size(); ++i)
  {
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> currscatra = scatravec_[i];
    Teuchos::RCP<DRT::Discretization> currdis = currscatra->ScaTraField().Discretization();
    Teuchos::RCP<LINALG::MultiMapExtractor> mapex = rcp(new LINALG::MultiMapExtractor());
    DRT::UTILS::MultiConditionSelector mcs;
    mcs.AddSelector(rcp(new DRT::UTILS::NDimConditionSelector(*currdis,"ScaTraCoupling",0,ndim)));
    mcs.SetupExtractor(*currdis,*currdis->DofRowMap(),*mapex);
    scatrafieldexvec_.push_back(mapex);
  }

  scatracoup_->SetupConditionCoupling(*(scatravec_[0]->ScaTraField().Discretization()),
                                     scatrafieldexvec_[0]->Map(1),
                                     *(scatravec_[1]->ScaTraField().Discretization()),
                                     scatrafieldexvec_[1]->Map(1),
                                     "ScaTraCoupling",
                                     1);

  // create map extractor for coupled scatra fields
  // the second field (currently structure) is always split
  std::vector<Teuchos::RCP<const Epetra_Map> > maps;

  // In the limiting case of an infinite permeability of the interface between
  // different scatra fields, the concentrations on both sides of the interface are
  // constrained to be equal. In this case, we keep the fluid scatra dofs at the
  // interface as unknowns in the overall system, whereas the structure scatra
  // dofs are condensed (cf. "structuresplit" in a monolithic FSI
  // system). Otherwise, both concentrations are kept in the overall system
  // and the equality of fluxes is considered explicitly.
  if (infperm_)
  {
    maps.push_back(scatrafieldexvec_[0]->FullMap());
    maps.push_back(scatrafieldexvec_[1]->Map(0));
  }
  else
  {
    maps.push_back(scatrafieldexvec_[0]->FullMap());
    maps.push_back(scatrafieldexvec_[1]->FullMap());
  }
  Teuchos::RCP<Epetra_Map> fullmap = LINALG::MultiMapExtractor::MergeMaps(maps);
  scatraglobalex_->Setup(*fullmap,maps);

  // create coupling vectors and matrices (only needed for finite surface permeabilities)
  if (!infperm_)
  {
    for (unsigned i=0; i<scatravec_.size(); ++i)
    {
      Teuchos::RCP<Epetra_Vector> scatracoupforce =
      Teuchos::rcp(new Epetra_Vector(*(scatraglobalex_->Map(i)),true));
      scatracoupforce_.push_back(scatracoupforce);

      Teuchos::RCP<LINALG::SparseMatrix> scatracoupmat =
        Teuchos::rcp(new LINALG::SparseMatrix(*(scatraglobalex_->Map(i)),27,false,true));
      scatracoupmat_.push_back(scatracoupmat);

      Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra = scatravec_[i];
      const Epetra_Map* dofrowmap = scatra->ScaTraField().Discretization()->DofRowMap();
      Teuchos::RCP<Epetra_Vector> zeros = LINALG::CreateVector(*dofrowmap,true);
      scatrazeros_.push_back(zeros);
    }
  }

  // create scatra block matrix
  scatrasystemmatrix_ =
    Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(*scatraglobalex_,
                                                                                   *scatraglobalex_,
                                                                                   27,
                                                                                   false,
                                                                                   true));

  // create scatra rhs vector
  scatrarhs_ = rcp(new Epetra_Vector(*scatraglobalex_->FullMap(),true));

  // create scatra increment vector
  scatraincrement_ = rcp(new Epetra_Vector(*scatraglobalex_->FullMap(),true));

  // check whether potential Dirichlet conditions at scatra interface are
  // defined for both discretizations
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
  const Teuchos::ParameterList& fs3icontrol = DRT::Problem::Instance()->FS3IControlParams();
  // get solver number used for fs3i
  const int linsolvernumber = fs3icontrol.get<int>("COUPLED_LINEAR_SOLVER");
  // check if LOMA solvers has a valid number
  if (linsolvernumber == (-1))
    dserror("no linear solver defined for FS3I problems. Please set COUPLED_LINEAR_SOLVER in FS3I CONTROL to a valid number!");

  const Teuchos::ParameterList& coupledscatrasolvparams =
    DRT::Problem::Instance()->SolverParams(linsolvernumber);
  const int solvertype =
    DRT::INPUT::IntegralValue<INPAR::SOLVER::SolverType>(coupledscatrasolvparams,"SOLVER");
  if (solvertype != INPAR::SOLVER::aztec_msr)
    dserror("aztec solver expected");
  const int azprectype =
    DRT::INPUT::IntegralValue<INPAR::SOLVER::AzPrecType>(coupledscatrasolvparams,"AZPREC");
  if (azprectype != INPAR::SOLVER::azprec_BGS2x2)
    dserror("Block Gauss-Seidel preconditioner expected");

  // use coupled scatra solver object
  scatrasolver_ = rcp(new LINALG::Solver(coupledscatrasolvparams,
                                         firstscatradis->Comm(),
                                         DRT::Problem::Instance()->ErrorFile()->Handle()));

  // get the solver number used for structural ScalarTransport solver
  const int linsolver1number = fs3icontrol.get<int>("LINEAR_SOLVER1");
  // get the solver number used for structural ScalarTransport solver
    const int linsolver2number = fs3icontrol.get<int>("LINEAR_SOLVER2");
  // check if the TSI solver has a valid solver number
  if (linsolver1number == (-1))
    dserror("no linear solver defined for fluid ScalarTransport solver. Please set LINEAR_SOLVER2 in FS3I CONTROL to a valid number!");
  if (linsolver2number == (-1))
    dserror("no linear solver defined for structural ScalarTransport solver. Please set LINEAR_SOLVER2 in FS3I CONTROL to a valid number!");

  scatrasolver_->PutSolverParamsToSubParams("Inverse1",DRT::Problem::Instance()->SolverParams(linsolver1number));
  scatrasolver_->PutSolverParamsToSubParams("Inverse2",DRT::Problem::Instance()->SolverParams(linsolver2number));

  (scatravec_[0])->ScaTraField().Discretization()->ComputeNullSpaceIfNecessary(scatrasolver_->Params().sublist("Inverse1"));
  (scatravec_[1])->ScaTraField().Discretization()->ComputeNullSpaceIfNecessary(scatrasolver_->Params().sublist("Inverse2"));
#endif
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I::TestResults(const Epetra_Comm& comm)
{
  DRT::Problem::Instance()->AddFieldTest(fsi_->FluidField().CreateFieldTest());
  DRT::Problem::Instance()->AddFieldTest(fsi_->StructureField()->CreateFieldTest());

  for (unsigned i=0; i<scatravec_.size(); ++i)
  {
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra = scatravec_[i];
    DRT::Problem::Instance()->AddFieldTest(scatra->CreateScaTraFieldTest());
  }
  DRT::Problem::Instance()->TestAll(comm);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I::SetFSISolution()
{
  SetMeshDisp();
  SetVelocityFields();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I::ScatraEvaluateSolveIterUpdate()
{
  EvaluateScatraFields();
  SetupCoupledScatraSystem();
  LinearSolveScatra();
  ScatraIterUpdate();
  // in case of later use of generalized-alpha time integration, a
  // routine for computing intermediate values is required at this point;
  // for the time being, this merely serves as a reminder for this
  // required inclusion
  //ComputeIntermediateValues();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I::EvaluateScatraFields()
{
  for (unsigned i=0; i<scatravec_.size(); ++i)
  {
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra_adap = scatravec_[i];
    SCATRA::ScaTraTimIntImpl& scatra = scatra_adap->ScaTraField();
    scatra.PrepareLinearSolve();

    // add contributions due to finite interface permeability
    if (!infperm_)
    {
      Teuchos::RCP<Epetra_Vector> coupforce = scatracoupforce_[i];
      Teuchos::RCP<LINALG::SparseMatrix> coupmat = scatracoupmat_[i];

      coupforce->PutScalar(0.0);
      coupmat->Zero();

      // evaluate interface flux condition
      scatra.SurfacePermeability(coupmat,coupforce);

      // apply Dirichlet boundary conditions to coupling matrix and vector
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
void FS3I::PartFS3I::SetupCoupledScatraSystem()
{
  // set up scatra rhs
  SetupCoupledScatraRHS();

  // set up scatra system matrix
  SetupCoupledScatraMatrix();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I::SetupCoupledScatraRHS()
{
  Teuchos::RCP<const Epetra_Vector> scatra1 = scatravec_[0]->ScaTraField().Residual();
  Teuchos::RCP<const Epetra_Vector> scatra2 = scatravec_[1]->ScaTraField().Residual();
  SetupCoupledScatraVector(scatrarhs_,scatra1,scatra2);

  // additional contributions in case of finite interface permeability
  if (!infperm_)
  {
    Teuchos::RCP<Epetra_Vector> coup1 = scatracoupforce_[0];
    Teuchos::RCP<Epetra_Vector> coup2 = scatracoupforce_[1];

    // contribution of the same field
    scatraglobalex_->AddVector(*coup1,0,*scatrarhs_,1.0);
    scatraglobalex_->AddVector(*coup2,1,*scatrarhs_,1.0);

    // contribution of the respective other field
    Teuchos::RCP<Epetra_Vector> coup1_boundary = scatrafieldexvec_[0]->ExtractVector(coup1,1);
    Teuchos::RCP<Epetra_Vector> temp = scatrafieldexvec_[1]->InsertVector(Scatra1ToScatra2(coup1_boundary),1);
    temp->Scale(-1.0);
    scatraglobalex_->AddVector(*temp,1,*scatrarhs_);

    Teuchos::RCP<Epetra_Vector> coup2_boundary = scatrafieldexvec_[1]->ExtractVector(coup2,1);
    temp = scatrafieldexvec_[0]->InsertVector(Scatra2ToScatra1(coup2_boundary),1);
    temp->Scale(-1.0);
    scatraglobalex_->AddVector(*temp,0,*scatrarhs_);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I::SetupCoupledScatraVector(Teuchos::RCP<Epetra_Vector>  globalvec,
                                              Teuchos::RCP<const Epetra_Vector>& vec1,
                                              Teuchos::RCP<const Epetra_Vector>& vec2)
{
  if (infperm_)
  {
    // concentrations are assumed to be equal at the interface
    // extract the inner (uncoupled) dofs from second field
    Teuchos::RCP<Epetra_Vector> vec2_other = scatrafieldexvec_[1]->ExtractVector(vec2,0);

    Teuchos::RCP<Epetra_Vector> vec2_boundary = scatrafieldexvec_[1]->ExtractVector(vec2,1);
    Teuchos::RCP<Epetra_Vector> temp = scatrafieldexvec_[0]->InsertVector(Scatra2ToScatra1(vec2_boundary),1);
    temp->Update(1.0,*vec1,1.0);

    scatraglobalex_->InsertVector(*temp,0,*globalvec);
    scatraglobalex_->InsertVector(*vec2_other,1,*globalvec);
  }
  else
  {
    scatraglobalex_->InsertVector(*vec1,0,*globalvec);
    scatraglobalex_->InsertVector(*vec2,1,*globalvec);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I::SetupCoupledScatraMatrix()
{
  Teuchos::RCP<LINALG::SparseMatrix> scatra1 = scatravec_[0]->ScaTraField().SystemMatrix();
  Teuchos::RCP<LINALG::SparseMatrix> scatra2 = scatravec_[1]->ScaTraField().SystemMatrix();

  if (scatra1==Teuchos::null)
    dserror("expect fluid scatra block matrix");
  if (scatra2==Teuchos::null)
    dserror("expect structure scatra block matrix");

  if (infperm_)
  {
    // Uncomplete system matrix to be able to deal with slightly defective
    // interface meshes.
    scatra1->UnComplete();

    // structure scatra
    // first split the matrix into 2x2 blocks (boundary vs. inner dofs)
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> blockscatra2 =
      scatra2->Split<LINALG::DefaultBlockMatrixStrategy>(*(scatrafieldexvec_[1]),*(scatrafieldexvec_[1]));
    blockscatra2->Complete();

    scatrasystemmatrix_->Assign(1,1,View,blockscatra2->Matrix(0,0));

    (*sibtransform_)(blockscatra2->FullRowMap(),
                     blockscatra2->FullColMap(),
                     blockscatra2->Matrix(0,1),
                     1.0,
                     ADAPTER::CouplingSlaveConverter(*scatracoup_),
                     scatrasystemmatrix_->Matrix(1,0));
    (*sbitransform_)(blockscatra2->Matrix(1,0),
                     1.0,
                     ADAPTER::CouplingSlaveConverter(*scatracoup_),
                     scatrasystemmatrix_->Matrix(0,1));
    (*sbbtransform_)(blockscatra2->Matrix(1,1),
                     1.0,
                     ADAPTER::CouplingSlaveConverter(*scatracoup_),
                     ADAPTER::CouplingSlaveConverter(*scatracoup_),
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
      = coup1->Split<LINALG::DefaultBlockMatrixStrategy>(*(scatrafieldexvec_[0]),*(scatrafieldexvec_[0]));
    coupblock1->Complete();
    (*fbitransform_)(coupblock1->Matrix(1,1),
                     -1.0,
                     ADAPTER::CouplingMasterConverter(*scatracoup_),
                     scatrasystemmatrix_->Matrix(1,0));

    Teuchos::RCP<LINALG::BlockSparseMatrixBase> coupblock2
      = coup2->Split<LINALG::DefaultBlockMatrixStrategy>(*(scatrafieldexvec_[1]),*(scatrafieldexvec_[1]));
    coupblock2->Complete();
    (*sbitransform_)(coupblock2->Matrix(1,1),
                     -1.0,
                     ADAPTER::CouplingSlaveConverter(*scatracoup_),
                     scatrasystemmatrix_->Matrix(0,1));
  }

  scatrasystemmatrix_->Complete();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I::LinearSolveScatra()
{
  scatraincrement_->PutScalar(0.0);

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


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I::ScatraIterUpdate()
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
void FS3I::PartFS3I::UpdateScatraFields()
{
  for (unsigned i=0; i<scatravec_.size(); ++i)
  {
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra = scatravec_[i];
    scatra->ScaTraField().Update();
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I::ScatraOutput()
{
  for (unsigned i=0; i<scatravec_.size(); ++i)
  {
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra = scatravec_[i];
    scatra->ScaTraField().Output();
  }
}


#ifdef CCADISCRET

#include <Teuchos_TimeMonitor.hpp>

#include "fsi_dyn.H"
#include "fsi_monolithicoverlap.H"
#include "fsi_monolithiclagrange.H"
#include "fsi_monolithicstructuresplit.H"
#include "fsi_utils.H"
#include "../drt_lib/drt_condition_selector.H"
#include "../drt_lib/drt_condition_utils.H"

#include "fs_monolithic.H"

#include "fsi_nox_aitken.H"
#include "fsi_nox_group.H"
#include "fsi_nox_newton.H"
#include "fsi_statustest.H"

#include <NOX_Epetra_Interface_Jacobian.H>
#include <NOX_Epetra_Interface_Preconditioner.H>
#include <NOX_Direction_UserDefinedFactory.H>

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

#include "fsi_lung_scatra.H"

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;


#define SCATRABLOCKMATRIXMERGE


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::LungScatra::LungScatra(Teuchos::RCP<FSI::Monolithic> fsi):
  fsi_(fsi)
{
  // access the problem-specific parameter lists
  const Teuchos::ParameterList& scatradyn = DRT::Problem::Instance()->ScalarTransportDynamicParams();
  const Teuchos::ParameterList& structdyn = DRT::Problem::Instance()->StructuralDynamicParams();

  // check time integration algo -> currently only one-step-theta scheme supported
  INPAR::SCATRA::TimeIntegrationScheme scatratimealgo = DRT::INPUT::IntegralValue<INPAR::SCATRA::TimeIntegrationScheme>(scatradyn,"TIMEINTEGR");
  INPAR::FLUID::TimeIntegrationScheme fluidtimealgo = fsi_->FluidAdapter().TimIntScheme();
  INPAR::STR::DynamicType structtimealgo = DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(structdyn,"DYNAMICTYP");

  if (scatratimealgo != INPAR::SCATRA::timeint_one_step_theta or
      fluidtimealgo != INPAR::FLUID::timeint_one_step_theta or
      structtimealgo != INPAR::STR::dyna_onesteptheta)
    dserror("lung gas exchange is limited in functionality (only one-step-theta scheme possible)");

  // create one-way coupling algorithm instances
  Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> fluidscatra = Teuchos::rcp(new ADAPTER::ScaTraBaseAlgorithm(scatradyn,true,0));
  Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> structscatra = Teuchos::rcp(new ADAPTER::ScaTraBaseAlgorithm(scatradyn,true,1));

  scatravec_.push_back(fluidscatra);
  scatravec_.push_back(structscatra);

  // check solver type -> it must be incremental, otherwise residual and
  // stiffness matrix determined by the scatra fields do not match the
  // formulation implemented below
  if (scatravec_[0]->ScaTraField().Incremental() == false)
    dserror("Incremental formulation needed for coupled lung scatra simulations");

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

  // if there are more than 2 scatra fields coupled, this needs to be done in
  // a loop
  scatracoup_.SetupConditionCoupling(*(scatravec_[0]->ScaTraField().Discretization()),
                                     scatrafieldexvec_[0].Map(1),
                                     *(scatravec_[1]->ScaTraField().Discretization()),
                                     scatrafieldexvec_[1].Map(1),
                                     "ScaTraCoupling",
                                     1);

  // create map extractor for coupled scatra fields
  // the second field (currently structure) is always split
  std::vector<Teuchos::RCP<const Epetra_Map> > maps;
  maps.push_back(scatrafieldexvec_[0].FullMap());
  maps.push_back(scatrafieldexvec_[1].Map(0));
  Teuchos::RCP<Epetra_Map> fullmap = LINALG::MultiMapExtractor::MergeMaps(maps);
  scatraglobalex_.Setup(*fullmap,maps);

  // create scatra block matrix
  scatrasystemmatrix_ = rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(scatraglobalex_,scatraglobalex_,81,false,true));

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
  dserror("currently only direct solution of merged scatra system possible");
#endif
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::LungScatra::Scatra2ToScatra1(Teuchos::RCP<Epetra_Vector> iv)
{
  return scatracoup_.SlaveToMaster(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::LungScatra::Scatra1ToScatra2(Teuchos::RCP<Epetra_Vector> iv)
{
  return scatracoup_.MasterToSlave(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungScatra::ReadRestart()
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
void FSI::LungScatra::SetupFSISystem()
{
  // now do the coupling setup and create the combined dofmap
  fsi_->SetupSystem();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungScatra::Timeloop()
{
  fsi_->PrepareTimeloop();

  while (fsi_->NotFinished())
  {
    DoFsiStep();
    DoScatraStep();
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungScatra::DoFsiStep()
{
  fsi_->PrepareTimeStep();
  fsi_->TimeStep(fsi_);
  fsi_->Update();
  fsi_->Output();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungScatra::DoScatraStep()
{
#ifdef PARALLEL
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  if (comm.MyPID()==0)
  {
    cout<<"\n***********************\n GAS TRANSPORT SOLVER \n***********************\n";
  }

  // first scatra field is associated with fluid, second scatra field is
  // associated with structure

  bool stopnonliniter=false;
  int itnum = 0;

  SetMeshDisp();
  SetVelocityFields();
  PrepareTimeStep();

  while (stopnonliniter==false)
  {
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

  UpdateAndOutput();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungScatra::PrepareTimeStep()
{
  SetVelocityFields();

  for (unsigned i=0; i<scatravec_.size(); ++i)
  {
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra = scatravec_[i];
    scatra->ScaTraField().PrepareTimeStep();
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungScatra::SetVelocityFields()
{
  std::vector<Teuchos::RCP<const Epetra_Vector> > vel;
  ExtractVel(vel);

  std::vector<Teuchos::RCP<DRT::Discretization> > discret;

  discret.push_back(fsi_->FluidAdapter().Discretization());
  discret.push_back(fsi_->StructureAdapter().Discretization());

  for (unsigned i=0; i<scatravec_.size(); ++i)
  {
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra = scatravec_[i];
    scatra->ScaTraField().SetVelocityField(vel[i],
                                           Teuchos::null,
                                           Teuchos::null,
                                           discret[i]);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungScatra::SetMeshDisp()
{
  // fluid field
  Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> fluidscatra = scatravec_[0];
  ADAPTER::Fluid& fluidadapter = fsi_->FluidAdapter();
  fluidscatra->ScaTraField().ApplyMeshMovement(fluidadapter.Dispnp(),
                                               fluidadapter.Discretization());

  // structure field
  Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> structscatra = scatravec_[1];
  ADAPTER::Structure& structadapter = fsi_->StructureAdapter();
  structscatra->ScaTraField().ApplyMeshMovement(structadapter.Dispnp(),
                                                structadapter.Discretization());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungScatra::EvaluateScatraFields()
{
  for (unsigned i=0; i<scatravec_.size(); ++i)
  {
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra = scatravec_[i];

    // this is the evaluation of the scatra field
    scatra->ScaTraField().PrepareLinearSolve();
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool FSI::LungScatra::AbortScatraNonlinIter(const int itnum)
{
#ifdef PARALLEL
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  int myrank = comm.MyPID();

  // some input parameters for the scatra fields
  const Teuchos::ParameterList& scatradyn = DRT::Problem::Instance()->ScalarTransportDynamicParams();
  const int itemax = scatradyn.sublist("NONLINEAR").get<int>("ITEMAX");
  const double ittol = scatradyn.sublist("NONLINEAR").get<double>("CONVTOL");
  const double abstolres = scatradyn.sublist("NONLINEAR").get<double>("ABSTOLRES");

  //----------------------------------------------------- compute norms
  double incconnorm_L2(0.0);
  scatraincrement_->Norm2(&incconnorm_L2);
  double conresnorm(0.0);
  scatrarhs_->Norm2(&conresnorm);
  // set up vector of absolute concentrations
  double connorm_L2(0.0);
  Teuchos::RCP<Epetra_Vector> con = rcp(new Epetra_Vector(scatraincrement_->Map()));
  Teuchos::RCP<const Epetra_Vector> scatra1 = scatravec_[0]->ScaTraField().Phinp();
  Teuchos::RCP<const Epetra_Vector> scatra2 = scatravec_[1]->ScaTraField().Phinp();
  SetupCoupledScatraVector(con,scatra1,scatra2);
  con->Norm2(&connorm_L2);

  // care for the case that nothing really happens in the concentration field
  if (connorm_L2 < 1e-5)
  {
    connorm_L2 = 1.0;
  }

  // absolute tolerance for deciding if residual is (already) zero
  // prevents additional solver calls that will not improve the residual anymore

  //-------------------------------------------------- output to screen
  /* special case of very first iteration step:
      - solution increment is not yet available
      - do not perform a solver call when the initial residuals are < EPS14*/
  if (itnum == 0)
  {
    if (myrank == 0)
    {
      printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E   |      --      |\n",itnum,itemax,ittol,conresnorm);
    }
    // abort iteration, when there's nothing more to do
    if (conresnorm < abstolres)
    {
      // print 'finish line'
      if (myrank == 0)
      {
        printf("+------------+-------------------+--------------+--------------+\n");
      }
      return true;
    }
    else
      return false;
  }
  /* ordinary case later iteration steps:
     - solution increment can be printed
     - convergence check should be done*/
  else
  {
    // print the screen info
    if (myrank == 0)
    {
      printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E   | %10.3E   |\n",itnum,itemax,ittol,conresnorm,incconnorm_L2/connorm_L2);
    }

    // this is the convergence check
    // We always require at least one solve. We test the L_2-norm of the
    // current residual. Norm of residual is just printed for information
    if (conresnorm <= ittol and incconnorm_L2/connorm_L2 <= ittol)
    {
      if (myrank == 0)
      {
        // print 'finish line'
        printf("+------------+-------------------+--------------+--------------+\n");
      }
      return true;
    }

    // abort iteration, when there's nothing more to do! -> more robustness
    else if (conresnorm < abstolres)
    {
      // print 'finish line'
      if (myrank == 0)
      {
        printf("+------------+-------------------+--------------+--------------+\n");
      }
      return true;
    }

    // warn if itemax is reached without convergence, but proceed to
    // next timestep...
    else if (itnum == itemax)
    {
      if (myrank == 0)
      {
        printf("+--------------------------------------------------------------+\n");
        printf("|        >>>>>> scatra not converged in itemax steps!          |\n");
        printf("+--------------------------------------------------------------+\n");
      }
      // yes, we stop the iteration
      return true;
    }
    else
      return false;
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungScatra::SetupCoupledScatraSystem()
{
  // set up scatra rhs
  Teuchos::RCP<const Epetra_Vector> scatra1 = scatravec_[0]->ScaTraField().Residual();
  Teuchos::RCP<const Epetra_Vector> scatra2 = scatravec_[1]->ScaTraField().Residual();
  SetupCoupledScatraVector(scatrarhs_,scatra1,scatra2);

  // set up scatra system matrix
  SetupCoupledScatraMatrix();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungScatra::SetupCoupledScatraVector(Teuchos::RCP<Epetra_Vector> globalvec,
                                               const Teuchos::RCP<const Epetra_Vector> vec1,
                                               const Teuchos::RCP<const Epetra_Vector> vec2)
{
  // extract the inner (uncoupled) dofs from second field
  Teuchos::RCP<Epetra_Vector> vec2_other = scatrafieldexvec_[1].ExtractVector(vec2,0);

  // add boundary dofs to first field (for the time being, concentrations are
  // assumed to be equal at the interface)
  Teuchos::RCP<Epetra_Vector> vec2_boundary = scatrafieldexvec_[1].ExtractVector(vec2,1);
  Teuchos::RCP<Epetra_Vector> temp = scatrafieldexvec_[0].InsertVector(Scatra2ToScatra1(vec2_boundary),1);
  temp->Update(1.0,*vec1,1.0);

  scatraglobalex_.InsertVector(*temp,0,*globalvec);
  scatraglobalex_.InsertVector(*vec2_other,1,*globalvec);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungScatra::ExtractScatraFieldVectors(Teuchos::RCP<const Epetra_Vector> globalvec,
                                                Teuchos::RCP<const Epetra_Vector>& vec1,
                                                Teuchos::RCP<const Epetra_Vector>& vec2)
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


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungScatra::SetupCoupledScatraMatrix()
{
  Teuchos::RCP<LINALG::SparseMatrix> scatra1 = scatravec_[0]->ScaTraField().SystemMatrix();
  Teuchos::RCP<LINALG::SparseMatrix> scatra2 = scatravec_[1]->ScaTraField().SystemMatrix();

  if (scatra1==Teuchos::null)
    dserror("expect fluid scatra block matrix");
  if (scatra2==Teuchos::null)
    dserror("expect structure scatra block matrix");

  // fluid scatra
  scatrasystemmatrix_->Assign(0,0,View,*scatra1);

  // structure scatra
  // first split the matrix into 2x2 blocks (boundary vs. inner dofs)
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> blockscatra2 = scatra2->Split<LINALG::DefaultBlockMatrixStrategy>(scatrafieldexvec_[1],scatrafieldexvec_[1]);
  blockscatra2->Complete();

  scatrasystemmatrix_->Assign(1,1,View,blockscatra2->Matrix(0,0));

  sibtransform_(blockscatra2->FullRowMap(),
                blockscatra2->FullColMap(),
                blockscatra2->Matrix(0,1),
                1.0,
                ADAPTER::Coupling::SlaveConverter(scatracoup_),
                scatrasystemmatrix_->Matrix(1,0));
  sbbtransform_(blockscatra2->Matrix(1,1),
                1.0,
                ADAPTER::Coupling::SlaveConverter(scatracoup_),
                ADAPTER::Coupling::SlaveConverter(scatracoup_),
                scatrasystemmatrix_->Matrix(0,0),
                true,
                true);
  sbitransform_(blockscatra2->Matrix(1,0),
                1.0,
                ADAPTER::Coupling::SlaveConverter(scatracoup_),
                scatrasystemmatrix_->Matrix(0,1));

  scatrasystemmatrix_->Complete();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungScatra::LinearSolveScatra()
{
  scatraincrement_->PutScalar(0.0);
  CoupledScatraSolve();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungScatra::FieldUpdateIter()
{
  Teuchos::RCP<const Epetra_Vector> inc1;
  Teuchos::RCP<const Epetra_Vector> inc2;

  ExtractScatraFieldVectors(scatraincrement_,inc1,inc2);

  scatravec_[0]->ScaTraField().UpdateIter(inc1);
  scatravec_[1]->ScaTraField().UpdateIter(inc2);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungScatra::UpdateAndOutput()
{
  for (unsigned i=0; i<scatravec_.size(); ++i)
  {
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra = scatravec_[i];
    scatra->ScaTraField().Update();
    scatra->ScaTraField().Output();
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungScatra::CoupledScatraSolve()
{
#ifdef SCATRABLOCKMATRIXMERGE
  Teuchos::RCP<LINALG::SparseMatrix> sparse = scatrasystemmatrix_->Merge();

  scatrasolver_->Solve(sparse->EpetraMatrix(),
                       scatraincrement_,
                       scatrarhs_,
                       true);
#else
  dserror("currently only direct solution of merged scatra system possible");
#endif
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungScatra::ExtractVel(std::vector<Teuchos::RCP<const Epetra_Vector> >& vel)
{
  // extract fluid velocities and accelerations

  ADAPTER::Fluid& fluid = fsi_->FluidAdapter();

  vel.push_back(fluid.ConvectiveVel());

  // extract structure velocities and accelerations

  ADAPTER::Structure& structure = fsi_->StructureAdapter();

  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();

  // major switch to different time integrators
  switch (DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(sdyn,"DYNAMICTYP"))
  {
  case INPAR::STR::dyna_centr_diff :
    dserror("no central differences in DRT");
    break;
  case INPAR::STR::dyna_gen_alfa :
  case INPAR::STR::dyna_genalpha :
  {
    Teuchos::RCP<Epetra_Vector> convel = rcp(new Epetra_Vector(*(structure.ExtractVelaf())));
    convel->Scale(-1.0);
    vel.push_back(convel);
    break;
  }
  case INPAR::STR::dyna_onesteptheta :
  {
    Teuchos::RCP<Epetra_Vector> convel = rcp(new Epetra_Vector(*(structure.ExtractVelnp())));
    convel->Scale(-1.0);
    vel.push_back(convel);
    break;
  }
  case INPAR::STR::dyna_Gen_EMM :
  case INPAR::STR::dyna_statics :
  case INPAR::STR::dyna_gen_alfa_statics :
  case INPAR::STR::dyna_gemm :
  case INPAR::STR::dyna_ab2:
  case INPAR::STR::dyna_euma :
  case INPAR::STR::dyna_euimsto :
  default :
  {
    dserror("structure time integration scheme not supported");
    break;
  }
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungScatra::CheckInterfaceDirichletBC()
{
#ifdef PARALLEL
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  Teuchos::RCP<DRT::Discretization> masterdis = scatravec_[0]->ScaTraField().Discretization();
  Teuchos::RCP<DRT::Discretization> slavedis = scatravec_[1]->ScaTraField().Discretization();

  Teuchos::RCP<const Epetra_Map> mastermap = scatracoup_.MasterDofMap();
  Teuchos::RCP<const Epetra_Map> permmastermap = scatracoup_.PermMasterDofMap();
  Teuchos::RCP<const Epetra_Map> slavemap = scatracoup_.SlaveDofMap();
  Teuchos::RCP<const Epetra_Map> permslavemap = scatracoup_.PermSlaveDofMap();

  const Teuchos::RCP<const LINALG::MapExtractor> masterdirichmapex = scatravec_[0]->ScaTraField().DirichMaps();
  const Teuchos::RCP<const Epetra_Map> masterdirichmap = masterdirichmapex->CondMap();

  // filter out master dirichlet dofs associated with the interface
  Teuchos::RCP<Epetra_Vector> masterifdirich = rcp(new Epetra_Vector(*mastermap,true));
  for (int i=0; i<mastermap->NumMyElements(); ++i)
  {
    int gid = mastermap->GID(i);
    if (masterdirichmap->MyGID(gid))
    {
      (*masterifdirich)[i] = 1.0;
    }
  }
  Teuchos::RCP<Epetra_Vector> test_slaveifdirich = scatracoup_.MasterToSlave(masterifdirich);

  const Teuchos::RCP<const LINALG::MapExtractor> slavedirichmapex = scatravec_[1]->ScaTraField().DirichMaps();
  const Teuchos::RCP<const Epetra_Map> slavedirichmap = slavedirichmapex->CondMap();

  // filter out slave dirichlet dofs associated with the interface
  Teuchos::RCP<Epetra_Vector> slaveifdirich = rcp(new Epetra_Vector(*slavemap,true));
  for (int i=0; i<slavemap->NumMyElements(); ++i)
  {
    int gid = slavemap->GID(i);
    if (slavedirichmap->MyGID(gid))
    {
      (*slaveifdirich)[i] = 1.0;
    }
  }
  Teuchos::RCP<Epetra_Vector> test_masterifdirich = scatracoup_.SlaveToMaster(slaveifdirich);

  // check if the locations of non-zero entries do not match
  for (int i=0; i<slavedis->DofRowMap()->NumMyElements(); ++i)
  {
    int gid = slavedis->DofRowMap()->GID(i);
    if (slavemap->MyGID(gid)) // in this case, the current dof is part of the interface
    {
      if ((*test_slaveifdirich)[slavemap->LID(gid)] == 1.0 and (*slaveifdirich)[slavemap->LID(gid)] != 1.0)
      {
        dserror("Dirichlet boundary conditions not matching at the interface");
      }
    }
  }

  for (int i=0; i<masterdis->DofRowMap()->NumMyElements(); ++i)
  {
    int gid = masterdis->DofRowMap()->GID(i);
    if (mastermap->MyGID(gid)) // in this case, the current dof is part of the interface
    {
      if ((*test_masterifdirich)[mastermap->LID(gid)] == 1.0 and (*masterifdirich)[mastermap->LID(gid)] != 1.0)
      {
        dserror("Dirichlet boundary conditions not matching at the interface");
      }
    }
  }
}

#endif

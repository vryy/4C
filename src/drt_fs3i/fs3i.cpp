/*!----------------------------------------------------------------------
\file fs3i.cpp
\brief cpp-file associated with general algorithmic routines for
       partitioned solution approaches to fluid-structure-scalar-scalar
       interaction (FS3I) and fluid-porous-structure-scalar-scalar
       interaction (FPS3I).

<pre>
Maintainers: Moritz Thon & Andre Hemmler
             thon@mhpc.mw.tum.de
             http://www.mhpc.mw.tum.de
             089 - 289-10364
</pre>

*----------------------------------------------------------------------*/


#include <Teuchos_TimeMonitor.hpp>

#include "../drt_fsi/fsi_dyn.H"
#include "../drt_fsi/fs_monolithic.H"
#include "../drt_fsi/fsi_monolithicfluidsplit.H"
#include "../drt_fsi/fsi_monolithicstructuresplit.H"
#include "../drt_fsi/fsi_matrixtransform.H"
#include "../drt_lib/drt_condition_selector.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_solver.H"
#include "../drt_fsi/fsi_utils.H"

#include "../drt_lib/drt_condition_selector.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/drt_validparameters.H"

#include "../drt_adapter/adapter_coupling.H"
#include "../drt_adapter/ad_str_fsiwrapper.H"

#include "../drt_scatra/scatra_algorithm.H"

#include "../drt_fluid/fluidimplicitintegration.H"
#include "../drt_fluid/fluid_utils.H"
#include "../drt_fluid/fluidresulttest.H"

#include "fs3i.H"



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FS3I::FS3I_Base::FS3I_Base()
{
  DRT::Problem* problem = DRT::Problem::Instance();
  //---------------------------------------------------------------------
  // read in and set private members of FS3I problem
  //---------------------------------------------------------------------
  const Teuchos::ParameterList& fs3icontrol = problem->FS3IControlParams();
  dt_      = fs3icontrol.get<double>("TIMESTEP");
  numstep_ = fs3icontrol.get<int>("NUMSTEP");
  timemax_ = fs3icontrol.get<double>("MAXTIME");

  infperm_ = DRT::INPUT::IntegralValue<int>(fs3icontrol,"INF_PERM");


  //---------------------------------------------------------------------
  // initialize step and time
  //---------------------------------------------------------------------
  step_ = 0;
  time_ = 0.;


  scatracoup_ = Teuchos::rcp(new ADAPTER::Coupling());
  scatraglobalex_ = Teuchos::rcp(new LINALG::MultiMapExtractor());
  sbbtransform_ = Teuchos::rcp(new FSI::UTILS::MatrixRowColTransform());
  sbitransform_ = Teuchos::rcp(new FSI::UTILS::MatrixRowTransform());
  sibtransform_ = Teuchos::rcp(new FSI::UTILS::MatrixColTransform());
  fbitransform_ = Teuchos::rcp(new FSI::UTILS::MatrixRowTransform());

}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::FS3I_Base::CheckInterfaceDirichletBC()
{
  Teuchos::RCP<DRT::Discretization> masterdis = scatravec_[0]->ScaTraField()->Discretization();
  Teuchos::RCP<DRT::Discretization> slavedis = scatravec_[1]->ScaTraField()->Discretization();

  Teuchos::RCP<const Epetra_Map> mastermap = scatracoup_->MasterDofMap();
  Teuchos::RCP<const Epetra_Map> permmastermap = scatracoup_->PermMasterDofMap();
  Teuchos::RCP<const Epetra_Map> slavemap = scatracoup_->SlaveDofMap();
  Teuchos::RCP<const Epetra_Map> permslavemap = scatracoup_->PermSlaveDofMap();

  const Teuchos::RCP<const LINALG::MapExtractor> masterdirichmapex = scatravec_[0]->ScaTraField()->DirichMaps();
  const Teuchos::RCP<const Epetra_Map> masterdirichmap = masterdirichmapex->CondMap();

  // filter out master dirichlet dofs associated with the interface
  Teuchos::RCP<Epetra_Vector> masterifdirich = Teuchos::rcp(new Epetra_Vector(*mastermap,true));
  for (int i=0; i<mastermap->NumMyElements(); ++i)
  {
    int gid = mastermap->GID(i);
    if (masterdirichmap->MyGID(gid))
    {
      (*masterifdirich)[i] = 1.0;
    }
  }
  Teuchos::RCP<Epetra_Vector> test_slaveifdirich = scatracoup_->MasterToSlave(masterifdirich);

  const Teuchos::RCP<const LINALG::MapExtractor> slavedirichmapex = scatravec_[1]->ScaTraField()->DirichMaps();
  const Teuchos::RCP<const Epetra_Map> slavedirichmap = slavedirichmapex->CondMap();

  // filter out slave dirichlet dofs associated with the interface
  Teuchos::RCP<Epetra_Vector> slaveifdirich = Teuchos::rcp(new Epetra_Vector(*slavemap,true));
  for (int i=0; i<slavemap->NumMyElements(); ++i)
  {
    int gid = slavemap->GID(i);
    if (slavedirichmap->MyGID(gid))
    {
      (*slaveifdirich)[i] = 1.0;
    }
  }
  Teuchos::RCP<Epetra_Vector> test_masterifdirich = scatracoup_->SlaveToMaster(slaveifdirich);

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
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::FS3I_Base::ScatraOutput()
{
  for (unsigned i=0; i<scatravec_.size(); ++i)
  {
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra = scatravec_[i];
    scatra->ScaTraField()->Output(i);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::FS3I_Base::IncrementTimeAndStep()
{
  step_ += 1;
  time_ += dt_;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::FS3I_Base::UpdateScatraFields()
{
  for (unsigned i=0; i<scatravec_.size(); ++i)
  {
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra = scatravec_[i];
    scatra->ScaTraField()->Update(i);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::FS3I_Base::ScatraEvaluateSolveIterUpdate()
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
void FS3I::FS3I_Base::EvaluateScatraFields()
{
  for (unsigned i=0; i<scatravec_.size(); ++i)
  {
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra_adap = scatravec_[i];
    Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra = scatra_adap->ScaTraField();
    scatra->PrepareLinearSolve();
    // add contributions due to finite interface permeability
    if (!infperm_)
    {

      Teuchos::RCP<Epetra_Vector> coupforce = scatracoupforce_[i];
      Teuchos::RCP<LINALG::SparseMatrix> coupmat = scatracoupmat_[i];

      coupforce->PutScalar(0.0);
      coupmat->Zero();

      scatra->SurfacePermeability(coupmat,coupforce);

      // apply Dirichlet boundary conditions to coupling matrix and vector
      Teuchos::RCP<Epetra_Vector> zeros = scatrazeros_[i];
      const Teuchos::RCP<const LINALG::MapExtractor> dbcmapex = scatra->DirichMaps();
      const Teuchos::RCP< const Epetra_Map > dbcmap = dbcmapex->CondMap();
      coupmat->ApplyDirichlet(*dbcmap,false);
      LINALG::ApplyDirichlettoSystem(coupforce,zeros,*dbcmap);

    }
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::FS3I_Base::SetupCoupledScatraSystem()
{
  // set up scatra rhs
  SetupCoupledScatraRHS();

  // set up scatra system matrix
  SetupCoupledScatraMatrix();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::FS3I_Base::SetupCoupledScatraRHS()
{
  Teuchos::RCP<const Epetra_Vector> scatra1 = scatravec_[0]->ScaTraField()->Residual();
  Teuchos::RCP<const Epetra_Vector> scatra2 = scatravec_[1]->ScaTraField()->Residual();
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
void FS3I::FS3I_Base::SetupCoupledScatraVector(Teuchos::RCP<Epetra_Vector>  globalvec,
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
void FS3I::FS3I_Base::SetupCoupledScatraMatrix()
{
  Teuchos::RCP<LINALG::SparseMatrix> scatra1 = scatravec_[0]->ScaTraField()->SystemMatrix();
  Teuchos::RCP<LINALG::SparseMatrix> scatra2 = scatravec_[1]->ScaTraField()->SystemMatrix();

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
Teuchos::RCP<Epetra_Vector> FS3I::FS3I_Base::Scatra2ToScatra1(Teuchos::RCP<Epetra_Vector> iv)
{
  return scatracoup_->SlaveToMaster(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FS3I::FS3I_Base::Scatra1ToScatra2(Teuchos::RCP<Epetra_Vector> iv)
{
  return scatracoup_->MasterToSlave(iv);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::FS3I_Base::LinearSolveScatra()
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
void FS3I::FS3I_Base::ScatraIterUpdate()
{
  // define incremental vectors for fluid- and structure-based scatra
  // fields and extract respective vectors
  Teuchos::RCP<const Epetra_Vector> inc1;
  Teuchos::RCP<const Epetra_Vector> inc2;
  ExtractScatraFieldVectors(scatraincrement_,inc1,inc2);

  // update both fluid- and structure-based solution vectors
  scatravec_[0]->ScaTraField()->UpdateIter(inc1);
  scatravec_[1]->ScaTraField()->UpdateIter(inc2);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::FS3I_Base::ExtractScatraFieldVectors(
  Teuchos::RCP<const Epetra_Vector>  globalvec,
  Teuchos::RCP<const Epetra_Vector>& vec1,
  Teuchos::RCP<const Epetra_Vector>& vec2
)
{
  if (infperm_)
  {
    // process fluid scatra unknowns
    vec1 = scatraglobalex_->ExtractVector(globalvec,0);

    // process structure scatra unknowns at the boundary
    Teuchos::RCP<Epetra_Vector> vec1_boundary = scatrafieldexvec_[0]->ExtractVector(vec1,1);
    Teuchos::RCP<const Epetra_Vector> vec2_inner = scatraglobalex_->ExtractVector(globalvec,1);
    Teuchos::RCP<Epetra_Vector> vec2_boundary = Scatra1ToScatra2(vec1_boundary);

    Teuchos::RCP<Epetra_Vector> vec2_temp = scatrafieldexvec_[1]->InsertVector(vec2_inner,0);
    scatrafieldexvec_[1]->InsertVector(vec2_boundary,1,vec2_temp);
    vec2 = vec2_temp;
  }
  else
  {
    vec1 = scatraglobalex_->ExtractVector(globalvec,0);
    vec2 = scatraglobalex_->ExtractVector(globalvec,1);
  }
}


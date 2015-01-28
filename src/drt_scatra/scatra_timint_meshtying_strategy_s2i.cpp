/*!----------------------------------------------------------------------
\file scatra_timint_meshtying_strategy_s2i.cpp

\brief Scatra-scatra interface coupling strategy for standard scalar transport problems

<pre>
Maintainer: Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15251
</pre>

*----------------------------------------------------------------------*/

#include "../drt_adapter/adapter_coupling.H"

#include "../drt_fluid/fluid_utils.H"

#include "../drt_fsi/fsi_matrixtransform.H"

#include "../drt_scatra/scatra_timint_implicit.H"

#include "../drt_scatra_ele/scatra_ele_action.H"

#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_utils.H"

#include "scatra_timint_meshtying_strategy_s2i.H"

/*----------------------------------------------------------------------*
 | constructor                                               fang 12/14 |
 *----------------------------------------------------------------------*/
SCATRA::MeshtyingStrategyS2I::MeshtyingStrategyS2I(
    SCATRA::ScaTraTimIntImpl* scatratimint
    ) :
MeshtyingStrategyBase(scatratimint),
maps_(Teuchos::null),
icoup_(Teuchos::null),
islavematrix_(Teuchos::null),
imastermatrix_(Teuchos::null),
islavetomastercoltransform_(Teuchos::null),
islavetomasterrowtransform_(Teuchos::null),
islavetomasterrowcoltransform_(Teuchos::null),
islaveresidual_(Teuchos::null),
imasterphinp_(Teuchos::null)
{
  return;
} // SCATRA::MeshtyingStrategyS2I::MeshtyingStrategyS2I


/*-----------------------------------------------------------------------*
 | evaluate scatra-scatra interface coupling conditions       fang 10/14 |
 *-----------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2I::EvaluateMeshtying() const
{
  // time measurement: evaluate condition 'S2ICoupling'
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:       + evaluate condition 'S2ICoupling'");

  // check matrix
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> blocksystemmatrix = scatratimint_->BlockSystemMatrix();
  if(blocksystemmatrix == Teuchos::null)
    dserror("System matrix is not a block matrix!");

  // create parameter list
  Teuchos::ParameterList condparams;

  // action for elements
  condparams.set<int>("action",SCATRA::bd_calc_s2icoupling);

  // further parameters
  condparams.set<int>("scatratype",scatratimint_->ScaTraType());
  condparams.set("isale",scatratimint_->IsALE());

  // set global and interface state vectors according to time-integration scheme
  scatratimint_->Discretization()->ClearState();
  scatratimint_->AddTimeIntegrationSpecificVectors();

  // fill interface state vector imasterphinp_ with transformed master dof values and add to discretization
  maps_->InsertVector(icoup_->MasterToSlave(maps_->ExtractVector(*(scatratimint_->Phiafnp()),2)),1,imasterphinp_);
  scatratimint_->Discretization()->SetState("imasterphinp",imasterphinp_);

  // evaluate scatra-scatra interface coupling at time t_{n+1} or t_{n+alpha_F}
  islavematrix_->Zero();
  imastermatrix_->Zero();
  islaveresidual_->PutScalar(0.);
  scatratimint_->Discretization()->EvaluateCondition(condparams,islavematrix_,imastermatrix_,islaveresidual_,Teuchos::null,Teuchos::null,"S2ICouplingSlave");
  scatratimint_->Discretization()->ClearState();

  // assemble linearizations of slave fluxes w.r.t. slave dofs into global system matrix
  islavematrix_->Complete();
  blocksystemmatrix->Matrix(1,1).Add(*islavematrix_,false,1.,1.);

  // transform linearizations of slave fluxes w.r.t. master dofs and assemble into global system matrix
  imastermatrix_->Complete();
  (*islavetomastercoltransform_)(imastermatrix_->RowMap(),imastermatrix_->ColMap(),*imastermatrix_,1.,
      ADAPTER::CouplingSlaveConverter(*icoup_),blocksystemmatrix->Matrix(1,2));

  // derive linearizations of master fluxes w.r.t. slave dofs and assemble into global system matrix
  (*islavetomasterrowtransform_)(*islavematrix_,-1.,ADAPTER::CouplingSlaveConverter(*icoup_),blocksystemmatrix->Matrix(2,1));

  // derive linearizations of master fluxes w.r.t. master dofs and assemble into global system matrix
  (*islavetomasterrowcoltransform_)(*imastermatrix_,-1.,ADAPTER::CouplingSlaveConverter(*icoup_),ADAPTER::CouplingSlaveConverter(*icoup_),blocksystemmatrix->Matrix(2,2),true,true);

  // assemble slave residuals into global residual vector
  maps_->AddVector(islaveresidual_,1,scatratimint_->Residual());

  // transform master residuals and assemble into global residual vector
  maps_->AddVector(icoup_->SlaveToMaster(islaveresidual_),2,scatratimint_->Residual(),-1.);

  return;
} // SCATRA::MeshtyingStrategyS2I::EvaluateMeshtying


/*----------------------------------------------------------------------*
 | perform setup of scatra-scatra interface coupling         fang 10/14 |
 *----------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2I::InitMeshtying()
{
  // extract scatra-scatra coupling conditions from discretization
  std::vector<DRT::Condition*> slaveconditions;
  scatratimint_->Discretization()->GetCondition("S2ICouplingSlave", slaveconditions);
  std::vector<DRT::Condition*> masterconditions;
  scatratimint_->Discretization()->GetCondition("S2ICouplingMaster", masterconditions);

  // initialize int vectors for global ids of slave and master interface nodes
  std::vector<int> islavenodegidvec;
  std::vector<int> imasternodegidvec;

  // fill vectors
  for (unsigned islavecondition=0; islavecondition<slaveconditions.size(); ++islavecondition)
  {
    const std::vector<int>* islavenodegids = slaveconditions[islavecondition]->Nodes();

    for (unsigned islavenode=0; islavenode<islavenodegids->size(); ++islavenode)
    {
      const int islavenodegid = (*islavenodegids)[islavenode];

      // insert global id of current node into associated vector only if node is owned by current processor
      // need to make sure that node is stored on current processor, otherwise cannot resolve "->Owner()"
      if(scatratimint_->Discretization()->HaveGlobalNode(islavenodegid) and scatratimint_->Discretization()->gNode(islavenodegid)->Owner() == scatratimint_->Discretization()->Comm().MyPID())
        islavenodegidvec.push_back(islavenodegid);
    }
  }
  for (unsigned imastercondition=0; imastercondition<masterconditions.size(); ++imastercondition)
  {
    const std::vector<int>* imasternodegids = masterconditions[imastercondition]->Nodes();

    for (unsigned imasternode=0; imasternode<imasternodegids->size(); ++imasternode)
    {
      const int imasternodegid = (*imasternodegids)[imasternode];

      // insert global id of current node into associated vector only if node is owned by current processor
      // need to make sure that node is stored on current processor, otherwise cannot resolve "->Owner()"
      if(scatratimint_->Discretization()->HaveGlobalNode(imasternodegid) and scatratimint_->Discretization()->gNode(imasternodegid)->Owner() == scatratimint_->Discretization()->Comm().MyPID())
        imasternodegidvec.push_back(imasternodegid);
    }
  }

  // remove potential duplicates from vectors
  std::sort(islavenodegidvec.begin(),islavenodegidvec.end());
  islavenodegidvec.erase(unique(islavenodegidvec.begin(),islavenodegidvec.end()),islavenodegidvec.end());
  std::sort(imasternodegidvec.begin(),imasternodegidvec.end());
  imasternodegidvec.erase(unique(imasternodegidvec.begin(),imasternodegidvec.end()),imasternodegidvec.end());

  // initialize coupling adapter
  if(scatratimint_->NumScal() < 1)
    dserror("Number of transported scalars not correctly set!");
  icoup_ = Teuchos::rcp(new ADAPTER::Coupling());
  icoup_->SetupCoupling(*(scatratimint_->Discretization()),*(scatratimint_->Discretization()),imasternodegidvec,islavenodegidvec,scatratimint_->NumScal());

  // generate interior and interface maps
  Teuchos::RCP<Epetra_Map> ifullmap = LINALG::MergeMap(icoup_->SlaveDofMap(),icoup_->MasterDofMap(),false);
  std::vector<Teuchos::RCP<const Epetra_Map> > maps;
  maps.push_back(LINALG::SplitMap(*(scatratimint_->Discretization()->DofRowMap()),*ifullmap));
  maps.push_back(icoup_->SlaveDofMap());
  maps.push_back(icoup_->MasterDofMap());

  // initialize global map extractor
  maps_ = Teuchos::rcp(new LINALG::MultiMapExtractor(*(scatratimint_->Discretization()->DofRowMap()),maps));
  maps_->CheckForValidMapExtractor();

  // initialize interface vector
  // Although the interface vector only contains the transformed master interface dofs, we still initialize it with
  // the full DofRowMap of the discretization to make it work for parallel computations.
  imasterphinp_ = LINALG::CreateVector(*(scatratimint_->Discretization()->DofRowMap()),false);

  // initialize auxiliary system matrices and associated transformation operators
  islavematrix_ = Teuchos::rcp(new LINALG::SparseMatrix(*(icoup_->SlaveDofMap()),81));
  imastermatrix_ = Teuchos::rcp(new LINALG::SparseMatrix(*(icoup_->SlaveDofMap()),81));
  islavetomastercoltransform_ = Teuchos::rcp(new FSI::UTILS::MatrixColTransform);
  islavetomasterrowtransform_ = Teuchos::rcp(new FSI::UTILS::MatrixRowTransform);
  islavetomasterrowcoltransform_ = Teuchos::rcp(new FSI::UTILS::MatrixRowColTransform);

  // initialize auxiliary residual vector
  islaveresidual_ = Teuchos::rcp(new Epetra_Vector(*(icoup_->SlaveDofMap())));

  return;
} // SCATRA::MeshtyingStrategyS2I::InitMeshtying


/*----------------------------------------------------------------------------*
 | initialize system matrix for scatra-scatra interface coupling   fang 10/14 |
 *----------------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseOperator> SCATRA::MeshtyingStrategyS2I::InitSystemMatrix() const
{
  // initialize system matrix and associated strategy for scatra-scatra interface coupling
  Teuchos::RCP<LINALG::SparseOperator> systemmatrix = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(*maps_,*maps_));

  return systemmatrix;
} // SCATRA::MeshtyingStrategyS2I::InitSystemMatrix


/*------------------------------------------------------------------------------------*
 | solve linear system of equations for scatra-scatra interface coupling   fang 12/14 |
 *------------------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2I::Solve(
    const Teuchos::RCP<LINALG::Solver>&            solver,         //! solver
    const Teuchos::RCP<LINALG::SparseOperator>&    systemmatrix,   //! system matrix
    const Teuchos::RCP<Epetra_Vector>&             increment,      //! increment vector
    const Teuchos::RCP<Epetra_Vector>&             residual,       //! residual vector
    const Teuchos::RCP<Epetra_Vector>&             phinp,          //! state vector at time n+1
    const int&                                     iteration,      //! number of current Newton-Raphson iteration
    const Teuchos::RCP<LINALG::KrylovProjector>&   projector       //! Krylov projector
    ) const
{
  solver->Solve(systemmatrix->EpetraOperator(),increment,residual,true,iteration==1,projector);

  return;
} // SCATRA::MeshtyingStrategyS2I::Solve

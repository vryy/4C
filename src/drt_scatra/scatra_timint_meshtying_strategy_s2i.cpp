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

#include "../drt_lib/drt_condition_utils.H"

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
imaps_(Teuchos::null),
icoup_(Teuchos::null),
auxmat_(Teuchos::null),
iphinp_(Teuchos::null),
imastertoslavetransform_(Teuchos::null),
islavetomastertransform_(Teuchos::null)
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
  condparams.set<int>("scatratype",scatratimint_->ScaTraType());
  condparams.set("isale",scatratimint_->IsALE());

  // set global and interface state vectors according to time-integration scheme
  scatratimint_->Discretization()->ClearState();
  scatratimint_->AddTimeIntegrationSpecificVectors();

  // set interface state vector iphinp_ with transformed dof values and add to element parameter list
  imaps_->InsertVector(icoup_->SlaveToMaster(maps_->ExtractVector(*(scatratimint_->Phiafnp()),1)),0,iphinp_);
  imaps_->InsertVector(icoup_->MasterToSlave(maps_->ExtractVector(*(scatratimint_->Phiafnp()),2)),1,iphinp_);
  condparams.set<Teuchos::RCP<const Epetra_Vector> >("iphinp",iphinp_);

  // evaluate scatra-scatra interface coupling at time t_{n+1} or t_{n+alpha_F}
  auxmat_->Zero();
  scatratimint_->Discretization()->EvaluateCondition(condparams,blocksystemmatrix,auxmat_,scatratimint_->Residual(),Teuchos::null,Teuchos::null,"S2ICoupling");
  scatratimint_->Discretization()->ClearState();

  // transform master and slave blocks of auxiliary system matrix and assemble into global system matrix
  auxmat_->Complete();
  (*imastertoslavetransform_)(auxmat_->FullRowMap(),auxmat_->FullColMap(),auxmat_->Matrix(0,0),1.,
      ADAPTER::CouplingMasterConverter(*icoup_),blocksystemmatrix->Matrix(2,1));
  (*islavetomastertransform_)(auxmat_->FullRowMap(),auxmat_->FullColMap(),auxmat_->Matrix(1,1),1.,
      ADAPTER::CouplingSlaveConverter(*icoup_),blocksystemmatrix->Matrix(1,2));

  return;
} // SCATRA::MeshtyingStrategyS2I::EvaluateMeshtying


/*----------------------------------------------------------------------*
 | perform setup of scatra-scatra interface coupling         fang 10/14 |
 *----------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2I::InitMeshtying()
{
  // extract scatra-scatra coupling conditions from discretization
  std::vector<DRT::Condition*> s2icouplingconditions;
  scatratimint_->Discretization()->GetCondition("S2ICoupling", s2icouplingconditions);

  // initialize int sets for global ids of interface nodes
  std::set<int> ithisnodegidset;
  std::set<int> iothernodegidset;

  // fill sets
  for (unsigned icond=0; icond<s2icouplingconditions.size(); ++icond)
  {
    const std::vector<int>* inodegids = s2icouplingconditions[icond]->Nodes();

    for (unsigned inode=0; inode<inodegids->size(); ++inode)
    {
      const int inodegid = (*inodegids)[inode];

      // insert global id of current node into associated set only if node is owned by current processor
      // need to make sure that node is stored on current processor, otherwise cannot resolve "->Owner()"
      if(scatratimint_->Discretization()->HaveGlobalNode(inodegid) and scatratimint_->Discretization()->gNode(inodegid)->Owner() == scatratimint_->Discretization()->Comm().MyPID())
      {
        // determine whether node is located on "This" or "Other" side of scatra-scatra interface
        if(*(s2icouplingconditions[icond]->Get<std::string>("Side")) == "This")
          ithisnodegidset.insert(inodegid);
        else if(*(s2icouplingconditions[icond]->Get<std::string>("Side")) == "Other")
          iothernodegidset.insert(inodegid);
        else
          dserror("Interface side must be either \"This\" or \"Other\"!");
      }
    }
  }

  // copy sets into vectors
  std::vector<int> ithisnodegidvec(ithisnodegidset.begin(),ithisnodegidset.end());
  std::vector<int> iothernodegidvec(iothernodegidset.begin(),iothernodegidset.end());

  // initialize coupling adapter
  if(scatratimint_->NumScal() < 1)
    dserror("Number of transported scalars not correctly set!");
  icoup_ = Teuchos::rcp(new ADAPTER::Coupling());
  icoup_->SetupCoupling(*(scatratimint_->Discretization()),*(scatratimint_->Discretization()),ithisnodegidvec,iothernodegidvec,scatratimint_->NumScal());

  // generate interior and interface maps
  Teuchos::RCP<Epetra_Map> ifullmap = LINALG::MergeMap(icoup_->MasterDofMap(),icoup_->SlaveDofMap(),false);
  std::vector<Teuchos::RCP<const Epetra_Map> > maps;
  maps.push_back(LINALG::SplitMap(*(scatratimint_->Discretization()->DofRowMap()),*ifullmap));
  maps.push_back(icoup_->SlaveDofMap());
  maps.push_back(icoup_->MasterDofMap());

  // initialize global and interface map extractors
  maps_ = Teuchos::rcp(new LINALG::MultiMapExtractor(*(scatratimint_->Discretization()->DofRowMap()),maps));
  maps_->CheckForValidMapExtractor();
  imaps_ = Teuchos::rcp(new LINALG::MapExtractor(*ifullmap,icoup_->SlaveDofMap(),icoup_->MasterDofMap()));
  imaps_->CheckForValidMapExtractor();

  // initialize interface vector
  iphinp_ = LINALG::CreateVector(*ifullmap,false);

  // initialize auxiliary system matrix and associated transformation operators
  auxmat_ = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(*imaps_,*imaps_));
  imastertoslavetransform_ = Teuchos::rcp(new FSI::UTILS::MatrixColTransform);
  islavetomastertransform_ = Teuchos::rcp(new FSI::UTILS::MatrixColTransform);

  return;
} // SCATRA::MeshtyingStrategyS2I::InitMeshtying


/*----------------------------------------------------------------------------*
 | initialize system matrix for scatra-scatra interface coupling   fang 10/14 |
 *----------------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseOperator> SCATRA::MeshtyingStrategyS2I::InitSystemMatrix() const
{
  // initialize system matrix and associated strategy for scatra-scatra interface coupling
  Teuchos::RCP<LINALG::SparseOperator> systemmatrix = Teuchos::rcp(new LINALG::BlockSparseMatrix<FLD::UTILS::InterfaceSplitStrategy>(*maps_,*maps_));
  Teuchos::rcp_static_cast<LINALG::BlockSparseMatrix<FLD::UTILS::InterfaceSplitStrategy> >(systemmatrix)->SetCondElements(DRT::UTILS::ConditionedElementMap(*(scatratimint_->Discretization()),"S2ICoupling"));

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

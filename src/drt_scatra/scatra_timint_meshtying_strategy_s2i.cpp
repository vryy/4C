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

#include "../drt_lib/drt_globalproblem.H"

#include "../drt_scatra/scatra_timint_implicit.H"

#include "../drt_scatra_ele/scatra_ele_action.H"

#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_utils.H"

#include "scatra_timint_meshtying_strategy_s2i.H"

/*----------------------------------------------------------------------*
 | constructor                                               fang 12/14 |
 *----------------------------------------------------------------------*/
SCATRA::MeshtyingStrategyS2I::MeshtyingStrategyS2I(
    SCATRA::ScaTraTimIntImpl*       scatratimint,   //! scalar transport time integrator
    const Teuchos::ParameterList&   parameters      //! input parameters for scatra-scatra interface coupling
    ) :
MeshtyingStrategyBase(scatratimint),
interfacemaps_(Teuchos::null),
conditionmaps_(Teuchos::null),
conditionmaps_slave_(Teuchos::null),
conditionmaps_master_(Teuchos::null),
icoup_(Teuchos::null),
islavematrix_(Teuchos::null),
imastermatrix_(Teuchos::null),
islavetomastercoltransform_(Teuchos::null),
islavetomasterrowtransform_(Teuchos::null),
islavetomasterrowcoltransform_(Teuchos::null),
islaveresidual_(Teuchos::null),
imasterphinp_(Teuchos::null),
invrowsums_(Teuchos::null),
invcolsums_(Teuchos::null),
parameters_(Teuchos::rcp(new Teuchos::ParameterList(parameters))),
rowequilibration_(
    DRT::INPUT::IntegralValue<int>(*parameters_,"EQUILIBRATION") == INPAR::S2I::equilibration_rows
    or
    DRT::INPUT::IntegralValue<int>(*parameters_,"EQUILIBRATION") == INPAR::S2I::equilibration_full
    ),
colequilibration_(
    DRT::INPUT::IntegralValue<int>(*parameters_,"EQUILIBRATION") == INPAR::S2I::equilibration_columns
    or
    DRT::INPUT::IntegralValue<int>(*parameters_,"EQUILIBRATION") == INPAR::S2I::equilibration_full
    ),
mortartype_(DRT::INPUT::IntegralValue<INPAR::S2I::MortarType>(*parameters_,"MORTARTYPE")),
matrixtype_(DRT::INPUT::IntegralValue<INPAR::S2I::MatrixType>(*parameters_,"MATRIXTYPE"))
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

  switch(mortartype_)
  {
  case INPAR::S2I::mortar_none:
  {
    // create parameter list
    Teuchos::ParameterList condparams;

    // action for elements
    condparams.set<int>("action",SCATRA::bd_calc_s2icoupling);

    // set global state vectors according to time-integration scheme
    scatratimint_->Discretization()->ClearState();
    scatratimint_->AddTimeIntegrationSpecificVectors();

    // fill interface state vector imasterphinp_ with transformed master dof values and add to discretization
    interfacemaps_->InsertVector(icoup_->MasterToSlave(interfacemaps_->ExtractVector(*(scatratimint_->Phiafnp()),2)),1,imasterphinp_);
    scatratimint_->Discretization()->SetState("imasterphinp",imasterphinp_);

    // evaluate scatra-scatra interface coupling at time t_{n+1} or t_{n+alpha_F}
    islavematrix_->Zero();
    imastermatrix_->Zero();
    islaveresidual_->PutScalar(0.);
    scatratimint_->Discretization()->EvaluateCondition(condparams,islavematrix_,imastermatrix_,islaveresidual_,Teuchos::null,Teuchos::null,"S2ICouplingSlave");
    scatratimint_->Discretization()->ClearState();

    // finalize interface matrices
    islavematrix_->Complete();
    imastermatrix_->Complete();

    // assemble global system matrix depending on matrix type
    switch(matrixtype_)
    {
      case INPAR::S2I::matrix_sparse:
      {
        // check matrix
        Teuchos::RCP<LINALG::SparseMatrix> systemmatrix = scatratimint_->SystemMatrix();
        if(systemmatrix == Teuchos::null)
          dserror("System matrix is not a sparse matrix!");

        // assemble linearizations of slave fluxes w.r.t. slave dofs into global system matrix
        systemmatrix->Add(*islavematrix_,false,1.,1.);

        // transform linearizations of slave fluxes w.r.t. master dofs and assemble into global system matrix
        (*islavetomastercoltransform_)(imastermatrix_->RowMap(),imastermatrix_->ColMap(),*imastermatrix_,1.,
            ADAPTER::CouplingSlaveConverter(*icoup_),*systemmatrix,true,true);

        // derive linearizations of master fluxes w.r.t. slave dofs and assemble into global system matrix
        (*islavetomasterrowtransform_)(*islavematrix_,-1.,ADAPTER::CouplingSlaveConverter(*icoup_),*systemmatrix,true);

        // derive linearizations of master fluxes w.r.t. master dofs and assemble into global system matrix
        (*islavetomasterrowcoltransform_)(*imastermatrix_,-1.,ADAPTER::CouplingSlaveConverter(*icoup_),ADAPTER::CouplingSlaveConverter(*icoup_),*systemmatrix,true,true);

        break;
      }

      case INPAR::S2I::matrix_block_geometry:
      {
        // check matrix
        Teuchos::RCP<LINALG::BlockSparseMatrixBase> blocksystemmatrix = scatratimint_->BlockSystemMatrix();
        if(blocksystemmatrix == Teuchos::null)
          dserror("System matrix is not a block matrix!");

        // assemble linearizations of slave fluxes w.r.t. slave dofs into global system matrix
        blocksystemmatrix->Matrix(1,1).Add(*islavematrix_,false,1.,1.);

        // transform linearizations of slave fluxes w.r.t. master dofs and assemble into global system matrix
        (*islavetomastercoltransform_)(imastermatrix_->RowMap(),imastermatrix_->ColMap(),*imastermatrix_,1.,
            ADAPTER::CouplingSlaveConverter(*icoup_),blocksystemmatrix->Matrix(1,2));

        // derive linearizations of master fluxes w.r.t. slave dofs and assemble into global system matrix
        (*islavetomasterrowtransform_)(*islavematrix_,-1.,ADAPTER::CouplingSlaveConverter(*icoup_),blocksystemmatrix->Matrix(2,1));

        // derive linearizations of master fluxes w.r.t. master dofs and assemble into global system matrix
        (*islavetomasterrowcoltransform_)(*imastermatrix_,-1.,ADAPTER::CouplingSlaveConverter(*icoup_),ADAPTER::CouplingSlaveConverter(*icoup_),blocksystemmatrix->Matrix(2,2),true,true);

        break;
      }

      case INPAR::S2I::matrix_block_condition:
      {
        // check matrix
        Teuchos::RCP<LINALG::BlockSparseMatrixBase> blocksystemmatrix = scatratimint_->BlockSystemMatrix();
        if(blocksystemmatrix == Teuchos::null)
          dserror("System matrix is not a block matrix!");

        Teuchos::RCP<LINALG::SparseMatrix> ksm(Teuchos::rcp(new LINALG::SparseMatrix(*icoup_->SlaveDofMap(),81,false)));
        Teuchos::RCP<LINALG::SparseMatrix> kms(Teuchos::rcp(new LINALG::SparseMatrix(*icoup_->MasterDofMap(),81,false)));
        Teuchos::RCP<LINALG::SparseMatrix> kmm(Teuchos::rcp(new LINALG::SparseMatrix(*icoup_->MasterDofMap(),81,false)));

        // transform linearizations of slave fluxes w.r.t. master dofs
        (*islavetomastercoltransform_)(imastermatrix_->RowMap(),imastermatrix_->ColMap(),*imastermatrix_,1.,ADAPTER::CouplingSlaveConverter(*icoup_),*ksm);
        ksm->Complete(*icoup_->MasterDofMap(),*icoup_->SlaveDofMap());

        // derive linearizations of master fluxes w.r.t. slave dofs
        (*islavetomasterrowtransform_)(*islavematrix_,-1.,ADAPTER::CouplingSlaveConverter(*icoup_),*kms);
        kms->Complete(*icoup_->SlaveDofMap(),*icoup_->MasterDofMap());

        // derive linearizations of master fluxes w.r.t. master dofs
        (*islavetomasterrowcoltransform_)(*imastermatrix_,-1.,ADAPTER::CouplingSlaveConverter(*icoup_),ADAPTER::CouplingSlaveConverter(*icoup_),*kmm);
        kmm->Complete();

        Teuchos::RCP<LINALG::BlockSparseMatrixBase> blockkss(islavematrix_->Split<LINALG::DefaultBlockMatrixStrategy>(*conditionmaps_slave_,*conditionmaps_slave_));
        blockkss->Complete();
        Teuchos::RCP<LINALG::BlockSparseMatrixBase> blockksm(ksm->Split<LINALG::DefaultBlockMatrixStrategy>(*conditionmaps_master_,*conditionmaps_slave_));
        blockksm->Complete();
        Teuchos::RCP<LINALG::BlockSparseMatrixBase> blockkms(kms->Split<LINALG::DefaultBlockMatrixStrategy>(*conditionmaps_slave_,*conditionmaps_master_));
        blockkms->Complete();
        Teuchos::RCP<LINALG::BlockSparseMatrixBase> blockkmm(kmm->Split<LINALG::DefaultBlockMatrixStrategy>(*conditionmaps_master_,*conditionmaps_master_));
        blockkmm->Complete();

        // assemble interface block matrices into global block system matrix
        blocksystemmatrix->Add(*blockkss,false,1.,1.);
        blocksystemmatrix->Add(*blockksm,false,1.,1.);
        blocksystemmatrix->Add(*blockkms,false,1.,1.);
        blocksystemmatrix->Add(*blockkmm,false,1.,1.);

        break;
      }

      default:
      {
        dserror("Type of global system matrix for scatra-scatra interface coupling not recognized!");
        break;
      }
    }

    // assemble slave residuals into global residual vector
    interfacemaps_->AddVector(islaveresidual_,1,scatratimint_->Residual());

    // transform master residuals and assemble into global residual vector
    interfacemaps_->AddVector(icoup_->SlaveToMaster(islaveresidual_),2,scatratimint_->Residual(),-1.);

    break;
  }

  case INPAR::S2I::mortar_standard:
  {
    dserror("Not yet implemented!");
    break;
  }

  case INPAR::S2I::mortar_saddlepoint:
  {
    dserror("Not yet implemented!");
    break;
  }

  case INPAR::S2I::mortar_condensed:
  {
    dserror("Not yet implemented!");
    break;
  }

  default:
  {
    dserror("Not yet implemented!");
    break;
  }
  }

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

  // determine type of mortar meshtying
  switch(mortartype_)
  {
  // setup scatra-scatra interface coupling for interfaces with pairwise overlapping interface nodes
  case INPAR::S2I::mortar_none:
  {
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

    // initialize non-mortar coupling adapter
    if(scatratimint_->NumScal() < 1)
      dserror("Number of transported scalars not correctly set!");
    icoup_ = Teuchos::rcp(new ADAPTER::Coupling());
    icoup_->SetupCoupling(*(scatratimint_->Discretization()),*(scatratimint_->Discretization()),imasternodegidvec,islavenodegidvec,scatratimint_->NumDofPerNode(),true,1.e-8);

    // generate interior and interface maps
    Teuchos::RCP<Epetra_Map> ifullmap = LINALG::MergeMap(icoup_->SlaveDofMap(),icoup_->MasterDofMap(),false);
    std::vector<Teuchos::RCP<const Epetra_Map> > imaps;
    imaps.push_back(LINALG::SplitMap(*(scatratimint_->Discretization()->DofRowMap()),*ifullmap));
    imaps.push_back(icoup_->SlaveDofMap());
    imaps.push_back(icoup_->MasterDofMap());

    // initialize global map extractor
    interfacemaps_ = Teuchos::rcp(new LINALG::MultiMapExtractor(*(scatratimint_->Discretization()->DofRowMap()),imaps));
    interfacemaps_->CheckForValidMapExtractor();

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

    if(matrixtype_ == INPAR::S2I::matrix_block_condition)
    {
      std::vector<Teuchos::RCP<DRT::Condition> > partitioningconditions;
      scatratimint_->Discretization()->GetCondition("S2ICouplingPartitioning",partitioningconditions);
      const unsigned ncond = partitioningconditions.size();
      if(!ncond)
        dserror("For block preconditioning based on domain partitioning, at least one associated condition needs to be specified in the input file!");
      std::vector<std::set<int> > dofids(ncond);
      std::vector<Teuchos::RCP<const Epetra_Map> > conditionmaps(ncond);
      for(unsigned icond=0; icond<ncond; ++icond)
      {
        const std::vector<int>* nodegids = partitioningconditions[icond]->Nodes();

        for (unsigned inode=0; inode<nodegids->size(); ++inode)
        {
          const int nodegid = (*nodegids)[inode];

          // consider current node only if node is owned by current processor
          // need to make sure that node is stored on current processor, otherwise cannot resolve "->Owner()"
          if(scatratimint_->Discretization()->HaveGlobalNode(nodegid) and scatratimint_->Discretization()->gNode(nodegid)->Owner() == scatratimint_->Discretization()->Comm().MyPID())
          {
            const std::vector<int> nodedofs = scatratimint_->Discretization()->Dof(scatratimint_->Discretization()->gNode(nodegid));
            std::copy(nodedofs.begin(),nodedofs.end(),std::inserter(dofids[icond],dofids[icond].end()));
          }
        }

        int nummyelements(0);
        int* myglobalelements(NULL);
        std::vector<int> dofidvec;
        if(dofids[icond].size() > 0)
        {
          dofidvec.reserve(dofids[icond].size());
          dofidvec.assign(dofids[icond].begin(),dofids[icond].end());
          nummyelements = dofidvec.size();
          myglobalelements = &(dofidvec[0]);
        }
        conditionmaps[icond] = Teuchos::rcp(new Epetra_Map(-1,nummyelements,myglobalelements,scatratimint_->DofRowMap()->IndexBase(),scatratimint_->DofRowMap()->Comm()));
      }

      // initialize condition map extractors
      conditionmaps_ = Teuchos::rcp(new LINALG::MultiMapExtractor(*(scatratimint_->Discretization()->DofRowMap()),conditionmaps));
      conditionmaps_->CheckForValidMapExtractor();

      std::vector<Teuchos::RCP<const Epetra_Map> > conditionmaps_slave(ncond);
      std::vector<Teuchos::RCP<const Epetra_Map> > conditionmaps_master(ncond);
      for(unsigned icond=0; icond<ncond; ++icond)
      {
        std::vector<Teuchos::RCP<const Epetra_Map> > maps(2);
        maps[0] = conditionmaps[icond];
        maps[1] = icoup_->SlaveDofMap();
        conditionmaps_slave[icond] = LINALG::MultiMapExtractor::IntersectMaps(maps);
        maps[1] = icoup_->MasterDofMap();
        conditionmaps_master[icond] = LINALG::MultiMapExtractor::IntersectMaps(maps);
      }
      conditionmaps_slave_ = Teuchos::rcp(new LINALG::MultiMapExtractor(*icoup_->SlaveDofMap(),conditionmaps_slave));
      conditionmaps_slave_->CheckForValidMapExtractor();
      conditionmaps_master_ = Teuchos::rcp(new LINALG::MultiMapExtractor(*icoup_->MasterDofMap(),conditionmaps_master));
      conditionmaps_master_->CheckForValidMapExtractor();

      // extract number of linear solver for each block of global block system matrix
      const int blocksolver = parameters_->get<int>("BLOCKSOLVER");

      // safety check
      if(blocksolver == -1)
        dserror("Invalid number of linear solver for each block of global block system matrix!");

      // feed block preconditioner with solver and nullspace information for each block of global block system matrix
      const Teuchos::ParameterList& blocksolverparams = DRT::Problem::Instance()->SolverParams(blocksolver);
      for(unsigned icond=0; icond<ncond; ++icond)
      {
        std::ostringstream icondstr;
        icondstr << icond+1;
        scatratimint_->Solver()->PutSolverParamsToSubParams("Inverse"+icondstr.str(),blocksolverparams);
        scatratimint_->Discretization()->ComputeNullSpaceIfNecessary(scatratimint_->Solver()->Params().sublist("Inverse"+icondstr.str()));
      }
    }

    break;
  }

  // setup scatra-scatra interface coupling for interfaces with non-overlapping interface nodes
  case INPAR::S2I::mortar_standard:
  case INPAR::S2I::mortar_saddlepoint:
  case INPAR::S2I::mortar_condensed:
  {
    dserror("Not yet implemented!");
    break;
  }

  default:
  {
    dserror("Type of mortar meshtying for scatra-scatra interface coupling not recognized!");
    break;
  }
  }

  // initialize vectors for row and column sums of global system matrix if necessary
  if(rowequilibration_)
    invrowsums_ = Teuchos::rcp(new Epetra_Vector(*scatratimint_->Discretization()->DofRowMap(),false));
  if(colequilibration_)
    invcolsums_ = Teuchos::rcp(new Epetra_Vector(*scatratimint_->Discretization()->DofRowMap(),false));

  return;
} // SCATRA::MeshtyingStrategyS2I::InitMeshtying


/*----------------------------------------------------------------------------*
 | initialize system matrix for scatra-scatra interface coupling   fang 10/14 |
 *----------------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseOperator> SCATRA::MeshtyingStrategyS2I::InitSystemMatrix() const
{
  Teuchos::RCP<LINALG::SparseOperator> systemmatrix(Teuchos::null);

  switch(matrixtype_)
  {
  case INPAR::S2I::matrix_sparse:
  {
    // initialize system matrix
    systemmatrix = Teuchos::rcp(new LINALG::SparseMatrix(*(scatratimint_->Discretization()->DofRowMap()),27,false,true));
    break;
  }

  case INPAR::S2I::matrix_block_geometry:
  {
    // initialize system matrix and associated strategy
    systemmatrix = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(*interfacemaps_,*interfacemaps_));
    break;
  }

  case INPAR::S2I::matrix_block_condition:
  {
    // initialize system matrix and associated strategy
    systemmatrix = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(*conditionmaps_,*conditionmaps_));

    break;
  }

  default:
  {
    dserror("Type of global system matrix for scatra-scatra interface coupling not recognized!");
    break;
  }
  }

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
  if(mortartype_ != INPAR::S2I::mortar_saddlepoint)
  {
    // equilibrate global system of equations if necessary
    EquilibrateSystem(systemmatrix,residual);

    // solve global system of equations
    solver->Solve(systemmatrix->EpetraOperator(),increment,residual,true,iteration==1,projector);

    // unequilibrate global increment vector if necessary
    UnequilibrateIncrement(increment);
  }

  else
    dserror("Not yet implemented!");

  return;
} // SCATRA::MeshtyingStrategyS2I::Solve


/*----------------------------------------------------------------------*
 | equilibrate global system of equations if necessary       fang 05/15 |
 *----------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2I::EquilibrateSystem(
    const Teuchos::RCP<LINALG::SparseOperator>&   systemmatrix,   //! system matrix
    const Teuchos::RCP<Epetra_Vector>&            residual        //! residual vector
    ) const
{
  if(rowequilibration_ or colequilibration_)
  {
    // perform equilibration depending on type of global system matrix
    switch(matrixtype_)
    {
    case INPAR::S2I::matrix_sparse:
    {
      // check matrix
      Teuchos::RCP<LINALG::SparseMatrix> sparsematrix = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(systemmatrix);
      if(sparsematrix == Teuchos::null)
        dserror("System matrix is not a sparse matrix!");

      // perform row equilibration
      if(rowequilibration_)
      {
        // compute inverse row sums of global system matrix
        ComputeInvRowSums(*sparsematrix,invrowsums_);

        // perform row equilibration of global system matrix
        EquilibrateMatrixRows(*sparsematrix,invrowsums_);
      }

      // perform column equilibration
      if(colequilibration_)
      {
        // compute inverse column sums of global system matrix
        ComputeInvColSums(*sparsematrix,invcolsums_);

        // perform column equilibration of global system matrix
        EquilibrateMatrixColumns(*sparsematrix,invcolsums_);
      }

      break;
    }

    case INPAR::S2I::matrix_block_condition:
    {
      // check matrix
      Teuchos::RCP<LINALG::BlockSparseMatrixBase> blocksparsematrix = Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(systemmatrix);
      if(blocksparsematrix == Teuchos::null)
        dserror("System matrix is not a block sparse matrix!");

      // perform row equilibration
      if(rowequilibration_)
      {
        for(int i=0; i<blocksparsematrix->Rows(); ++i)
        {
          // compute inverse row sums of current main diagonal matrix block
          Teuchos::RCP<Epetra_Vector> invrowsums(Teuchos::rcp(new Epetra_Vector(blocksparsematrix->Matrix(i,i).RangeMap())));
          ComputeInvRowSums(blocksparsematrix->Matrix(i,i),invrowsums);

          // perform row equilibration of matrix blocks in current row block of global system matrix
          for(int j=0; j<blocksparsematrix->Cols(); ++j)
            EquilibrateMatrixRows(blocksparsematrix->Matrix(i,j),invrowsums);

          // insert inverse row sums of current main diagonal matrix block into global vector
          conditionmaps_->InsertVector(invrowsums,i,invrowsums_);
        }
      }

      // perform column equilibration
      if(colequilibration_)
      {
        for(int j=0; j<blocksparsematrix->Cols(); ++j)
        {
          // compute inverse column sums of current main diagonal matrix block
          Teuchos::RCP<Epetra_Vector> invcolsums(Teuchos::rcp(new Epetra_Vector(blocksparsematrix->Matrix(j,j).DomainMap())));
          ComputeInvColSums(blocksparsematrix->Matrix(j,j),invcolsums);

          // perform column equilibration of matrix blocks in current column block of global system matrix
          for(int i=0; i<blocksparsematrix->Rows(); ++i)
            EquilibrateMatrixColumns(blocksparsematrix->Matrix(i,j),invcolsums);

          // insert inverse column sums of current main diagonal matrix block into global vector
          conditionmaps_->InsertVector(invcolsums,j,invcolsums_);
        }
      }

      break;
    }

    default:
    {
      dserror("Equilibration of global system of equations for scatra-scatra interface coupling is not implemented for chosen type of global system matrix!");
      break;
    }
    }

    // perform equilibration of global residual vector
    if(rowequilibration_)
      if(residual->Multiply(1.,*invrowsums_,*residual,0.))
        dserror("Equilibration of global residual vector failed!");
  }

  return;
} // SCATRA::MeshtyingStrategyS2I::EquilibrateSystem


/*----------------------------------------------------------------------------*
 | compute inverse sums of absolute values of matrix row entries   fang 06/15 |
 *----------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2I::ComputeInvRowSums(
    const LINALG::SparseMatrix&          matrix,      //! matrix
    const Teuchos::RCP<Epetra_Vector>&   invrowsums   //! inverse sums of absolute values of row entries in matrix
    ) const
{
  // compute inverse row sums of matrix
  if(matrix.EpetraMatrix()->InvRowSums(*invrowsums))
    dserror("Inverse row sums of matrix could not be successfully computed!");

  // take square root of inverse row sums if matrix is scaled from left and right
  if(colequilibration_)
    for(int i=0; i<invrowsums->MyLength(); ++i)
      (*invrowsums)[i] = sqrt((*invrowsums)[i]);

  return;
} // SCATRA::MeshtyingStrategyS2I::ComputeInvRowSums


/*-------------------------------------------------------------------------------*
 | compute inverse sums of absolute values of matrix column entries   fang 06/15 |
 *-------------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2I::ComputeInvColSums(
    const LINALG::SparseMatrix&          matrix,      //! matrix
    const Teuchos::RCP<Epetra_Vector>&   invcolsums   //! inverse sums of absolute values of column entries in matrix
    ) const
{
  // compute inverse column sums of matrix
  if(matrix.EpetraMatrix()->InvColSums(*invcolsums))
    dserror("Inverse column sums of matrix could not be successfully computed!");

  // take square root of inverse column sums if matrix is scaled from left and right
  if(rowequilibration_)
    for(int i=0; i<invcolsums->MyLength(); ++i)
      (*invcolsums)[i] = sqrt((*invcolsums)[i]);

  return;
} // SCATRA::MeshtyingStrategyS2I::ComputeInvColSums


/*----------------------------------------------------------------------*
 | equilibrate matrix rows                                   fang 06/15 |
 *----------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2I::EquilibrateMatrixRows(
    LINALG::SparseMatrix&                matrix,      //! matrix
    const Teuchos::RCP<Epetra_Vector>&   invrowsums   //! sums of absolute values of row entries in matrix
    ) const
{
  if(matrix.LeftScale(*invrowsums))
    dserror("Row equilibration of matrix failed!");

  return;
} // SCATRA::MeshtyingStrategyS2I::EquilibrateMatrixRows


/*----------------------------------------------------------------------*
 | equilibrate matrix columns                                fang 06/15 |
 *----------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2I::EquilibrateMatrixColumns(
    LINALG::SparseMatrix&                matrix,      //! matrix
    const Teuchos::RCP<Epetra_Vector>&   invcolsums   //! sums of absolute values of column entries in matrix
    ) const
{
  if(matrix.RightScale(*invcolsums))
    dserror("Column equilibration of matrix failed!");

  return;
} // SCATRA::MeshtyingStrategyS2I::EquilibrateMatrixColumns


/*----------------------------------------------------------------------*
 | unequilibrate global increment vector if necessary        fang 05/15 |
 *----------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2I::UnequilibrateIncrement(
    const Teuchos::RCP<Epetra_Vector>&   increment   //! increment vector
    ) const
{
  // unequilibrate global increment vector if necessary
  if(colequilibration_)
  {
    if(increment->Multiply(1.,*invcolsums_,*increment,0.))
      dserror("Unequilibration of global increment vector failed!");
  }

  return;
} // SCATRA::MeshtyingStrategyS2I::UnequilibrateIncrement

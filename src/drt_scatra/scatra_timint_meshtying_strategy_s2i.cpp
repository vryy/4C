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

#include "scatra_timint_meshtying_strategy_s2i.H"
#include "scatra_timint_implicit.H"

#include "../drt_adapter/adapter_coupling.H"
#include "../drt_adapter/adapter_coupling_mortar.H"

#include "../drt_fluid/fluid_utils.H"

#include "../drt_fsi/fsi_matrixtransform.H"

#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../drt_mortar/mortar_coupling3d_classes.H"
#include "../drt_mortar/mortar_interface.H"
#include "../drt_mortar/mortar_projector.H"
#include "../drt_mortar/mortar_shape_utils.H"

#include "../drt_scatra_ele/scatra_ele_action.H"
#include "../drt_scatra_ele/scatra_ele_calc_utils.H"

#include "../linalg/linalg_solver.H"

/*----------------------------------------------------------------------*
 | constructor                                               fang 12/14 |
 *----------------------------------------------------------------------*/
SCATRA::MeshtyingStrategyS2I::MeshtyingStrategyS2I(
    SCATRA::ScaTraTimIntImpl*       scatratimint,   //! scalar transport time integrator
    const Teuchos::ParameterList&   parameters      //! input parameters for scatra-scatra interface coupling
    ) :
MeshtyingStrategyBase(scatratimint),
interfacemaps_(Teuchos::null),
blockmaps_(Teuchos::null),
blockmaps_slave_(Teuchos::null),
blockmaps_master_(Teuchos::null),
icoup_(Teuchos::null),
icoupmortar_(),
imortarcells_(),
islavematrix_(Teuchos::null),
imastermatrix_(Teuchos::null),
imastertoslaverowtransform_(Teuchos::null),
islavetomastercoltransform_(Teuchos::null),
islavetomasterrowtransform_(Teuchos::null),
islavetomasterrowcoltransform_(Teuchos::null),
islaveresidual_(Teuchos::null),
imasterresidual_(Teuchos::null),
imasterphinp_(Teuchos::null),
invrowsums_(Teuchos::null),
invcolsums_(Teuchos::null),
parameters_(Teuchos::rcp(new Teuchos::ParameterList(parameters))),
matrixtype_(DRT::INPUT::IntegralValue<INPAR::S2I::MatrixType>(*parameters_,"MATRIXTYPE")),
slaveconditions_(),
rowequilibration_(
    DRT::INPUT::IntegralValue<INPAR::S2I::EquilibrationMethods>(*parameters_,"EQUILIBRATION") == INPAR::S2I::equilibration_rows
    or
    DRT::INPUT::IntegralValue<INPAR::S2I::EquilibrationMethods>(*parameters_,"EQUILIBRATION") == INPAR::S2I::equilibration_full
    ),
colequilibration_(
    DRT::INPUT::IntegralValue<INPAR::S2I::EquilibrationMethods>(*parameters_,"EQUILIBRATION") == INPAR::S2I::equilibration_columns
    or
    DRT::INPUT::IntegralValue<INPAR::S2I::EquilibrationMethods>(*parameters_,"EQUILIBRATION") == INPAR::S2I::equilibration_full
    ),
mortartype_(DRT::INPUT::IntegralValue<INPAR::S2I::MortarType>(*parameters_,"MORTARTYPE")),
slaveonly_(DRT::INPUT::IntegralValue<bool>(*parameters_,"SLAVEONLY"))
{
  return;
} // SCATRA::MeshtyingStrategyS2I::MeshtyingStrategyS2I


/*-----------------------------------------------------------------------*
 | evaluate scatra-scatra interface coupling conditions       fang 10/14 |
 *-----------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2I::EvaluateMeshtying()
{
  // time measurement: evaluate condition 'S2ICoupling'
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:       + evaluate condition 'S2ICoupling'");

  // extract scatra-scatra coupling conditions from discretization
  std::vector<DRT::Condition*> conditions;
  scatratimint_->Discretization()->GetCondition("S2ICoupling",conditions);

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
    if(not slaveonly_)
      imastermatrix_->Zero();
    islaveresidual_->PutScalar(0.);
    for(unsigned icondition=0; icondition<conditions.size(); ++icondition)
    {
      if(conditions[icondition]->GetInt("interface side") == INPAR::S2I::side_slave)
      {
        if(not slaveonly_)
          scatratimint_->Discretization()->EvaluateCondition(condparams,islavematrix_,imastermatrix_,islaveresidual_,Teuchos::null,Teuchos::null,"S2ICoupling",conditions[icondition]->GetInt("ConditionID"));
        else
          scatratimint_->Discretization()->EvaluateCondition(condparams,islavematrix_,Teuchos::null,islaveresidual_,Teuchos::null,Teuchos::null,"S2ICoupling",conditions[icondition]->GetInt("ConditionID"));
      }
    }
    scatratimint_->Discretization()->ClearState();

    // finalize interface matrices
    islavematrix_->Complete();
    if(not slaveonly_)
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

        if(not slaveonly_)
        {
          // transform linearizations of slave fluxes w.r.t. master dofs and assemble into global system matrix
          (*islavetomastercoltransform_)(imastermatrix_->RowMap(),imastermatrix_->ColMap(),*imastermatrix_,1.,
              ADAPTER::CouplingSlaveConverter(*icoup_),*systemmatrix,true,true);

          // derive linearizations of master fluxes w.r.t. slave dofs and assemble into global system matrix
          (*islavetomasterrowtransform_)(*islavematrix_,-1.,ADAPTER::CouplingSlaveConverter(*icoup_),*systemmatrix,true);

          // derive linearizations of master fluxes w.r.t. master dofs and assemble into global system matrix
          (*islavetomasterrowcoltransform_)(*imastermatrix_,-1.,ADAPTER::CouplingSlaveConverter(*icoup_),ADAPTER::CouplingSlaveConverter(*icoup_),*systemmatrix,true,true);
        }

        // In case the interface linearizations and residuals are evaluated on slave side only,
        // we now apply a standard meshtying algorithm to condense out the master-side degrees of freedom.
        else
        {
          // initialize temporary matrix for master-side rows of system matrix
          LINALG::SparseMatrix systemmatrixrowsmaster(*icoup_->MasterDofMap(),81);

          // loop over all master-side rows of system matrix
          for(int masterdoflid=0; masterdoflid<icoup_->MasterDofMap()->NumMyElements(); ++masterdoflid)
          {
            // determine global ID of current matrix row
            const int masterdofgid = icoup_->MasterDofMap()->GID(masterdoflid);
            if(masterdofgid < 0)
              dserror("Couldn't find local ID %d in map!",masterdoflid);

            // extract current matrix row from system matrix
            const int length = systemmatrix->EpetraMatrix()->NumGlobalEntries(masterdofgid);
            int numentries(0);
            std::vector<double> values(length,0.);
            std::vector<int> indices(length,0);
            if(systemmatrix->EpetraMatrix()->ExtractGlobalRowCopy(masterdofgid,length,numentries,&values[0],&indices[0]) != 0)
              dserror("Cannot extract matrix row with global ID %d from system matrix!",masterdofgid);

            // copy current matrix row of system matrix into temporary matrix
            if(systemmatrixrowsmaster.EpetraMatrix()->InsertGlobalValues(masterdofgid,numentries,&values[0],&indices[0]) < 0)
              dserror("Cannot insert matrix row with global ID %d into temporary matrix!",masterdofgid);
          }

          // zero out master-side rows of system matrix and put a one on the main diagonal
          systemmatrix->Complete();
          systemmatrix->ApplyDirichlet(*icoup_->MasterDofMap(),true);
          systemmatrix->UnComplete();

          // loop over all master-side rows of system matrix
          for(int masterdoflid=0; masterdoflid<icoup_->MasterDofMap()->NumMyElements(); ++masterdoflid)
          {
            // determine global ID of current matrix row
            const int masterdofgid = icoup_->MasterDofMap()->GID(masterdoflid);
            if(masterdofgid < 0)
              dserror("Couldn't find local ID %d in map!",masterdoflid);

            // determine global ID of associated slave-side matrix column
            const int slavedofgid = icoup_->PermSlaveDofMap()->GID(masterdoflid);
            if(slavedofgid < 0)
              dserror("Couldn't find local ID %d in permuted map!",masterdoflid);

            // insert value -1. into intersection of master-side row and slave-side column in system matrix
            // this effectively forces the master-side degree of freedom to assume the same value as the slave-side degree of freedom
            const double value(-1.);
            if(systemmatrix->EpetraMatrix()->InsertGlobalValues(masterdofgid,1,&value,&slavedofgid) < 0)
              dserror("Cannot insert value -1. into matrix row with global ID %d and matrix column with global ID %d!",masterdofgid,slavedofgid);

            // insert zero into intersection of master-side row and slave-side column in temporary matrix
            // this prevents the system matrix from changing its graph when calling this function again during the next Newton iteration
            const double zero(0.);
            if(systemmatrixrowsmaster.EpetraMatrix()->InsertGlobalValues(masterdofgid,1,&zero,&slavedofgid) < 0)
              dserror("Cannot insert zero into matrix row with global ID %d and matrix column with global ID %d!",masterdofgid,slavedofgid);
          }

          // finalize temporary matrix with master-side rows of system matrix
          systemmatrixrowsmaster.Complete(*scatratimint_->DofRowMap(),*icoup_->MasterDofMap());

          // add master-side rows of system matrix to corresponding slave-side rows to finalize matrix condensation of master-side degrees of freedom
          (*imastertoslaverowtransform_)(
              systemmatrixrowsmaster,
              1.,
              ADAPTER::CouplingMasterConverter(*icoup_),
              *systemmatrix,
              true
              );
        }

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

        if(not slaveonly_)
        {
          // transform linearizations of slave fluxes w.r.t. master dofs and assemble into global system matrix
          (*islavetomastercoltransform_)(imastermatrix_->RowMap(),imastermatrix_->ColMap(),*imastermatrix_,1.,
              ADAPTER::CouplingSlaveConverter(*icoup_),blocksystemmatrix->Matrix(1,2));

          // derive linearizations of master fluxes w.r.t. slave dofs and assemble into global system matrix
          (*islavetomasterrowtransform_)(*islavematrix_,-1.,ADAPTER::CouplingSlaveConverter(*icoup_),blocksystemmatrix->Matrix(2,1));

          // derive linearizations of master fluxes w.r.t. master dofs and assemble into global system matrix
          (*islavetomasterrowcoltransform_)(*imastermatrix_,-1.,ADAPTER::CouplingSlaveConverter(*icoup_),ADAPTER::CouplingSlaveConverter(*icoup_),blocksystemmatrix->Matrix(2,2),true,true);
        }

        // safety check
        else
          dserror("Scatra-scatra interface coupling with evaluation of interface linearizations and residuals on slave side only is not yet available for block system matrices!");

        break;
      }

      case INPAR::S2I::matrix_block_condition:
      case INPAR::S2I::matrix_block_condition_dof:
      {
        // check matrix
        Teuchos::RCP<LINALG::BlockSparseMatrixBase> blocksystemmatrix = scatratimint_->BlockSystemMatrix();
        if(blocksystemmatrix == Teuchos::null)
          dserror("System matrix is not a block matrix!");

        Teuchos::RCP<LINALG::BlockSparseMatrixBase> blockkss(islavematrix_->Split<LINALG::DefaultBlockMatrixStrategy>(*blockmaps_slave_,*blockmaps_slave_));
        blockkss->Complete();

        // assemble interface block matrix into global block system matrix
        blocksystemmatrix->Add(*blockkss,false,1.,1.);

        if(not slaveonly_)
        {
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

          Teuchos::RCP<LINALG::BlockSparseMatrixBase> blockksm(ksm->Split<LINALG::DefaultBlockMatrixStrategy>(*blockmaps_master_,*blockmaps_slave_));
          blockksm->Complete();
          Teuchos::RCP<LINALG::BlockSparseMatrixBase> blockkms(kms->Split<LINALG::DefaultBlockMatrixStrategy>(*blockmaps_slave_,*blockmaps_master_));
          blockkms->Complete();
          Teuchos::RCP<LINALG::BlockSparseMatrixBase> blockkmm(kmm->Split<LINALG::DefaultBlockMatrixStrategy>(*blockmaps_master_,*blockmaps_master_));
          blockkmm->Complete();

          // assemble interface block matrices into global block system matrix
          blocksystemmatrix->Add(*blockksm,false,1.,1.);
          blocksystemmatrix->Add(*blockkms,false,1.,1.);
          blocksystemmatrix->Add(*blockkmm,false,1.,1.);
        }

        // safety check
        else
          dserror("Scatra-scatra interface coupling with evaluation of interface linearizations and residuals on slave side only is not yet available for block system matrices!");

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

    if(not slaveonly_)
      // transform master residuals and assemble into global residual vector
      interfacemaps_->AddVector(icoup_->SlaveToMaster(islaveresidual_),2,scatratimint_->Residual(),-1.);

    // In case the interface linearizations and residuals are evaluated on slave side only,
    // we now apply a standard meshtying algorithm to condense out the master-side degrees of freedom.
    else
    {
      // initialize temporary vector for master-side entries of residual vector
      Teuchos::RCP<Epetra_Vector> residualmaster = Teuchos::rcp(new Epetra_Vector(*icoup_->MasterDofMap()));

      // loop over all master-side entries of residual vector
      for(int masterdoflid=0; masterdoflid<icoup_->MasterDofMap()->NumMyElements(); ++masterdoflid)
      {
        // determine global ID of current vector entry
        const int masterdofgid = icoup_->MasterDofMap()->GID(masterdoflid);
        if(masterdofgid < 0)
          dserror("Couldn't find local ID %d in map!",masterdoflid);

        // copy current vector entry into temporary vector
        if(residualmaster->ReplaceGlobalValue(masterdofgid,0,(*scatratimint_->Residual())[scatratimint_->DofRowMap()->LID(masterdofgid)]))
          dserror("Cannot insert residual vector entry with global ID %d into temporary vector!",masterdofgid);

        // zero out current vector entry
        if(scatratimint_->Residual()->ReplaceGlobalValue(masterdofgid,0,0.))
          dserror("Cannot insert zero into residual vector entry with global ID %d!",masterdofgid);
      }

      // add master-side entries of residual vector to corresponding slave-side entries to finalize vector condensation of master-side degrees of freedom
      interfacemaps_->AddVector(icoup_->MasterToSlave(residualmaster),1,scatratimint_->Residual());
    }

    break;
  }

  case INPAR::S2I::mortar_standard:
  {
    // initialize auxiliary system matrices and vectors
    islavematrix_->Zero();
    imastermatrix_->Zero();
    islaveresidual_->PutScalar(0.);
    imasterresidual_->PutScalar(0.);

    // loop over all scatra-scatra coupling interfaces
    for(std::map<const int,DRT::Condition* const>::iterator islavecondition=slaveconditions_.begin(); islavecondition!=slaveconditions_.end(); ++islavecondition)
    {
      // extract mortar interface discretization
      const DRT::Discretization& idiscret = icoupmortar_[islavecondition->first]->Interface()->Discret();

      // export global state vector to mortar interface
      Teuchos::RCP<Epetra_Vector> iphinp = Teuchos::rcp(new Epetra_Vector(*idiscret.DofColMap(),false));
      LINALG::Export(*scatratimint_->Phiafnp(),*iphinp);

      // extract current condition
      DRT::Condition& condition = *islavecondition->second;

      // loop over all mortar integration cells at current interface
      for(unsigned icell=0; icell<imortarcells_[islavecondition->first].size(); ++icell)
      {
        // safety check
        if(imortarcells_[islavecondition->first][icell] == Teuchos::null)
          dserror("Invalid mortar integration cell!");

        // evaluate current cell
        EvaluateMortarCell(*imortarcells_[islavecondition->first][icell],condition,*iphinp,idiscret);
      }
    }

    // finalize auxiliary system matrices
    islavematrix_->Complete(*interfacemaps_->FullMap(),*interfacemaps_->Map(1));
    imastermatrix_->Complete(*interfacemaps_->FullMap(),*interfacemaps_->Map(2));

    // finalize auxiliary residual vector for master residuals
    if(imasterresidual_->GlobalAssemble(Add,true))
      dserror("Assembly of auxiliary residual vector for master residuals not successful!");

    // assemble interface contributions into global system of equations
    const Teuchos::RCP<LINALG::SparseMatrix> systemmatrix = scatratimint_->SystemMatrix();
    if(systemmatrix == Teuchos::null)
      dserror("System matrix is not a sparse matrix!");
    systemmatrix->Add(*islavematrix_,false,1.,1.);
    systemmatrix->Add(*imastermatrix_,false,1.,1.);
    interfacemaps_->AddVector(islaveresidual_,1,scatratimint_->Residual());
    interfacemaps_->AddVector(*imasterresidual_,2,*scatratimint_->Residual());

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


/*--------------------------------------------------------------------------------------*
 | evaluate single mortar integration cell                                   fang 01/16 |
 *--------------------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2I::EvaluateMortarCell(
    MORTAR::IntCell&             cell,        //!< mortar integration cell
    DRT::Condition&              condition,   //!< scatra-scatra interface coupling condition
    const Epetra_Vector&         iphinp,      //!< interface state vector
    const DRT::Discretization&   idiscret     //!< interface discretization
    ) const
{
  // extract slave-side element associated with current cell
  MORTAR::MortarElement* slaveelement = dynamic_cast<MORTAR::MortarElement*>(idiscret.gElement(cell.GetSlaveId()));
  if(!slaveelement)
    dserror("Couldn't extract slave element from mortar interface discretization!");

  // extract master-side element associated with current cell
  MORTAR::MortarElement* masterelement = dynamic_cast<MORTAR::MortarElement*>(idiscret.gElement(cell.GetMasterId()));
  if(!masterelement)
    dserror("Couldn't extract master element from mortar interface discretization!");

  // safety check
  if(!slaveelement->IsSlave() or masterelement->IsSlave())
    dserror("Something is wrong with the slave-master element pairing!");

  // construct slave-side and master-side location arrays
  DRT::Element::LocationArray la_slave(idiscret.NumDofSets());
  slaveelement->LocationVector(idiscret,la_slave,false);
  DRT::Element::LocationArray la_master(idiscret.NumDofSets());
  masterelement->LocationVector(idiscret,la_master,false);

  // evaluate single mortar integration cell
  MortarCellCalc(
      *slaveelement,
      *masterelement,
      cell,
      condition,
      iphinp,
      idiscret,
      la_slave,
      la_master
      );

  return;
}


/*----------------------------------------------------------------------------------------*
 | evaluate single mortar integration cell                                     fang 01/16 |
 *----------------------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2I::MortarCellCalc(
    MORTAR::MortarElement&         slaveelement,    //!< slave-side mortar element
    MORTAR::MortarElement&         masterelement,   //!< master-side mortar element
    MORTAR::IntCell&               cell,            //!< mortar integration cell
    DRT::Condition&                condition,       //!< scatra-scatra interface coupling condition
    const Epetra_Vector&           iphinp,          //!< interface state vector
    const DRT::Discretization&     idiscret,        //!< interface discretization
    DRT::Element::LocationArray&   la_slave,        //!< slave-side location array
    DRT::Element::LocationArray&   la_master        //!< master-side location array
    ) const
{
  switch(slaveelement.Shape())
  {
    case DRT::Element::DiscretizationType::tri3:
    {
      MortarCellCalc<DRT::Element::DiscretizationType::tri3>(
          slaveelement,
          masterelement,
          cell,
          condition,
          iphinp,
          idiscret,
          la_slave,
          la_master
          );
      break;
    }

    default:
    {
      dserror("Invalid slave-side discretization type!");
      break;
    }
  }

  return;
}


/*---------------------------------------------------------------------------------------------------*
 | evaluate single mortar integration cell of particular slave-side discretization type   fang 01/16 |
 *---------------------------------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeS>
void SCATRA::MeshtyingStrategyS2I::MortarCellCalc(
    MORTAR::MortarElement&         slaveelement,    //!< slave-side mortar element
    MORTAR::MortarElement&         masterelement,   //!< master-side mortar element
    MORTAR::IntCell&               cell,            //!< mortar integration cell
    DRT::Condition&                condition,       //!< scatra-scatra interface coupling condition
    const Epetra_Vector&           iphinp,          //!< interface state vector
    const DRT::Discretization&     idiscret,        //!< interface discretization
    DRT::Element::LocationArray&   la_slave,        //!< slave-side location array
    DRT::Element::LocationArray&   la_master        //!< master-side location array
    ) const
{
  switch(masterelement.Shape())
  {
    case DRT::Element::DiscretizationType::tri3:
    {
      SCATRA::MortarCellCalc<distypeS,DRT::Element::DiscretizationType::tri3>::Instance()->Evaluate(
          *islavematrix_,
          *imastermatrix_,
          *islaveresidual_,
          *imasterresidual_,
          iphinp,
          idiscret,
          condition,
          slaveelement,
          masterelement,
          cell,
          la_slave,
          la_master
          );

      break;
    }

    default:
    {
      dserror("Invalid master-side discretization type!");
      break;
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 | perform setup of scatra-scatra interface coupling         fang 10/14 |
 *----------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2I::InitMeshtying()
{
  // extract scatra-scatra coupling conditions from discretization
  std::vector<DRT::Condition*> conditions(0,NULL);
  scatratimint_->Discretization()->GetCondition("S2ICoupling",conditions);
  slaveconditions_.clear();
  std::map<const int,DRT::Condition* const> masterconditions;
  for(unsigned icondition=0; icondition<conditions.size(); ++icondition)
  {
    DRT::Condition* const condition = conditions[icondition];
    const int condid = condition->GetInt("ConditionID");

    switch(condition->GetInt("interface side"))
    {
      case INPAR::S2I::side_slave:
      {
        if(slaveconditions_.find(condid) == slaveconditions_.end())
          slaveconditions_.insert(std::pair<const int,DRT::Condition* const>(condid,condition));
        else
          dserror("Cannot have multiple slave-side scatra-scatra interface coupling conditions with the same ID!");
        break;
      }

      case INPAR::S2I::side_master:
      {
        if(masterconditions.find(condid) == masterconditions.end())
          masterconditions.insert(std::pair<const int,DRT::Condition* const>(condid,condition));
        else
          dserror("Cannot have multiple master-side scatra-scatra interface coupling conditions with the same ID!");
        break;
      }

      default:
      {
        dserror("Invalid scatra-scatra interface coupling condition!");
        break;
      }
    }
  }

  // determine type of mortar meshtying
  switch(mortartype_)
  {
  // setup scatra-scatra interface coupling for interfaces with pairwise overlapping interface nodes
  case INPAR::S2I::mortar_none:
  {
    // overwrite IDs of master-side scatra-scatra interface coupling conditions with the value -1
    // to prevent them from being evaluated when calling EvaluateCondition on the discretization
    for(std::map<const int,DRT::Condition* const>::iterator imastercondition=masterconditions.begin(); imastercondition!=masterconditions.end(); ++imastercondition)
      imastercondition->second->Add("ConditionID",-1);

    // initialize int vectors for global ids of slave and master interface nodes
    std::vector<int> islavenodegidvec;
    std::vector<int> imasternodegidvec;

    // fill vectors
    for (std::map<const int,DRT::Condition* const>::iterator islavecondition=slaveconditions_.begin(); islavecondition!=slaveconditions_.end(); ++islavecondition)
    {
      const std::vector<int>* islavenodegids = islavecondition->second->Nodes();

      for (unsigned islavenode=0; islavenode<islavenodegids->size(); ++islavenode)
      {
        const int islavenodegid = (*islavenodegids)[islavenode];

        // insert global id of current node into associated vector only if node is owned by current processor
        // need to make sure that node is stored on current processor, otherwise cannot resolve "->Owner()"
        if(scatratimint_->Discretization()->HaveGlobalNode(islavenodegid) and scatratimint_->Discretization()->gNode(islavenodegid)->Owner() == scatratimint_->Discretization()->Comm().MyPID())
          islavenodegidvec.push_back(islavenodegid);
      }
    }
    for (std::map<const int,DRT::Condition* const>::iterator imastercondition=masterconditions.begin(); imastercondition!=masterconditions.end(); ++imastercondition)
    {
      const std::vector<int>* imasternodegids = imastercondition->second->Nodes();

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
    if(not slaveonly_)
    {
      islavetomastercoltransform_ = Teuchos::rcp(new FSI::UTILS::MatrixColTransform);
      islavetomasterrowtransform_ = Teuchos::rcp(new FSI::UTILS::MatrixRowTransform);
      islavetomasterrowcoltransform_ = Teuchos::rcp(new FSI::UTILS::MatrixRowColTransform);
    }
    else
      imastertoslaverowtransform_ = Teuchos::rcp(new FSI::UTILS::MatrixRowTransform);

    // initialize auxiliary residual vector
    islaveresidual_ = Teuchos::rcp(new Epetra_Vector(*(icoup_->SlaveDofMap())));

    // further initializations depending on type of global system matrix
    switch(matrixtype_)
    {
    case INPAR::S2I::matrix_sparse:
      // nothing needs to be done in this case
      break;

    case INPAR::S2I::matrix_block_condition:
    case INPAR::S2I::matrix_block_condition_dof:
    case INPAR::S2I::matrix_block_geometry:
    {
      // safety check
      if(!scatratimint_->Solver()->Params().isSublist("AMGnxn Parameters"))
        dserror("Global system matrix with block structure requires AMGnxn block preconditioner!");

      // initialize map extractors associated with blocks of global system matrix
      BuildBlockMapExtractors();

      // feed AMGnxn block preconditioner with null space information for each block of global block system matrix
      BuildBlockNullSpaces();

      break;
    }
    }

    break;
  }

  // setup scatra-scatra interface coupling for interfaces with non-overlapping interface nodes
  case INPAR::S2I::mortar_standard:
  case INPAR::S2I::mortar_saddlepoint:
  case INPAR::S2I::mortar_condensed:
  {
    // safety checks
    if(DRT::INPUT::IntegralValue<INPAR::MORTAR::ParRedist>(DRT::Problem::Instance()->MortarCouplingParams(),"PARALLEL_REDIST") != INPAR::MORTAR::parredist_none)
      dserror("Parallel redistribution not yet implemented for scatra-scatra interface coupling!");
    if(DRT::INPUT::IntegralValue<INPAR::MORTAR::MeshRelocation>(DRT::Problem::Instance()->MortarCouplingParams(),"MESH_RELOCATION") != INPAR::MORTAR::relocation_none)
      dserror("Mesh relocation not yet implemented for scatra-scatra interface coupling!");

    // initialize empty interface maps
    Teuchos::RCP<Epetra_Map> imastermap = Teuchos::rcp(new Epetra_Map(0,0,scatratimint_->Discretization()->Comm()));
    Teuchos::RCP<Epetra_Map> islavemap = Teuchos::rcp(new Epetra_Map(0,0,scatratimint_->Discretization()->Comm()));
    Teuchos::RCP<Epetra_Map> ifullmap = Teuchos::rcp(new Epetra_Map(0,0,scatratimint_->Discretization()->Comm()));

    // loop over all slave-side scatra-scatra interface coupling conditions
    for(std::map<const int,DRT::Condition* const>::iterator islavecondition=slaveconditions_.begin(); islavecondition!=slaveconditions_.end(); ++islavecondition)
    {
      // initialize maps for row nodes associated with current condition
      std::map<int,DRT::Node*> masternodes;
      std::map<int,DRT::Node*> slavenodes;

      // initialize maps for column nodes associated with current condition
      std::map<int,DRT::Node*> mastergnodes;
      std::map<int,DRT::Node*> slavegnodes;

      // initialize maps for elements associated with current condition
      std::map<int,Teuchos::RCP<DRT::Element> > masterelements;
      std::map<int,Teuchos::RCP<DRT::Element> > slaveelements;

      // extract current slave-side and associated master-side scatra-scatra interface coupling conditions
      std::vector<DRT::Condition*> mastercondition(1,masterconditions[islavecondition->first]);
      std::vector<DRT::Condition*> slavecondition(1,islavecondition->second);

      // fill maps
      DRT::UTILS::FindConditionObjects(*scatratimint_->Discretization(),masternodes,mastergnodes,masterelements,mastercondition);
      DRT::UTILS::FindConditionObjects(*scatratimint_->Discretization(),slavenodes,slavegnodes,slaveelements,slavecondition);

      // initialize mortar coupling adapter
      icoupmortar_[islavecondition->first] = Teuchos::rcp(new ADAPTER::CouplingMortar());
      std::vector<int> coupleddof(scatratimint_->NumDofPerNode(),1);
      icoupmortar_[islavecondition->first]->SetupInterface(
          scatratimint_->Discretization(),
          scatratimint_->Discretization(),
          Teuchos::null,
          coupleddof,
          mastergnodes,
          slavegnodes,
          masterelements,
          slaveelements,
          scatratimint_->Discretization()->Comm()
          );

      // generate mortar integration cells
      imortarcells_[islavecondition->first] = icoupmortar_[islavecondition->first]->EvaluateGeometry();

      // build interface maps
      imastermap = LINALG::MergeMap(imastermap,icoupmortar_[islavecondition->first]->MasterDofRowMap(),false);
      islavemap = LINALG::MergeMap(islavemap,icoupmortar_[islavecondition->first]->SlaveDofRowMap(),false);
      ifullmap = LINALG::MergeMap(ifullmap,LINALG::MergeMap(icoupmortar_[islavecondition->first]->MasterDofRowMap(),icoupmortar_[islavecondition->first]->SlaveDofRowMap(),false),false);
    }

    // generate interior and interface maps
    std::vector<Teuchos::RCP<const Epetra_Map> > imaps;
    imaps.push_back(LINALG::SplitMap(*(scatratimint_->Discretization()->DofRowMap()),*ifullmap));
    imaps.push_back(islavemap);
    imaps.push_back(imastermap);

    // initialize global map extractor
    interfacemaps_ = Teuchos::rcp(new LINALG::MultiMapExtractor(*(scatratimint_->Discretization()->DofRowMap()),imaps));
    interfacemaps_->CheckForValidMapExtractor();

    // initialize auxiliary system matrices
    islavematrix_ = Teuchos::rcp(new LINALG::SparseMatrix(*interfacemaps_->Map(1),81));
    imastermatrix_ = Teuchos::rcp(new LINALG::SparseMatrix(*interfacemaps_->Map(2),81,true,false,LINALG::SparseMatrix::FE_MATRIX));

    // initialize auxiliary residual vectors
    islaveresidual_ = Teuchos::rcp(new Epetra_Vector(*interfacemaps_->Map(1)));
    imasterresidual_ = Teuchos::rcp(new Epetra_FEVector(*interfacemaps_->Map(2)));

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


/*----------------------------------------------------------------------------------*
 | build map extractors associated with blocks of global system matrix   fang 07/15 |
 *----------------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2I::BuildBlockMapExtractors()
{
  if(matrixtype_ == INPAR::S2I::matrix_block_condition or matrixtype_ == INPAR::S2I::matrix_block_condition_dof)
  {
    // extract domain partitioning conditions from discretization
    std::vector<Teuchos::RCP<DRT::Condition> > partitioningconditions;
    scatratimint_->Discretization()->GetCondition("S2ICouplingPartitioning",partitioningconditions);

    // safety check
    if(!partitioningconditions.size())
      dserror("For block preconditioning based on domain partitioning, at least one associated condition needs to be specified in the input file!");

    // build maps associated with blocks of global system matrix
    std::vector<Teuchos::RCP<const Epetra_Map> > blockmaps;
    BuildBlockMaps(partitioningconditions,blockmaps);

    // initialize full map extractor associated with blocks of global system matrix
    blockmaps_ = Teuchos::rcp(new LINALG::MultiMapExtractor(*(scatratimint_->Discretization()->DofRowMap()),blockmaps));
    blockmaps_->CheckForValidMapExtractor();

    // initialize reduced interface map extractors associated with blocks of global system matrix
    const unsigned nblocks = blockmaps.size();
    std::vector<Teuchos::RCP<const Epetra_Map> > blockmaps_slave(nblocks);
    std::vector<Teuchos::RCP<const Epetra_Map> > blockmaps_master(nblocks);
    for(unsigned iblock=0; iblock<nblocks; ++iblock)
    {
      std::vector<Teuchos::RCP<const Epetra_Map> > maps(2);
      maps[0] = blockmaps[iblock];
      maps[1] = icoup_->SlaveDofMap();
      blockmaps_slave[iblock] = LINALG::MultiMapExtractor::IntersectMaps(maps);
      maps[1] = icoup_->MasterDofMap();
      blockmaps_master[iblock] = LINALG::MultiMapExtractor::IntersectMaps(maps);
    }
    blockmaps_slave_ = Teuchos::rcp(new LINALG::MultiMapExtractor(*icoup_->SlaveDofMap(),blockmaps_slave));
    blockmaps_slave_->CheckForValidMapExtractor();
    blockmaps_master_ = Teuchos::rcp(new LINALG::MultiMapExtractor(*icoup_->MasterDofMap(),blockmaps_master));
    blockmaps_master_->CheckForValidMapExtractor();
  }

  else if(matrixtype_ == INPAR::S2I::matrix_block_geometry)
    // matrix block map extractor equals interface map extractor in this case
    blockmaps_ = interfacemaps_;

  return;
} // SCATRA::MeshtyingStrategyS2I::BuildBlockMapExtractors


/*----------------------------------------------------------------------------*
 | build maps associated with blocks of global system matrix       fang 06/15 |
 *----------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2I::BuildBlockMaps(
    const std::vector<Teuchos::RCP<DRT::Condition> >&   partitioningconditions,   //! domain partitioning conditions
    std::vector<Teuchos::RCP<const Epetra_Map> >&       blockmaps                 //! empty vector for maps to be built
    ) const
{
  if(matrixtype_ == INPAR::S2I::matrix_block_condition)
  {
    // extract number of domain partitioning conditions
    const unsigned ncond = partitioningconditions.size();

    // prepare vector for maps to be built
    blockmaps.resize(ncond,Teuchos::null);

    // loop over all domain partitioning conditions
    for(unsigned icond=0; icond<ncond; ++icond)
    {
      // initialize set for dof IDs associated with current partitioning condition
      std::set<int> dofids;

      // extract nodes associated with current domain partitioning condition
      const std::vector<int>* nodegids = partitioningconditions[icond]->Nodes();

      // loop over all nodes associated with current domain partitioning condition
      for (unsigned inode=0; inode<nodegids->size(); ++inode)
      {
        // extract global ID of current node
        const int nodegid = (*nodegids)[inode];

        // consider current node only if node is owned by current processor
        // need to make sure that node is stored on current processor, otherwise cannot resolve "->Owner()"
        if(scatratimint_->Discretization()->HaveGlobalNode(nodegid) and scatratimint_->Discretization()->gNode(nodegid)->Owner() == scatratimint_->Discretization()->Comm().MyPID())
        {
          // add dof IDs associated with current node to corresponding set
          const std::vector<int> nodedofs = scatratimint_->Discretization()->Dof(scatratimint_->Discretization()->gNode(nodegid));
          std::copy(nodedofs.begin(),nodedofs.end(),std::inserter(dofids,dofids.end()));
        }
      }

      // transform set for dof IDs into vector and then into Epetra map
      int nummyelements(0);
      int* myglobalelements(NULL);
      std::vector<int> dofidvec;
      if(dofids.size() > 0)
      {
        dofidvec.reserve(dofids.size());
        dofidvec.assign(dofids.begin(),dofids.end());
        nummyelements = dofidvec.size();
        myglobalelements = &(dofidvec[0]);
      }
      blockmaps[icond] = Teuchos::rcp(new Epetra_Map(-1,nummyelements,myglobalelements,scatratimint_->DofRowMap()->IndexBase(),scatratimint_->DofRowMap()->Comm()));
    }
  }

  // safety check
  else
    dserror("Invalid type of global system matrix!");

  return;
} // SCATRA::MeshtyingStrategyS2I::BuildBlockMaps


/*-------------------------------------------------------------------------------*
 | build null spaces associated with blocks of global system matrix   fang 07/15 |
 *-------------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2I::BuildBlockNullSpaces() const
{
  // loop over blocks of global system matrix
  for(int iblock=0; iblock<blockmaps_->NumMaps(); ++iblock)
  {
    // store number of current block as string, starting from 1
    std::ostringstream iblockstr;
    iblockstr << iblock+1;

    // equip smoother for current matrix block with empty parameter sublists to trigger null space computation
    Teuchos::ParameterList& blocksmootherparams = scatratimint_->Solver()->Params().sublist("Inverse"+iblockstr.str());
    blocksmootherparams.sublist("Aztec Parameters");
    blocksmootherparams.sublist("MueLu Parameters");

    // equip smoother for current matrix block with null space associated with all degrees of freedom on discretization
    scatratimint_->Discretization()->ComputeNullSpaceIfNecessary(blocksmootherparams);

    // reduce full null space to match degrees of freedom associated with current matrix block
    LINALG::Solver::FixMLNullspace("Block "+iblockstr.str(),*scatratimint_->Discretization()->DofRowMap(),*blockmaps_->Map(iblock),blocksmootherparams);
  }

  return;
} // SCATRA::MeshtyingStrategyS2I::BuildBlockNullSpaces


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
  case INPAR::S2I::matrix_block_condition_dof:
  {
    // initialize system matrix and associated strategy
    systemmatrix = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(*blockmaps_,*blockmaps_));

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
    case INPAR::S2I::matrix_block_condition_dof:
    case INPAR::S2I::matrix_block_geometry:
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
          blockmaps_->InsertVector(invrowsums,i,invrowsums_);
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
          blockmaps_->InsertVector(invcolsums,j,invcolsums_);
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


/*----------------------------------------------------------------------*
 | singleton access method                                   fang 01/16 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeS,DRT::Element::DiscretizationType distypeM>
SCATRA::MortarCellCalc<distypeS,distypeM>* SCATRA::MortarCellCalc<distypeS,distypeM>::Instance(bool create)
{
  static MortarCellCalc<distypeS,distypeM>* instance;

  if(create)
  {
    if(instance == NULL)
      instance = new MortarCellCalc<distypeS,distypeM>();
  }

  else if(instance != NULL)
  {
    delete instance;
    instance = NULL;
  }

  return instance;
}


/*----------------------------------------------------------------------*
 | singleton destruction                                     fang 01/16 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeS,DRT::Element::DiscretizationType distypeM>
void SCATRA::MortarCellCalc<distypeS,distypeM>::Done()
{
  // delete singleton
  Instance(false);

  return;
}


/*--------------------------------------------------------------------------------------------------------------------*
 | evaluate single mortar integration cell of particular slave-side and master-side discretization types   fang 01/16 |
 *--------------------------------------------------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeS,DRT::Element::DiscretizationType distypeM>
void SCATRA::MortarCellCalc<distypeS,distypeM>::Evaluate(
    LINALG::SparseMatrix&          islavematrix,      //!< linearizations of slave-side residuals
    LINALG::SparseMatrix&          imastermatrix,     //!< linearizations of master-side residuals
    Epetra_Vector&                 islaveresidual,    //!< slave-side residual vector
    Epetra_FEVector&               imasterresidual,   //!< master-side residual vector
    const Epetra_Vector&           iphinp,            //!< interface state vector
    const DRT::Discretization&     idiscret,          //!< interface discretization
    DRT::Condition&                condition,         //!< scatra-scatra interface coupling condition
    MORTAR::MortarElement&         slaveelement,      //!< slave-side mortar element
    MORTAR::MortarElement&         masterelement,     //!< master-side mortar element
    MORTAR::IntCell&               cell,              //!< mortar integration cell
    DRT::Element::LocationArray&   la_slave,          //!< slave-side location array
    DRT::Element::LocationArray&   la_master          //!< master-side location array
    ) const
{
  // extract nodal state variables associated with slave and master elements
  static const int nen_slave = DRT::UTILS::DisTypeToNumNodePerEle<distypeS>::numNodePerElement;
  static const int nen_master = DRT::UTILS::DisTypeToNumNodePerEle<distypeM>::numNodePerElement;
  std::vector<LINALG::Matrix<nen_slave,1> > ephinp_slave(slaveelement.NumDofPerNode(*slaveelement.Nodes()[0]),LINALG::Matrix<nen_slave,1>(true));
  std::vector<LINALG::Matrix<nen_master,1> > ephinp_master(masterelement.NumDofPerNode(*masterelement.Nodes()[0]),LINALG::Matrix<nen_master,1>(true));
  ExtractNodeValues(ephinp_slave,ephinp_master,iphinp,la_slave,la_master);

  // evaluate and assemble interface linearizations and residuals
  CalcMatAndRhs(
      islavematrix,
      imastermatrix,
      islaveresidual,
      imasterresidual,
      idiscret,
      condition,
      slaveelement,
      masterelement,
      cell,
      ephinp_slave,
      ephinp_master,
      la_slave,
      la_master
      );

  return;
}


/*--------------------------------------------------------------------------------------*
 | protected constructor for singletons                                      fang 01/16 |
 *--------------------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeS,DRT::Element::DiscretizationType distypeM>
SCATRA::MortarCellCalc<distypeS,distypeM>::MortarCellCalc()
{
  // safety check
  if(DRT::UTILS::DisTypeToDim<distypeS>::dim != 2 or DRT::UTILS::DisTypeToDim<distypeM>::dim != 2)
    dserror("Scatra-scatra interface coupling with non-matching interface discretization currently only implemented for two-dimensional interface manifolds!");

  return;
}


/*--------------------------------------------------------------------------------------*
 | extract nodal state variables associated with slave and master elements   fang 01/16 |
 *--------------------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeS,DRT::Element::DiscretizationType distypeM>
void SCATRA::MortarCellCalc<distypeS,distypeM>::ExtractNodeValues(
    std::vector<LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distypeS>::numNodePerElement,1> >&   ephinp_slave,    //!< state variables at slave-side nodes
    std::vector<LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distypeM>::numNodePerElement,1> >&   ephinp_master,   //!< state variables at master-side nodes
    const Epetra_Vector&                                                                               iphinp,          //!< interface state vector
    DRT::Element::LocationArray&                                                                       la_slave,        //!< slave-side location array
    DRT::Element::LocationArray&                                                                       la_master        //!< master-side location array
    ) const
{
  // extract nodal state variables associated with slave and master elements
  DRT::UTILS::ExtractMyValues<LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distypeS>::numNodePerElement,1> >(iphinp,ephinp_slave,la_slave[0].lm_);
  DRT::UTILS::ExtractMyValues<LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distypeM>::numNodePerElement,1> >(iphinp,ephinp_master,la_master[0].lm_);

  return;
}


/*------------------------------------------------------------------------------------------*
 | evaluate shape functions and domain integration factor at integration point   fang 01/16 |
 *------------------------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeS,DRT::Element::DiscretizationType distypeM>
const double SCATRA::MortarCellCalc<distypeS,distypeM>::EvalShapeFuncAndDomIntFacAtIntPoint(
    LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distypeS>::numNodePerElement,1>&   funct_slave,     //!< slave-side shape function values
    LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distypeM>::numNodePerElement,1>&   funct_master,    //!< master-side shape function values
    MORTAR::MortarElement&                                                               slaveelement,    //!< slave-side mortar element
    MORTAR::MortarElement&                                                               masterelement,   //!< master-side mortar element
    MORTAR::IntCell&                                                                     cell,            //!< mortar integration cell
    const DRT::UTILS::IntPointsAndWeights<DRT::UTILS::DisTypeToDim<distypeS>::dim>&      intpoints,       //!< quadrature rule
    const int                                                                            iquad            //!< ID of integration point
    ) const
{
  // number of reference dimensions
  static const int dim_ref = DRT::UTILS::DisTypeToDim<distypeS>::dim;

  // reference coordinates of integration point
  double coordinates_ref[dim_ref] = {};
  for(int idim=0; idim<dim_ref; ++idim)
    coordinates_ref[idim] = intpoints.IP().qxg[iquad][idim];

  // global coordinates of integration point
  double coordinates_global[dim_ref+1] = {};
  cell.LocalToGlobal(coordinates_ref,coordinates_global,0);

  // project integration point onto slave and master elements
  double coordinates_slave[dim_ref] = {};
  double coordinates_master[dim_ref] = {};
  double dummy(0.);
  MORTAR::MortarProjector::Impl(slaveelement)->ProjectGaussPointAuxn3D(coordinates_global,cell.Auxn(),slaveelement,coordinates_slave,dummy);
  MORTAR::MortarProjector::Impl(masterelement)->ProjectGaussPointAuxn3D(coordinates_global,cell.Auxn(),masterelement,coordinates_master,dummy);

  // evaluate shape functions at current integration point on slave and master elements
  MORTAR::UTILS::EvaluateShape_Displ(coordinates_slave,funct_slave,slaveelement,false);
  MORTAR::UTILS::EvaluateShape_Displ(coordinates_master,funct_master,masterelement,false);

  // integration weight
  const double weight = intpoints.IP().qwgt[iquad];

  // Jacobian determinant
  const double jacobian = cell.Jacobian(coordinates_ref);

  // domain integration factor
  return jacobian*weight;
}


/*---------------------------------------------------------------------------*
 | evaluate and assemble interface linearizations and residuals   fang 01/16 |
 *---------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeS,DRT::Element::DiscretizationType distypeM>
void SCATRA::MortarCellCalc<distypeS,distypeM>::CalcMatAndRhs(
    LINALG::SparseMatrix&                                                                              islavematrix,      //!< linearizations of slave-side residuals
    LINALG::SparseMatrix&                                                                              imastermatrix,     //!< linearizations of master-side residuals
    Epetra_Vector&                                                                                     islaveresidual,    //!< slave-side residual vector
    Epetra_FEVector&                                                                                   imasterresidual,   //!< master-side residual vector
    const DRT::Discretization&                                                                         idiscret,          //!< interface discretization
    DRT::Condition&                                                                                    condition,         //!< scatra-scatra interface coupling condition
    MORTAR::MortarElement&                                                                             slaveelement,      //!< slave-side mortar element
    MORTAR::MortarElement&                                                                             masterelement,     //!< master-side mortar element
    MORTAR::IntCell&                                                                                   cell,              //!< mortar integration cell
    std::vector<LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distypeS>::numNodePerElement,1> >&   ephinp_slave,      //!< state variables at slave-side nodes
    std::vector<LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distypeM>::numNodePerElement,1> >&   ephinp_master,     //!< state variables at master-side nodes
    DRT::Element::LocationArray&                                                                       la_slave,          //!< slave-side location array
    DRT::Element::LocationArray&                                                                       la_master          //!< master-side location array
    ) const
{
  // TODO: implement and activate!
  dserror("Not yet implemented!");

  return;
}


// forward declaration
template class SCATRA::MortarCellCalc<DRT::Element::DiscretizationType::tri3,DRT::Element::DiscretizationType::tri3>;

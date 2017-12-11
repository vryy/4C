/*--------------------------------------------------------------------------*/
/*!
\file ssi_monolithic.cpp

\brief monolithic scalar-structure interaction

\level 2

<pre>
\maintainer Rui Fang & Christoph Schmidt
            {fang,schmidt}@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>
*/
/*--------------------------------------------------------------------------*/
#include "ssi_monolithic.H"
#include "ssi_coupling.H"
#include "ssi_monolithic_convcheck_strategies.H"
#include "ssi_monolithic_resulttest.H"

#include <Epetra_Time.h>

#include "../drt_adapter/ad_str_ssiwrapper.H"
#include "../drt_adapter/ad_str_structure_new.H"
#include "../drt_adapter/adapter_coupling.H"
#include "../drt_adapter/adapter_scatra_base_algorithm.H"

#include "../drt_fsi/fsi_matrixtransform.H"

#include "../drt_inpar/inpar_scatra.H"
#include "../drt_inpar/inpar_ssi.H"

#include "../drt_io/io_control.H"

#include "../drt_lib/drt_assemblestrategy.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_locsys.H"

#include "../drt_scatra/scatra_timint_implicit.H"
#include "../drt_scatra/scatra_timint_meshtying_strategy_s2i.H"

#include "../drt_scatra_ele/scatra_ele_action.H"

#include "../linalg/linalg_mapextractor.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_utils.H"

/*--------------------------------------------------------------------------*
 | constructor                                                   fang 08/17 |
 *--------------------------------------------------------------------------*/
SSI::SSI_Mono::SSI_Mono(
    const Epetra_Comm&              comm,              //!< communicator
    const Teuchos::ParameterList&   globaltimeparams   //!< parameter list for time integration
    )
    // call base class constructor
  : SSI_Base(comm, globaltimeparams),

    // initialize member variables
    dtele_(0.),
    dtsolve_(0.),
    icoup_structure_(Teuchos::null),
    iter_(0),
    map_structure_(Teuchos::null),
    maps_(Teuchos::null),
    maps_structure_(Teuchos::null),
    matrixtype_(DRT::INPUT::IntegralValue<INPAR::SSI::MatrixType>(globaltimeparams.sublist("MONOLITHIC"),"MATRIXTYPE")),
    residual_(Teuchos::null),
    scatrastructureblock_(Teuchos::null),
    solver_(Teuchos::rcp(new LINALG::Solver(
        DRT::Problem::Instance()->SolverParams(globaltimeparams.sublist("MONOLITHIC").get<int>("LINEAR_SOLVER")),
        comm,
        DRT::Problem::Instance()->ErrorFile()->Handle()
        ))),
    strategy_convcheck_(Teuchos::null),
    strategy_scatra_(Teuchos::null),
    structurescatrablock_(Teuchos::null),
    systemmatrix_(Teuchos::null),
    timer_(Teuchos::rcp(new Epetra_Time(comm)))
{
  return;
}


/*--------------------------------------------------------------------------*
 | assemble global system of equations                           fang 08/17 |
 *--------------------------------------------------------------------------*/
void SSI::SSI_Mono::AssembleMatAndRHS()
{
  // pass scalar transport degrees of freedom to structural discretization
  SetScatraSolution(scatra_->ScaTraField()->Phinp());

  // build system matrix and residual for structure field
  if(iter_ == 1)
    structure_->Evaluate();
  else
  {
    // set up structural increment vector
    const Teuchos::RCP<Epetra_Vector> increment_structure = maps_->ExtractVector(increment_,1);

    // consider structural meshtying if applicable
    if(scatra_->ScaTraField()->S2ICoupling())
    {
      maps_structure_->InsertVector(icoup_structure_->MasterToSlave(maps_structure_->ExtractVector(structure_->Dispnp(),2)),1,structure_->WriteAccessDispnp());
      structure_->SetState(structure_->WriteAccessDispnp());
      maps_structure_->InsertVector(icoup_structure_->MasterToSlave(maps_structure_->ExtractVector(increment_structure,2)),1,increment_structure);
    }

    // evaluate structural field
    structure_->Evaluate(increment_structure);
  }

  // pass structural degrees of freedom to scalar transport discretization
  SetStructSolution(structure_->Dispnp(),structure_->Velnp());

  // build system matrix and residual for scalar transport field
  scatra_->ScaTraField()->PrepareLinearSolve();

  // assemble off-diagonal scatra-structure block of global system matrix
  AssembleODBlockScatraStructure();

  //! assemble off-diagonal structure-scatra block of global system matrix
  AssembleODBlockStructureScatra();

  // build global system matrix
  switch(matrixtype_)
  {
    case INPAR::SSI::matrix_sparse:
    {
      // check global system matrix
      Teuchos::RCP<LINALG::SparseMatrix> systemmatrix = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(systemmatrix_);
      if(systemmatrix == Teuchos::null)
        dserror("System matrix is not a sparse matrix!");

      // add scalar transport system matrix to global system matrix
      systemmatrix->Add(*scatra_->ScaTraField()->SystemMatrix(),false,1.,0.);

      // perform structural meshtying before adding remaining matrix blocks
      if(scatra_->ScaTraField()->S2ICoupling())
      {
        // assemble interior and master-side columns of scatra-structure block into global system matrix
        const LINALG::SparseMatrix& scatrastructureblock = *Teuchos::rcp_dynamic_cast<const LINALG::SparseMatrix>(scatrastructureblock_);
        FSI::UTILS::MatrixLogicalSplitAndTransform()(
            scatrastructureblock,
            scatrastructureblock.RangeMap(),
            *map_structure_,
            1.,
            NULL,
            NULL,
            *systemmatrix,
            true,
            true
            );

        // transform and assemble slave-side columns of scatra-structure block into global system matrix
        ADAPTER::CouplingSlaveConverter converter(*icoup_structure_);
        FSI::UTILS::MatrixLogicalSplitAndTransform()(
            scatrastructureblock,
            scatrastructureblock.RangeMap(),
            *icoup_structure_->SlaveDofMap(),
            1.,
            NULL,
            &converter,
            *systemmatrix,
            true,
            true
            );

        // assemble interior and master-side rows of structure-scatra block into global system matrix
        const LINALG::SparseMatrix& structurescatrablock = *Teuchos::rcp_dynamic_cast<const LINALG::SparseMatrix>(structurescatrablock_);
        FSI::UTILS::MatrixLogicalSplitAndTransform()(
            structurescatrablock,
            *map_structure_,
            structurescatrablock.DomainMap(),
            1.,
            NULL,
            NULL,
            *systemmatrix,
            true,
            true
            );

        // transform and assemble slave-side rows of structure-scatra block into global system matrix
        FSI::UTILS::MatrixLogicalSplitAndTransform()(
            structurescatrablock,
            *icoup_structure_->SlaveDofMap(),
            structurescatrablock.DomainMap(),
            1.,
            &converter,
            NULL,
            *systemmatrix,
            true,
            true
            );

        // assemble interior and master-side rows and columns of structural system matrix into global system matrix
        FSI::UTILS::MatrixLogicalSplitAndTransform()(
            *structure_->SystemMatrix(),
            *map_structure_,
            *map_structure_,
            1.,
            NULL,
            NULL,
            *systemmatrix,
            true,
            true
            );

        // transform and assemble slave-side rows of structural system matrix into global system matrix
        FSI::UTILS::MatrixLogicalSplitAndTransform()(
            *structure_->SystemMatrix(),
            *icoup_structure_->SlaveDofMap(),
            *map_structure_,
            1.,
            &converter,
            NULL,
            *systemmatrix,
            true,
            true
            );

        // transform and assemble slave-side columns of structural system matrix into global system matrix
        FSI::UTILS::MatrixLogicalSplitAndTransform()(
            *structure_->SystemMatrix(),
            *map_structure_,
            *icoup_structure_->SlaveDofMap(),
            1.,
            NULL,
            &converter,
            *systemmatrix,
            true,
            true
            );

        // transform and assemble slave-side rows and columns of structural system matrix into global system matrix
        FSI::UTILS::MatrixLogicalSplitAndTransform()(
            *structure_->SystemMatrix(),
            *icoup_structure_->SlaveDofMap(),
            *icoup_structure_->SlaveDofMap(),
            1.,
            &converter,
            &converter,
            *systemmatrix,
            true,
            true
            );

        // subject slave-side rows of structural system matrix to pseudo Dirichlet conditions to finalize structural meshtying
        const double value(1.);
        for(int doflid_slave=0; doflid_slave<icoup_structure_->SlaveDofMap()->NumMyElements(); ++doflid_slave)
        {
          // extract global ID of current slave-side row
          const int dofgid_slave = icoup_structure_->SlaveDofMap()->GID(doflid_slave);
          if(dofgid_slave < 0)
            dserror("Local ID not found!");

          // apply pseudo Dirichlet conditions to filled matrix, i.e., to local row and column indices
          if(systemmatrix->Filled())
          {
            const int rowlid_slave = systemmatrix->RowMap().LID(dofgid_slave);
            if(rowlid_slave < 0)
              dserror("Global ID not found!");
            if(systemmatrix->EpetraMatrix()->ReplaceMyValues(rowlid_slave,1,&value,&rowlid_slave))
              dserror("ReplaceMyValues failed!");
          }

          // apply pseudo Dirichlet conditions to unfilled matrix, i.e., to global row and column indices
          else if(systemmatrix->EpetraMatrix()->InsertGlobalValues(dofgid_slave,1,&value,&dofgid_slave))
            dserror("InsertGlobalValues failed!");
        }
      }

      // add remaining matrix blocks directly
      else
      {
        systemmatrix->Add(*scatrastructureblock_,false,1.,1.);
        systemmatrix->Add(*structurescatrablock_,false,1.,1.);
        systemmatrix->Add(*structure_->SystemMatrix(),false,1.,1.);
      }

      break;
    }

    default:
    {
      dserror("Type of global system matrix for scalar-structure interaction not recognized!");
      break;
    }
  }

  // finalize global system matrix
  systemmatrix_->Complete();

  // apply structural Dirichlet conditions
  if(structure_->LocsysManager() == Teuchos::null)
    systemmatrix_->ApplyDirichlet(*structure_->GetDBCMapExtractor()->CondMap());
  else
  {
    switch(matrixtype_)
    {
      case INPAR::SSI::matrix_sparse:
      {
        // check global system matrix
        Teuchos::RCP<LINALG::SparseMatrix> systemmatrix = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(systemmatrix_);
        if(systemmatrix == Teuchos::null)
          dserror("System matrix is not a sparse matrix!");

        // extract structural rows of global system matrix
        const Teuchos::RCP<LINALG::SparseMatrix> systemmatrix_structure = Teuchos::rcp(new LINALG::SparseMatrix(
            *structure_->DofRowMap(),
            27,
            false,
            true
            ));
        FSI::UTILS::MatrixLogicalSplitAndTransform()(
            *systemmatrix,
            *structure_->DofRowMap(),
            systemmatrix->DomainMap(),
            1.,
            NULL,
            NULL,
            *systemmatrix_structure
            );
        systemmatrix_structure->Complete(systemmatrix->DomainMap(),*structure_->DofRowMap());

        // apply structural Dirichlet conditions
        structure_->LocsysManager()->RotateGlobalToLocal(systemmatrix_structure);
        systemmatrix_structure->ApplyDirichletWithTrafo(structure_->LocsysManager()->Trafo(),*structure_->GetDBCMapExtractor()->CondMap());
        structure_->LocsysManager()->RotateLocalToGlobal(systemmatrix_structure);

        // assemble structural rows of global system matrix back into global system matrix
        systemmatrix->Put(*systemmatrix_structure,1.,structure_->DofRowMap());

        break;
      }

      default:
      {
        dserror("Type of global system matrix for scalar-structure interaction not recognized!");
        break;
      }
    }
  }

  // initialize monolithic right-hand side vector
  residual_->PutScalar(0.);

  // assemble scalar transport right-hand side vector into monolithic right-hand side vector
  maps_->InsertVector(scatra_->ScaTraField()->Residual(),0,residual_);

  // perform structural meshtying before assembling structural right-hand side vector into monolithic right-hand side vector
  if(scatra_->ScaTraField()->S2ICoupling())
  {
    // make copy of structural right-hand side vector
    Epetra_Vector residual_structure(*structure_->RHS());

    // transform slave-side part of structural right-hand side vector to master side
    Teuchos::RCP<Epetra_Vector> slavetomaster = maps_structure_->InsertVector(icoup_structure_->SlaveToMaster(maps_structure_->ExtractVector(residual_structure,1)),2);

    // apply Dirichlet conditions to transformed slave-side part of structural right-hand side vector
    const Teuchos::RCP<const Epetra_Vector> zeros = Teuchos::rcp(new Epetra_Vector(slavetomaster->Map()));
    if(structure_->LocsysManager() != Teuchos::null)
      structure_->LocsysManager()->RotateGlobalToLocal(slavetomaster);
    LINALG::ApplyDirichlettoSystem(slavetomaster,zeros,*structure_->GetDBCMapExtractor()->CondMap());
    if(structure_->LocsysManager() != Teuchos::null)
      structure_->LocsysManager()->RotateLocalToGlobal(slavetomaster);

    // assemble transformed slave-side part of structural right-hand side vector
    residual_structure.Update(1.,*slavetomaster,1.);

    // zero out slave-side part of structural right-hand side vector
    maps_structure_->PutScalar(residual_structure,1,0.);

    // assemble final structural right-hand side vector into monolithic right-hand side vector
    maps_->AddVector(residual_structure,1,*residual_,-1.);
  }

  // directly assemble structural right-hand side vector into monolithic right-hand side vector
  else
    maps_->AddVector(structure_->RHS(),1,residual_,-1.);

  return;
}


/*-----------------------------------------------------------------------------------*
 | assemble off-diagonal scatra-structure block of global system matrix   fang 08/17 |
 *-----------------------------------------------------------------------------------*/
void SSI::SSI_Mono::AssembleODBlockScatraStructure() const
{
  // initialize scatra-structure matrix block
  scatrastructureblock_->Zero();

  // create parameter list for element evaluation
  Teuchos::ParameterList eleparams;

  // action for elements
  eleparams.set<int>("action",SCATRA::calc_scatra_mono_odblock_mesh);

  // number of dofset associated with displacement-related dofs on scalar transport discretization
  eleparams.set<int>("ndsdisp",1);

  // number of dofset associated with velocity-related dofs on scalar transport discretization
  eleparams.set<int>("ndsvel",1);

  // remove state vectors from scalar transport discretization
  scatra_->ScaTraField()->Discretization()->ClearState();

  // add state vectors to scalar transport discretization
  scatra_->ScaTraField()->AddTimeIntegrationSpecificVectors();

  // create strategy for assembly of scatra-structure matrix block
  DRT::AssembleStrategy strategyscatrastructure(
      0,                       // row assembly based on number of dofset associated with scalar transport dofs on scalar transport discretization
      1,                       // column assembly based on number of dofset associated with structural dofs on scalar transport discretization
      scatrastructureblock_,   // scatra-structure matrix block
      Teuchos::null,           // no additional matrices or vectors
      Teuchos::null,
      Teuchos::null,
      Teuchos::null
      );

  // assemble scatra-structure matrix block
  scatra_->ScaTraField()->Discretization()->Evaluate(eleparams,strategyscatrastructure);

  // provide scatra-structure matrix block with contributions from scatra-scatra interface coupling if applicable
  if(scatra_->ScaTraField()->S2ICoupling())
  {
    // create parameter list for element evaluation
    Teuchos::ParameterList condparams;

    // action for elements
    condparams.set<int>("action",SCATRA::bd_calc_s2icoupling_od);

    // number of dofset associated with displacement-related dofs on scalar transport discretization
    condparams.set<int>("ndsdisp",1);

    // add state vector containing master-side scatra degrees of freedom to scatra discretization
    scatra_->ScaTraField()->Discretization()->SetState("imasterphinp",strategy_scatra_->MasterPhinp());

    switch(matrixtype_)
    {
      case INPAR::SSI::matrix_sparse:
      {
        // initialize auxiliary system matrix for linearizations of slave-side scatra fluxes w.r.t. slave-side structural dofs
        strategy_scatra_->SlaveMatrix()->Zero();

        // create strategy for assembly of auxiliary system matrix
        DRT::AssembleStrategy strategyscatrastructures2i(
            0,                                 // row assembly based on number of dofset associated with scalar transport dofs on scalar transport discretization
            1,                                 // column assembly based on number of dofset associated with structural dofs on structural discretization
            strategy_scatra_->SlaveMatrix(),   // auxiliary system matrix
            Teuchos::null,                     // no additional matrices of vectors
            Teuchos::null,
            Teuchos::null,
            Teuchos::null
            );

        // evaluate scatra-scatra interface coupling
        std::vector<DRT::Condition*> conditions;
        scatra_->ScaTraField()->Discretization()->GetCondition("S2ICoupling",conditions);
        for(unsigned icondition=0; icondition<conditions.size(); ++icondition)
          if(conditions[icondition]->GetInt("interface side") == INPAR::S2I::side_slave)
            scatra_->ScaTraField()->Discretization()->EvaluateCondition(condparams,strategyscatrastructures2i,"S2ICoupling",conditions[icondition]->GetInt("ConditionID"));

        // finalize auxiliary system matrix
        strategy_scatra_->SlaveMatrix()->Complete(*maps_->Map(1),*maps_->Map(0));

        // assemble linearizations of slave-side scatra fluxes w.r.t. slave-side structural dofs into scatra-structure matrix block
        scatrastructureblock_->Add(*strategy_scatra_->SlaveMatrix(),false,1.,1.);

        // derive linearizations of master-side scatra fluxes w.r.t. master-side structural dofs and assemble into scatra-structure matrix block
        FSI::UTILS::MatrixRowColTransform()(
            *strategy_scatra_->SlaveMatrix(),
            -1.,
            ADAPTER::CouplingSlaveConverter(*strategy_scatra_->CouplingAdapter()),
            ADAPTER::CouplingSlaveConverter(*icoup_structure_),
            *Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(scatrastructureblock_),
            true,
            true
            );

        break;
      }

      default:
      {
        dserror("Type of global system matrix for scalar-structure interaction not recognized!");
        break;
      }
    }
  }

  // finalize scatra-structure matrix block
  scatrastructureblock_->Complete(*maps_->Map(1),*maps_->Map(0));

  // apply Dirichlet boundary conditions to scatra-structure matrix block
  scatrastructureblock_->ApplyDirichlet(*scatra_->ScaTraField()->DirichMaps()->CondMap(),false);

  // remove state vectors from scalar transport discretization
  scatra_->ScaTraField()->Discretization()->ClearState();

  return;
}


/*-----------------------------------------------------------------------------------*
 | assemble off-diagonal structure-scatra block of global system matrix   fang 08/17 |
 *-----------------------------------------------------------------------------------*/
void SSI::SSI_Mono::AssembleODBlockStructureScatra() const
{
  // initialize structure-scatra matrix block
  structurescatrablock_->Zero();

  // create parameter list for element evaluation and fill it
  Teuchos::ParameterList eleparams;
  // set action
  eleparams.set("action","calc_struct_stiffscalar");
  // set time
  eleparams.set<double>("total time",Time());

  // remove state vectors from structure discretization
  structure_->Discretization()->ClearState();

  // and set the current displacement and acceleration state vector
  structure_->Discretization()->SetState("displacement",structure_->Dispnp());

  // create strategy for assembly of structure-scatra matrix block
  DRT::AssembleStrategy strategystructurescatra(
      0,                       // row assembly based on number of dofset associated with structure dofs on structural discretization
      1,                       // column assembly based on number of dofset associated with scalar transport dofs on structural discretization
      structurescatrablock_,   // structure-scatra matrix block
      Teuchos::null,           // no additional matrices or vectors needed
      Teuchos::null,
      Teuchos::null,
      Teuchos::null
      );

  // assemble structure-scatra matrix block
  structure_->Discretization()->Evaluate(eleparams,strategystructurescatra);

  // need to scale structurescatrablock_ with 'timefac' (e.g. with theta for OST-scheme) to get correct implementation
  const double timeintparam = structure_->TimIntParam();
  // scale with theta
  structurescatrablock_->Scale(1.0 - timeintparam);

  // complete the structure scatra off diagonal block
  structurescatrablock_->Complete(*maps_->Map(0),*maps_->Map(1));

  return;
}


/*--------------------------------------------------------------------------*
 | return global map of degrees of freedom                       fang 08/17 |
 *--------------------------------------------------------------------------*/
const Teuchos::RCP<const Epetra_Map>& SSI::SSI_Mono::DofRowMap() const
{
  return maps_->FullMap();
}


/*--------------------------------------------------------------------------*
 | finite difference check for global system matrix              fang 08/17 |
 *--------------------------------------------------------------------------*/
void SSI::SSI_Mono::FDCheck()
{
  // initial screen output
  if(Comm().MyPID() == 0)
    std::cout << std::endl << "FINITE DIFFERENCE CHECK FOR SSI SYSTEM MATRIX" << std::endl;

  // create global state vector
  Teuchos::RCP<Epetra_Vector> statenp(LINALG::CreateVector(*DofRowMap(),true));
  maps_->InsertVector(scatra_->ScaTraField()->Phinp(),0,statenp);
  maps_->InsertVector(structure_->Dispnp(),1,statenp);

  // make a copy of global state vector to undo perturbations later
  Teuchos::RCP<Epetra_Vector> statenp_original = Teuchos::rcp(new Epetra_Vector(*statenp));

  // make a copy of system matrix as Epetra_CrsMatrix
  Teuchos::RCP<Epetra_CrsMatrix> sysmat_original = Teuchos::null;
  if(Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(systemmatrix_) != Teuchos::null)
    sysmat_original = (new LINALG::SparseMatrix(*(Teuchos::rcp_static_cast<LINALG::SparseMatrix>(systemmatrix_))))->EpetraMatrix();
  else if(Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(systemmatrix_) != Teuchos::null)
    sysmat_original = (new LINALG::SparseMatrix(*(Teuchos::rcp_static_cast<LINALG::BlockSparseMatrixBase>(systemmatrix_)->Merge())))->EpetraMatrix();
  else
    dserror("Type of system matrix unknown!");
  sysmat_original->FillComplete();

  // make a copy of system right-hand side vector
  Teuchos::RCP<Epetra_Vector> rhs_original = Teuchos::rcp(new Epetra_Vector(*residual_));

  // copy and zero out system increment vector if necessary
  Teuchos::RCP<Epetra_Vector> increment_original(Teuchos::null);
  if(iter_ != 1)
  {
    increment_original = Teuchos::rcp(new Epetra_Vector(*increment_));
    increment_->PutScalar(0.);
  }

  // initialize counter for system matrix entries with failing finite difference check
  int counter(0);

  // initialize tracking variable for maximum absolute and relative errors
  double maxabserr(0.);
  double maxrelerr(0.);

  for (int colgid=0; colgid<=sysmat_original->ColMap().MaxAllGID(); ++colgid)
  {
    // continue loop if current column index is not a valid global column index
    int collid(sysmat_original->ColMap().LID(colgid));
    int maxcollid(-1);
    Comm().MaxAll(&collid,&maxcollid,1);
    if(maxcollid < 0)
      continue;

    // continue loop if current column index is associated with slave side of structural meshtying interface
    if(scatra_->ScaTraField()->S2ICoupling())
    {
      collid = icoup_structure_->SlaveDofMap()->LID(colgid);
      Comm().MaxAll(&collid,&maxcollid,1);
      if(maxcollid >= 0)
        continue;
    }

    // fill global state vector with original state variables
    statenp->Update(1.,*statenp_original,0.);

    // impose perturbation
    if(statenp->Map().MyGID(colgid))
      if(statenp->SumIntoGlobalValue(colgid,0,scatra_->ScaTraField()->FDCheckEps()))
        dserror("Perturbation could not be imposed on state vector for finite difference check!");
    if(scatra_->ScaTraField()->S2ICoupling() and icoup_structure_->PermMasterDofMap()->MyGID(colgid))
      if(statenp->SumIntoGlobalValue(icoup_structure_->SlaveDofMap()->GID(icoup_structure_->PermMasterDofMap()->LID(colgid)),0,scatra_->ScaTraField()->FDCheckEps()))
        dserror("Perturbation could not be imposed on state vector for finite difference check!");
    scatra_->ScaTraField()->Phinp()->Update(1.,*maps_->ExtractVector(statenp,0),0.);
    structure_->SetState(maps_->ExtractVector(statenp,1));

    // calculate element right-hand side vector for perturbed state
    AssembleMatAndRHS();

    // Now we compare the difference between the current entries in the system matrix
    // and their finite difference approximations according to
    // entries ?= (residual_perturbed - residual_original) / epsilon

    // Note that the residual_ vector actually denotes the right-hand side of the linear
    // system of equations, i.e., the negative system residual.
    // To account for errors due to numerical cancellation, we additionally consider
    // entries + residual_original / epsilon ?= residual_perturbed / epsilon

    // Note that we still need to evaluate the first comparison as well. For small entries in the system
    // matrix, the second comparison might yield good agreement in spite of the entries being wrong!
    for(int rowlid=0; rowlid<DofRowMap()->NumMyElements(); ++rowlid)
    {
      // get global index of current matrix row
      const int rowgid = sysmat_original->RowMap().GID(rowlid);
      if(rowgid < 0)
        dserror("Invalid global ID of matrix row!");

      // skip matrix rows associated with Dirichlet boundary conditions and slave side of structural meshtying interface
      if(scatra_->ScaTraField()->DirichMaps()->CondMap()->MyGID(rowgid) or structure_->GetDBCMapExtractor()->CondMap()->MyGID(rowgid) or (scatra_->ScaTraField()->S2ICoupling() and icoup_structure_->SlaveDofMap()->MyGID(rowgid)))
        continue;

      // get relevant entry in current row of original system matrix
      double entry(0.);
      int length = sysmat_original->NumMyEntries(rowlid);
      int numentries;
      std::vector<double> values(length);
      std::vector<int> indices(length);
      sysmat_original->ExtractMyRowCopy(rowlid,length,numentries,&values[0],&indices[0]);
      for(int ientry=0; ientry<length; ++ientry)
      {
        if(sysmat_original->ColMap().GID(indices[ientry]) == colgid)
        {
          entry = values[ientry];
          break;
        }
      }

      // finite difference suggestion (first divide by epsilon and then add for better conditioning)
      const double fdval = -(*residual_)[rowlid] / scatra_->ScaTraField()->FDCheckEps() + (*rhs_original)[rowlid] / scatra_->ScaTraField()->FDCheckEps();

      // confirm accuracy of first comparison
      if(abs(fdval) > 1.e-20 and abs(fdval) < 1.e-15)
      {
        // output warning
        std::cout << "WARNING: Finite difference check involves values very close to numerical zero!" << std::endl;

        // skip comparison if current entry is very small
        if(abs(entry) < 1.e-15)
          continue;
      }

      // absolute and relative errors in first comparison
      const double abserr1 = entry - fdval;
      if(abs(abserr1) > maxabserr)
        maxabserr = abs(abserr1);
      double relerr1(0.);
      if(abs(entry) > 1.e-17)
        relerr1 = abserr1 / abs(entry);
      else if(abs(fdval) > 1.e-17)
        relerr1 = abserr1 / abs(fdval);
      if(abs(relerr1) > maxrelerr)
        maxrelerr = abs(relerr1);

      // evaluate first comparison
      if(abs(relerr1) > scatra_->ScaTraField()->FDCheckTol())
      {
        std::cout << "sysmat[" << rowgid << "," << colgid << "]:  " << entry << "   ";
        std::cout << "finite difference suggestion:  " << fdval << "   ";
        std::cout << "absolute error:  " << abserr1 << "   ";
        std::cout << "relative error:  " << relerr1 << std::endl;

        counter++;
      }

      // first comparison OK
      else
      {
        // left-hand side in second comparison
        const double left  = entry - (*rhs_original)[rowlid] / scatra_->ScaTraField()->FDCheckEps();

        // right-hand side in second comparison
        const double right = -(*residual_)[rowlid] / scatra_->ScaTraField()->FDCheckEps();

        // confirm accuracy of second comparison
        if(abs(right) > 1.e-20 and abs(right) < 1.e-15)
        {
          // output warning
          std::cout << "WARNING: Finite difference check involves values very close to numerical zero!" << std::endl;

          // skip comparison if current left-hand side is very small
          if(abs(left) < 1.e-15)
            continue;
        }

        // absolute and relative errors in second comparison
        const double abserr2 = left - right;
        if(abs(abserr2) > maxabserr)
          maxabserr = abs(abserr2);
        double relerr2(0.);
        if(abs(left) > 1.e-17)
          relerr2 = abserr2 / abs(left);
        else if(abs(right) > 1.e-17)
          relerr2 = abserr2 / abs(right);
        if(abs(relerr2) > maxrelerr)
          maxrelerr = abs(relerr2);

        // evaluate second comparison
        if(abs(relerr2) > scatra_->ScaTraField()->FDCheckTol())
        {
          std::cout << "sysmat[" << rowgid << "," << colgid << "]-rhs[" << rowgid << "]/eps:  " << left << "   ";
          std::cout << "-rhs_perturbed[" << rowgid << "]/eps:  " << right << "   ";
          std::cout << "absolute error:  " << abserr2 << "   ";
          std::cout << "relative error:  " << relerr2 << std::endl;

          counter++;
        }
      }
    }
  }

  // communicate tracking variables
  int counterglobal(0);
  Comm().SumAll(&counter,&counterglobal,1);
  double maxabserrglobal(0.);
  Comm().MaxAll(&maxabserr,&maxabserrglobal,1);
  double maxrelerrglobal(0.);
  Comm().MaxAll(&maxrelerr,&maxrelerrglobal,1);

  // final screen output
  if(Comm().MyPID() == 0)
  {
    if(counterglobal)
    {
      printf("--> FAILED AS LISTED ABOVE WITH %d CRITICAL MATRIX ENTRIES IN TOTAL\n\n",counterglobal);
      printf("--> FAILED WITH MAXIMUM ABSOLUTE ERROR %+12.5e AND MAXIMUM RELATIVE ERROR %+12.5e\n\n",maxabserrglobal,maxrelerrglobal);
      dserror("Finite difference check failed for SSI system matrix!");
    }
    else
      printf("--> PASSED WITH MAXIMUM ABSOLUTE ERROR %+12.5e AND MAXIMUM RELATIVE ERROR %+12.5e\n\n",maxabserrglobal,maxrelerrglobal);
  }

  // undo perturbations of state variables
  scatra_->ScaTraField()->Phinp()->Update(1.,*maps_->ExtractVector(statenp_original,0),0.);
  structure_->SetState(maps_->ExtractVector(statenp_original,1));

  // recompute system matrix and right-hand side vector based on original state variables
  AssembleMatAndRHS();

  // restore system increment vector if necessary
  if(increment_original != Teuchos::null)
    increment_ = increment_original;

  return;
}


/*--------------------------------------------------------------------------*
 | initialize monolithic algorithm                               fang 08/17 |
 *--------------------------------------------------------------------------*/
int SSI::SSI_Mono::Init(
    const Epetra_Comm&              comm,               //!< communicator
    const Teuchos::ParameterList&   globaltimeparams,   //!< parameter list for time integration
    const Teuchos::ParameterList&   scatraparams,       //!< parameter list for scalar transport
    const Teuchos::ParameterList&   structparams,       //!< parameter list for structure
    const std::string               struct_disname,     //!< name of structural discretization
    const std::string               scatra_disname,     //!< name of scalar transport discretization
    bool                            isAle               //!< flag for ALE
    )
{
  // check input parameters for scalar transport field
  if(DRT::INPUT::IntegralValue<INPAR::SCATRA::VelocityField>(scatraparams,"VELOCITYFIELD") != INPAR::SCATRA::velocity_Navier_Stokes)
    dserror("Invalid type of velocity field for scalar-structure interaction!");

  // initialize strategy for Newton-Raphson convergence check
  switch(DRT::INPUT::IntegralValue<INPAR::SSI::ScaTraTimIntType>(globaltimeparams,"SCATRATIMINTTYPE"))
  {
    case INPAR::SSI::scatratiminttype_elch:
    {
      strategy_convcheck_ = Teuchos::rcp(new SSI::SSI_Mono::ConvCheckStrategyElch(globaltimeparams));
      break;
    }

    case INPAR::SSI::scatratiminttype_standard:
    {
      strategy_convcheck_ = Teuchos::rcp(new SSI::SSI_Mono::ConvCheckStrategyStd(globaltimeparams));
      break;
    }

    default:
    {
      dserror("Type of scalar transport time integrator currently not supported!");
      break;
    }
  }

  // call base class routine
  return SSI_Base::Init(comm,globaltimeparams,scatraparams,structparams,struct_disname,scatra_disname,isAle);
}


/*--------------------------------------------------------------------------*
 | output solution to screen and files                           fang 08/17 |
 *--------------------------------------------------------------------------*/
void SSI::SSI_Mono::Output()
{
  // output scalar transport field
  scatra_->ScaTraField()->Output();

  // output structure field
  structure_->Output();

  return;
}


/*--------------------------------------------------------------------------*
 | prepare time step                                             fang 08/17 |
 *--------------------------------------------------------------------------*/
void SSI::SSI_Mono::PrepareTimeStep()
{
  // update time and time step
  IncrementTimeAndStep();

  // pass structural degrees of freedom to scalar transport discretization
  SetStructSolution(structure_->Dispnp(),structure_->Velnp());

  // prepare time step for scalar transport field
  scatra_->ScaTraField()->PrepareTimeStep();

  // pass scalar transport degrees of freedom to structural discretization
  // has to be called AFTER scatra_->ScaTraField()->PrepareTimeStep(), because in the Timestep 0
  // the potential field is calculated consistently from defined initial and boundary conditions within this call
  SetScatraSolution(scatra_->ScaTraField()->Phinp());

  // prepare time step for structural field
  structure_->PrepareTimeStep();

  // print time step information to screen
  scatra_->ScaTraField()->PrintTimeStepInfo();

  return;
}


/*--------------------------------------------------------------------------*
 | setup monolithic algorithm                                    fang 08/17 |
 *--------------------------------------------------------------------------*/
void SSI::SSI_Mono::Setup()
{
  // check initialization
  CheckIsInit();

  // set up scalar transport field
  scatra_->ScaTraField()->Setup();

  // pass initial scalar field to structural discretization to correctly compute initial accelerations
  DRT::Problem::Instance()->GetDis("structure")->SetState(1,"temperature",scatra_->ScaTraField()->Phinp());

  // set up structural base algorithm
  struct_adapterbase_ptr_->Setup();

  // extract and cast wrapper
  structure_ = Teuchos::rcp_dynamic_cast<::ADAPTER::SSIStructureWrapper>(struct_adapterbase_ptr_->StructureField());
  if(structure_ == Teuchos::null)
    dserror("Invalid structural wrapper!");

  // set up helper class for field coupling
  ssicoupling_->Setup();

  // set up materials
  ssicoupling_->AssignMaterialPointers(struct_adapterbase_ptr_->StructureField()->Discretization(),scatra_->ScaTraField()->Discretization());

  // check maps from scalar transport and structure discretizations
  if(scatra_->ScaTraField()->DofRowMap()->NumGlobalElements() == 0)
    dserror("Scalar transport discretization does not have any degrees of freedom!");
  if(structure_->DofRowMap()->NumGlobalElements() == 0)
    dserror("Structure discretization does not have any degrees of freedom!");

  // additional safety checks
  if(!scatra_->ScaTraField()->IsIncremental())
    dserror("Must have incremental solution approach for scalar-structure interaction!");
  if(scatra_->ScaTraField()->SystemMatrix() == Teuchos::null)
    dserror("System matrix associated with scalar transport field must be a sparse matrix!");

  // set up scatra-scatra interface coupling
  if(scatra_->ScaTraField()->S2ICoupling())
  {
    // extract meshtying strategy for scatra-scatra interface coupling on scatra discretization
    strategy_scatra_ = Teuchos::rcp_dynamic_cast<const SCATRA::MeshtyingStrategyS2I>(scatra_->ScaTraField()->Strategy());

    // safety checks
    if(strategy_scatra_ == Teuchos::null)
      dserror("Invalid scatra-scatra interface coupling strategy!");
    if(strategy_scatra_->CouplingType() != INPAR::S2I::coupling_matching_nodes)
      dserror("Monolithic scalar-structure interaction only implemented for scatra-scatra interface coupling with matching interface nodes!");

    // set up scatra-scatra interface coupling adapter for structure field
    SetupCouplingAdapterStructure();

    // set up map for interior and master-side structural degrees of freedom
    map_structure_ = LINALG::SplitMap(*structure_->Discretization()->DofRowMap(),*icoup_structure_->SlaveDofMap());

    // set up structural map extractor holding interior and interface maps of degrees of freedom
    std::vector<Teuchos::RCP<const Epetra_Map> > maps(0,Teuchos::null);
    maps.push_back(LINALG::SplitMap(*map_structure_,*icoup_structure_->MasterDofMap()));
    maps.push_back(icoup_structure_->SlaveDofMap());
    maps.push_back(icoup_structure_->MasterDofMap());
    maps_structure_ = Teuchos::rcp(new LINALG::MultiMapExtractor(*structure_->Discretization()->DofRowMap(),maps));
    maps_structure_->CheckForValidMapExtractor();
  }

  // construct vector of zeroes
  zeros_ = LINALG::CreateVector(*structure_->DofRowMap());

  // set flag
  SetIsSetup(true);

  return;
}


/*--------------------------------------------------------------------------*
 | setup global system of equations                              fang 08/17 |
 *--------------------------------------------------------------------------*/
void SSI::SSI_Mono::SetupSystem()
{
  if(scatra_->ScaTraField()->S2ICoupling())
  {
    // check whether slave-side degrees of freedom are Dirichlet-free
    std::vector<Teuchos::RCP<const Epetra_Map> > maps(2,Teuchos::null);
    maps[0] = icoup_structure_->SlaveDofMap();
    maps[1] = structure_->GetDBCMapExtractor()->CondMap();
    if(LINALG::MultiMapExtractor::IntersectMaps(maps)->NumGlobalElements() > 0)
      dserror("Must not apply Dirichlet conditions to slave-side structural displacements!");
  }

  // initialize global map extractor
  maps_ = Teuchos::rcp(new LINALG::MapExtractor(
      *LINALG::MergeMap(
          *scatra_->ScaTraField()->DofRowMap(),
          *structure_->DofRowMap(),
          false
          ),
      structure_->DofRowMap(),
      scatra_->ScaTraField()->DofRowMap()
      ));

  // check global map extractor
  maps_->CheckForValidMapExtractor();

  // initialize global increment vector for Newton-Raphson iteration
  increment_ = LINALG::CreateVector(
      *DofRowMap(),
      true
      );

  // initialize global residual vector
  residual_ = LINALG::CreateVector(
      *DofRowMap(),
      true
      );

  // perform initializations associated with global system matrix
  switch(matrixtype_)
  {
    case INPAR::SSI::matrix_sparse:
    {
      // safety check
      if(scatra_->ScaTraField()->SystemMatrix() == Teuchos::null)
        dserror("Incompatible matrix type associated with scalar transport field!");

      // initialize global system matrix
      systemmatrix_ = Teuchos::rcp(new LINALG::SparseMatrix(*DofRowMap(),27,false,true));

      break;
    }

    default:
    {
      dserror("Type of global system matrix for scalar-structure interaction not recognized!");
      break;
    }
  }

  // initialize scatra-structure block
  scatrastructureblock_ = Teuchos::rcp(new LINALG::SparseMatrix(
      *scatra_->ScaTraField()->DofRowMap(),
      27,
      false,
      true
      ));

  // initialize structure-scatra block
  structurescatrablock_ = Teuchos::rcp(new LINALG::SparseMatrix(
      *structure_->DofRowMap(),
      27,
      false,
      true
      ));

  return;
}


/*---------------------------------------------------------------------------------*
 | setup scatra-scatra interface coupling adapter for structure field   fang 08/17 |
 *---------------------------------------------------------------------------------*/
void SSI::SSI_Mono::SetupCouplingAdapterStructure()
{
  // initialize integer vectors for global IDs of master-side and slave-side interface nodes on structure discretization
  std::vector<int> inodegidvec_master;
  std::vector<int> inodegidvec_slave;

  // extract scatra-scatra interface coupling conditions from structure discretization
  std::vector<DRT::Condition*> conditions(0,NULL);
  structure_->Discretization()->GetCondition("S2ICoupling",conditions);

  // loop over all conditions
  for(unsigned icondition=0; icondition<conditions.size(); ++icondition)
  {
    // extract current condition
    DRT::Condition* const condition = conditions[icondition];

    // extract interface side associated with current condition
    const int side = condition->GetInt("interface side");

    // extract nodes associated with current condition
    const std::vector<int>* const inodegids = condition->Nodes();

    // loop over all nodes
    for(unsigned inode=0; inode<inodegids->size(); ++inode)
    {
      // extract ID of current node
      const int inodegid = (*inodegids)[inode];

      // insert global ID of current node into associated vector only if node is owned by current processor
      // need to make sure that node is stored on current processor, otherwise cannot resolve "->Owner()"
      if(structure_->Discretization()->HaveGlobalNode(inodegid) and structure_->Discretization()->gNode(inodegid)->Owner() == structure_->Discretization()->Comm().MyPID())
        side == INPAR::S2I::side_master ? inodegidvec_master.push_back(inodegid) : inodegidvec_slave.push_back(inodegid);
    }
  }

  // remove potential duplicates from vectors
  std::sort(inodegidvec_master.begin(),inodegidvec_master.end());
  inodegidvec_master.erase(unique(inodegidvec_master.begin(),inodegidvec_master.end()),inodegidvec_master.end());
  std::sort(inodegidvec_slave.begin(),inodegidvec_slave.end());
  inodegidvec_slave.erase(unique(inodegidvec_slave.begin(),inodegidvec_slave.end()),inodegidvec_slave.end());

  // setup scatra-scatra interface coupling adapter for structure field
  icoup_structure_ = Teuchos::rcp(new ADAPTER::Coupling());
  icoup_structure_->SetupCoupling(
      *structure_->Discretization(),
      *structure_->Discretization(),
      inodegidvec_master,
      inodegidvec_slave,
      DRT::Problem::Instance()->NDim(),
      true,
      1.e-8
      );

  return;
}


/*--------------------------------------------------------------------------*
 | evaluate time step using Newton-Raphson iteration             fang 08/17 |
 *--------------------------------------------------------------------------*/
void SSI::SSI_Mono::Solve()
{
  // initialize counter for Newton-Raphson iteration
  iter_ = 0;

  // start Newton-Raphson iteration
  while(true)
  {
    // update iteration counter
    ++iter_;

    // reset timer
    timer_->ResetStartTime();

    // store time before evaluating elements and assembling global system of equations
    double time = timer_->WallTime();

    // assemble global system of equations
    AssembleMatAndRHS();

    // determine time needed for evaluating elements and assembling global system of equations,
    // and take maximum over all processors via communication
    double mydtele = timer_->WallTime()-time;
    Comm().MaxAll(&mydtele,&dtele_,1);

    // safety check
    if(!systemmatrix_->Filled())
      dserror("Complete() has not been called on global system matrix yet!");

    // perform finite difference check on time integrator level
    if((scatra_->ScaTraField()->FDCheckType() == INPAR::SCATRA::fdcheck_global) and (Step() > 1))
      FDCheck();

    // check termination criterion for Newton-Raphson iteration
    if(strategy_convcheck_->ExitNewtonRaphson(*this))
      break;

    // initialize global increment vector
    increment_->PutScalar(0.);

    // store time before solving global system of equations
    time = timer_->WallTime();

    // solve global system of equations
    // Dirichlet boundary conditions have already been applied to global system of equations
    solver_->Solve(
        systemmatrix_->EpetraOperator(),
        increment_,
        residual_,
        true,
        iter_==1
        );

    // determine time needed for solving global system of equations,
    // and take maximum over all processors via communication
    double mydtsolve = timer_->WallTime()-time;
    Comm().MaxAll(&mydtsolve,&dtsolve_,1);

    // update scalar transport field
    scatra_->ScaTraField()->UpdateIter(maps_->ExtractVector(increment_,0));
    scatra_->ScaTraField()->ComputeIntermediateValues();

    // structure field is updated during the next Newton-Raphson iteration step
  } // Newton-Raphson iteration

  return;
}


/*--------------------------------------------------------------------------*
 | test results                                                  fang 11/17 |
 *--------------------------------------------------------------------------*/
void SSI::SSI_Mono::TestResults(const Epetra_Comm& comm) const
{
  // activate result testing for monolithic scalar-structure interaction
  DRT::Problem::Instance()->AddFieldTest(Teuchos::rcp(new SSI::SSI_Mono_ResultTest(Teuchos::rcp(this,false))));

  // call base class routine
  SSI_Base::TestResults(comm);

  return;
}


/*--------------------------------------------------------------------------*
 | time loop                                                     fang 08/17 |
 *--------------------------------------------------------------------------*/
void SSI::SSI_Mono::Timeloop()
{
  // output initial solution to screen and files
  if(Step() == 0)
  {
    SetStructSolution(structure_->Dispnp(),structure_->Velnp());
    Output();
  }

  // time loop
  while(NotFinished())
  {
    // prepare time step
    PrepareTimeStep();

    // evaluate time step
    Solve();

    // update scalar transport and structure fields
    Update();

    // output solution to screen and files
    Output();
  } // while(NotFinished())

  return;
}


/*--------------------------------------------------------------------------------------*
 | update scalar transport and structure fields after time step evaluation   fang 08/17 |
 *--------------------------------------------------------------------------------------*/
void SSI::SSI_Mono::Update()
{
  // update scalar transport field
  scatra_->ScaTraField()->Update();

  // update structure field
  structure_->Update();

  return;
}

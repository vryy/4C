/*----------------------------------------------------------------------*/
/*!
\file sti_algorithm.cpp

\brief monolithic algorithm for scatra-thermo interaction

\level 2

<pre>
\maintainer Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>
*/
/*----------------------------------------------------------------------*/
#include "sti_algorithm.H"

#include <Epetra_Time.h>

#include "../drt_adapter/adapter_coupling.H"
#include "../drt_adapter/adapter_scatra_base_algorithm.H"

#include "../drt_fsi/fsi_matrixtransform.H"

#include "../drt_io/io_control.H"

#include "../drt_lib/drt_assemblestrategy.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../drt_scatra/scatra_timint_implicit.H"
#include "../drt_scatra/scatra_timint_meshtying_strategy_s2i.H"

#include "../drt_scatra_ele/scatra_ele_action.H"

#include "../linalg/linalg_multiply.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_utils.H"

/*--------------------------------------------------------------------------------*
 | constructor                                                         fang 04/15 |
 *--------------------------------------------------------------------------------*/
STI::Algorithm::Algorithm(
    const Epetra_Comm&              comm,                  //! communicator
    const Teuchos::ParameterList&   stidyn,                //! parameter list for scatra-thermo interaction
    const Teuchos::ParameterList&   scatradyn,             //! scalar transport parameter list for scatra and thermo fields
    const Teuchos::ParameterList&   solverparams,          //! solver parameter list for scatra-thermo interaction
    const Teuchos::ParameterList&   solverparams_scatra,   //! solver parameter list for scatra field
    const Teuchos::ParameterList&   solverparams_thermo    //! solver parameter list for thermo field
    ) :
    // instantiate base class
    AlgorithmBase(comm,scatradyn),

    scatra_(Teuchos::null),
    thermo_(Teuchos::null),
    strategyscatra_(Teuchos::null),
    strategythermo_(Teuchos::null),
    stiparameters_(Teuchos::rcp(new Teuchos::ParameterList(stidyn))),
    fieldparameters_(Teuchos::rcp(new Teuchos::ParameterList(scatradyn))),
    iter_(0),
    itermax_(fieldparameters_->sublist("NONLINEAR").get<int>("ITEMAX")),
    itertol_(fieldparameters_->sublist("NONLINEAR").get<double>("CONVTOL")),
    restol_(fieldparameters_->sublist("NONLINEAR").get<double>("ABSTOLRES")),
    maps_(Teuchos::null),
    condensationthermo_(DRT::INPUT::IntegralValue<bool>(stidyn,"THERMO_CONDENSATION")),
    systemmatrix_(Teuchos::null),
    matrixtype_(DRT::INPUT::IntegralValue<INPAR::STI::MatrixType>(stidyn,"MATRIXTYPE")),
    scatrathermoblock_(Teuchos::null),
    thermoscatrablock_(Teuchos::null),
    blockmaps_(Teuchos::null),
    blockmapthermo_(Teuchos::null),
    increment_(Teuchos::null),

    // initialize timer for Newton-Raphson iteration
    timer_(Teuchos::rcp(new Epetra_Time(comm))),

    residual_(Teuchos::null),
    dtele_(0.),
    dtsolve_(0.),

    // initialize algebraic solver for global system of equations
    solver_(Teuchos::rcp(new LINALG::Solver(
        solverparams,
        comm,
        DRT::Problem::Instance()->ErrorFile()->Handle()
        ))),

    equilibration_(DRT::INPUT::IntegralValue<INPAR::S2I::EquilibrationMethods>(fieldparameters_->sublist("S2I COUPLING"),"EQUILIBRATION")),
    invrowsums_(Teuchos::null),

    // initialize L2 norms for Newton-Raphson convergence check
    scatradofnorm_(0.),
    scatraresnorm_(0.),
    scatraincnorm_(0.),
    thermodofnorm_(0.),
    thermoresnorm_(0.),
    thermoincnorm_(0.),

    icoupscatra_(Teuchos::null),
    icoupthermo_(Teuchos::null),
    islavetomasterrowtransformscatraod_(Teuchos::null),
    islavetomastercoltransformthermood_(Teuchos::null),
    islavetomasterrowtransformthermood_(Teuchos::null)
{
  // check input parameters for scatra and thermo fields
  if(DRT::INPUT::IntegralValue<INPAR::SCATRA::VelocityField>(*fieldparameters_,"VELOCITYFIELD") != INPAR::SCATRA::velocity_zero)
    dserror("Scatra-thermo interaction with convection not yet implemented!");

  // initialize scatra time integrator
  Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm > basealgo =
      Teuchos::rcp(new ADAPTER::ScaTraBaseAlgorithm());
  basealgo->Init(*fieldparameters_,*fieldparameters_,solverparams_scatra);
  basealgo->Setup();

  // get pointer to time integrator
  scatra_ = basealgo->ScaTraField();

  // modify field parameters for thermo field
  ModifyFieldParametersForThermoField();

  // initialize thermo time integrator
  Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm > thermo_basealgo =
      Teuchos::rcp(new ADAPTER::ScaTraBaseAlgorithm());
  thermo_basealgo->Init(*fieldparameters_,*fieldparameters_,solverparams_thermo,"thermo");
  thermo_basealgo->Setup();

  // get pointer to time integrator
  thermo_ = thermo_basealgo->ScaTraField();

  // check maps from scatra and thermo discretizations
  if(scatra_->Discretization()->DofRowMap()->NumGlobalElements() == 0)
    dserror("Scatra discretization does not have any degrees of freedom!");
  if(thermo_->Discretization()->DofRowMap()->NumGlobalElements() == 0)
    dserror("Thermo discretization does not have any degrees of freedom!");

  // additional safety checks
  if(!scatra_->IsIncremental())
    dserror("Must have incremental solution approach for scatra-thermo interaction!");
  if(thermo_->NumScal() != 1)
    dserror("Thermo field must involve exactly one transported scalar!");
  if(thermo_->SystemMatrix() == Teuchos::null)
    dserror("System matrix associated with temperature field must be a sparse matrix!");

  // perform initializations associated with scatra-scatra interface coupling
  if(scatra_->S2ICoupling())
  {
    // safety check
    if(!thermo_->S2ICoupling())
      dserror("Can't evaluate scatra-scatra interface coupling in scatra field, but not in thermo field!");

    // extract meshtying strategies for scatra-scatra interface coupling from scatra and thermo time integrators
    strategyscatra_ = Teuchos::rcp_dynamic_cast<SCATRA::MeshtyingStrategyS2I>(scatra_->Strategy());
    strategythermo_ = Teuchos::rcp_dynamic_cast<SCATRA::MeshtyingStrategyS2I>(thermo_->Strategy());

    // perform initializations depending on type of meshtying method
    switch(strategyscatra_->CouplingType())
    {
      case INPAR::S2I::coupling_matching_nodes:
      {
        // safety check
        if(strategythermo_->CouplingType() != INPAR::S2I::coupling_matching_nodes)
          dserror("Must have matching nodes at scatra-scatra coupling interfaces in both the scatra and the thermo fields!");

        // extract coupling adapters
        icoupscatra_ = strategyscatra_->CouplingAdapter();
        icoupthermo_ = strategythermo_->CouplingAdapter();

        break;
      }

      case INPAR::S2I::coupling_mortar_standard:
      {
        // safety check
        if(strategythermo_->CouplingType() != INPAR::S2I::coupling_mortar_condensed_bubnov)
          dserror("Invalid type of scatra-scatra interface coupling for thermo field!");

        // extract scatra-scatra interface coupling conditions
        std::vector<DRT::Condition*> conditions;
        scatra_->Discretization()->GetCondition("S2ICoupling",conditions);

        // loop over all conditions
        for(unsigned icondition=0; icondition<conditions.size(); ++icondition)
          // consider conditions for slave side only
          if(conditions[icondition]->GetInt("interface side") == INPAR::S2I::side_slave)
          {
            // extract ID of current condition
            const int condid = conditions[icondition]->GetInt("ConditionID");
            if(condid < 0)
              dserror("Invalid condition ID!");

            // extract mortar discretizations associated with current condition
            DRT::Discretization& scatradis = strategyscatra_->MortarDiscretization(condid);
            DRT::Discretization& thermodis = strategythermo_->MortarDiscretization(condid);

            // exchange dofsets between discretizations
            scatradis.AddDofSet(thermodis.GetDofSetProxy());
            thermodis.AddDofSet(scatradis.GetDofSetProxy());
          }

        break;
      }

      default:
      {
        dserror("Invalid type of scatra-scatra interface coupling!");
        break;
      }
    }
  }

  // initialize map associated with single thermo block of global system matrix
  Teuchos::RCP<const Epetra_Map> mapthermo(Teuchos::null);
  if(condensationthermo_)
    mapthermo = LINALG::MergeMap(*strategythermo_->InterfaceMaps()->Map(0),*strategythermo_->InterfaceMaps()->Map(2));
  else
    mapthermo = thermo_->DofRowMap();

  // initialize global map extractor
  maps_ = Teuchos::rcp(new LINALG::MapExtractor(
      *LINALG::MergeMap(
          *scatra_->Discretization()->DofRowMap(),
          *mapthermo,
          false
          ),
      mapthermo,
      scatra_->DofRowMap()
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

  // initialize transformation operators
  islavetomasterrowtransformscatraod_ = Teuchos::rcp(new FSI::UTILS::MatrixRowTransform);
  islavetomastercoltransformthermood_ = Teuchos::rcp(new FSI::UTILS::MatrixColTransform);
  islavetomasterrowtransformthermood_ = Teuchos::rcp(new FSI::UTILS::MatrixRowTransform);

  // initialize map extractors associated with blocks of global system matrix
  switch(strategyscatra_->MatrixType())
  {
    // one single main-diagonal matrix block associated with scalar transport field
    case INPAR::S2I::matrix_sparse:
    {
      blockmaps_ = maps_;
      break;
    }

    // several main-diagonal matrix blocks associated with scalar transport field
    case INPAR::S2I::matrix_block_condition:
    {
      // extract maps underlying main-diagonal matrix blocks associated with scalar transport field
      const unsigned nblockmapsscatra = strategyscatra_->BlockMaps().NumMaps();
      std::vector<Teuchos::RCP<const Epetra_Map> > blockmaps(nblockmapsscatra+1);
      for(unsigned iblockmap=0; iblockmap<nblockmapsscatra; ++iblockmap)
        blockmaps[iblockmap] = strategyscatra_->BlockMaps().Map(iblockmap);

      // extract map underlying single main-diagonal matrix block associated with temperature field
      blockmaps[nblockmapsscatra] = mapthermo;

      // initialize map extractor associated with blocks of global system matrix
      blockmaps_ = Teuchos::rcp(new LINALG::MultiMapExtractor(
          *DofRowMap(),
          blockmaps
          ));

      // initialize map extractor associated with all degrees of freedom inside temperature field
      blockmapthermo_ = Teuchos::rcp(new LINALG::MultiMapExtractor(
          *thermo_->Discretization()->DofRowMap(),
          std::vector<Teuchos::RCP<const Epetra_Map> >(1,thermo_->DofRowMap())
          ));

      // safety check
      blockmapthermo_->CheckForValidMapExtractor();

      break;
    }

    default:
    {
      dserror("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }

  // safety check
  blockmaps_->CheckForValidMapExtractor();

  // perform initializations associated with global system matrix
  switch(matrixtype_)
  {
    case INPAR::STI::matrix_block:
    {
      // safety check
      if(!solver_->Params().isSublist("AMGnxn Parameters"))
        dserror("Global system matrix with block structure requires AMGnxn block preconditioner!");

      // initialize global system matrix
      systemmatrix_ = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
          *blockmaps_,
          *blockmaps_,
          81,
          false,
          true
          ));

      // feed AMGnxn block preconditioner with null space information for each block of global block system matrix
      BuildBlockNullSpaces();

      break;
    }

    case INPAR::STI::matrix_sparse:
    {
      // safety check
      if(strategyscatra_->MatrixType() != INPAR::S2I::matrix_sparse)
        dserror("Incompatible matrix type associated with scalar transport field!");

      // initialize global system matrix
      systemmatrix_ = Teuchos::rcp(new LINALG::SparseMatrix(*DofRowMap(),27,false,true));

      // feed AMG preconditioner with null space information associated with global system matrix if applicable
      ComputeNullSpaceIfNecessary(solver_->Params());

      break;
    }

    default:
    {
      dserror("Type of global system matrix for scatra-thermo interaction not recognized!");
      break;
    }
  }

  // initialize scatra-thermo block and thermo-scatra block of global system matrix
  switch(strategyscatra_->MatrixType())
  {
    case INPAR::S2I::matrix_block_condition:
    {
      // initialize scatra-thermo block
      scatrathermoblock_ = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
          *blockmapthermo_,
          strategyscatra_->BlockMaps(),
          81,
          false,
          true
          ));

      // initialize thermo-scatra block
      thermoscatrablock_ = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
          strategyscatra_->BlockMaps(),
          *blockmapthermo_,
          81,
          false,
          true
          ));

      break;
    }

    case INPAR::S2I::matrix_sparse:
    {
      // initialize scatra-thermo block
      scatrathermoblock_ = Teuchos::rcp(new LINALG::SparseMatrix(
          *scatra_->Discretization()->DofRowMap(),
          27,
          false,
          true
          ));

      // initialize thermo-scatra block
      thermoscatrablock_ = Teuchos::rcp(new LINALG::SparseMatrix(
          *thermo_->Discretization()->DofRowMap(),
          27,
          false,
          true
          ));

      break;
    }

    default:
    {
      dserror("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }

  // perform initialization associated with equilibration of global system of equations
  switch(equilibration_)
  {
    case INPAR::S2I::equilibration_none:
    {
      // do nothing
      break;
    }

    case INPAR::S2I::equilibration_rows_full:
    case INPAR::S2I::equilibration_rows_maindiag:
    {
      // initialize vector for row sums of global system matrix if necessary
      invrowsums_ = Teuchos::rcp(new Epetra_Vector(*DofRowMap(),false));
      break;
    }

    default:
    {
      dserror("Equilibration method not yet implemented!");
      break;
    }
  }

  return;
} // STI::Algorithm::Algorithm


/*----------------------------------------------------------------------*
 | assemble global system of equations                       fang 07/15 |
 *----------------------------------------------------------------------*/
void STI::Algorithm::AssembleMatAndRHS()
{
  // pass scatra degrees of freedom to thermo discretization and vice versa
  ExchangeStateVectors();

  // build system matrix and residual for scatra field
  scatra_->PrepareLinearSolve();

  // pass master-side scatra degrees of freedom to thermo discretization for evaluation of scatra-scatra interface coupling
  if(thermo_->S2ICoupling() and strategythermo_->CouplingType() == INPAR::S2I::coupling_matching_nodes)
    thermo_->Discretization()->SetState(2,"imasterscatra",strategyscatra_->MasterPhinp());

  // build system matrix and residual for thermo field
  thermo_->PrepareLinearSolve();

  // assemble off-diagonal scatra-thermo block of global system matrix (derivatives of scatra residuals w.r.t. thermo degrees of freedom)
  AssembleODBlockScatraThermo();

  // assemble off-diagonal thermo-scatra block of global system matrix (derivatives of thermo residuals w.r.t. scatra degrees of freedom)
  AssembleODBlockThermoScatra();

  // build global system matrix
  switch(matrixtype_)
  {
    case INPAR::STI::matrix_block:
    {
      // check global system matrix
      Teuchos::RCP<LINALG::BlockSparseMatrixBase> blocksystemmatrix = Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(systemmatrix_);
      if(blocksystemmatrix == Teuchos::null)
        dserror("System matrix is not a block matrix!");

      switch(strategyscatra_->MatrixType())
      {
        case INPAR::S2I::matrix_block_condition:
        {
          // extract number of matrix row or column blocks associated with scalar transport field
          const unsigned nblockmapsscatra = strategyscatra_->BlockMaps().NumMaps();

          // construct global system matrix by assigning matrix blocks
          for(unsigned iblock=0; iblock<nblockmapsscatra; ++iblock)
          {
            for(unsigned jblock=0; jblock<nblockmapsscatra; ++jblock)
              blocksystemmatrix->Assign(
                  iblock,
                  jblock,
                  LINALG::View,
                  scatra_->BlockSystemMatrix()->Matrix(iblock,jblock)
                  );

            // perform second condensation before assigning matrix blocks
            if(condensationthermo_)
            {
              const LINALG::SparseMatrix& scatrathermoblock = Teuchos::rcp_dynamic_cast<const LINALG::BlockSparseMatrixBase>(scatrathermoblock_)->Matrix(iblock,0);
              FSI::UTILS::MatrixLogicalSplitAndTransform()(
                  scatrathermoblock,
                  scatrathermoblock.RangeMap(),
                  *maps_->Map(1),
                  1.,
                  NULL,
                  NULL,
                  blocksystemmatrix->Matrix(iblock,nblockmapsscatra)
                  );

              switch(strategyscatra_->CouplingType())
              {
                case INPAR::S2I::coupling_matching_nodes:
                {
                  ADAPTER::CouplingSlaveConverter converter(*icoupthermo_);
                  FSI::UTILS::MatrixLogicalSplitAndTransform()(
                      scatrathermoblock,
                      scatrathermoblock.RangeMap(),
                      *strategythermo_->InterfaceMaps()->Map(1),
                      1.,
                      NULL,
                      &converter,
                      blocksystemmatrix->Matrix(iblock,nblockmapsscatra),
                      true,
                      true
                      );
                  break;
                }
                case INPAR::S2I::coupling_mortar_standard:
                {
                  // initialize temporary matrix for slave-side columns of current matrix block
                  LINALG::SparseMatrix scatrathermocolsslave(scatrathermoblock.RowMap(),81);

                  // fill temporary matrix for slave-side columns of current matrix block
                  FSI::UTILS::MatrixLogicalSplitAndTransform()(
                      scatrathermoblock,
                      scatrathermoblock.RangeMap(),
                      *strategythermo_->InterfaceMaps()->Map(1),
                      1.,
                      NULL,
                      NULL,
                      scatrathermocolsslave
                      );

                  // finalize temporary matrix for slave-side columns of current matrix block
                  scatrathermocolsslave.Complete(*strategythermo_->InterfaceMaps()->Map(1),scatrathermoblock.RowMap());

                  // transform and assemble temporary matrix for slave-side columns of current matrix block
                  blocksystemmatrix->Matrix(iblock,nblockmapsscatra).Add(*LINALG::MLMultiply(scatrathermocolsslave,*strategythermo_->P(),true),false,1.,1.);

                  break;
                }
                default:
                {
                  dserror("Invalid type of scatra-scatra interface coupling!");
                  break;
                }
              }

              const LINALG::SparseMatrix& thermoscatrablock = Teuchos::rcp_dynamic_cast<const LINALG::BlockSparseMatrixBase>(thermoscatrablock_)->Matrix(0,iblock);
              FSI::UTILS::MatrixLogicalSplitAndTransform()(
                  thermoscatrablock,
                  *maps_->Map(1),
                  thermoscatrablock.DomainMap(),
                  1.,
                  NULL,
                  NULL,
                  blocksystemmatrix->Matrix(nblockmapsscatra,iblock)
                  );
            }

            // assign matrix blocks directly
            else
            {
              blocksystemmatrix->Assign(
                  iblock,
                  nblockmapsscatra,
                  LINALG::View,
                  Teuchos::rcp_dynamic_cast<const LINALG::BlockSparseMatrixBase>(scatrathermoblock_)->Matrix(iblock,0)
                  );

              blocksystemmatrix->Assign(
                  nblockmapsscatra,
                  iblock,
                  LINALG::View,
                  Teuchos::rcp_dynamic_cast<const LINALG::BlockSparseMatrixBase>(thermoscatrablock_)->Matrix(0,iblock)
                  );
            }
          }

          // perform second condensation before assigning thermo-thermo matrix block
          if(condensationthermo_)
          {
            FSI::UTILS::MatrixLogicalSplitAndTransform()(
                *thermo_->SystemMatrix(),
                *maps_->Map(1),
                *maps_->Map(1),
                1.,
                NULL,
                NULL,
                blocksystemmatrix->Matrix(nblockmapsscatra,nblockmapsscatra)
                );

            switch(strategyscatra_->CouplingType())
            {
              case INPAR::S2I::coupling_matching_nodes:
              {
                ADAPTER::CouplingSlaveConverter converter(*icoupthermo_);
                FSI::UTILS::MatrixLogicalSplitAndTransform()(
                    *thermo_->SystemMatrix(),
                    *maps_->Map(1),
                    *strategythermo_->InterfaceMaps()->Map(1),
                    1.,
                    NULL,
                    &converter,
                    blocksystemmatrix->Matrix(nblockmapsscatra,nblockmapsscatra),
                    true,
                    true
                    );
                break;
              }
              case INPAR::S2I::coupling_mortar_standard:
              {
                // initialize temporary matrix for slave-side columns of thermo-thermo matrix block
                LINALG::SparseMatrix thermothermocolsslave(*maps_->Map(1),81);

                // fill temporary matrix for slave-side columns of thermo-thermo matrix block
                FSI::UTILS::MatrixLogicalSplitAndTransform()(
                    *thermo_->SystemMatrix(),
                    *maps_->Map(1),
                    *strategythermo_->InterfaceMaps()->Map(1),
                    1.,
                    NULL,
                    NULL,
                    thermothermocolsslave
                    );

                // finalize temporary matrix for slave-side columns of thermo-thermo matrix block
                thermothermocolsslave.Complete(*strategythermo_->InterfaceMaps()->Map(1),*maps_->Map(1));

                // transform and assemble temporary matrix for slave-side columns of thermo-thermo matrix block
                blocksystemmatrix->Matrix(nblockmapsscatra,nblockmapsscatra).Add(*LINALG::MLMultiply(thermothermocolsslave,*strategythermo_->P(),true),false,1.,1.);

                break;
              }
              default:
              {
                dserror("Invalid type of scatra-scatra interface coupling!");
                break;
              }
            }
          }

          // assign thermo-thermo matrix block directly
          else
            blocksystemmatrix->Assign(nblockmapsscatra,nblockmapsscatra,LINALG::View,*thermo_->SystemMatrix());

          break;
        }

        case INPAR::S2I::matrix_sparse:
        {
          // construct global system matrix by assigning matrix blocks
          blocksystemmatrix->Assign(0,0,LINALG::View,*scatra_->SystemMatrix());

          // perform second condensation before assigning matrix blocks
          if(condensationthermo_)
          {
            const LINALG::SparseMatrix& scatrathermoblock = *Teuchos::rcp_dynamic_cast<const LINALG::SparseMatrix>(scatrathermoblock_);
            FSI::UTILS::MatrixLogicalSplitAndTransform()(
                scatrathermoblock,
                scatrathermoblock.RangeMap(),
                *maps_->Map(1),
                1.,
                NULL,
                NULL,
                blocksystemmatrix->Matrix(0,1)
                );

            FSI::UTILS::MatrixLogicalSplitAndTransform()(
                *Teuchos::rcp_dynamic_cast<const LINALG::SparseMatrix>(thermoscatrablock_),
                *maps_->Map(1),
                thermoscatrablock_->DomainMap(),
                1.,
                NULL,
                NULL,
                blocksystemmatrix->Matrix(1,0)
                );

            FSI::UTILS::MatrixLogicalSplitAndTransform()(
                *thermo_->SystemMatrix(),
                *maps_->Map(1),
                *maps_->Map(1),
                1.,
                NULL,
                NULL,
                blocksystemmatrix->Matrix(1,1)
                );

            switch(strategyscatra_->CouplingType())
            {
              case INPAR::S2I::coupling_matching_nodes:
              {
                ADAPTER::CouplingSlaveConverter converter(*icoupthermo_);
                FSI::UTILS::MatrixLogicalSplitAndTransform()(
                    scatrathermoblock,
                    scatrathermoblock.RangeMap(),
                    *strategythermo_->InterfaceMaps()->Map(1),
                    1.,
                    NULL,
                    &converter,
                    blocksystemmatrix->Matrix(0,1),
                    true,
                    true
                    );

                FSI::UTILS::MatrixLogicalSplitAndTransform()(
                    *thermo_->SystemMatrix(),
                    *maps_->Map(1),
                    *strategythermo_->InterfaceMaps()->Map(1),
                    1.,
                    NULL,
                    &converter,
                    blocksystemmatrix->Matrix(1,1),
                    true,
                    true
                    );

                break;
              }

              case INPAR::S2I::coupling_mortar_standard:
              {
                // initialize temporary matrix for slave-side columns of scatra-thermo matrix block
                LINALG::SparseMatrix scatrathermocolsslave(*scatra_->Discretization()->DofRowMap(),81);

                // fill temporary matrix for slave-side columns of scatra-thermo matrix block
                FSI::UTILS::MatrixLogicalSplitAndTransform()(
                    scatrathermoblock,
                    scatrathermoblock.RangeMap(),
                    *strategythermo_->InterfaceMaps()->Map(1),
                    1.,
                    NULL,
                    NULL,
                    scatrathermocolsslave
                    );

                // finalize temporary matrix for slave-side columns of scatra-thermo matrix block
                scatrathermocolsslave.Complete(*strategythermo_->InterfaceMaps()->Map(1),*scatra_->Discretization()->DofRowMap());

                // transform and assemble temporary matrix for slave-side columns of scatra-thermo matrix block
                blocksystemmatrix->Matrix(0,1).Add(*LINALG::MLMultiply(scatrathermocolsslave,*strategythermo_->P(),true),false,1.,1.);

                // initialize temporary matrix for slave-side columns of thermo-thermo matrix block
                LINALG::SparseMatrix thermothermocolsslave(*maps_->Map(1),81);

                // fill temporary matrix for slave-side columns of thermo-thermo matrix block
                FSI::UTILS::MatrixLogicalSplitAndTransform()(
                    *thermo_->SystemMatrix(),
                    *maps_->Map(1),
                    *strategythermo_->InterfaceMaps()->Map(1),
                    1.,
                    NULL,
                    NULL,
                    thermothermocolsslave
                    );

                // finalize temporary matrix for slave-side columns of thermo-thermo matrix block
                thermothermocolsslave.Complete(*strategythermo_->InterfaceMaps()->Map(1),*maps_->Map(1));

                // transform and assemble temporary matrix for slave-side columns of thermo-thermo matrix block
                blocksystemmatrix->Matrix(1,1).Add(*LINALG::MLMultiply(thermothermocolsslave,*strategythermo_->P(),true),false,1.,1.);

                break;
              }

              default:
              {
                dserror("Invalid type of scatra-scatra interface coupling!");
                break;
              }
            }
          }

          // assign matrix blocks directly
          else
          {
            blocksystemmatrix->Assign(0,1,LINALG::View,*Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(scatrathermoblock_));
            blocksystemmatrix->Assign(1,0,LINALG::View,*Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(thermoscatrablock_));
            blocksystemmatrix->Assign(1,1,LINALG::View,*thermo_->SystemMatrix());
          }

          break;
        }

        default:
        {
          dserror("Invalid matrix type associated with scalar transport field!");
          break;
        }
      }

      break;
    }

    case INPAR::STI::matrix_sparse:
    {
      // check global system matrix
      Teuchos::RCP<LINALG::SparseMatrix> systemmatrix = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(systemmatrix_);
      if(systemmatrix == Teuchos::null)
        dserror("System matrix is not a sparse matrix!");

      // construct global system matrix by adding matrix blocks
      systemmatrix->Add(*scatra_->SystemMatrix(),false,1.,0.);

      // perform second condensation before adding matrix blocks
      if(condensationthermo_)
      {
        const LINALG::SparseMatrix& scatrathermoblock = *Teuchos::rcp_dynamic_cast<const LINALG::SparseMatrix>(scatrathermoblock_);
        FSI::UTILS::MatrixLogicalSplitAndTransform()(
            scatrathermoblock,
            scatrathermoblock.RangeMap(),
            *maps_->Map(1),
            1.,
            NULL,
            NULL,
            *systemmatrix,
            true,
            true
            );

        FSI::UTILS::MatrixLogicalSplitAndTransform()(
            *Teuchos::rcp_dynamic_cast<const LINALG::SparseMatrix>(thermoscatrablock_),
            *maps_->Map(1),
            thermoscatrablock_->DomainMap(),
            1.,
            NULL,
            NULL,
            *systemmatrix,
            true,
            true
            );

        FSI::UTILS::MatrixLogicalSplitAndTransform()(
            *thermo_->SystemMatrix(),
            *maps_->Map(1),
            *maps_->Map(1),
            1.,
            NULL,
            NULL,
            *systemmatrix,
            true,
            true
            );

        switch(strategyscatra_->CouplingType())
        {
          case INPAR::S2I::coupling_matching_nodes:
          {
            ADAPTER::CouplingSlaveConverter converter(*icoupthermo_);
            FSI::UTILS::MatrixLogicalSplitAndTransform()(
                scatrathermoblock,
                scatrathermoblock.RangeMap(),
                *strategythermo_->InterfaceMaps()->Map(1),
                1.,
                NULL,
                &converter,
                *systemmatrix,
                true,
                true
                );

            FSI::UTILS::MatrixLogicalSplitAndTransform()(
                *thermo_->SystemMatrix(),
                *maps_->Map(1),
                *strategythermo_->InterfaceMaps()->Map(1),
                1.,
                NULL,
                &converter,
                *systemmatrix,
                true,
                true
                );

            break;
          }

          case INPAR::S2I::coupling_mortar_standard:
          {
            // initialize temporary matrix for slave-side columns of global system matrix
            LINALG::SparseMatrix systemmatrixcolsslave(*DofRowMap(),81);

            // fill temporary matrix for slave-side columns of global system matrix
            FSI::UTILS::MatrixLogicalSplitAndTransform()(
                scatrathermoblock,
                scatrathermoblock.RangeMap(),
                *strategythermo_->InterfaceMaps()->Map(1),
                1.,
                NULL,
                NULL,
                systemmatrixcolsslave
                );

            FSI::UTILS::MatrixLogicalSplitAndTransform()(
                *thermo_->SystemMatrix(),
                *maps_->Map(1),
                *strategythermo_->InterfaceMaps()->Map(1),
                1.,
                NULL,
                NULL,
                systemmatrixcolsslave,
                true,
                true
                );

            // finalize temporary matrix for slave-side columns of global system matrix
            systemmatrixcolsslave.Complete(*strategythermo_->InterfaceMaps()->Map(1),*DofRowMap());

            // transform and assemble temporary matrix for slave-side columns of global system matrix
            systemmatrix->Add(*LINALG::MLMultiply(systemmatrixcolsslave,*strategythermo_->P(),true),false,1.,1.);

            break;
          }

          default:
          {
            dserror("Invalid type of scatra-scatra interface coupling!");
            break;
          }
        }
      }

      // add matrix blocks directly
      else
      {
        systemmatrix->Add(*scatrathermoblock_,false,1.,1.);
        systemmatrix->Add(*thermoscatrablock_,false,1.,1.);
        systemmatrix->Add(*thermo_->SystemMatrix(),false,1.,1.);
      }

      break;
    }

    default:
    {
      dserror("Type of global system matrix for scatra-thermo interaction not recognized!");
      break;
    }
  }

  // finalize global system matrix
  systemmatrix_->Complete();

  // create full monolithic right-hand side vector
  maps_->InsertVector(scatra_->Residual(),0,residual_);
  Teuchos::RCP<Epetra_Vector> thermoresidual(Teuchos::null);
  if(condensationthermo_)
  {
    thermoresidual = Teuchos::rcp(new Epetra_Vector(*maps_->Map(1)));
    LINALG::Export(*thermo_->Residual(),*thermoresidual);
  }
  else
    thermoresidual = thermo_->Residual();
  maps_->InsertVector(thermoresidual,1,residual_);

  return;
} // STI::Algorithm::AssembleMatAndRHS()


/*-------------------------------------------------------------------------------*
 | build null spaces associated with blocks of global system matrix   fang 02/17 |
 *-------------------------------------------------------------------------------*/
void STI::Algorithm::BuildBlockNullSpaces() const
{
  switch(strategyscatra_->MatrixType())
  {
    case INPAR::S2I::matrix_block_condition:
    {
      // loop over block(s) of global system matrix associated with scalar transport field
      for(int iblock=0; iblock<strategyscatra_->BlockMaps().NumMaps(); ++iblock)
      {
        // store number of current block as string, starting from 1
        std::stringstream iblockstr;
        iblockstr << iblock+1;

        // equip smoother for current matrix block with empty parameter sublists to trigger null space computation
        Teuchos::ParameterList& blocksmootherparams = solver_->Params().sublist("Inverse"+iblockstr.str());
        blocksmootherparams.sublist("Aztec Parameters");
        blocksmootherparams.sublist("MueLu Parameters");

        // equip smoother for current matrix block with null space associated with all degrees of freedom on scalar transport discretization
        scatra_->Discretization()->ComputeNullSpaceIfNecessary(blocksmootherparams);

        // reduce full null space to match degrees of freedom associated with current matrix block
        LINALG::Solver::FixMLNullspace("Block "+iblockstr.str(),*scatra_->Discretization()->DofRowMap(),*strategyscatra_->BlockMaps().Map(iblock),blocksmootherparams);
      }

      break;
    }

    case INPAR::S2I::matrix_sparse:
    {
      // equip smoother for scatra matrix block with empty parameter sublists to trigger null space computation
      Teuchos::ParameterList& blocksmootherparams = solver_->Params().sublist("Inverse1");
      blocksmootherparams.sublist("Aztec Parameters");
      blocksmootherparams.sublist("MueLu Parameters");

      // equip smoother for scatra matrix block with null space associated with all degrees of freedom on scatra discretization
      scatra_->Discretization()->ComputeNullSpaceIfNecessary(blocksmootherparams);

      break;
    }

    default:
    {
      dserror("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }

  // store number of matrix block associated with temperature field as string
  std::stringstream iblockstr;
  iblockstr << blockmaps_->NumMaps();

  // equip smoother for thermo matrix block with empty parameter sublists to trigger null space computation
  Teuchos::ParameterList& blocksmootherparams = solver_->Params().sublist("Inverse"+iblockstr.str());
  blocksmootherparams.sublist("Aztec Parameters");
  blocksmootherparams.sublist("MueLu Parameters");

  // equip smoother for thermo matrix block with null space associated with all degrees of freedom on thermo discretization
  thermo_->Discretization()->ComputeNullSpaceIfNecessary(blocksmootherparams);

  // reduce full null space to match degrees of freedom associated with thermo matrix block if necessary
  if(condensationthermo_)
    LINALG::Solver::FixMLNullspace("Block "+iblockstr.str(),*thermo_->Discretization()->DofRowMap(),*maps_->Map(1),blocksmootherparams);

  return;
} // STI::Algorithm::BuildBlockNullSpaces


/*----------------------------------------------------------------------------*
 | compute inverse sums of absolute values of matrix row entries   fang 02/17 |
 *----------------------------------------------------------------------------*/
void STI::Algorithm::ComputeInvRowSums(
    const LINALG::SparseMatrix&          matrix,      //!< matrix
    const Teuchos::RCP<Epetra_Vector>&   invrowsums   //!< inverse sums of absolute values of row entries in matrix
    ) const
{
  // compute inverse row sums of matrix
  if(matrix.EpetraMatrix()->InvRowSums(*invrowsums))
    dserror("Inverse row sums of matrix could not be successfully computed!");

  return;
} // STI::Algorithm::ComputeInvRowSums


/*------------------------------------------------------------------------------------------------*
 | compute null space information associated with global system matrix if applicable   fang 06/17 |
 *------------------------------------------------------------------------------------------------*/
void STI::Algorithm::ComputeNullSpaceIfNecessary(
    Teuchos::ParameterList&   solverparams   //! solver parameter list for scatra-thermo interaction
    ) const
{
  // compute vector-based null space information for ML preconditioner
  if(solverparams.isSublist("ML Parameters"))
  {
    // extract parameter list for ML preconditioner
    Teuchos::ParameterList& mllist = solverparams.sublist("ML Parameters",true);

    // determine null space dimension
    const int numdofpernode_scatra = scatra_->NumDofPerNode();
    const int numdofpernode_thermo = thermo_->NumDofPerNode();
    const int dimns = numdofpernode_scatra+numdofpernode_thermo;

    // allocate vector for null space information
    const Teuchos::RCP<std::vector<double> > ns = Teuchos::rcp(new std::vector<double>(dimns*DofRowMap()->NumMyElements(),0.));

    // compute null space modes associated with scatra field
    const DRT::Discretization& scatradis = *scatra_->Discretization();
    double* modes_scatra[numdofpernode_scatra];
    for(int i=0; i<numdofpernode_scatra; ++i)
      modes_scatra[i] = &((*ns)[i*DofRowMap()->NumMyElements()]);
    for(int i=0; i<scatradis.NumMyRowNodes(); ++i)
    {
      const int lid = DofRowMap()->LID(scatradis.Dof(0,scatradis.lRowNode(i),0));
      if(lid < 0)
        dserror("Cannot find scatra degree of freedom!");
      for(int j=0; j<numdofpernode_scatra; ++j)
        modes_scatra[j][lid+j] = 1.;
    }

    // compute null space modes associated with thermo field
    const DRT::Discretization& thermodis = *thermo_->Discretization();
    double* modes_thermo[numdofpernode_thermo];
    for(int i=0; i<numdofpernode_thermo; ++i)
      modes_thermo[i] = &((*ns)[(numdofpernode_scatra+i)*DofRowMap()->NumMyElements()]);
    for(int i=0; i<thermodis.NumMyRowNodes(); ++i)
    {
      const int lid = DofRowMap()->LID(thermodis.Dof(0,thermodis.lRowNode(i),0));
      if(lid < 0)
        dserror("Cannot find thermo degree of freedom!");
      for(int j=0; j<numdofpernode_thermo; ++j)
        modes_thermo[j][lid+j] = 1.;
    }

    // fill parameter list
    mllist.set("PDE equations",dimns);
    mllist.set("null space: dimension",dimns);
    mllist.set("null space: type","pre-computed");
    mllist.set("null space: add default vectors",false);
    mllist.set<Teuchos::RCP<std::vector<double> > >("nullspace",ns);
    mllist.set("null space: vectors",&((*ns)[0]));
    mllist.set<bool>("ML validate parameter list",false);
  }

  // compute point-based null space information for MueLu preconditioner
  else if(solverparams.isSublist("MueLu Parameters"))
  {
    // extract and fill parameter list for MueLu preconditioner
    Teuchos::ParameterList& mllist = solverparams.sublist("MueLu Parameters",true);
    mllist.set("PDE equations",1);
    mllist.set("null space: dimension",1);
    mllist.set("null space: type","pre-computed");
    mllist.set("null space: add default vectors",false);
    const Teuchos::RCP<std::vector<double> > ns = Teuchos::rcp(new std::vector<double>(DofRowMap()->NumMyElements(),1.));
    mllist.set<Teuchos::RCP<std::vector<double> > >("nullspace",ns);
    mllist.set("null space: vectors",&((*ns)[0]));
    mllist.set<bool>("ML validate parameter list",false);
  }

  return;
} // STI::Algorithm::ComputeNullSpaceIfNecessary


/*----------------------------------------------------------------------*
 | global map of degrees of freedom                          fang 04/15 |
 *----------------------------------------------------------------------*/
const Teuchos::RCP<const Epetra_Map>& STI::Algorithm::DofRowMap() const
{
  return maps_->FullMap();
} // STI::Algorithm::DofRowMap()


/*----------------------------------------------------------------------*
 | equilibrate matrix rows                                   fang 02/17 |
 *----------------------------------------------------------------------*/
void STI::Algorithm::EquilibrateMatrixRows(
    LINALG::SparseMatrix&                matrix,      //!< matrix
    const Teuchos::RCP<Epetra_Vector>&   invrowsums   //!< sums of absolute values of row entries in matrix
    ) const
{
  if(matrix.LeftScale(*invrowsums))
    dserror("Row equilibration of matrix failed!");

  return;
} // STI::Algorithm::EquilibrateMatrixRows


/*----------------------------------------------------------------------*
 | equilibrate global system of equations if necessary       fang 02/17 |
 *----------------------------------------------------------------------*/
void STI::Algorithm::EquilibrateSystem(
    const Teuchos::RCP<LINALG::SparseOperator>&   systemmatrix,   //!< system matrix
    const Teuchos::RCP<Epetra_Vector>&            residual        //!< residual vector
    ) const
{
  switch(equilibration_)
  {
    case INPAR::S2I::equilibration_none:
    {
      // do nothing
      break;
    }

    case INPAR::S2I::equilibration_rows_full:
    case INPAR::S2I::equilibration_rows_maindiag:
    {
      switch(matrixtype_)
      {
        case INPAR::STI::matrix_block:
        {
          // check matrix
          const Teuchos::RCP<LINALG::BlockSparseMatrixBase> blocksparsematrix = Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(systemmatrix);
          if(blocksparsematrix == Teuchos::null)
            dserror("System matrix is not a block sparse matrix!");

          // perform row equilibration
          for(int i=0; i<blocksparsematrix->Rows(); ++i)
          {
            // initialize vector for inverse row sums
            const Teuchos::RCP<Epetra_Vector> invrowsums(Teuchos::rcp(new Epetra_Vector(blocksparsematrix->Matrix(i,i).RowMap())));

            // compute inverse row sums of current main diagonal matrix block
            if(equilibration_ == INPAR::S2I::equilibration_rows_maindiag)
              ComputeInvRowSums(blocksparsematrix->Matrix(i,i),invrowsums);

            // compute inverse row sums of current row block of global system matrix
            else
            {
              // loop over all column blocks of global system matrix
              for(int j=0; j<blocksparsematrix->Cols(); ++j)
              {
                // extract current block of global system matrix
                const LINALG::SparseMatrix& matrix = blocksparsematrix->Matrix(i,j);

                // loop over all rows of current matrix block
                for(int irow=0; irow<matrix.RowMap().NumMyElements(); ++irow)
                {
                  // determine length of current matrix row
                  const int length = matrix.EpetraMatrix()->NumMyEntries(irow);

                  if(length > 0)
                  {
                    // extract current matrix row from matrix block
                    int numentries(0);
                    std::vector<double> values(length,0.);
                    if(matrix.EpetraMatrix()->ExtractMyRowCopy(irow,length,numentries,&values[0]))
                      dserror("Cannot extract matrix row with local ID %d from matrix block!",irow);

                    // compute and store current row sum
                    double rowsum(0.);
                    for(int ientry=0; ientry<numentries; ++ientry)
                      rowsum += std::abs(values[ientry]);
                    (*invrowsums)[irow] += rowsum;
                  }
                }
              }

              // invert row sums
              if(invrowsums->Reciprocal(*invrowsums))
                dserror("Vector could not be inverted!");
            }

            // perform row equilibration of matrix blocks in current row block of global system matrix
            for(int j=0; j<blocksparsematrix->Cols(); ++j)
              EquilibrateMatrixRows(blocksparsematrix->Matrix(i,j),invrowsums);

            // insert inverse row sums of current main diagonal matrix block into global vector
            blockmaps_->InsertVector(invrowsums,i,invrowsums_);
          }

          break;
        }

        case INPAR::STI::matrix_sparse:
        {
          // check matrix
          const Teuchos::RCP<LINALG::SparseMatrix> sparsematrix = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(systemmatrix_);
          if(sparsematrix == Teuchos::null)
            dserror("System matrix is not a sparse matrix!");

          // compute inverse row sums of global system matrix
          ComputeInvRowSums(*sparsematrix,invrowsums_);

          // perform row equilibration of global system matrix
          EquilibrateMatrixRows(*sparsematrix,invrowsums_);

          break;
        }
      }

      // perform equilibration of global residual vector
      if(residual->Multiply(1.,*invrowsums_,*residual,0.))
         dserror("Equilibration of global residual vector failed!");

      break;
    }

    default:
    {
      dserror("Equilibration method not yet implemented!");
      break;
    }
  }

  return;
} // STI::Algorithm::EquilibrateSystem


/*-------------------------------------------------------------------------------------*
 | pass scatra degrees of freedom to thermo discretization and vice versa   fang 04/15 |
 *-------------------------------------------------------------------------------------*/
void STI::Algorithm::ExchangeStateVectors()
{
  // pass scatra degrees of freedom to thermo discretization and vice versa
  scatra_->Discretization()->SetState(2,"thermo",thermo_->Phiafnp());
  thermo_->Discretization()->SetState(2,"scatra",scatra_->Phiafnp());

  // exchange state vectors for evaluation of scatra-scatra interface coupling
  if(scatra_->S2ICoupling() and strategyscatra_->CouplingType() == INPAR::S2I::coupling_mortar_standard)
  {
    // extract scatra-scatra interface coupling conditions
    std::vector<DRT::Condition*> conditions;
    scatra_->Discretization()->GetCondition("S2ICoupling",conditions);

    // loop over all conditions
    for(unsigned icondition=0; icondition<conditions.size(); ++icondition)
      // consider conditions for slave side only
      if(conditions[icondition]->GetInt("interface side") == INPAR::S2I::side_slave)
      {
        // extract ID of current condition
        const int condid = conditions[icondition]->GetInt("ConditionID");
        if(condid < 0)
          dserror("Invalid condition ID!");

        // extract mortar discretizations associated with current condition
        DRT::Discretization& scatradis = strategyscatra_->MortarDiscretization(condid);
        DRT::Discretization& thermodis = strategythermo_->MortarDiscretization(condid);

        // pass interfacial scatra degrees of freedom to thermo discretization and vice versa
        const Teuchos::RCP<Epetra_Vector> scatra = Teuchos::rcp(new Epetra_Vector(*thermodis.DofRowMap(1)));
        LINALG::Export(*scatra_->Phiafnp(),*scatra);
        thermodis.SetState(1,"scatra",scatra);
        const Teuchos::RCP<Epetra_Vector> thermo = Teuchos::rcp(new Epetra_Vector(*scatradis.DofRowMap(1)));
        LINALG::Export(*thermo_->Phiafnp(),*thermo);
        scatradis.SetState(1,"thermo",thermo);
      }
  }

  return;
} // STI::Algorithm::ExchangeStateVectors()


/*-----------------------------------------------------------------------*
 | check termination criterion for Newton-Raphson iteration   fang 04/15 |
 *-----------------------------------------------------------------------*/
bool STI::Algorithm::ExitNewtonRaphson()
{
  // initialize exit flag
  bool exit(false);

  // compute vector norms for convergence check
  scatra_->Phinp()->Norm2(&scatradofnorm_);
  maps_->ExtractVector(residual_,0)->Norm2(&scatraresnorm_);
  maps_->ExtractVector(increment_,0)->Norm2(&scatraincnorm_);
  thermo_->Phinp()->Norm2(&thermodofnorm_);
  maps_->ExtractVector(residual_,1)->Norm2(&thermoresnorm_);
  maps_->ExtractVector(increment_,1)->Norm2(&thermoincnorm_);

  // safety checks
  if(std::isnan(scatradofnorm_) or
     std::isnan(scatraresnorm_) or
     std::isnan(scatraincnorm_) or
     std::isnan(thermodofnorm_) or
     std::isnan(thermoresnorm_) or
     std::isnan(thermoincnorm_))
    dserror("Vector norm is not a number!");
  if(std::isinf(scatradofnorm_) or
     std::isinf(scatraresnorm_) or
     std::isinf(scatraincnorm_) or
     std::isinf(thermodofnorm_) or
     std::isinf(thermoresnorm_) or
     std::isinf(thermoincnorm_))
    dserror("Vector norm is infinity!");

  // prevent division by zero
  if(scatradofnorm_ < 1.e-5)
    scatradofnorm_ = 1.e-5;
  if(thermodofnorm_ < 1.e-5)
    scatradofnorm_ = 1.e-5;

  // first Newton-Raphson iteration
  if(iter_ == 1)
  {
    // print first line of convergence table to screen
    // solution increment not yet available during first Newton-Raphson iteration
    if(Comm().MyPID() == 0)
      std::cout << "|  " << std::setw(3) << iter_ << "/" << std::setw(3) << itermax_ << "   | "
                << std::setw(10) << std::setprecision(3) << std::scientific << itertol_ << "[L_2 ]  | "
                << std::setw(10) << std::setprecision(3) << std::scientific << scatraresnorm_
                << "   |      --      | "
                << std::setw(10) << std::setprecision(3) << std::scientific << thermoresnorm_
                << "   |      --      | "
                << "(       --      , te = "
                << std::setw(10) << std::setprecision(3) << dtele_ << ")" << std::endl;
  }

  // subsequent Newton-Raphson iterations
  else
  {
    // print current line of convergence table to screen
    if(Comm().MyPID() == 0)
      std::cout << "|  " << std::setw(3) << iter_ << "/" << std::setw(3) << itermax_ << "   | "
                << std::setw(10) << std::setprecision(3) << std::scientific << itertol_ << "[L_2 ]  | "
                << std::setw(10) << std::setprecision(3) << std::scientific << scatraresnorm_ << "   | "
                << std::setw(10) << std::setprecision(3) << std::scientific << scatraincnorm_/scatradofnorm_ << "   | "
                << std::setw(10) << std::setprecision(3) << std::scientific << thermoresnorm_ << "   | "
                << std::setw(10) << std::setprecision(3) << std::scientific << thermoincnorm_/thermodofnorm_ << "   | (ts = "
                << std::setw(10) << std::setprecision(3) << dtsolve_ << ", te = "
                << std::setw(10) << std::setprecision(3) << dtele_ << ")" << std::endl;

    // convergence check
    if(scatraresnorm_ <= itertol_ and
       thermoresnorm_ <= itertol_ and
       scatraincnorm_/scatradofnorm_ <= itertol_ and
       thermoincnorm_/thermodofnorm_ <= itertol_)
      // exit Newton-Raphson iteration upon convergence
      exit = true;
  }

  // exit Newton-Raphson iteration when residuals are small enough to prevent unnecessary additional solver calls
  if(scatraresnorm_ < restol_ and thermoresnorm_ < restol_)
    exit = true;

  // print warning to screen if maximum number of Newton-Raphson iterations is reached without convergence
  if(iter_ == itermax_)
  {
    if(Comm().MyPID() == 0)
    {
      std::cout << "+------------+-------------------+--------------+--------------+--------------+--------------+" << std::endl;
      std::cout << "|      Newton-Raphson method has not converged after a maximum number of " << std::setw(2) << itermax_ << " iterations!      |" << std::endl;
    }

    // proceed to next time step
    exit = true;
  }

  // print finish line of convergence table to screen
  if(exit and Comm().MyPID() == 0)
    std::cout << "+------------+-------------------+--------------+--------------+--------------+--------------+" << std::endl << std::endl;

  return exit;
} // STI::Algorithm::ExitNewtonRaphson()


/*----------------------------------------------------------------------*
 | modify field parameters for thermo field                  fang 06/15 |
 *----------------------------------------------------------------------*/
void STI::Algorithm::ModifyFieldParametersForThermoField()
{
  // extract parameters for initial temperature field from parameter list for scatra-thermo interaction
  // and overwrite corresponding parameters in parameter list for thermo field
  if(!fieldparameters_->isParameter("INITIALFIELD") or !fieldparameters_->isParameter("INITFUNCNO"))
    dserror("Initial field parameters not properly set in input file section SCALAR TRANSPORT DYNAMIC!");
  if(!stiparameters_->isParameter("THERMO_INITIALFIELD") or !stiparameters_->isParameter("THERMO_INITFUNCNO"))
    dserror("Initial field parameters not properly set in input file section SCALAR TRANSPORT DYNAMIC!");
  fieldparameters_->set<std::string>("INITIALFIELD",stiparameters_->get<std::string>("THERMO_INITIALFIELD"));
  fieldparameters_->set<int>("INITFUNCNO",stiparameters_->get<int>("THERMO_INITFUNCNO"));

  // perform additional manipulations associated with scatra-scatra interface coupling
  if(scatra_->S2ICoupling())
  {
    // set flag for matrix type associated with thermo field
    fieldparameters_->sublist("S2I COUPLING").set<std::string>("MATRIXTYPE","sparse");

    // set flag in thermo meshtying strategy for evaluation of interface linearizations and residuals on slave side only
    fieldparameters_->sublist("S2I COUPLING").set<std::string>("SLAVEONLY","Yes");

    // adapt type of meshtying method for thermo field
    if(fieldparameters_->sublist("S2I COUPLING").get<std::string>("COUPLINGTYPE") == "StandardMortar")
      fieldparameters_->sublist("S2I COUPLING").set<std::string>("COUPLINGTYPE","CondensedMortar_Bubnov");
    else if(fieldparameters_->sublist("S2I COUPLING").get<std::string>("COUPLINGTYPE") != "MatchingNodes")
      dserror("Invalid type of scatra-scatra interface coupling!");

    // make sure that interface side underlying Lagrange multiplier definition is slave side
    fieldparameters_->sublist("S2I COUPLING").set<std::string>("LMSIDE","slave");
  }

  return;
} // STI::Algorithm::ModifyFieldParametersForThermoField()


/*----------------------------------------------------------------------*
 | output solution to screen and files                       fang 04/15 |
 *----------------------------------------------------------------------*/
void STI::Algorithm::Output()
{
  // output scatra field
  scatra_->Output();

  // output thermo field
  thermo_->Output();

  return;
}


/*---------------------------------------------------------------------------------------------------------------------------------------------------------------*
 | output matrix to *.csv file for debugging purposes, with global row and column IDs of matrix components in ascending order across all processors   fang 01/17 |
 *---------------------------------------------------------------------------------------------------------------------------------------------------------------*/
void STI::Algorithm::OutputMatrixToFile(
    const Teuchos::RCP<const LINALG::SparseOperator>   sparseoperator,   //!< sparse or block sparse matrix to be output
    const int                                          precision,        //!< output precision
    const double                                       tolerance         //!< output omission tolerance
    )
{
  // safety check
  if(!sparseoperator->Filled())
    dserror("Sparse operator must be filled for output!");

  // extract communicator
  const Epetra_Comm& comm = sparseoperator->Comm();

  // determine whether sparse matrix or block sparse matrix should be output
  const Teuchos::RCP<const LINALG::SparseMatrix> sparsematrix = Teuchos::rcp_dynamic_cast<const LINALG::SparseMatrix>(sparseoperator);
  const Teuchos::RCP<const LINALG::BlockSparseMatrixBase> blocksparsematrix = Teuchos::rcp_dynamic_cast<const LINALG::BlockSparseMatrixBase>(sparseoperator);
  if(sparsematrix == Teuchos::null and blocksparsematrix == Teuchos::null)
    dserror("Unknown type of sparse operator!");

  // extract row map
  const Epetra_Map& rowmap = sparsematrix != Teuchos::null ? sparsematrix->RowMap() : blocksparsematrix->FullRowMap();

  // safety check
  if(!rowmap.UniqueGIDs())
    dserror("Row map of matrix must be non-overlapping!");

  // copy global IDs of matrix rows stored on current processor into vector
  std::vector<int> myrowgids(rowmap.NumMyElements(),0);
  int* myglobalelements = rowmap.MyGlobalElements();
  std::copy(myglobalelements,myglobalelements+rowmap.NumMyElements(),&myrowgids[0]);

  // communicate global IDs
  std::vector<int> rowgids(0,0);
  LINALG::AllreduceVector(myrowgids,rowgids,comm);

  // retain communicated global IDs only on processor with ID 0
  if(comm.MyPID())
    rowgids.clear();

  // create full row map on processor with ID 0
  const Epetra_Map fullrowmap(-1,rowgids.size(),rowgids.size() ? &rowgids[0] : NULL,0,comm);

  // import matrix to processor with ID 0
  Epetra_CrsMatrix crsmatrix(Copy,fullrowmap,0);
  if(sparsematrix != Teuchos::null)
  {
    if(crsmatrix.Import(*sparsematrix->EpetraMatrix(),Epetra_Import(fullrowmap,rowmap),Insert))
      dserror("Matrix import failed!");
  }
  else
  {
    for(int i=0; i<blocksparsematrix->Rows(); ++i)
      for(int j=0; j<blocksparsematrix->Cols(); ++j)
        if(crsmatrix.Import(*blocksparsematrix->Matrix(i,j).EpetraMatrix(),Epetra_Import(fullrowmap,blocksparsematrix->RangeMap(i)),Insert))
          dserror("Matrix import failed!");
  }

  // let processor with ID 0 output matrix to file
  if(comm.MyPID() == 0)
  {
    // set file name
    std::ostringstream nproc;
    nproc << comm.NumProc();
    const std::string filename(DRT::Problem::Instance()->OutputControlFile()->FileName()+".matrix_"+nproc.str()+"proc.csv");

    // open file and write header at beginning
    std::ofstream file;
    file.open(filename.c_str(),std::fstream::trunc);
    file << std::setprecision(precision) << std::scientific << "RowGIDs,ColumnGIDs,Values" << std::endl;

    // write matrix to file
    for(int rowlid=0; rowlid<crsmatrix.NumMyRows(); ++rowlid)
    {
      // extract global ID of current matrix row
      const int rowgid = fullrowmap.GID(rowlid);

      // extract current matrix row
      int numentries;
      double* values;
      int* indices;
      if(crsmatrix.ExtractGlobalRowView(rowgid,numentries,values,indices))
        dserror("Cannot extract matrix row with global ID %d!",rowgid);

      // sort entries in current matrix row in ascending order of column global ID via map
      std::map<int,double> entries;

      // loop over all entries in current matrix row
      for(int j=0; j<numentries; ++j)
        // add current matrix entry to map
        entries[indices[j]] = values[j];

      // loop over all sorted entries in current matrix row
      for(std::map<int,double>::iterator j=entries.begin(); j!=entries.end(); ++j)
        // write current matrix entry to file
        if(std::abs(j->second) > tolerance)
          file << rowgid << "," << j->first << "," << j->second << std::endl;
    }

    // close file
    file.close();
  }

  // wait until output is complete
  comm.Barrier();

  // throw error to abort simulation for debugging
  dserror("Matrix was output to *.csv file!");

  return;
}


/*------------------------------------------------------------------------------------------------------------------------------------------------*
 | output vector to *.csv file for debugging purposes, with global IDs of vector components in ascending order across all processors   fang 01/17 |
 *------------------------------------------------------------------------------------------------------------------------------------------------*/
void STI::Algorithm::OutputVectorToFile(
    const Epetra_MultiVector&   vector,      //!< vector to be output
    const int                   precision,   //!< output precision
    const double                tolerance    //!< output omission tolerance
    )
{
  // extract communicator
  const Epetra_Comm& comm = vector.Comm();

  // extract vector map
  const Epetra_BlockMap& map = vector.Map();

  // safety check
  if(!map.UniqueGIDs())
    dserror("Vector output to *.csv file currently only works for non-overlapping vector maps!");

  // copy global IDs of vector components stored on current processor into vector
  std::vector<int> mygids(map.NumMyElements(),0);
  int* myglobalelements = map.MyGlobalElements();
  std::copy(myglobalelements,myglobalelements+map.NumMyElements(),&mygids[0]);

  // communicate global IDs
  std::vector<int> gids(0,0);
  LINALG::AllreduceVector(mygids,gids,comm);

  // retain communicated global IDs only on processor with ID 0
  if(comm.MyPID())
    gids.clear();

  // create full vector map on processor with ID 0
  const Epetra_Map fullmap(-1,gids.size(),gids.size() ? &gids[0] : NULL,0,comm);

  // export vector to processor with ID 0
  Epetra_MultiVector fullvector(fullmap,vector.NumVectors(),true);
  LINALG::Export(vector,fullvector);

  // let processor with ID 0 output vector to file
  if(comm.MyPID() == 0)
  {
    // set file name
    std::ostringstream nproc;
    nproc << comm.NumProc();
    const std::string filename(DRT::Problem::Instance()->OutputControlFile()->FileName()+".vector_"+nproc.str()+"proc.csv");

    // open file and write header at beginning
    std::ofstream file;
    file.open(filename.c_str(),std::fstream::trunc);
    file << std::setprecision(precision) << std::scientific << "GIDs,Values" << std::endl;

    // write vector to file
    for(int lid=0; lid<fullvector.MyLength(); ++lid)
    {
      // inner loop index
      int j(-1);

      // check output omission tolerance
      for(j=0; j<fullvector.NumVectors(); ++j)
        if(std::abs(fullvector[j][lid]) > tolerance)
          break;

      // perform output if applicable
      if(j<fullvector.NumVectors())
      {
        // write global ID of current vector component
        file << fullmap.GID(lid);

        // loop over all subvectors
        for(j=0; j<fullvector.NumVectors(); ++j)
          // write current vector component to file
          file << "," << fullvector[j][lid];

        // output line break
        file << std::endl;
      }
    }

    // close file
    file.close();
  }

  // wait until output is complete
  comm.Barrier();

  // throw error to abort simulation for debugging
  dserror("Vector was output to *.csv file!");

  return;
}


/*----------------------------------------------------------------------*
 | prepare time step                                         fang 04/15 |
 *----------------------------------------------------------------------*/
void STI::Algorithm::PrepareTimeStep()
{
  // update time and time step
  IncrementTimeAndStep();

  // provide scatra and thermo fields with velocities
  scatra_->SetVelocityField(1);
  thermo_->SetVelocityField(1);

  // pass scatra degrees of freedom to thermo discretization and vice versa
  ExchangeStateVectors();

  // prepare time step for scatra field
  scatra_->PrepareTimeStep();

  // pass scatra degrees of freedom to thermo discretization and vice versa
  // this only needs to be done for the first time step, when the initial values of the electric potential state variables are computed
  if(Step() == 1)
  {
    ExchangeStateVectors();

    // pass master-side scatra degrees of freedom to thermo discretization for evaluation of scatra-scatra interface coupling
    if(thermo_->S2ICoupling() and strategythermo_->CouplingType() == INPAR::S2I::coupling_matching_nodes)
      thermo_->Discretization()->SetState(2,"imasterscatra",strategyscatra_->MasterPhinp());
  }

  // prepare time step for thermo field
  thermo_->PrepareTimeStep();

  // print time step information to screen
  scatra_->PrintTimeStepInfo();

  return;
} // STI::Algorithm::PrepareTimeStep()


/*----------------------------------------------------------------------*
 | read restart data                                         fang 04/15 |
 *----------------------------------------------------------------------*/
void STI::Algorithm::ReadRestart(
    int step   //! time step for restart
    )
{
  // read scatra and thermo restart variables
  scatra_->ReadRestart(step);
  thermo_->ReadRestart(step);

  // pass scatra degrees of freedom to thermo discretization and vice versa
  ExchangeStateVectors();

  // set time and time step
  SetTimeStep(scatra_->Time(),step);

  return;
} // STI::Algorithm::ReadRestart


/*----------------------------------------------------------------------*
 | evaluate time step using Newton-Raphson iteration         fang 04/15 |
 *----------------------------------------------------------------------*/
void STI::Algorithm::Solve()
{
  // initialize counter for Newton-Raphson iterations
  iter_ = 0;

  // print header of convergence table to screen
  if(Comm().MyPID() == 0)
  {
    std::cout << "+------------+-------------------+--------------+--------------+--------------+--------------+" << std::endl;
    std::cout << "|- step/max -|- tolerance[norm] -|- scatra-res -|- scatra-inc -|- thermo-res -|- thermo-inc -|" << std::endl;
  }

  // start Newton-Raphson iteration
  while(true)
  {
    // update iteration counter
    iter_ += 1;

    // reset timer
    timer_->ResetStartTime();

    // store time before evaluating elements and assembling global system of equations
    double time = timer_->WallTime();

    // assemble global system of equations
    AssembleMatAndRHS();

    // determine time needed for evaluating elements and assembling global system of equations,
    // and take average over all processors via communication
    double mydtele = timer_->WallTime()-time;
    Comm().MaxAll(&mydtele,&dtele_,1);

    // safety check
    if(!systemmatrix_->Filled())
      dserror("Complete() has not been called on global system matrix yet!");

    // perform finite difference check on time integrator level
    if(scatra_->FDCheckType() == INPAR::SCATRA::fdcheck_global)
      FDCheck();

    // check termination criterion for Newton-Raphson iteration
    if(ExitNewtonRaphson())
      break;

    // initialize global increment vector
    increment_->PutScalar(0.);

    // store time before solving global system of equations
    time = timer_->WallTime();

    // equilibrate global system of equations if necessary
    EquilibrateSystem(systemmatrix_,residual_);

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
    // and take average over all processors via communication
    double mydtsolve = timer_->WallTime()-time;
    Comm().MaxAll(&mydtsolve,&dtsolve_,1);

    // output performance statistics associated with linear solver into text file if applicable
    if(DRT::INPUT::IntegralValue<int>(*fieldparameters_,"OUTPUTSOLVERSTATS"))
      scatra_->OutputSolverStats(*solver_,dtsolve_,Step(),iter_,residual_->Map().NumGlobalElements());

    // update scatra field
    scatra_->UpdateIter(maps_->ExtractVector(increment_,0));
    scatra_->ComputeIntermediateValues();

    // update thermo field
    Teuchos::RCP<Epetra_Vector> thermoincrement(Teuchos::null);
    if(condensationthermo_)
    {
      thermoincrement = Teuchos::rcp(new Epetra_Vector(*thermo_->Discretization()->DofRowMap()));
      LINALG::Export(*maps_->ExtractVector(increment_,1),*thermoincrement);
      const Teuchos::RCP<const Epetra_Vector> masterincrement = strategythermo_->InterfaceMaps()->ExtractVector(*thermoincrement,2);
      const Teuchos::RCP<Epetra_Vector> slaveincrement = LINALG::CreateVector(*strategythermo_->InterfaceMaps()->Map(1));
      switch(strategyscatra_->CouplingType())
      {
        case INPAR::S2I::coupling_matching_nodes:
        {
          icoupthermo_->MasterToSlave(masterincrement,slaveincrement);
          break;
        }
        case INPAR::S2I::coupling_mortar_standard:
        {
          strategythermo_->P()->Multiply(false,*masterincrement,*slaveincrement);
          break;
        }
        default:
        {
          dserror("Invalid type of scatra-scatra interface coupling!");
          break;
        }
      }
      strategythermo_->InterfaceMaps()->InsertVector(slaveincrement,1,thermoincrement);
    }
    else
      thermoincrement = maps_->ExtractVector(increment_,1);
    thermo_->UpdateIter(thermoincrement);
    thermo_->ComputeIntermediateValues();
  } // Newton-Raphson iteration

  return;
} // STI::Algorithm::Solve


/*--------------------------------------------------------------------------------*
 | assemble off-diagonal scatra-thermo block of global system matrix   fang 12/15 |
 *--------------------------------------------------------------------------------*/
void STI::Algorithm::AssembleODBlockScatraThermo()
{
  // initialize scatra-thermo matrix block
  scatrathermoblock_->Zero();

  // create parameter list for element evaluation
  Teuchos::ParameterList eleparams;

  // action for elements
  eleparams.set<int>("action",SCATRA::calc_scatra_mono_odblock_scatrathermo);

  // number of dofset associated with velocity-related dofs on scatra discretization
  eleparams.set<int>("ndsvel",1);

  // remove state vectors from scatra discretization
  scatra_->Discretization()->ClearState();

  // add state vectors to scatra discretization
  scatra_->AddTimeIntegrationSpecificVectors();

  // create strategy for assembly of scatra-thermo matrix block
  DRT::AssembleStrategy strategyscatrathermo(
      0,                    // row assembly based on number of dofset associated with scatra dofs on scatra discretization
      2,                    // column assembly based on number of dofset associated with thermo dofs on scatra discretization
      scatrathermoblock_,   // scatra-thermo matrix block
      Teuchos::null,        // no additional matrices or vectors
      Teuchos::null,
      Teuchos::null,
      Teuchos::null
      );

  // assemble scatra-thermo matrix block
  scatra_->Discretization()->Evaluate(eleparams,strategyscatrathermo);

  // provide scatra-thermo matrix block with contributions from scatra-scatra interface coupling if applicable
  if(scatra_->S2ICoupling())
  {
    // differentiate between different meshtying methods
    switch(strategyscatra_->CouplingType())
    {
      case INPAR::S2I::coupling_matching_nodes:
      {
        // create parameter list for element evaluation
        Teuchos::ParameterList condparams;

        // action for elements
        condparams.set<int>("action",SCATRA::bd_calc_s2icoupling_od);

        // add state vector containing master-side scatra degrees of freedom to scatra discretization
        scatra_->Discretization()->SetState("imasterphinp",strategyscatra_->MasterPhinp());

        switch(strategyscatra_->MatrixType())
        {
          case INPAR::S2I::matrix_block_condition:
          {
            // initialize auxiliary system matrix for linearizations of slave-side scatra fluxes w.r.t. slave-side thermo dofs
            Teuchos::RCP<LINALG::BlockSparseMatrixBase> blockslavematrix = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
                *blockmapthermo_,
                strategyscatra_->BlockMapsSlave(),
                81,
                false,
                true
                ));

            // create strategy for assembly of auxiliary system matrix
            DRT::AssembleStrategy strategyscatrathermos2i(
                0,                  // row assembly based on number of dofset associated with scatra dofs on scatra discretization
                2,                  // column assembly based on number of dofset associated with thermo dofs on scatra discretization
                blockslavematrix,   // auxiliary system matrix
                Teuchos::null,      // no additional matrices of vectors
                Teuchos::null,
                Teuchos::null,
                Teuchos::null
                );

            // evaluate scatra-scatra interface coupling
            std::vector<DRT::Condition*> conditions;
            scatra_->Discretization()->GetCondition("S2ICoupling",conditions);
            for(unsigned icondition=0; icondition<conditions.size(); ++icondition)
              if(conditions[icondition]->GetInt("interface side") == INPAR::S2I::side_slave)
                scatra_->Discretization()->EvaluateCondition(condparams,strategyscatrathermos2i,"S2ICoupling",conditions[icondition]->GetInt("ConditionID"));

            // finalize auxiliary system matrix
            blockslavematrix->Complete();

            // assemble linearizations of slave-side scatra fluxes w.r.t. slave-side thermo dofs into scatra-thermo matrix block
            scatrathermoblock_->Add(*blockslavematrix,false,1.,1.);

            // initialize auxiliary system matrix for linearizations of master-side scatra fluxes w.r.t. slave-side thermo dofs
            LINALG::SparseMatrix mastermatrix(*icoupscatra_->MasterDofMap(),27,false,true);

            // derive linearizations of master-side scatra fluxes w.r.t. slave-side thermo dofs and assemble into auxiliary system matrix
            for(int iblock=0; iblock<strategyscatra_->BlockMapsSlave().NumMaps(); ++iblock)
              FSI::UTILS::MatrixRowTransform()(
                  blockslavematrix->Matrix(iblock,0),
                  -1.,
                  ADAPTER::CouplingSlaveConverter(*icoupscatra_),
                  mastermatrix,
                  true
                  );

            // finalize auxiliary system matrix
            mastermatrix.Complete(*icoupthermo_->SlaveDofMap(),*icoupscatra_->MasterDofMap());

            // split auxiliary system matrix and assemble into scatra-thermo matrix block
            const Teuchos::RCP<LINALG::BlockSparseMatrixBase> blockmastermatrix = mastermatrix.Split<LINALG::DefaultBlockMatrixStrategy>(*blockmapthermo_,strategyscatra_->BlockMaps());
            blockmastermatrix->Complete();
            scatrathermoblock_->Add(*blockmastermatrix,false,1.,1.);

            // linearizations of scatra fluxes w.r.t. master-side thermo dofs are not needed, since these dofs will be condensed out later

            break;
          }

          case INPAR::S2I::matrix_sparse:
          {
            // initialize auxiliary system matrix for linearizations of slave-side scatra fluxes w.r.t. slave-side thermo dofs
            strategyscatra_->SlaveMatrix()->Zero();

            // create strategy for assembly of auxiliary system matrix
            DRT::AssembleStrategy strategyscatrathermos2i(
                0,                                // row assembly based on number of dofset associated with scatra dofs on scatra discretization
                2,                                // column assembly based on number of dofset associated with thermo dofs on scatra discretization
                strategyscatra_->SlaveMatrix(),   // auxiliary system matrix
                Teuchos::null,                    // no additional matrices of vectors
                Teuchos::null,
                Teuchos::null,
                Teuchos::null
                );

            // evaluate scatra-scatra interface coupling
            std::vector<DRT::Condition*> conditions;
            scatra_->Discretization()->GetCondition("S2ICoupling",conditions);
            for(unsigned icondition=0; icondition<conditions.size(); ++icondition)
              if(conditions[icondition]->GetInt("interface side") == INPAR::S2I::side_slave)
                scatra_->Discretization()->EvaluateCondition(condparams,strategyscatrathermos2i,"S2ICoupling",conditions[icondition]->GetInt("ConditionID"));

            // finalize auxiliary system matrix
            strategyscatra_->SlaveMatrix()->Complete(*thermo_->Discretization()->DofRowMap(),*maps_->Map(0));

            // assemble linearizations of slave-side scatra fluxes w.r.t. slave-side thermo dofs into scatra-thermo matrix block
            scatrathermoblock_->Add(*strategyscatra_->SlaveMatrix(),false,1.,1.);

            // derive linearizations of master-side scatra fluxes w.r.t. slave-side thermo dofs and assemble into scatra-thermo matrix block
            (*islavetomasterrowtransformscatraod_)(
                *strategyscatra_->SlaveMatrix(),
                -1.,
                ADAPTER::CouplingSlaveConverter(*icoupscatra_),
                *Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(scatrathermoblock_),
                true
                );

            // linearizations of scatra fluxes w.r.t. master-side thermo dofs are not needed, since these dofs will be condensed out later

            break;
          }

          default:
          {
            dserror("Invalid matrix type associated with scalar transport field!");
            break;
          }
        }

        break;
      }

      case INPAR::S2I::coupling_mortar_standard:
      {
        // initialize auxiliary system matrices for linearizations of slave-side and master-side scatra fluxes w.r.t. slave-side thermo dofs
        Teuchos::RCP<LINALG::SparseOperator> slavematrix(Teuchos::null);
        switch(strategyscatra_->MatrixType())
        {
          case INPAR::S2I::matrix_block_condition:
          {
            slavematrix = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
                *blockmapthermo_,
                strategyscatra_->BlockMapsSlave(),
                81,
                false,
                true
                ));
            break;
          }

          case INPAR::S2I::matrix_sparse:
          {
            slavematrix = strategyscatra_->SlaveMatrix();
            slavematrix->Zero();
            break;
          }

          default:
          {
            dserror("Invalid matrix type associated with scalar transport field!");
            break;
          }
        }
        strategyscatra_->MasterMatrix()->Zero();

        // create parameter list for element evaluation
        Teuchos::ParameterList condparams;

        // action for elements
        condparams.set<int>("action",INPAR::S2I::evaluate_condition_od);

        // create strategy for assembly of auxiliary system matrices
        SCATRA::MortarCellAssemblyStrategy strategyscatrathermos2i(
            slavematrix,
            INPAR::S2I::side_slave,
            INPAR::S2I::side_slave,
            Teuchos::null,
            INPAR::S2I::side_undefined,
            INPAR::S2I::side_undefined,
            strategyscatra_->MasterMatrix(),
            INPAR::S2I::side_master,
            INPAR::S2I::side_slave,
            Teuchos::null,
            INPAR::S2I::side_undefined,
            INPAR::S2I::side_undefined,
            Teuchos::null,
            INPAR::S2I::side_undefined,
            Teuchos::null,
            INPAR::S2I::side_undefined,
            0,
            1
            );

        // extract scatra-scatra interface coupling conditions
        std::vector<DRT::Condition*> conditions;
        scatra_->Discretization()->GetCondition("S2ICoupling",conditions);

        // loop over all conditions
        for(unsigned icondition=0; icondition<conditions.size(); ++icondition)
        {
          // extract current condition
          DRT::Condition& condition = *conditions[icondition];

          // consider conditions for slave side only
          if(condition.GetInt("interface side") == INPAR::S2I::side_slave)
          {
            // add condition to parameter list
            condparams.set<DRT::Condition*>("condition",&condition);

            // evaluate mortar integration cells
            strategyscatra_->EvaluateMortarCells(strategyscatra_->MortarDiscretization(condition.GetInt("ConditionID")),condparams,strategyscatrathermos2i);
          }
        }

        // finalize auxiliary system matrices
        strategyscatra_->MasterMatrix()->Complete(*thermo_->Discretization()->DofRowMap(),*maps_->Map(0));
        Teuchos::RCP<LINALG::SparseOperator> mastermatrix(Teuchos::null);
        switch(strategyscatra_->MatrixType())
        {
          case INPAR::S2I::matrix_block_condition:
          {
            slavematrix->Complete();
            mastermatrix = strategyscatra_->MasterMatrix()->Split<LINALG::DefaultBlockMatrixStrategy>(*blockmapthermo_,strategyscatra_->BlockMapsMaster());
            mastermatrix->Complete();

            break;
          }

          case INPAR::S2I::matrix_sparse:
          {
            slavematrix->Complete(*thermo_->Discretization()->DofRowMap(),*maps_->Map(0));
            mastermatrix = strategyscatra_->MasterMatrix();

            break;
          }

          default:
          {
            dserror("Invalid matrix type associated with scalar transport field!");
            break;
          }
        }

        // assemble linearizations of slave-side and master-side scatra fluxes w.r.t. slave-side thermo dofs into scatra-thermo matrix block
        scatrathermoblock_->Add(*slavematrix,false,1.,1.);
        scatrathermoblock_->Add(*mastermatrix,false,1.,1.);

        // linearizations of scatra fluxes w.r.t. master-side thermo dofs are not needed, since these dofs will be condensed out later

        break;
      }

      default:
      {
        dserror("Invalid type of scatra-scatra interface coupling!");
        break;
      }
    }
  }

  // finalize scatra-thermo matrix block
  switch(strategyscatra_->MatrixType())
  {
    case INPAR::S2I::matrix_block_condition:
    {
      scatrathermoblock_->Complete();
      break;
    }

    case INPAR::S2I::matrix_sparse:
    {
      scatrathermoblock_->Complete(*thermo_->Discretization()->DofRowMap(),*maps_->Map(0));
      break;
    }

    default:
    {
      dserror("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }

  // apply Dirichlet boundary conditions to scatra-thermo matrix block
  scatrathermoblock_->ApplyDirichlet(*scatra_->DirichMaps()->CondMap(),false);

  // remove state vectors from scatra discretization
  scatra_->Discretization()->ClearState();

  return;
} // STI::Algorithm::AssembleODBlockScatraThermo()


/*--------------------------------------------------------------------------------*
 | assemble off-diagonal thermo-scatra block of global system matrix   fang 12/15 |
 *--------------------------------------------------------------------------------*/
void STI::Algorithm::AssembleODBlockThermoScatra()
{
  // initialize thermo-scatra matrix block
  thermoscatrablock_->Zero();

  // create parameter list for element evaluation
  Teuchos::ParameterList eleparams;

  // action for elements
  eleparams.set<int>("action",SCATRA::calc_scatra_mono_odblock_thermoscatra);

  // number of dofset associated with velocity-related dofs on thermo discretization
  eleparams.set<int>("ndsvel",1);

  // remove state vectors from thermo discretization
  thermo_->Discretization()->ClearState();

  // add state vectors to thermo discretization
  thermo_->AddTimeIntegrationSpecificVectors();

  // create strategy for assembly of thermo-scatra matrix block
  DRT::AssembleStrategy strategythermoscatra(
      0,                    // row assembly based on number of dofset associated with thermo dofs on thermo discretization
      2,                    // column assembly based on number of dofset associated with scatra dofs on thermo discretization
      thermoscatrablock_,   // thermo-scatra matrix block
      Teuchos::null,        // no additional matrices or vectors
      Teuchos::null,
      Teuchos::null,
      Teuchos::null
      );

  // assemble thermo-scatra matrix block
  thermo_->Discretization()->Evaluate(eleparams,strategythermoscatra);

  // provide thermo-scatra matrix block with contributions from scatra-scatra interface coupling if applicable
  if(thermo_->S2ICoupling())
  {
    // differentiate between different meshtying methods
    switch(strategythermo_->CouplingType())
    {
      case INPAR::S2I::coupling_matching_nodes:
      {
        // initialize auxiliary system matrix for linearizations of slave-side thermo fluxes w.r.t. master-side scatra dofs
        strategythermo_->MasterMatrix()->Zero();

        // create parameter list for element evaluation
        Teuchos::ParameterList condparams;

        // action for elements
        condparams.set<int>("action",SCATRA::bd_calc_s2icoupling_od);

        switch(strategyscatra_->MatrixType())
        {
          case INPAR::S2I::matrix_block_condition:
          {
            // initialize auxiliary system matrix for linearizations of slave-side thermo fluxes w.r.t. slave-side scatra dofs
            Teuchos::RCP<LINALG::BlockSparseMatrixBase> slavematrix = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
                strategyscatra_->BlockMapsSlave(),
                *blockmapthermo_,
                81,
                false,
                true
                ));

            // create strategy for assembly of auxiliary system matrices
            DRT::AssembleStrategy strategythermoscatras2i(
                0,                                 // row assembly based on number of dofset associated with thermo dofs on thermo discretization
                2,                                 // column assembly based on number of dofset associated with scatra dofs on thermo discretization
                slavematrix,                       // auxiliary system matrix for slave side
                strategythermo_->MasterMatrix(),   // auxiliary system matrix for master side
                Teuchos::null,                     // no additional matrices of vectors
                Teuchos::null,
                Teuchos::null
                );

            // evaluate scatra-scatra interface coupling
            std::vector<DRT::Condition*> conditions;
            thermo_->Discretization()->GetCondition("S2ICoupling",conditions);
            for(unsigned icondition=0; icondition<conditions.size(); ++icondition)
              if(conditions[icondition]->GetInt("interface side") == INPAR::S2I::side_slave)
                thermo_->Discretization()->EvaluateCondition(condparams,strategythermoscatras2i,"S2ICoupling",conditions[icondition]->GetInt("ConditionID"));

            // finalize auxiliary system matrices
            slavematrix->Complete();
            strategythermo_->MasterMatrix()->Complete(*icoupscatra_->SlaveDofMap(),*icoupthermo_->SlaveDofMap());

            // assemble linearizations of slave-side thermo fluxes w.r.t. slave-side scatra dofs into thermo-scatra matrix block
            thermoscatrablock_->Add(*slavematrix,false,1.,1.);

            // initialize temporary matrix
            LINALG::SparseMatrix ksm(*icoupthermo_->SlaveDofMap(),27,false,true);

            // transform linearizations of slave-side thermo fluxes w.r.t. master-side scatra dofs
            (*islavetomastercoltransformthermood_)(
                strategythermo_->MasterMatrix()->RowMap(),
                strategythermo_->MasterMatrix()->ColMap(),
                *strategythermo_->MasterMatrix(),
                1.,
                ADAPTER::CouplingSlaveConverter(*icoupscatra_),
                ksm,
                true,
                false
                );

            // finalize temporary matrix
            ksm.Complete(*icoupscatra_->MasterDofMap(),*icoupthermo_->SlaveDofMap());

            // split temporary matrix and assemble into thermo-scatra matrix block
            const Teuchos::RCP<LINALG::BlockSparseMatrixBase> blockksm(ksm.Split<LINALG::DefaultBlockMatrixStrategy>(strategyscatra_->BlockMaps(),*blockmapthermo_));
            blockksm->Complete();
            thermoscatrablock_->Add(*blockksm,false,1.,1.);

            break;
          }

          case INPAR::S2I::matrix_sparse:
          {
            // initialize auxiliary system matrix for linearizations of slave-side thermo fluxes w.r.t. slave-side scatra dofs
            strategythermo_->SlaveMatrix()->Zero();

            // create strategy for assembly of auxiliary system matrices
            DRT::AssembleStrategy strategythermoscatras2i(
                0,                                 // row assembly based on number of dofset associated with thermo dofs on thermo discretization
                2,                                 // column assembly based on number of dofset associated with scatra dofs on thermo discretization
                strategythermo_->SlaveMatrix(),    // auxiliary system matrix for slave side
                strategythermo_->MasterMatrix(),   // auxiliary system matrix for master side
                Teuchos::null,                     // no additional matrices of vectors
                Teuchos::null,
                Teuchos::null
                );

            // evaluate scatra-scatra interface coupling
            std::vector<DRT::Condition*> conditions;
            thermo_->Discretization()->GetCondition("S2ICoupling",conditions);
            for(unsigned icondition=0; icondition<conditions.size(); ++icondition)
              if(conditions[icondition]->GetInt("interface side") == INPAR::S2I::side_slave)
                thermo_->Discretization()->EvaluateCondition(condparams,strategythermoscatras2i,"S2ICoupling",conditions[icondition]->GetInt("ConditionID"));

            // finalize auxiliary system matrices
            strategythermo_->SlaveMatrix()->Complete(*icoupscatra_->SlaveDofMap(),*icoupthermo_->SlaveDofMap());
            strategythermo_->MasterMatrix()->Complete(*icoupscatra_->SlaveDofMap(),*icoupthermo_->SlaveDofMap());

            // assemble linearizations of slave-side thermo fluxes w.r.t. slave-side scatra dofs into thermo-scatra matrix block
            thermoscatrablock_->Add(*strategythermo_->SlaveMatrix(),false,1.,1.);

            // derive linearizations of slave-side thermo fluxes w.r.t. master-side scatra dofs and assemble into thermo-scatra matrix block
            (*islavetomastercoltransformthermood_)(
                strategythermo_->MasterMatrix()->RowMap(),
                strategythermo_->MasterMatrix()->ColMap(),
                *strategythermo_->MasterMatrix(),
                1.,
                ADAPTER::CouplingSlaveConverter(*icoupscatra_),
                *Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(thermoscatrablock_),
                true,
                true
                );

            break;
          }

          default:
          {
            dserror("Invalid matrix type associated with scalar transport field!");
            break;
          }
        }

        // linearizations of master-side thermo fluxes w.r.t. scatra dofs are not needed, since thermo fluxes are source terms and thus only evaluated once on slave side

        // standard meshtying algorithm with Lagrange multipliers condensed out
        if(!thermo_->Discretization()->GetCondition("PointCoupling"))
        {
          // during the very first run of the following code, Complete() has not yet been called on the thermo-scatra block of the global system matrix
          // experiments have shown that the ExtractMatrixRows routine called in the following does not properly work in this case if the thermo-scatra block exhibits a block structure
          // in particular, some matrix entries are simply omitted during extraction, although the underlying routine ExtractGlobalRowCopy returns a zero error code
          // as a consequence, the final thermo-scatra block lacks some entries, and so does the final matrix graph after calling Complete()
          // this leads to errors involving unknown global indices of matrix columns during the second run of the following code, with the thermo-scatra block now being Complete()
          // to circumvent this problem, we call Complete() on the thermo-scatra block before the very first run of the following code
          if(strategyscatra_->MatrixType() == INPAR::S2I::matrix_block_condition and Step() == 1 and iter_ == 1)
            thermoscatrablock_->Complete();

          // loop over all thermo-scatra matrix blocks
          for(int iblock=0; iblock<blockmaps_->NumMaps()-1; ++iblock)
          {
            // initialize temporary matrix for slave-side rows of current thermo-scatra matrix block
            LINALG::SparseMatrix thermoscatrarowsslave(*icoupthermo_->SlaveDofMap(),27,false,true);

            // extract current thermo-scatra matrix block
            LINALG::SparseMatrix& thermoscatrablock = strategyscatra_->MatrixType() == INPAR::S2I::matrix_block_condition ? Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(thermoscatrablock_)->Matrix(0,iblock) : *Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(thermoscatrablock_);

            // extract slave-side rows of thermo-scatra matrix block into temporary matrix
            SCATRA::MeshtyingStrategyS2I::ExtractMatrixRows(thermoscatrablock,thermoscatrarowsslave,*icoupthermo_->SlaveDofMap());

            // finalize temporary matrix with slave-side rows of thermo-scatra matrix block
            thermoscatrarowsslave.Complete(*maps_->Map(0),*icoupthermo_->SlaveDofMap());

            // undo Complete() from above before performing subsequent matrix row transformation
            if(strategyscatra_->MatrixType() == INPAR::S2I::matrix_block_condition and Step() == 1 and iter_ == 1)
              thermoscatrablock.UnComplete();

            // add slave-side rows of thermo-scatra matrix block to corresponding slave-side rows
            const Teuchos::RCP<FSI::UTILS::MatrixRowTransform> islavetomasterrowtransformthermood = strategyscatra_->MatrixType() == INPAR::S2I::matrix_block_condition ? Teuchos::rcp(new FSI::UTILS::MatrixRowTransform()) : islavetomasterrowtransformthermood_;
            (*islavetomasterrowtransformthermood)(
                thermoscatrarowsslave,
                1.,
                ADAPTER::CouplingSlaveConverter(*icoupthermo_),
                thermoscatrablock,
                true
                );
          }
        }

        break;
      }

      case INPAR::S2I::coupling_mortar_condensed_bubnov:
      {
        // initialize auxiliary system matrix for linearizations of slave-side thermo fluxes w.r.t. slave-side and master-side scatra dofs
        Teuchos::RCP<LINALG::SparseOperator> slavematrix(Teuchos::null);
        switch(strategyscatra_->MatrixType())
        {
          case INPAR::S2I::matrix_block_condition:
          {
            slavematrix = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
                strategyscatra_->BlockMaps(),
                *blockmapthermo_,
                81,
                false,
                true
                ));
            break;
          }

          case INPAR::S2I::matrix_sparse:
          {
            slavematrix = strategythermo_->SlaveMatrix();
            slavematrix->Zero();
            break;
          }

          default:
          {
            dserror("Invalid matrix type associated with scalar transport field!");
            break;
          }
        }

        // create parameter list for element evaluation
        Teuchos::ParameterList condparams;

        // action for elements
        condparams.set<int>("action",INPAR::S2I::evaluate_condition_od);

        // create strategy for assembly of auxiliary system matrix
        SCATRA::MortarCellAssemblyStrategy strategythermoscatras2i(
            slavematrix,
            INPAR::S2I::side_slave,
            INPAR::S2I::side_slave,
            slavematrix,
            INPAR::S2I::side_slave,
            INPAR::S2I::side_master,
            Teuchos::null,
            INPAR::S2I::side_undefined,
            INPAR::S2I::side_undefined,
            Teuchos::null,
            INPAR::S2I::side_undefined,
            INPAR::S2I::side_undefined,
            Teuchos::null,
            INPAR::S2I::side_undefined,
            Teuchos::null,
            INPAR::S2I::side_undefined,
            0,
            1
            );

        // extract scatra-scatra interface coupling conditions
        std::vector<DRT::Condition*> conditions;
        thermo_->Discretization()->GetCondition("S2ICoupling",conditions);

        // loop over all conditions
        for(unsigned icondition=0; icondition<conditions.size(); ++icondition)
        {
          // extract current condition
          DRT::Condition& condition = *conditions[icondition];

          // consider conditions for slave side only
          if(condition.GetInt("interface side") == INPAR::S2I::side_slave)
          {
            // add condition to parameter list
            condparams.set<DRT::Condition*>("condition",&condition);

            // evaluate mortar integration cells
            strategythermo_->EvaluateMortarCells(strategythermo_->MortarDiscretization(condition.GetInt("ConditionID")),condparams,strategythermoscatras2i);
          }
        }

        // finalize auxiliary system matrix
        switch(strategyscatra_->MatrixType())
        {
          case INPAR::S2I::matrix_block_condition:
          {
            slavematrix->Complete();
            break;
          }

          case INPAR::S2I::matrix_sparse:
          {
            slavematrix->Complete(*maps_->Map(0),*maps_->Map(1));
            break;
          }

          default:
          {
            dserror("Invalid matrix type associated with scalar transport field!");
            break;
          }
        }

        // assemble linearizations of slave-side thermo fluxes w.r.t. slave-side and master-side scatra dofs into thermo-scatra matrix block
        thermoscatrablock_->Add(*slavematrix,false,1.,1.);

        // linearizations of master-side thermo fluxes w.r.t. scatra dofs are not needed, since thermo fluxes are source terms and thus only evaluated once on slave side

        // standard meshtying algorithm with Lagrange multipliers condensed out
        // loop over all thermo-scatra matrix blocks
        for(int iblock=0; iblock<blockmaps_->NumMaps()-1; ++iblock)
        {
          // initialize temporary matrix for slave-side rows of current thermo-scatra matrix block
          LINALG::SparseMatrix thermoscatrarowsslave(*strategythermo_->InterfaceMaps()->Map(1),27,false,true);

          // extract current thermo-scatra matrix block
          LINALG::SparseMatrix& thermoscatrablock = strategyscatra_->MatrixType() == INPAR::S2I::matrix_block_condition ? Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(thermoscatrablock_)->Matrix(0,iblock) : *Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(thermoscatrablock_);

          // extract slave-side rows of thermo-scatra matrix block into temporary matrix
          SCATRA::MeshtyingStrategyS2I::ExtractMatrixRows(thermoscatrablock,thermoscatrarowsslave,*strategythermo_->InterfaceMaps()->Map(1));

          // finalize temporary matrix with slave-side rows of thermo-scatra matrix block
          thermoscatrarowsslave.Complete(*maps_->Map(0),*strategythermo_->InterfaceMaps()->Map(1));

          // add projected slave-side rows of thermo-scatra matrix block to corresponding master-side rows
          thermoscatrablock.Add(*LINALG::MLMultiply(*strategythermo_->P(),true,thermoscatrarowsslave,false,false,false,true),false,1.,1.);
        }

        break;
      }

      default:
      {
        dserror("Invalid type of scatra-scatra interface coupling!");
        break;
      }
    }
  }

  // finalize thermo-scatra matrix block
  switch(strategyscatra_->MatrixType())
  {
    case INPAR::S2I::matrix_block_condition:
    {
      thermoscatrablock_->Complete();
      break;
    }

    case INPAR::S2I::matrix_sparse:
    {
      thermoscatrablock_->Complete(*maps_->Map(0),*maps_->Map(1));
      break;
    }

    default:
    {
      dserror("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }

  // apply Dirichlet boundary conditions to scatra-thermo matrix block
  thermoscatrablock_->ApplyDirichlet(*thermo_->DirichMaps()->CondMap(),false);

  // zero out slave-side rows of thermo-scatra matrix block after having added them to the
  // corresponding master-side rows to finalize condensation of slave-side thermo dofs
  if(thermo_->S2ICoupling())
  {
    switch(strategythermo_->CouplingType())
    {
      case INPAR::S2I::coupling_matching_nodes:
      {
        if(!thermo_->Discretization()->GetCondition("PointCoupling"))
          thermoscatrablock_->ApplyDirichlet(*icoupthermo_->SlaveDofMap(),false);
        break;
      }

      case INPAR::S2I::coupling_mortar_condensed_bubnov:
      {
        thermoscatrablock_->ApplyDirichlet(*(strategythermo_->InterfaceMaps()->Map(1)),false);
        break;
      }

      default:
      {
        dserror("Invalid type of scatra-scatra interface coupling!");
        break;
      }
    }
  }

  // remove state vectors from thermo discretization
  thermo_->Discretization()->ClearState();

  return;
} // STI::Algorithm::AssembleODBlockThermoScatra()


/*----------------------------------------------------------------------*
 | time loop                                                 fang 04/15 |
 *----------------------------------------------------------------------*/
void STI::Algorithm::TimeLoop()
{
  // output initial solution to screen and files
  if(Step() == 0)
    Output();

  // time loop
  while(NotFinished())
  {
    // prepare time step
    PrepareTimeStep();

    // evaluate time step
    Solve();

    // update scatra and thermo fields
    Update();

    // output solution to screen and files
    Output();
  } // while(NotFinished())

  return;
} // STI::Algorithm::TimeLoop()


/*-------------------------------------------------------------------------*
 | update scatra and thermo fields after time step evaluation   fang 04/15 |
 *-------------------------------------------------------------------------*/
void STI::Algorithm::Update()
{
  // update scatra field
  scatra_->Update();

  // compare scatra field to analytical solution if applicable
  scatra_->EvaluateErrorComparedToAnalyticalSol();

  // update thermo field
  thermo_->Update();

  // compare thermo field to analytical solution if applicable
  thermo_->EvaluateErrorComparedToAnalyticalSol();

  return;
} // STI::Algorithm::Update()


/*---------------------------------------------------------------------------------------------*
 | finite difference check for global system matrix (for debugging only)            fang 07/15 |
 *---------------------------------------------------------------------------------------------*/
void STI::Algorithm::FDCheck()
{
  // initial screen output
  if(Comm().MyPID() == 0)
    std::cout << std::endl << "FINITE DIFFERENCE CHECK FOR STI SYSTEM MATRIX" << std::endl;

  // create global state vector
  Teuchos::RCP<Epetra_Vector> statenp(LINALG::CreateVector(*DofRowMap(),true));
  maps_->InsertVector(scatra_->Phinp(),0,statenp);
  maps_->InsertVector(thermo_->Phinp(),1,statenp);

  // make a copy of global state vector to undo perturbations later
  Teuchos::RCP<Epetra_Vector> statenp_original = Teuchos::rcp(new Epetra_Vector(*statenp));

  // make a copy of system matrix as Epetra_CrsMatrix
  Teuchos::RCP<Epetra_CrsMatrix> sysmat_original = Teuchos::null;
  if(Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(systemmatrix_) != Teuchos::null)
    sysmat_original = (new LINALG::SparseMatrix(*(Teuchos::rcp_static_cast<LINALG::BlockSparseMatrixBase>(systemmatrix_)->Merge())))->EpetraMatrix();
  else
    dserror("Global system matrix must be a block sparse matrix!");
  sysmat_original->FillComplete();

  // make a copy of system right-hand side vector
  Teuchos::RCP<Epetra_Vector> rhs_original = Teuchos::rcp(new Epetra_Vector(*residual_));

  // initialize counter for system matrix entries with failing finite difference check
  int counter(0);

  // initialize tracking variable for maximum absolute and relative errors
  double maxabserr(0.);
  double maxrelerr(0.);

  for (int colgid=0; colgid<=sysmat_original->ColMap().MaxAllGID(); ++colgid)
  {
    // check whether current column index is a valid global column index and continue loop if not
    int collid(sysmat_original->ColMap().LID(colgid));
    int maxcollid(-1);
    Comm().MaxAll(&collid,&maxcollid,1);
    if(maxcollid < 0)
      continue;

    // fill global state vector with original state variables
    statenp->Update(1.,*statenp_original,0.);

    // impose perturbation
    if(statenp->Map().MyGID(colgid))
      if(statenp->SumIntoGlobalValue(colgid,0,scatra_->FDCheckEps()))
        dserror("Perturbation could not be imposed on state vector for finite difference check!");
    scatra_->Phinp()->Update(1.,*maps_->ExtractVector(statenp,0),0.);
    thermo_->Phinp()->Update(1.,*maps_->ExtractVector(statenp,1),0.);

    // carry perturbation over to state vectors at intermediate time stages if necessary
    scatra_->ComputeIntermediateValues();
    thermo_->ComputeIntermediateValues();

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
      const double fdval = -(*residual_)[rowlid] / scatra_->FDCheckEps() + (*rhs_original)[rowlid] / scatra_->FDCheckEps();

      // confirm accuracy of first comparison
      if(abs(fdval) > 1.e-17 and abs(fdval) < 1.e-15)
        dserror("Finite difference check involves values too close to numerical zero!");

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
      if(abs(relerr1) > scatra_->FDCheckTol())
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
        const double left  = entry - (*rhs_original)[rowlid] / scatra_->FDCheckEps();

        // right-hand side in second comparison
        const double right = -(*residual_)[rowlid] / scatra_->FDCheckEps();

        // confirm accuracy of second comparison
        if(abs(right) > 1.e-17 and abs(right) < 1.e-15)
          dserror("Finite difference check involves values too close to numerical zero!");

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
        if(abs(relerr2) > scatra_->FDCheckTol())
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
      dserror("Finite difference check failed for STI system matrix!");
    }
    else
      printf("--> PASSED WITH MAXIMUM ABSOLUTE ERROR %+12.5e AND MAXIMUM RELATIVE ERROR %+12.5e\n\n",maxabserrglobal,maxrelerrglobal);
  }

  // undo perturbations of state variables
  scatra_->Phinp()->Update(1.,*maps_->ExtractVector(statenp_original,0),0.);
  scatra_->ComputeIntermediateValues();
  thermo_->Phinp()->Update(1.,*maps_->ExtractVector(statenp_original,1),0.);
  thermo_->ComputeIntermediateValues();

  // recompute system matrix and right-hand side vector based on original state variables
  AssembleMatAndRHS();

  return;
}

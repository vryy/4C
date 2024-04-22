/*----------------------------------------------------------------------*/
/*! \file

 \brief  monolithic structure split poroelasticity algorithm

\level 2

 *------------------------------------------------------------------------------------------------*/

#include "baci_poroelast_monolithicstructuresplit.hpp"

#include "baci_adapter_fld_poro.hpp"
#include "baci_adapter_str_fpsiwrapper.hpp"
#include "baci_coupling_adapter_converter.hpp"
#include "baci_fluid_utils_mapextractor.hpp"
#include "baci_linalg_matrixtransform.hpp"
#include "baci_linalg_utils_sparse_algebra_assemble.hpp"
#include "baci_linalg_utils_sparse_algebra_create.hpp"
#include "baci_structure_aux.hpp"

#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN


POROELAST::MonolithicStructureSplit::MonolithicStructureSplit(const Epetra_Comm& comm,
    const Teuchos::ParameterList& timeparams,
    Teuchos::RCP<CORE::LINALG::MapExtractor> porosity_splitter)
    : MonolithicSplit(comm, timeparams, porosity_splitter)
{
  sggtransform_ = Teuchos::rcp(new CORE::LINALG::MatrixRowColTransform);
  sgitransform_ = Teuchos::rcp(new CORE::LINALG::MatrixRowTransform);
  sigtransform_ = Teuchos::rcp(new CORE::LINALG::MatrixColTransform);
  csggtransform_ = Teuchos::rcp(new CORE::LINALG::MatrixRowTransform);
  cfggtransform_ = Teuchos::rcp(new CORE::LINALG::MatrixColTransform);
  csgitransform_ = Teuchos::rcp(new CORE::LINALG::MatrixRowTransform);
  cfigtransform_ = Teuchos::rcp(new CORE::LINALG::MatrixColTransform);

  // Recovering of Lagrange multiplier happens on structure field
  lambda_ = Teuchos::rcp(new Epetra_Vector(*StructureField()->Interface()->FSICondMap()));
}

void POROELAST::MonolithicStructureSplit::SetupSystem()
{
  {
    // create combined map
    std::vector<Teuchos::RCP<const Epetra_Map>> vecSpaces;

    vecSpaces.push_back(StructureField()->Interface()->OtherMap());
    vecSpaces.push_back(FluidField()->DofRowMap());

    if (vecSpaces[0]->NumGlobalElements() == 0) FOUR_C_THROW("No structure equation. Panic.");
    if (vecSpaces[1]->NumGlobalElements() == 0) FOUR_C_THROW("No fluid equation. Panic.");

    // full Poroelasticity-map
    fullmap_ = CORE::LINALG::MultiMapExtractor::MergeMaps(vecSpaces);
    // full Poroelasticity-blockmap
    blockrowdofmap_->Setup(*fullmap_, vecSpaces);
  }

  // Use splitted structure matrix
  StructureField()->UseBlockMatrix();

  SetupCouplingAndMatrices();

  BuildCombinedDBCMap();

  SetupEquilibration();
}

void POROELAST::MonolithicStructureSplit::SetupRHS(bool firstcall)
{
  TEUCHOS_FUNC_TIME_MONITOR("POROELAST::MonolithicStructureSplit::SetupRHS");

  // create full monolithic rhs vector
  rhs_ = Teuchos::rcp(new Epetra_Vector(*DofRowMap(), true));

  SetupVector(*rhs_, StructureField()->RHS(), FluidField()->RHS(), FluidField()->ResidualScaling());

  // store interface force onto the structure to know it in the next time step as previous force
  // in order to recover the Lagrange multiplier
  // fgpre_ = fgcur_;
  fgcur_ = StructureField()->Interface()->ExtractFSICondVector(StructureField()->RHS());
}

void POROELAST::MonolithicStructureSplit::SetupSystemMatrix(
    CORE::LINALG::BlockSparseMatrixBase& mat)
{
  TEUCHOS_FUNC_TIME_MONITOR("POROELAST::MonolithicStructureSplit::SetupSystemMatrix");

  Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> s = StructureField()->BlockSystemMatrix();
  if (s == Teuchos::null) FOUR_C_THROW("expect structure block matrix");
  Teuchos::RCP<CORE::LINALG::SparseMatrix> f = FluidField()->SystemMatrix();
  if (f == Teuchos::null) FOUR_C_THROW("expect fluid matrix");

  // just to play it safe ...
  mat.Zero();

  /*----------------------------------------------------------------------*/

  // build block matrix
  // The maps of the block matrix have to match the maps of the blocks we
  // insert here.

  /*----------------------------------------------------------------------*/
  // structural part k_sf (3nxn)
  // build mechanical-fluid block

  // create empty matrix
  Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> k_sf = StructFluidCouplingBlockMatrix();
  if (k_sf == Teuchos::null) FOUR_C_THROW("expect coupling block matrix");

  // call the element and calculate the matrix block
  ApplyStrCouplMatrix(k_sf);

  /*----------------------------------------------------------------------*/
  // fluid part k_fs ( (3n+1)x3n )
  // build fluid-mechanical block

  // create empty matrix
  Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> k_fs = FluidStructCouplingBlockMatrix();
  if (k_fs == Teuchos::null) FOUR_C_THROW("expect coupling block matrix");

  // call the element and calculate the matrix block
  ApplyFluidCouplMatrix(k_fs);

  /*----------------------------------------------------------------------*/

  k_fs->Complete();
  k_sf->Complete();

  f->UnComplete();

  /*----------------------------------------------------------------------*/
  if (evaluateinterface_)
  {
    double scale = FluidField()->ResidualScaling();
    double timescale = FluidField()->TimeScaling();

    // get time integration parameters of structure an fluid time integrators
    // to enable consistent time integration among the fields
    double stiparam = StructureField()->TimIntParam();
    double ftiparam = FluidField()->TimIntParam();

    (*sigtransform_)(s->FullRowMap(), s->FullColMap(), s->Matrix(0, 1), 1. / timescale,
        CORE::ADAPTER::CouplingMasterConverter(*icoupfs_),
        k_sf->Matrix(0, 1),  // k_sf->Matrix(0,1),mat.Matrix(0,1)
        true, true);

    (*sggtransform_)(s->Matrix(1, 1), (1.0 - ftiparam) / ((1.0 - stiparam) * scale * timescale),
        CORE::ADAPTER::CouplingMasterConverter(*icoupfs_),
        CORE::ADAPTER::CouplingMasterConverter(*icoupfs_), *f, true, true);

    (*sgitransform_)(s->Matrix(1, 0), (1.0 - ftiparam) / ((1.0 - stiparam) * scale),
        CORE::ADAPTER::CouplingMasterConverter(*icoupfs_),
        k_fs->Matrix(1, 0),  // k_fs->Matrix(1,0), mat.Matrix(1,0)
        true);

    (*cfggtransform_)(s->FullRowMap(),  // k_fs->FullRowMap(),
        s->FullColMap(),                // k_fs->FullColMap(),
        k_fs->Matrix(1, 1), 1. / timescale, CORE::ADAPTER::CouplingMasterConverter(*icoupfs_), *f,
        true, true);

    (*cfigtransform_)(s->FullRowMap(),  // k_fs->FullRowMap(),
        s->FullColMap(),                // k_fs->FullColMap(),
        k_fs->Matrix(0, 1), 1. / timescale, CORE::ADAPTER::CouplingMasterConverter(*icoupfs_), *f,
        true, true);

    (*csggtransform_)(k_sf->Matrix(1, 1), (1.0 - ftiparam) / ((1.0 - stiparam) * scale),
        CORE::ADAPTER::CouplingMasterConverter(*icoupfs_), *f, true);

    (*csgitransform_)(k_sf->Matrix(1, 0), (1.0 - ftiparam) / ((1.0 - stiparam) * scale),
        CORE::ADAPTER::CouplingMasterConverter(*icoupfs_), *f, true);
  }

  /*----------------------------------------------------------------------*/
  // pure structural part
  mat.Matrix(0, 0).Add(s->Matrix(0, 0), false, 1., 0.0);

  // structure coupling part
  mat.Matrix(0, 1).Add(k_sf->Matrix(0, 0), false, 1.0, 0.0);
  mat.Matrix(0, 1).Add(k_sf->Matrix(0, 1), false, 1.0, 1.0);

  // pure fluid part
  mat.Assign(1, 1, CORE::LINALG::View, *f);

  // fluid coupling part
  mat.Matrix(1, 0).Add(k_fs->Matrix(0, 0), false, 1.0, 0.0);
  mat.Matrix(1, 0).Add(k_fs->Matrix(1, 0), false, 1.0, 1.0);
  /*----------------------------------------------------------------------*/
  // done. make sure all blocks are filled.
  mat.Complete();

  sgicur_ = Teuchos::rcp(new CORE::LINALG::SparseMatrix(s->Matrix(1, 0)));
  sggcur_ = Teuchos::rcp(new CORE::LINALG::SparseMatrix(s->Matrix(1, 1)));
  cgicur_ = Teuchos::rcp(new CORE::LINALG::SparseMatrix(k_sf->Matrix(1, 0)));
  cggcur_ = Teuchos::rcp(new CORE::LINALG::SparseMatrix(k_sf->Matrix(1, 1)));
}

void POROELAST::MonolithicStructureSplit::SetupVector(Epetra_Vector& f,
    Teuchos::RCP<const Epetra_Vector> sv, Teuchos::RCP<const Epetra_Vector> fv, double fluidscale)
{
  // extract the inner and boundary dofs of all three fields

  Teuchos::RCP<Epetra_Vector> sov = StructureField()->Interface()->ExtractOtherVector(sv);

  if (fluidscale != 0.0)
  {
    // get time integration parameters of structure an fluid time integrators
    // to enable consistent time integration among the fields
    double stiparam = StructureField()->TimIntParam();
    double ftiparam = FluidField()->TimIntParam();

    // add fluid interface values to structure vector
    Teuchos::RCP<Epetra_Vector> scv = StructureField()->Interface()->ExtractFSICondVector(sv);
    Teuchos::RCP<Epetra_Vector> modfv =
        FluidField()->Interface()->InsertFSICondVector(StructureToFluidAtInterface(scv));
    modfv->Update(1.0, *fv, (1.0 - ftiparam) / ((1.0 - stiparam) * fluidscale));

    // add contribution of Lagrange multiplier from previous time step
    if (lambda_ != Teuchos::null)
      modfv->Update(-ftiparam + stiparam * (1.0 - ftiparam) / (1.0 - stiparam),
          *StructureToFluidAtInterface(lambda_), 1.0);

    Extractor()->InsertVector(*modfv, 1, f);
  }
  else
  {
    Extractor()->InsertVector(*fv, 1, f);
  }

  Extractor()->InsertVector(*sov, 0, f);
}

void POROELAST::MonolithicStructureSplit::ExtractFieldVectors(Teuchos::RCP<const Epetra_Vector> x,
    Teuchos::RCP<const Epetra_Vector>& sx, Teuchos::RCP<const Epetra_Vector>& fx, bool firstcall)
{
  TEUCHOS_FUNC_TIME_MONITOR("POROELAST::MonolithicStructureSplit::ExtractFieldVectors");

  // process fluid unknowns of the second field
  fx = Extractor()->ExtractVector(x, 1);

  // process structure unknowns
  if (evaluateinterface_)
  {
    Teuchos::RCP<Epetra_Vector> fcx = FluidField()->Interface()->ExtractFSICondVector(fx);

    {
      double timescale = 1. / FluidField()->TimeScaling();
      fcx->Scale(timescale);
    }

    Teuchos::RCP<Epetra_Vector> scx = FluidToStructureAtInterface(fcx);
    Teuchos::RCP<const Epetra_Vector> sox = Extractor()->ExtractVector(x, 0);

    Teuchos::RCP<Epetra_Vector> s = StructureField()->Interface()->InsertOtherVector(sox);
    StructureField()->Interface()->InsertFSICondVector(scx, s);

    auto zeros = Teuchos::rcp(new const Epetra_Vector(s->Map(), true));
    CORE::LINALG::ApplyDirichletToSystem(
        *s, *zeros, *(StructureField()->GetDBCMapExtractor()->CondMap()));

    sx = s;

    Teuchos::RCP<Epetra_Vector> fox = FluidField()->Interface()->ExtractOtherVector(fx);

    // Store field vectors to know them later on as previous quantities
    if (solipre_ != Teuchos::null)
      ddiinc_->Update(1.0, *sox, -1.0, *solipre_, 0.0);  // compute current iteration increment
    else
      ddiinc_ = Teuchos::rcp(new Epetra_Vector(*sox));  // first iteration increment

    solipre_ = sox;  // store current step increment

    if (solgpre_ != Teuchos::null)
      ddginc_->Update(1.0, *scx, -1.0, *solgpre_, 0.0);  // compute current iteration increment
    else
      ddginc_ = Teuchos::rcp(new Epetra_Vector(*scx));  // first iteration increment

    solgpre_ = scx;  // store current step increment

    if (solivelpre_ != Teuchos::null)
      duiinc_->Update(1.0, *fox, -1.0, *solivelpre_, 0.0);  // compute current iteration increment
    else
      duiinc_ = Teuchos::rcp(new Epetra_Vector(*fox));  // first iteration increment

    solivelpre_ = fox;  // store current step increment
  }
  else
    sx = Extractor()->ExtractVector(x, 0);
}

void POROELAST::MonolithicStructureSplit::RecoverLagrangeMultiplierAfterTimeStep()
{
  TEUCHOS_FUNC_TIME_MONITOR(
      "POROELAST::MonolithicStructureSplit::RecoverLagrangeMultiplierAfterTimeStep");
  if (evaluateinterface_)
  {
    // get time integration parameter of structural time integrator
    // to enable consistent time integration among the fields
    double stiparam = StructureField()->TimIntParam();
    double timescale = FluidField()->TimeScaling();

    // store the product S_{\GammaI} \Delta d_I^{n+1} in here
    Teuchos::RCP<Epetra_Vector> sgiddi =
        CORE::LINALG::CreateVector(*StructureField()->Interface()->FSICondMap(), true);
    // compute the above mentioned product
    sgicur_->Multiply(false, *ddiinc_, *sgiddi);

    // store the product C_{\GammaI} \Delta u_I^{n+1} in here
    Teuchos::RCP<Epetra_Vector> fgiddi =
        CORE::LINALG::CreateVector(*StructureField()->Interface()->FSICondMap(), true);
    // compute the above mentioned product
    cgicur_->Multiply(false, *duiinc_, *fgiddi);

    // store the product S_{\Gamma\Gamma} \Delta d_\Gamma^{n+1} in here
    Teuchos::RCP<Epetra_Vector> sggddg =
        CORE::LINALG::CreateVector(*StructureField()->Interface()->FSICondMap(), true);
    // compute the above mentioned product
    sggcur_->Multiply(false, *ddginc_, *sggddg);

    // store the prodcut C_{\Gamma\Gamma} \Delta u_\Gamma^{n+1} in here
    Teuchos::RCP<Epetra_Vector> cggddg =
        CORE::LINALG::CreateVector(*StructureField()->Interface()->FSICondMap(), true);
    // compute the above mentioned product
    cggcur_->Multiply(false, *ddginc_, *cggddg);
    cggddg->Scale(timescale);

    // Update the Lagrange multiplier:
    /* \lambda^{n+1} =  1/b * [ - a*\lambda^n - f_\Gamma^S
     *                          - S_{\Gamma I} \Delta d_I
     *                          - C_{\Gamma I} \Delta u_I
     *                          - S_{\Gamma\Gamma} \Delta d_\Gamma]
     *                          - C_{\Gamma\Gamma} * 2 / \Delta t * \Delta d_\Gamma]
     */
    // lambda_->Update(1.0, *fgpre_, -stiparam);
    lambda_->Update(1.0, *fgcur_, -stiparam);
    lambda_->Update(-1.0, *sgiddi, -1.0, *sggddg, 1.0);
    lambda_->Update(-1.0, *fgiddi, -1.0, *cggddg, 1.0);
    lambda_->Scale(1 / (1.0 - stiparam));  // entire Lagrange multiplier is divided by (1.-stiparam)
  }
}

FOUR_C_NAMESPACE_CLOSE
